import os
import re
import time
import torch
import esm
import sys
import platform
import pandas as pd
import tempfile
import shutil
import csv
import json
from Bio import SeqIO
from torch.utils.data import DataLoader, Dataset
from pathlib import Path
import argparse
from huggingface_hub import hf_hub_download

start_time = time.time()

# ----------------------------- Classifier ----------------------------- #
class ProteinClassifier(torch.nn.Module):
    def __init__(self, esm_model, embedding_dim=2560, num_classes=2):  #t36=2560, t33=1280
        super().__init__()
        self.esm_model = esm_model
        self.fc = torch.nn.Linear(embedding_dim, num_classes)

    def forward(self, tokens):
        with torch.no_grad():
            results = self.esm_model(tokens, repr_layers=[36])
        embeddings = results["representations"][36].mean(1)
        logits = self.fc(embeddings)
        return logits, embeddings

# --------------------------- ENV Variables ---------------------------- #
def parse_envfile(envfile_path='envfile'):
    if not os.path.exists(envfile_path):
        print(f"Error: The environment file '{envfile_path}' is missing.", file=sys.stderr)
        sys.exit(1)

    paths = {}
    with open(envfile_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#') and '=' in line:
                key, value = line.split('=', 1)
                paths[key.strip()] = value.strip()
    return paths

def get_os_specific_key(base_key):
    os_name = platform.system().lower()
    if 'linux' in os_name:
        return f"{base_key}_ubuntu"
    elif 'windows' in os_name:
        return f"{base_key}_windows"
    elif 'darwin' in os_name:  # macOS
        return f"{base_key}_macos"
    else:
        print(f"Unsupported OS: {os_name}", file=sys.stderr)
        sys.exit(1)

# Load environment variables
env_paths = parse_envfile()

# Determine correct BLAST path based on OS
blast_key = get_os_specific_key('BLAST')
blastp_path = env_paths.get(blast_key)
if blastp_path is None:
    print(f"Error: Could not find path for {blast_key} in envfile.", file=sys.stderr)
    sys.exit(1)

# Get common paths
blast_db_path = env_paths.get('BLAST_database')
merci_script_path = env_paths.get('MERCI')
merci_motif_file = env_paths.get('MERCI_motif_file')

# Verify critical paths
if not blast_db_path or not merci_script_path or not merci_motif_file:
    print("Error: Missing one or more required paths (BLAST_database, MERCI, MERCI_motif_file) in envfile.", file=sys.stderr)
    sys.exit(1)

# Device setup
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# ---------------------------- Dataset ---------------------------- #
class ProteinDataset(Dataset):
    def __init__(self, fasta_file, batch_converter):
        self.data = [(r.id, str(r.seq)) for r in SeqIO.parse(fasta_file, "fasta")]
        if not self.data:
            raise ValueError(f"No sequences found in {fasta_file}")
        self.batch_converter = batch_converter

    def __len__(self): return len(self.data)
    def __getitem__(self, idx): return self.data[idx]
    def collate_fn(self, batch):
        _, _, tokens = self.batch_converter(batch)
        return tokens.to(device), [sid for sid, _ in batch]

# ---------------------------- MERCI ---------------------------- #
def run_merci(val_fasta: Path, output_dir: Path):
    out_locate = output_dir / f"val_pos.locate"
    if merci_script_path and merci_motif_file:
        merci_command = f"{merci_script_path} -p {val_fasta} -i {merci_motif_file} -o {out_locate} -c KOOLMAN-ROHM"
        print(f"Running MERCI:\n{merci_command}")
        os.system(merci_command)
    else:
        print("Error: MERCI paths are not properly configured in the envfile.", file=sys.stderr)
        sys.exit(1)
    return out_locate

# ---------------------------- BLAST ---------------------------- #
def run_blast(val_fasta: Path, output_dir: Path):
    out_file = output_dir / f"blast_out.csv"
    if blastp_path and blast_db_path:
        blast_command = f"{blastp_path} -db {blast_db_path} -query {val_fasta} -out {out_file} -outfmt 6 -max_target_seqs 1 -num_threads 8 -evalue 0.001 -subject_besthit"
        print(f"Running BLAST:\n{blast_command}")
        os.system(blast_command)
        print(f"BLAST run complete. Output saved to {out_file}")
    else:
        print("Error: BLAST paths are not properly configured in the envfile.", file=sys.stderr)
        sys.exit(1)
    return out_file

# ---------------------------- Adjust ---------------------------- #
def parse_coverage_section(file_content):
    hits = []
    match = re.search(r'COVERAGE\s*\n\*+\n(.*)', file_content, re.DOTALL)
    if match:
        lines = match.group(1).strip().splitlines()
        for line in lines:
            m = re.match(r"(\S+)\s+\((\d+)\s+motifs match\)", line.strip())
            if m:
                hits.append(m.group(1))
    return set(hits)

def adjust_with_blast_and_motif(df, blast_file, motif_file):
    df["blast_adjustment"] = 0.0
    df["motif_adjustment"] = 0.0

    blast_hits = set()
    motif_hits = set()

    if blast_file and Path(blast_file).exists() and Path(blast_file).stat().st_size > 0:
        try:
            blast_df = pd.read_csv(blast_file, sep="\t", header=None)
            for _, row in blast_df.iterrows():
                qid, sid = row[0], row[1]
                if qid in df["ID"].values:
                    if sid.endswith("_1"):
                        df.loc[df["ID"] == qid, "blast_adjustment"] = 0.5
                        blast_hits.add(qid)
                    elif sid.endswith("_0"):
                        df.loc[df["ID"] == qid, "blast_adjustment"] = -0.5
                        blast_hits.add(qid)
        except pd.errors.EmptyDataError:
            print("Warning: BLAST output is empty. Skipping BLAST adjustment.")

    if motif_file and Path(motif_file).exists() and Path(motif_file).stat().st_size > 0:
        try:
            content = Path(motif_file).read_text()
            motif_hits = parse_coverage_section(content)
            df.loc[df["ID"].isin(motif_hits), "motif_adjustment"] = 0.5
        except Exception as e:
            print(f"Warning: Could not parse MERCI output: {e}")

    df["combined"] = df["probability"] + df["blast_adjustment"] + df["motif_adjustment"]
    df["combined"] = df["combined"].clip(0, 1)

    # Logging
    total = len(df)
    blast_only = len(blast_hits - motif_hits)
    motif_only = len(motif_hits - blast_hits)
    both = len(blast_hits & motif_hits)
    none = total - len(blast_hits | motif_hits)

    print(f"Total sequences: {total}")
    print(f"Adjusted by BLAST only: {blast_only}")
    print(f"Adjusted by Motif only: {motif_only}")
    print(f"Adjusted by both BLAST and Motif: {both}")
    print(f"No adjustment (no hits): {none}")

    return df

# ------------------------ Main Predict Function ------------------------ #
def main():


    parser = argparse.ArgumentParser(description="Predict using fine-tuned ESM2 + BLAST + MERCI")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("--outdir", default=None, help="Optional output directory (default: same as FASTA)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Prediction threshold")
    parser.add_argument("--no-cleanup", action="store_true", help="Do not delete intermediate files")
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir) if args.outdir else None
    cleanup = not args.no_cleanup
    threshold=0.5
    
    combo_name = fasta_path.stem
    output_dir = outdir if outdir else fasta_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / f"esm2_predictions.csv"
       
    # --- Setup rejected log ---
    valid_fasta = output_dir / f"{combo_name}_valid.fasta"
    rejected_log_path = output_csv.parent / "rejected_log.txt"
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    valid_sequences = []

    # --- Filter sequences and write to valid FASTA ---
    with open(rejected_log_path, "w") as rejlog, open(valid_fasta, "w") as vfile:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            if len(seq) < 50:
                rejlog.write(f"{record.id}: Rejected (length < 50)\n")
                continue
            if len(seq) > 3000:
                rejlog.write(f"{record.id}: Rejected (length > 3000)\n")
                continue
            if any(aa not in valid_aas for aa in seq):
                rejlog.write(f"{record.id}: Rejected (non-standard amino acids)\n")
                continue
            SeqIO.write(record, vfile, "fasta")
            valid_sequences.append(record)

    if not valid_sequences:
        raise ValueError("❌ No valid sequences found. Check rejected_log.txt for details.")

    
    # Step 1: Run MERCI and BLAST
    motif_file = run_merci(valid_fasta, output_dir)
    blast_file = run_blast(valid_fasta, output_dir)
    
    
    # --- Load model and alphabet ---
    print("⏬ Downloading model files from Hugging Face (if not cached)...")
    repo_id = "raghavagps-group/antifp2"
    config_path = hf_hub_download(repo_id=repo_id, filename="config.json")
    weights_path = hf_hub_download(repo_id=repo_id, filename="pytorch_model.bin")
    alphabet_path = hf_hub_download(repo_id=repo_id, filename="alphabet.bin")

    with open(config_path, "r") as f:
        config = json.load(f)

    embedding_dim = config["embedding_dim"]
    num_classes = config["num_classes"]

    torch.serialization.add_safe_globals([esm.data.Alphabet])
    alphabet = torch.load(alphabet_path, map_location="cpu", weights_only=False)
    batch_converter = alphabet.get_batch_converter()

    esm_model, _ = esm.pretrained.esm2_t36_3B_UR50D()
    esm_model = esm_model.to(device)

    classifier = ProteinClassifier(esm_model, embedding_dim, num_classes).to(device)
    classifier.load_state_dict(torch.load(weights_path, map_location=device))
    classifier.eval()

    print("✅ Model loaded. Beginning predictions...")

    # --- Predict and write output ---
    tmp_dir = tempfile.mkdtemp()

    try:
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["ID", "probability"])  # header

            for record in valid_sequences:
                seq = str(record.seq).upper()
                batch = [(record.id, seq)]
                _, _, tokens = batch_converter(batch)
                tokens = tokens.to(device)

                with torch.no_grad():
                    logits, _ = classifier(tokens)
                    prob = torch.softmax(logits, dim=1)[0][1].item()

                writer.writerow([record.id, prob])
                csvfile.flush()
                print(f"✅ {record.id}: prob={prob:.4f}")
                
        
        df_preds= pd.read_csv(output_csv)
        
        df_final = adjust_with_blast_and_motif(df_preds, blast_file, motif_file)
        df_final["prediction"] = (df_final["combined"] >= threshold).astype(int)
        
        final_output = output_dir / f"{combo_name}_predictions.csv"
        df_final.to_csv(final_output, index=False)
        
        end_time = time.time()
        elapsed = end_time - start_time
        print(f"✅ Prediction complete in {elapsed:.2f} seconds.")
        print(f"\n✅ Predictions written to {output_csv}")
        print(f"📄 Rejected sequences logged in: {rejected_log_path}")

    finally:
        if cleanup:
            try:
                os.remove(valid_fasta)
                os.remove(blast_file)
                os.remove(motif_file)
                os.remove(output_csv)
                print("🧹 Intermediate files removed.")
            except Exception as e:
                print(f"Warning during cleanup: {e}")

# ------------------------ CLI ------------------------ #
if __name__ == "__main__":
    main()

