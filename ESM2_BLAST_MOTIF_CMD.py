import os
import re
import time
import torch
import esm
import pandas as pd
from Bio import SeqIO
from torch.utils.data import DataLoader, Dataset
from pathlib import Path
import argparse
import sys
from huggingface_hub import hf_hub_download

# Load paths from envfile
def parse_envfile(envfile_path='envfile'):
    if not os.path.exists(envfile_path):
        print(f"Error: The environment file '{envfile_path}' is missing.", file=sys.stderr)
        sys.exit(1)

    paths = {}
    with open(envfile_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#') and ':' in line:
                key, value = line.split(':', 1)
                paths[key.strip()] = value.strip()
    return paths

env_paths = parse_envfile()
blastp_path = env_paths.get('BLAST')
blast_db_path = env_paths.get('BLAST database')
merci_script_path = env_paths.get('MERCI')
merci_motif_file = env_paths.get('MERCI motif file')

# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


# Download the model files from Hugging Face
repo_id = "raghavagps-group/antifp2"
classifier_path = hf_hub_download(repo_id=repo_id, filename="fine_tuned_classifier.pth")
alphabet_path = hf_hub_download(repo_id=repo_id, filename="esm_alphabet.pth")

torch.serialization.add_safe_globals([esm.data.Alphabet])
alphabet = torch.load(alphabet_path, map_location="cpu", weights_only=False)
batch_converter = alphabet.get_batch_converter()
model, _ = esm.pretrained.esm2_t36_3B_UR50D()

class ProteinClassifier(torch.nn.Module):
    def __init__(self, esm_model, embedding_dim=2560, num_classes=2):
        super().__init__()
        self.esm_model = esm_model
        self.fc = torch.nn.Linear(embedding_dim, num_classes)

    def forward(self, tokens):
        with torch.no_grad():
            results = self.esm_model(tokens, repr_layers=[36])
        embeddings = results["representations"][36].mean(1)
        logits = self.fc(embeddings)
        return logits, embeddings

torch.serialization.add_safe_globals({'ProteinClassifier': ProteinClassifier})
classifier = torch.load(classifier_path, map_location=device, weights_only=False)
classifier = classifier.to(device)
classifier.eval()
print("Model and classifier loaded")

# Dataset
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

# MERCI
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

# BLAST
def run_blast(val_fasta: Path, output_dir: Path):
    out_file = output_dir / f"blast_out.csv"
    if blastp_path and blast_db_path:
        blast_command = f"{blastp_path} -db {blast_db_path} -query {val_fasta} -out {out_file} -outfmt 6 -max_target_seqs 1 -num_threads 8 -evalue 0.001 -subject_besthit"
        print(f"Running BLAST:\n{blast_command}")
        os.system(blast_command)
    else:
        print("Error: BLAST paths are not properly configured in the envfile.", file=sys.stderr)
        sys.exit(1)
    return out_file

# Parse MERCI coverage section
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

# Adjust probabilities
def adjust_with_blast_and_motif(df, blast_file, motif_file):
    df["blast_adjustment"] = 0.0
    df["motif_adjustment"] = 0.0

    blast_hits = set()
    motif_hits = set()

    # Process BLAST hits
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
    else:
        print("BLAST file missing or empty — no BLAST adjustment applied.")

    # Process MERCI motif hits
    if motif_file and Path(motif_file).exists() and Path(motif_file).stat().st_size > 0:
        try:
            content = Path(motif_file).read_text()
            motif_hits = parse_coverage_section(content)
            df.loc[df["ID"].isin(motif_hits), "motif_adjustment"] = 0.5
        except Exception as e:
            print(f"Warning: Could not parse MERCI output: {e}")
    else:
        print("Motif file missing or empty — no motif adjustment applied.")

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


# Prediction pipeline
def predict_all(fasta_path: Path, threshold=0.5, outdir: Path = None, batch_size=4):
    combo_name = fasta_path.stem
    output_dir = outdir if outdir else fasta_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    motif_file = run_merci(fasta_path, output_dir)
    blast_file = run_blast(fasta_path, output_dir)

    dataset = ProteinDataset(fasta_path, alphabet.get_batch_converter())
    loader = DataLoader(dataset, batch_size=batch_size, collate_fn=dataset.collate_fn)

    ids, probs = [], []
    with torch.no_grad():
        for tokens, sids in loader:
            logits, _ = classifier(tokens)
            p = torch.softmax(logits, dim=1)[:, 1].cpu().numpy()
            ids.extend(sids)
            probs.extend(p)

    df = pd.DataFrame({"ID": ids, "probability": probs})
    df = adjust_with_blast_and_motif(df, blast_file, motif_file)
    df["prediction"] = (df["combined"] >= threshold).astype(int)

    output_csv = output_dir / f"{combo_name}_predictions.csv"
    df.to_csv(output_csv, index=False)
    print(f"\nPredictions saved to: {output_csv}")
    return df

# CLI
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict using fine-tuned ESM2 + BLAST + MERCI")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("--outdir", default=None, help="Optional output directory (default: same as FASTA)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Prediction threshold")
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir) if args.outdir else None

    predict_all(fasta_path, threshold=args.threshold, outdir=outdir)

