#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full prediction pipeline:
 - Validates sequences (length, amino acids)
 - Loads fine-tuned ESM2-t36 model from Hugging Face cache
 - Runs predictions
 - Adjusts with BLAST and MERCI
 - Writes final CSV
"""

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
from torch.utils.data import Dataset
from pathlib import Path
import argparse
from huggingface_hub import hf_hub_download

start_time = time.time()

# ----------------------------- Classifier ----------------------------- #
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
    elif 'darwin' in os_name:
        return f"{base_key}_macos"
    else:
        print(f"Unsupported OS: {os_name}", file=sys.stderr)
        sys.exit(1)

# ---------------------------- MERCI ---------------------------- #
def run_merci(val_fasta: Path, output_dir: Path, merci_script_path, merci_motif_file):
    out_locate = output_dir / f"val_pos.locate"
    if merci_script_path and merci_motif_file:
        merci_command = f"{merci_script_path} -p {val_fasta} -i {merci_motif_file} -o {out_locate} -c KOOLMAN-ROHM"
        print(f"Running MERCI:\n{merci_command}")
        os.system(merci_command)
    else:
        print("Error: MERCI paths not configured in envfile.", file=sys.stderr)
        sys.exit(1)
    return out_locate

# ---------------------------- BLAST ---------------------------- #
def run_blast(val_fasta: Path, output_dir: Path, blastp_path, blast_db_path):
    out_file = output_dir / f"blast_out.csv"
    if blastp_path and blast_db_path:
        blast_command = f"{blastp_path} -db {blast_db_path} -query {val_fasta} -out {out_file} -outfmt 6 -max_target_seqs 1 -num_threads 8 -evalue 0.001 -subject_besthit"
        print(f"Running BLAST:\n{blast_command}")
        os.system(blast_command)
        print(f"BLAST run complete. Output saved to {out_file}")
    else:
        print("Error: BLAST paths not configured in envfile.", file=sys.stderr)
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
            print("Warning: BLAST output empty. Skipping adjustment.")

    if motif_file and Path(motif_file).exists() and Path(motif_file).stat().st_size > 0:
        try:
            content = Path(motif_file).read_text()
            motif_hits = parse_coverage_section(content)
            df.loc[df["ID"].isin(motif_hits), "motif_adjustment"] = 0.5
        except Exception as e:
            print(f"Warning: Could not parse MERCI output: {e}")

    df["combined"] = df["probability"] + df["blast_adjustment"] + df["motif_adjustment"]
    df["combined"] = df["combined"].clip(0, 1)

    return df

# ------------------------ Main Predict Function ------------------------ #
def main():
    parser = argparse.ArgumentParser(description="Predict using fine-tuned ESM2 + BLAST + MERCI")
    parser.add_argument("--input", required=True, help="Path to input FASTA file")
    parser.add_argument("--outdir", default=None, help="Optional output directory (default: same as FASTA)")
    parser.add_argument("--threshold", type=float, default=0.5, help="Prediction threshold")
    parser.add_argument("--no-cleanup", action="store_true", help="Do not delete intermediate files")
    args = parser.parse_args()

    fasta_path = Path(args.input)
    output_dir = Path(args.outdir) if args.outdir else fasta_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    combo_name = fasta_path.stem
    output_csv = output_dir / f"{combo_name}_esm2_predictions.csv"

    # ---------------- Validate sequences ---------------- #
    rejected_log_path = output_dir / f"{combo_name}_rejected_log.txt"
    valid_fasta = output_dir / f"{combo_name}_valid.fasta"
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    valid_sequences = []

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
        raise ValueError("❌ No valid sequences found. See rejected_log.txt")

    # ---------------- BLAST + MERCI ---------------- #
    env_paths = parse_envfile()
    blast_key = get_os_specific_key('BLAST')
    blastp_path = env_paths.get(blast_key)
    blast_db_path = env_paths.get('BLAST_database')
    merci_script_path = env_paths.get('MERCI')
    merci_motif_file = env_paths.get('MERCI_motif_file')

    motif_file = run_merci(valid_fasta, output_dir, merci_script_path, merci_motif_file)
    blast_file = run_blast(valid_fasta, output_dir, blastp_path, blast_db_path)

    # ---------------- Load fine-tuned model ---------------- #
    repo_id = "raghavagps-group/antifp2"
    cache_dir = Path("~/.cache/huggingface/antifp2").expanduser()
    print(f"📦 Loading fine-tuned model from Hugging Face repo: {repo_id}")

    config_path = hf_hub_download(repo_id=repo_id, filename="config.json", cache_dir=cache_dir)
    weights_path = hf_hub_download(repo_id=repo_id, filename="pytorch_model.bin", cache_dir=cache_dir)
    alphabet_path = hf_hub_download(repo_id=repo_id, filename="alphabet.bin", cache_dir=cache_dir)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    with open(config_path, "r") as f:
        config = json.load(f)
    embedding_dim, num_classes = config["embedding_dim"], config["num_classes"]

    torch.serialization.add_safe_globals([esm.data.Alphabet])
    alphabet = torch.load(alphabet_path, map_location="cpu", weights_only=False)
    batch_converter = alphabet.get_batch_converter()

    esm_model, _ = esm.pretrained.esm2_t36_3B_UR50D()
    esm_model = esm_model.to(device)

    classifier = ProteinClassifier(esm_model, embedding_dim, num_classes).to(device)
    classifier.load_state_dict(torch.load(weights_path, map_location=device, weights_only=True))
    classifier.eval()

    print("✅ Model loaded. Beginning predictions...")

    # ---------------- Run predictions ---------------- #
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "probability"])  

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

    df_preds = pd.read_csv(output_csv)
    df_final = adjust_with_blast_and_motif(df_preds, blast_file, motif_file)
    df_final["prediction"] = (df_final["combined"] >= args.threshold).astype(int)

    final_output = output_dir / f"{combo_name}_predictions.csv"
    df_final.to_csv(final_output, index=False)

    print(f"\n✅ Predictions written to {final_output}")
    print(f"📄 Rejected sequences logged in: {rejected_log_path}")
    print(f"⏱️ Time taken: {time.time() - start_time:.2f} seconds")

    # Cleanup
    if not args.no_cleanup:
        for f in [valid_fasta, blast_file, motif_file, output_csv]:
            try:
                os.remove(f)
            except Exception:
                pass
        print("🧹 Cleaned up intermediate files")

# ------------------------ CLI ------------------------ #
if __name__ == "__main__":
    main()

