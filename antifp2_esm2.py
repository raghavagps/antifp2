#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Predict each sequence one-by-one using the fine-tuned ESM model
and write results incrementally to the output CSV. Sequences that are
too short, too long, or contain non-standard amino acids are logged.
"""

import torch
import esm
from pathlib import Path
from Bio import SeqIO
import time
import tempfile
import shutil
import csv
import argparse
import json
from huggingface_hub import hf_hub_download


class ProteinClassifier(torch.nn.Module):
    def __init__(self, esm_model, embedding_dim, num_classes):
        super(ProteinClassifier, self).__init__()
        self.esm_model = esm_model
        self.fc = torch.nn.Linear(embedding_dim, num_classes)

    def forward(self, tokens):
        with torch.no_grad():
            results = self.esm_model(tokens, repr_layers=[36])
        embeddings = results["representations"][36].mean(1)
        logits = self.fc(embeddings)
        return logits, embeddings


def run_prediction(fasta_input_path, output_csv_path, threshold=0.5, cleanup=True):
    start_time = time.time()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"✅ Using device: {device}")

    fasta_input = Path(fasta_input_path)
    output_csv = Path(output_csv_path)

    # Setup rejected log
    rejected_log_path = output_csv.parent / "rejected_log.txt"
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    valid_sequences = []

    with open(rejected_log_path, "w") as rejlog:
        for record in SeqIO.parse(fasta_input, "fasta"):
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
            valid_sequences.append(record)

    if not valid_sequences:
        raise ValueError("❌ No valid sequences found. Check rejected_log.txt for details.")

    # -------------------------------
    # Load your fine-tuned files via cache
    # -------------------------------
    repo_id = "raghavagps-group/antifp2"
    cache_dir = Path("~/.cache/huggingface/antifp2").expanduser()

    print(f"📦 Loading model files via Hugging Face cache: {cache_dir}")

    config_path = hf_hub_download(repo_id=repo_id, filename="config.json", cache_dir=cache_dir)
    weights_path = hf_hub_download(repo_id=repo_id, filename="pytorch_model.bin", cache_dir=cache_dir)
    alphabet_path = hf_hub_download(repo_id=repo_id, filename="alphabet.bin", cache_dir=cache_dir)

    # Load config
    with open(config_path, "r") as f:
        config = json.load(f)

    embedding_dim = config["embedding_dim"]
    num_classes = config["num_classes"]

    # Load alphabet + batch converter
    torch.serialization.add_safe_globals([esm.data.Alphabet])
    alphabet = torch.load(alphabet_path, map_location="cpu", weights_only=False)
    batch_converter = alphabet.get_batch_converter()

    # Load pretrained backbone
    esm_model, _ = esm.pretrained.esm2_t36_3B_UR50D()
    esm_model = esm_model.to(device)

    # Build classifier and load fine-tuned head
    classifier = ProteinClassifier(esm_model, embedding_dim, num_classes).to(device)
    classifier.load_state_dict(torch.load(weights_path, map_location=device, weights_only=True))
    classifier.eval()

    print("✅ Model loaded. Beginning predictions...")

    tmp_dir = tempfile.mkdtemp()

    try:
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["ID", "probability", "prediction"])  # header

            for record in valid_sequences:
                seq = str(record.seq).upper()
                batch = [(record.id, seq)]
                _, _, tokens = batch_converter(batch)
                tokens = tokens.to(device)

                with torch.no_grad():
                    logits, _ = classifier(tokens)
                    prob = torch.softmax(logits, dim=1)[0][1].item()
                    pred = int(prob >= threshold)

                writer.writerow([record.id, prob, pred])
                csvfile.flush()
                print(f"✅ {record.id}: prob={prob:.4f}, pred={pred}")

        elapsed_time = time.time() - start_time
        print(f"\n✅ Predictions written to {output_csv}")
        print(f"📄 Rejected sequences logged in: {rejected_log_path}")
        print(f"⏱️ Time taken: {elapsed_time:.2f} seconds")

    finally:
        if cleanup:
            shutil.rmtree(tmp_dir)
            print("🧹 Cleaned up temporary files")
        else:
            print(f"⚠️ Temporary files retained at: {tmp_dir}")


def main():
    parser = argparse.ArgumentParser(description="Predict with fine-tuned ESM2-t36 model (one-by-one mode)")
    parser.add_argument("--input", type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--output", type=str, required=True, help="Path to output CSV file")
    parser.add_argument("--threshold", type=float, default=0.5, help="Prediction threshold (default: 0.5)")
    parser.add_argument("--no-cleanup", action="store_true", help="Do not delete intermediate files")
    args = parser.parse_args()

    run_prediction(
        fasta_input_path=args.input,
        output_csv_path=args.output,
        threshold=args.threshold,
        cleanup=not args.no_cleanup
    )


if __name__ == "__main__":
    main()

