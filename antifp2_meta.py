#!/usr/bin/env python3

import os
import argparse
import subprocess
import tempfile
import shutil
import time
from pathlib import Path
from Bio import SeqIO
import torch
import esm
import csv
import json
import multiprocessing
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


def run_prokka(contigs_path, prokka_dir, threads, metagenome=False):
    print(f"🚀 Running Prokka with {threads} threads...")
    cmd = [
        "prokka",
        "--outdir", str(prokka_dir),
        "--prefix", "annotation",
        "--force",
        "--kingdom", "Bacteria",
        "--norrna",
        "--notrna",
        "--quiet",
        "--cpus", str(threads),
        str(contigs_path)
    ]
    if metagenome:
        cmd.append("--metagenome")
    subprocess.run(cmd, check=True)
    faa_file = prokka_dir / "annotation.faa"
    if not faa_file.exists():
        raise FileNotFoundError(f"❌ Prokka failed: {faa_file} not found.")
    print(f"✅ Prokka done: {faa_file}")
    return faa_file


def predict_sequences(fasta_file, output_csv, threshold=0.5):
    start_time = time.time()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"✅ Using device: {device}")

    rejected_log_path = Path(output_csv).parent / "rejected_log.txt"
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    valid_sequences = []

    with open(rejected_log_path, "w") as rejlog:
        for record in SeqIO.parse(fasta_file, "fasta"):
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
        raise ValueError("❌ No valid sequences found after filtering.")

    print("⏬ Loading model from Hugging Face...")
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
    classifier.load_state_dict(torch.load(weights_path, map_location=device, weights_only=True))
    classifier.eval()

    print("🔬 Running predictions...")
    predicted_ids = []

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "probability", "prediction"])

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
            if pred == 1:
                predicted_ids.append(record.id)

            print(f"✅ {record.id}: prob={prob:.4f}, pred={pred}")

    print(f"✅ Results saved to {output_csv}")
    print(f"📄 Rejected log: {rejected_log_path}")
    print(f"⏱️ Done in {time.time() - start_time:.2f} seconds")

    return predicted_ids


def extract_positive_sequences(faa_file, positive_ids, output_fasta):
    id_file = f"{output_fasta}.ids.txt"
    with open(id_file, "w") as f:
        for pid in positive_ids:
            f.write(f"{pid}\n")
    cmd = ["seqkit", "grep", "-f", id_file, "-o", output_fasta, str(faa_file)]
    subprocess.run(cmd, check=True)
    os.remove(id_file)
    print(f"🧬 Saved {len(positive_ids)} antifungal sequences to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(description="MetaPipeline: Prokka + ESM2 Antifungal Prediction")
    parser.add_argument("--contigs", required=True, help="Input contigs FASTA")
    parser.add_argument("--outdir", required=True, help="Directory to save all outputs")
    parser.add_argument("--threshold", type=float, default=0.5, help="Classification threshold (default: 0.5)")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Threads for Prokka")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep intermediate files (Prokka)")
    parser.add_argument("--metagenome", action="store_true", help="Enable Prokka metagenome mode")

    contigs_path = Path(args.contigs)
    contig_basename = contigs_path.stem
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prokka_dir = outdir / "prokka"
    output_csv = outdir / f"{contig_basename}_metapred.csv"
    output_fasta = outdir / f"{contig_basename}_antifp2.fasta"

    temp_dir = tempfile.mkdtemp()
    try:
        faa_file = run_prokka(contigs_path, prokka_dir, threads=args.threads, metagenome=args.metagenome)
        positive_ids = predict_sequences(faa_file, output_csv, threshold=args.threshold)
        extract_positive_sequences(faa_file, positive_ids, output_fasta)

        if not args.no_cleanup:
            shutil.rmtree(prokka_dir, ignore_errors=True)
            print("🧹 Cleaned up intermediate Prokka output.")

    except Exception as e:
        print(f"❌ Error: {e}")
        exit(1)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
