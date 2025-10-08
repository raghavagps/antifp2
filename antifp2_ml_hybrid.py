#!/usr/bin/env python3
import os
import re
import argparse
import subprocess
import time
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import joblib   # sklearn model
import platform
import sys

start_time = time.time()

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

# ---------------------------- MERCI ----------------------------------- #
def run_merci(val_fasta: Path, output_dir: Path, merci_script_path: str, merci_motif_file: str):
    out_locate = output_dir / f"{val_fasta.stem}_pos.locate"
    if merci_script_path and merci_motif_file:
        merci_command = f"{merci_script_path} -p {val_fasta} -i {merci_motif_file} -o {out_locate} -c KOOLMAN-ROHM"
        print(f"🔎 Running MERCI:\n{merci_command}")
        os.system(merci_command)
    else:
        print("Error: MERCI paths are not properly configured in the envfile.", file=sys.stderr)
        sys.exit(1)
    return out_locate

# ---------------------------- BLAST ----------------------------------- #
def run_blast(val_fasta: Path, output_dir: Path, blastp_path: str, blast_db_path: str):
    out_file = output_dir / f"{val_fasta.stem}_blast_out.csv"
    if blastp_path and blast_db_path:
        blast_command = (
            f"{blastp_path} -db {blast_db_path} -query {val_fasta} -out {out_file} "
            f"-outfmt 6 -max_target_seqs 1 -num_threads 8 -evalue 0.001 -subject_besthit"
        )
        print(f"🧬 Running BLAST:\n{blast_command}")
        os.system(blast_command)
        print(f"BLAST run complete. Output saved to {out_file}")
    else:
        print("Error: BLAST paths are not properly configured in the envfile.", file=sys.stderr)
        sys.exit(1)
    return out_file

# ---------------------------- Adjust ---------------------------------- #
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
                    if str(sid).endswith("_1"):
                        df.loc[df["ID"] == qid, "blast_adjustment"] = 0.5
                        blast_hits.add(qid)
                    elif str(sid).endswith("_0"):
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

    df["combined"] = (df["probability"]
                      + df["blast_adjustment"]
                      + df["motif_adjustment"]).clip(0, 1)
    return df

# ------------------------- Feature Extraction ------------------------- #
def run_pfeature(fasta_file: Path, out_csv: Path):
    cmd = f"pfeature_comp -i {fasta_file} -o {out_csv} -j PAAC -s 4"
    print(f"🔎 Running pfeature:\n{cmd}")
    subprocess.run(cmd, shell=True, check=True)

    if not out_csv.exists():
        raise FileNotFoundError(f"❌ Feature extraction failed, {out_csv} not found")

    ids = [rec.id for rec in SeqIO.parse(fasta_file, "fasta")]
    df = pd.read_csv(out_csv)

    # Ensure SampleName column exists
    if "SampleName" not in df.columns:
        if len(ids) != len(df):
            raise ValueError(f"❌ Mismatch: {len(ids)} sequences vs {len(df)} feature rows")
        df.insert(0, "SampleName", ids)
        df.to_csv(out_csv, index=False)

    return out_csv

# ------------------------- Prediction Core ---------------------------- #
def predict_sequences_ml(fasta_file: Path, model_path: Path, outdir: Path):
    features_csv = outdir / f"{fasta_file.stem}_PAAC.csv"
    run_pfeature(fasta_file, features_csv)

    features_df = pd.read_csv(features_csv)
    if "SampleName" in features_df.columns:
        ids = features_df["SampleName"].tolist()
        X = features_df.drop(columns=["SampleName"])
    else:
        raise ValueError("❌ pfeature output missing SampleName column")

    print(f"📦 Loading ML model from: {model_path}")
    model = joblib.load(model_path)

    probs = model.predict_proba(X)[:, 1]

    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
    df = pd.DataFrame({"ID": ids, "sequence": [seqs[i] for i in ids], "probability": probs})
    return df

# ------------------------------- Main --------------------------------- #
def main():
    parser = argparse.ArgumentParser(description="Protein Prediction: ML (PAAC) + BLAST + MERCI Antifungal Prediction")
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument("--outdir", required=True, help="Directory to save outputs")
    parser.add_argument("--threshold", type=float, default=0.5, help="Classification threshold (default: 0.5)")
    parser.add_argument("--envfile", default="envfile", help="Path to envfile with BLAST/MERCI config (default: envfile)")
    args = parser.parse_args()

    proteins_fasta = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_output_csv = outdir / f"{proteins_fasta.stem}_raw.csv"
    final_output_csv = outdir / f"{proteins_fasta.stem}_combined.csv"
    model_path = Path("./antifp2_xgboost_PAAC.pkl")

    # Load envfile
    env_paths = parse_envfile(args.envfile)
    blast_key = get_os_specific_key('BLAST')
    blastp_path = env_paths.get(blast_key)
    if blastp_path is None:
        print(f"Error: Could not find path for {blast_key} in envfile.", file=sys.stderr)
        sys.exit(1)
    blast_db_path = env_paths.get('BLAST_database')
    merci_script_path = env_paths.get('MERCI')
    merci_motif_file = env_paths.get('MERCI_motif_file')

    if not blast_db_path or not merci_script_path or not merci_motif_file:
        print("Error: Missing one or more required paths (BLAST_database, MERCI, MERCI_motif_file) in envfile.", file=sys.stderr)
        sys.exit(1)

    try:
        # 1) Run MERCI & BLAST
        motif_file = run_merci(proteins_fasta, outdir, merci_script_path, merci_motif_file)
        blast_file = run_blast(proteins_fasta, outdir, blastp_path, blast_db_path)

        # 2) Predict with ML model
        df_preds = predict_sequences_ml(proteins_fasta, model_path, outdir)
        df_preds.to_csv(raw_output_csv, index=False)

        # 3) Adjust with BLAST + MERCI
        df_final = adjust_with_blast_and_motif(df_preds.copy(), blast_file, motif_file)
        df_final["prediction"] = (df_final["combined"] >= args.threshold).astype(int)
        df_final.to_csv(final_output_csv, index=False)

        print(f"✅ Raw predictions: {raw_output_csv}")
        print(f"✅ Final predictions (with adjustments): {final_output_csv}")
        print(f"⏱️ Done in {time.time() - start_time:.2f} seconds")
        
        # ---------------- Cleanup ---------------- #
        stem = proteins_fasta.stem
        cleanup_files = [
            outdir / f"{stem}_PAAC.csv",
            outdir / f"{stem}_blast_out.csv",
            outdir / f"{stem}_pos.locate",
            outdir / f"{stem}_raw.csv",
        ]
        for f in cleanup_files:
            try:
                if f.exists():
                    f.unlink()
                    print(f"🧹 Deleted intermediate file: {f}")
            except Exception as e:
                print(f"⚠️ Could not delete {f}: {e}")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

