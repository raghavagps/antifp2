#!/usr/bin/env python3
import os
import re
import argparse
import subprocess
import tempfile
import shutil
import time
from pathlib import Path
from Bio import SeqIO
import csv
import json
import multiprocessing
import platform
import sys
import pandas as pd
import joblib   # sklearn model

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

# ------------------------------ Prokka -------------------------------- #
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

# ---------------------------- Filtering ------------------------------- #
def filter_proteins(fasta_in: Path, fasta_out: Path, rejected_log_path: Path):
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    valid_records = []
    with open(rejected_log_path, "w") as rejlog, open(fasta_out, "w") as vfile:
        for record in SeqIO.parse(fasta_in, "fasta"):
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
            valid_records.append(record)
    if not valid_records:
        raise ValueError("❌ No valid sequences found after filtering.")
    return valid_records

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

    # Logging
    total = len(df)
    blast_only = len(blast_hits - motif_hits)
    motif_only = len(motif_hits - blast_hits)
    both = len(blast_hits & motif_hits)
    none = total - len(blast_hits | motif_hits)

    print(f"📊 Total sequences: {total}")
    print(f"   Adjusted by BLAST only: {blast_only}")
    print(f"   Adjusted by Motif only: {motif_only}")
    print(f"   Adjusted by both BLAST and Motif: {both}")
    print(f"   No adjustment (no hits): {none}")

    return df

# ------------------------- Feature Extraction ------------------------- #
def run_pfeature(fasta_file: Path, out_csv: Path):
    """
    Run pfeature_comp to extract PAAC features and ensure SampleName column is added.
    """
    seqkit_cmd1 = f"seqkit replace -p ' .*' -r '' {fasta_file} > new_input.fasta"
    seqkit_cmd2 = f"seqkit seq -w 0 new_input.fasta > input_singleline.fasta"
    cmd = f"pfeature_comp -i input_singleline.fasta -o {out_csv} -j PAAC -s 4"

    print(f"🔎 Running pfeature:\n{cmd}")
    subprocess.run(seqkit_cmd1, shell=True, check=True)
    subprocess.run(seqkit_cmd2, shell=True, check=True)
    subprocess.run(cmd, shell=True, check=True)

    if not out_csv.exists():
        raise FileNotFoundError(f"❌ Feature extraction failed, {out_csv} not found")

    # ✅ Fix missing SampleName column
    from Bio import SeqIO
    import pandas as pd

    ids = [rec.id for rec in SeqIO.parse("input_singleline.fasta", "fasta")]
    df = pd.read_csv(out_csv)

    # Only add if missing
    if "SampleName" not in df.columns:
        if len(ids) != len(df):
            raise ValueError(f"❌ Mismatch: {len(ids)} sequences vs {len(df)} feature rows")
        df.insert(0, "SampleName", ids)
        df.to_csv(out_csv, index=False)

    return out_csv

# ------------------------- Prediction Core ---------------------------- #
def predict_sequences_ml(fasta_file: Path, model_path: Path, outdir: Path, threshold=0.5):
    """
    Extract PAAC features with pfeature and predict with ML model.
    """
    features_csv = outdir / f"{fasta_file.stem}_PAAC.csv"
    run_pfeature(fasta_file, features_csv)

    # Load features
    features_df = pd.read_csv(features_csv)
    if "SampleName" in features_df.columns:
        ids = features_df["SampleName"].tolist()
        X = features_df.drop(columns=["SampleName"])
    else:
        raise ValueError("❌ pfeature output missing SampleName column")

    # Load model
    print(f"📦 Loading ML model from: {model_path}")
    model = joblib.load(model_path)

    probs = model.predict_proba(X)[:, 1]
    df = pd.DataFrame({"ID": ids, "probability": probs})
    return df

# ------------------------ Extract Positives ---------------------------- #
def extract_positive_sequences(faa_file: Path, positive_ids, output_fasta: Path):
    id_file = f"{output_fasta}.ids.txt"
    with open(id_file, "w") as f:
        for pid in positive_ids:
            f.write(f"{pid}\n")
    cmd = ["seqkit", "grep", "-f", id_file, "-o", str(output_fasta), str(faa_file)]
    subprocess.run(cmd, check=True)
    os.remove(id_file)
    print(f"🧬 Saved {len(positive_ids)} antifungal sequences to {output_fasta}")

# ------------------------------- Main --------------------------------- #
def main():
    parser = argparse.ArgumentParser(description="MetaPipeline: Prokka + ML (PAAC) + BLAST + MERCI Antifungal Prediction")
    parser.add_argument("--input", required=True, help="Input contigs FASTA")
    parser.add_argument("--outdir", required=True, help="Directory to save all outputs")
    parser.add_argument("--threshold", type=float, default=0.5, help="Classification threshold on combined score (default: 0.5)")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Threads for Prokka")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep intermediate files (Prokka, temps)")
    parser.add_argument("--metagenome", action="store_true", help="Enable Prokka metagenome mode")
    parser.add_argument("--envfile", default="envfile", help="Path to envfile with BLAST/MERCI config (default: envfile)")
    args = parser.parse_args()

    # Resolve paths
    contigs_path = Path(args.input)
    contig_basename = contigs_path.stem
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prokka_dir = outdir / "prokka"
    raw_output_csv = outdir / f"{contig_basename}_metapred_raw.csv"
    final_output_csv = outdir / f"{contig_basename}_meta_predictions.csv"
    output_fasta = outdir / f"{contig_basename}_antifp2.fasta"
    valid_fasta = outdir / f"{contig_basename}_valid.fasta"
    rejected_log_path = outdir / "rejected_log.txt"

    model_path = Path("./antifp2_xgboost_PAAC.pkl")

    # Load envfile for BLAST/MERCI
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

    temp_dir = tempfile.mkdtemp()
    to_cleanup = []

    try:
        # 1) Annotate contigs → proteins
        faa_file = run_prokka(contigs_path, prokka_dir, threads=args.threads, metagenome=args.metagenome)

        # 2) Filter proteins and write valid FASTA
        valid_records = filter_proteins(faa_file, valid_fasta, rejected_log_path)

        # 3) Run MERCI & BLAST on valid set
        motif_file = run_merci(valid_fasta, outdir, merci_script_path, merci_motif_file)
        blast_file = run_blast(valid_fasta, outdir, blastp_path, blast_db_path)
        to_cleanup.extend([motif_file, blast_file])

        # 4) Predict with ML model on PAAC features
        df_preds = predict_sequences_ml(valid_fasta, model_path, outdir, threshold=args.threshold)
        #df_preds.to_csv(raw_output_csv, index=False)

        # 5) Adjust with BLAST + MERCI
        df_final = adjust_with_blast_and_motif(df_preds.copy(), blast_file, motif_file)
        df_final["prediction"] = (df_final["combined"] >= args.threshold).astype(int)
        df_final.to_csv(final_output_csv, index=False)

        # 6) Extract positives based on combined score
        positive_ids = df_final.loc[df_final["prediction"] == 1, "ID"].tolist()
        if positive_ids:
            extract_positive_sequences(faa_file, positive_ids, output_fasta)
        else:
            print("ℹ️ No positive sequences detected at the given threshold.")

        print(f"✅ Raw predictions: {raw_output_csv}")
        print(f"✅ Final predictions (with adjustments): {final_output_csv}")
        print(f"📄 Rejected sequences log: {rejected_log_path}")
        print(f"⏱️ Done in {time.time() - start_time:.2f} seconds")

        # 7) Cleanup
        if not args.no_cleanup:
            try:
                # Remove Prokka directory
                shutil.rmtree(prokka_dir, ignore_errors=True)

                # Remove intermediate FASTAs
                if valid_fasta.exists():
                    os.remove(valid_fasta)

                # Remove MERCI / BLAST outputs
                for p in to_cleanup:
                    if p and Path(p).exists():
                        os.remove(p)

                # Remove raw predictions
                if raw_output_csv.exists():
                    os.remove(raw_output_csv)

                # Remove rejected log
                if rejected_log_path.exists():
                    os.remove(rejected_log_path)

                # Remove PAAC features
                paac_file = outdir / f"{contig_basename}_valid_PAAC.csv"
                if paac_file.exists():
                    os.remove(paac_file)

                # Remove seqkit temp FASTAs
                for fname in ["new_input.fasta", "input_singleline.fasta"]:
                    fpath = Path(fname)
                    if fpath.exists():
                        os.remove(fpath)

                print("🧹 Cleaned up intermediate files. Only final results kept.")

            except Exception as e:
                print(f"Warning during cleanup: {e}")


    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()

