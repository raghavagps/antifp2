#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

# ----------------- CONFIG ----------------- #
SCRIPT_MAP = {
    "esm2": {
        "script": "antifp2_esm2.py",
        "desc": "Predict with fine-tuned ESM2-t36 model (one-by-one mode)"
    },
    "esm2-hybrid": {
        "script": "antifp2_esm2_hybrid.py",
        "desc": "Predict using fine-tuned ESM2 + BLAST + MERCI"
    },
    "esm2-meta": {
        "script": "antifp2_esm2_meta.py",
        "desc": "MetaPipeline: Prokka + ESM2 Antifungal Prediction"
    },
    "ml-hybrid": {
        "script": "antifp2_ml_hybrid.py",
        "desc": "Protein Prediction: ML (PAAC) + BLAST + MERCI Antifungal Prediction"
    },
    "ml-hybrid-meta": {
        "script": "antifp2_ml_hybrid_meta.py",
        "desc": "MetaPipeline: Prokka + ML (PAAC) + BLAST + MERCI Antifungal Prediction"
    },
}

def list_methods():
    print("\nAvailable methods:\n")
    for key, val in SCRIPT_MAP.items():
        print(f"  {key:15} {val['desc']}")
    print("\nExample usage:")
    print("  python antifp2_master.py --method esm2 --fasta input.fasta --output result.csv")
    print("  python antifp2_master.py --method ml-hybrid --proteins proteins.fasta --outdir results/ --model model.pkl\n")

def main():
    parser = argparse.ArgumentParser(
        description="AntiFP2 Master Launcher (Unified Entry Point)",
        add_help=False   # disable default -h/--help
    )
    parser.add_argument(
        "--method", choices=SCRIPT_MAP.keys(),
        help="Which prediction method to run"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available methods and descriptions"
    )

    args, unknown = parser.parse_known_args()

    if args.list:
        list_methods()
        sys.exit(0)

    if not args.method:
        parser.print_help()
        sys.exit(1)

    # Resolve script path
    script_name = SCRIPT_MAP[args.method]["script"]
    script_path = Path(__file__).resolve().parent / script_name

    if not script_path.exists():
        sys.exit(f"❌ Error: Script not found: {script_path}")

    # Build command: forward all unknown args (including -h/--help)
    cmd = [sys.executable, str(script_path)] + unknown

    try:
        result = subprocess.run(cmd)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        sys.exit("\n❌ Interrupted by user")

if __name__ == "__main__":
    main()

