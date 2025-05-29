
# AntiFP2

**AntiFP2** is a tool for the prediction of antifungal proteins using a fine-tuned [ESM2](https://github.com/facebookresearch/esm) language model, optionally enhanced by post-prediction adjustment with **BLAST** and **MERCI** motif detection.

This pipeline combines deep learning-based embeddings with classical bioinformatics methods for improved reliability in antifungal protein prediction.

---

## 🚀 Features

- Fine-tuned [ESM2-t36_3B_UR50D](https://huggingface.co/models) model for antifungal prediction
- Post-prediction adjustment using:
  - **BLAST**: Sequence similarity matching to known antifungal/negative examples
  - **MERCI**: Motif Enrichment Recognition to enhance biological relevance
- One-by-one and batch prediction modes
- Rejection logging for low-quality or invalid sequences
- Hugging Face integration for model loading

---

## 📦 Installation

```bash
pip install git+https://github.com/patrik-ackerman/antifp2.git
````

> **Note**: Requires Python ≥ 3.12

Ensure that BLAST+ and MERCI binaries are properly configured via the `envfile` as shown below.

---

## 📁 Project Structure

```
antifp2/
│
├── python_scripts/
│   ├── antifp2_ESM2.py      # Main pipeline with ESM2 + BLAST + MERCI
│   ├── antifp2_BLAST.py     # ESM2-only one-by-one predictor
│   └── envfile              # Config file for paths to BLAST and MERCI tools
│
├── MERCI/                   # MERCI motif files
├── blast_db/                # Preformatted BLAST database
├── README.md
├── setup.py
└── ...
```

---

## 🧪 Usage

### 🔮 ESM2-Only Prediction (one-by-one)

```bash
antifp2_blast --fasta path/to/input.fasta --output results.csv
```

* Output will include:

  * `ID`, `probability`, `prediction` columns
* Logs invalid sequences to `rejected_log.txt`

### 🧬 Full Pipeline (ESM2 + BLAST + MERCI)

```bash
antifp2_esm --fasta path/to/input.fasta --output ./output_dir/
```

* Performs predictions
* Runs BLAST against provided database
* Executes MERCI with motif file
* Adjusts predictions and saves final output to:

  * `output_dir/<input>.adjusted.csv`

Optional flag:

```bash
--no-cleanup   # Retains intermediate files like raw BLAST output, logs, etc.
```

---

## 🔧 Configuring Environment

The tool reads environment-specific paths from a file named `envfile`. Example format:

```ini
# Path settings for different OS
BLAST_ubuntu=/usr/bin/blastp
BLAST_windows=C:/Program Files/NCBI/blastp.exe
BLAST_macos=/usr/local/bin/blastp

BLAST_database=antifp2/blast_db/antifungal_db
MERCI=antifp2/MERCI/merci
MERCI_motif_file=antifp2/MERCI/motifs.motif
```

Make sure this file is located in `antifp2/python_scripts/envfile`.

---

## 📋 Output Format

### `adjusted.csv` columns:

| Column            | Description                           |
| ----------------- | ------------------------------------- |
| ID                | Sequence ID from FASTA                |
| probability       | Raw ESM2-based antifungal probability |
| blast\_adjustment | Adjustment based on BLAST hit         |
| motif\_adjustment | Adjustment based on MERCI hit         |
| combined          | Final adjusted probability            |
| prediction        | 1 if `combined` ≥ 0.5, else 0         |

---

## 💾 Model Files

Downloaded automatically from Hugging Face:

* `config.json`
* `pytorch_model.bin`
* `alphabet.bin`

Repo: [raghavagps-group/antifp2](https://huggingface.co/raghavagps-group/antifp2)

---

## 📝 License

This project is licensed under the terms of the **MIT License**. See the `LICENSE.txt` file for details.

---

## 👨‍🔬 Author

**Pratik Shinde**
Indian Institute of Information Technology Delhi
[Email](mailto:pratiks@iiitd.ac.in)

---

## 🌐 Links

* 🔗 GitHub: [https://github.com/patrik-ackerman/antifp2](https://github.com/patrik-ackerman/antifp2)
* 🤗 Hugging Face Model: [raghavagps-group/antifp2](https://huggingface.co/raghavagps-group/antifp2)

```
