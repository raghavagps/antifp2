====================================
AntiFP2:Antifungal Protein Prediction Model
====================================

This toolkit uses a fine-tuned ESM2 model to classify protein sequences, with optional integration of BLAST and MERCI motif scanning. It provides two standalone prediction scripts:

1. **ESM2-only** classification using embeddings and a trained classifier
2. **ESM2 + BLAST + MERCI** for combined evidence-based prediction

|huggingface| |license|

.. |huggingface| image:: https://img.shields.io/badge/model-huggingface-blue.svg
    :alt: Model on Hugging Face
    :scale: 100%
    :target: https://huggingface.co/raghavagps-group/antifp2

.. |license| image:: https://img.shields.io/badge/license-MIT-green.svg
    :alt: License
    :scale: 100%
    :target: https://opensource.org/licenses/MIT

Overview
========

These scripts classify protein sequences into functional categories using a fine-tuned ESM2 model. Optionally, predictions can be enhanced by incorporating:

- **BLAST results** (from precomputed TSV)
- **MERCI motif hits** (from `.locate` files)

Installation
============

Install dependencies:

.. code-block:: sh

    pip install torch biopython pandas esm huggingface_hub

Clone this repository:

.. code-block:: sh

    git clone https://github.com/YOUR_USERNAME/esm2_classifier.git
    cd esm2_classifier

(Optional) Create a virtual environment:

.. code-block:: sh

    python -m venv venv
    source venv/bin/activate

Dependencies
------------

- Python ≥ 3.7
- PyTorch (with CUDA support recommended)
- Biopython
- pandas
- ESM (from Facebook Research)
- Hugging Face Transformers
- BLAST+ (optional)
- MERCI motif tool (optional)

Scripts and Usage
=================

1. ESM2-Only Prediction
-----------------------

.. code-block:: sh

    python predict_esm2.py --fasta input.fasta --output esm2_results.csv

Arguments:

- ``--fasta``: Input FASTA file
- ``--output``: Output CSV file
- ``--threshold``: (Optional) Classification threshold (default: 0.5)

Output: A CSV with ``sequence_id``, ``probability``, and binary ``prediction`` column.

2. ESM2 + BLAST + MERCI Prediction
----------------------------------

.. code-block:: sh

    python predict_esm2_blast_merci.py \
        --fasta input.fasta \
        --output combo_results.csv \
        --blast blast_results.tsv \
        --merci motif_results.locate \
        --threshold 0.5

Arguments:

- ``--fasta``: Input FASTA file
- ``--output``: Output CSV file
- ``--blast``: BLAST output in TSV format (qseqid sseqid ...)
- ``--merci``: MERCI motif `.locate` file
- ``--threshold``: (Optional) Base threshold for classifier (default: 0.5)

Output: A CSV with ESM2 score, adjusted score based on BLAST and MERCI hits, and final prediction.

Scoring Rules
=============

- **ESM2 probability** is first calculated using the fine-tuned classifier.
- **BLAST match** to a known positive adds +0.5 to the probability.
- **MERCI motif hit** adds +0.5 to the probability.
- Final prediction = 1 if adjusted probability ≥ threshold, else 0.

Model Details
=============

- Architecture: ESM2 (esm2_t36_3B_UR50D)
- Fine-tuned on a curated dataset for binary classification.
- Model available on Hugging Face:
  `raghavagps-group/antifp2 <https://huggingface.co/raghavagps-group/antifp2>`_

Citation
========

If you use this tool, please cite:

