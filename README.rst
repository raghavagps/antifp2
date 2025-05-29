====================================
AntiFP2: Antifungal Protein Prediction Toolkit
====================================

This toolkit utilizes a fine-tuned ESM2 model to classify protein sequences, with optional integration of BLAST and MERCI motif scanning. It offers three standalone prediction scripts:

1. **ESM2-only** classification using embeddings and a trained classifier.
2. **ESM2 + BLAST + MERCI** for combined evidence-based prediction.
3. **MetaPipeline**: Annotate contigs with Prokka and predict antifungal proteins using ESM2.

.. image:: https://img.shields.io/badge/model-huggingface-blue.svg
   :target: https://huggingface.co/raghavagps-group/antifp2
.. image:: https://img.shields.io/badge/license-GPLv3-green.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

Overview
--------

These scripts classify protein sequences into functional categories using a fine-tuned ESM2 model. Optionally, predictions can be enhanced by incorporating:

- **BLAST results** (from precomputed TSV files).
- **MERCI motif hits** (from `.locate` files).

Installation
------------

Install dependencies using pip:

.. code-block:: bash

    pip install torch==2.4.1 biopython==1.85 pandas fair-esm==2.0.0 huggingface_hub==0.30.2

Or using Conda:

.. code-block:: bash

    conda install -c conda-forge -c bioconda -c defaults torch==2.4.1 biopython==1.85 pandas fair-esm==2.0.0 huggingface_hub==0.30.2

Clone this repository:

.. code-block:: bash

    git clone https://github.com/patrik-ackerman/antifp2.git
    cd antifp2

(Optional) Create a virtual environment:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate antifp2

Scripts and Usage
-----------------

1. ESM2-Only Prediction
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python3 antifp2_ESM2.py --fasta input.fasta --output esm2_results.csv

**Arguments:**

- ``--fasta``: Input FASTA file containing protein sequences.
- ``--output``: Output CSV file to save predictions.
- ``--threshold``: (Optional) Classification threshold (default: 0.5).

**Output:**

A CSV file with the following columns:

- ``ID``: Sequence identifier.
- ``probability``: Predicted probability of being an antifungal protein.
- ``prediction``: Binary prediction (1 for antifungal, 0 for non-antifungal) based on the threshold.

**Help:**

.. code-block:: bash

    python3 antifp2_ESM2.py --help

2. ESM2 + BLAST + MERCI Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python3 antifp2_ESM_blast.py \
        --fasta input.fasta \
        --outdir output_directory \
        --threshold 0.5 \
        --no-cleanup

**Arguments:**

- ``--fasta``: Input FASTA file containing protein sequences.
- ``--outdir``: Output directory to save predictions and intermediate files (optional; defaults to input FASTA directory).
- ``--threshold``: Classification threshold (default: 0.5).
- ``--no-cleanup``: Flag to keep intermediate files like BLAST and MERCI outputs.

**Description:**

This script runs the fine-tuned ESM2 model to predict antifungal proteins, then adjusts predictions based on BLAST and MERCI motif hits:

- Filters sequences to keep lengths between 50 and 3000 and only standard amino acids.
- Runs MERCI motif scanning and BLASTp search against a configured database.
- Adjusts probabilities by adding +0.5 for BLAST matches to known positives and +0.5 for motif hits.
- Clips combined probabilities between 0 and 1.
- Outputs a CSV file with columns: ``ID``, ``probability``, ``blast_adjustment``, ``motif_adjustment``, ``combined``, ``prediction``.

**Outputs:**

- ``<prefix>_predictions.csv``: Final prediction CSV file with adjusted probabilities and binary predictions.
- ``rejected_log.txt``: Log of sequences rejected during filtering.

**Example:**

.. code-block:: bash

    python antifp2_ESM_blast.py --fasta proteins.fasta --outdir results --threshold 0.6


3. MetaPipeline: Prokka + ESM2 Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python3 meta_pipeline.py \
        --contigs AM09.contigs.fa \
        --outdir results_dir \
        --threshold 0.5 \
        --threads 8

**Arguments:**

- ``--contigs``: Input contigs FASTA file.
- ``--outdir``: Output directory to save prediction results and intermediate files.
- ``--threshold``: (Optional) Classification threshold (default: 0.5).
- ``--threads``: (Optional) Number of threads to use for Prokka (default: all available).
- ``--no-cleanup``: (Optional) Retain Prokka intermediate files.
- ``--metagenome``: (Optional) Enable Prokka's metagenome mode.

**Workflow:**

1. **Prokka Annotation**: Annotates input contigs to predict coding sequences.
2. **Sequence Filtering**: Filters out proteins:
   - Shorter than 50 or longer than 3000 amino acids.
   - Containing non-standard amino acids.
3. **Prediction**: Runs the fine-tuned ESM2 model on valid sequences.
4. **Extraction**: Saves predicted antifungal sequences to a FASTA file.

**Outputs (in ``--outdir``):**

- ``*_metapred.csv``: CSV file with prediction results:
  - ``ID``, ``probability``, ``prediction``.
- ``*_antifp2.fasta``: FASTA file of positively predicted antifungal proteins.
- ``rejected_log.txt``: Log of sequences excluded during filtering.

**Help:**

.. code-block:: bash

    python3 meta_pipeline.py --help


Citation
--------

If you use this tool, please cite the following resources:

- **ESM2 Model**:

  Lin, Z., Akin, H., Rao, R., et al. (2023). Language models of protein sequences at the scale of evolution enable accurate structure prediction. *Nature*, 601(7891), 277–284.  
  DOI: https://doi.org/10.1038/s41586-021-03819-2

- **Prokka**:

  Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(14), 2068–2069.  
  DOI: https://doi.org/10.1093/bioinformatics/btu153  
  PMID: https://pubmed.ncbi.nlm.nih.gov/24642063/

- **BLAST**:

  Camacho, C., Coulouris, G., Avagyan, V., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.  
  DOI: https://doi.org/10.1186/1471-2105-10-421

- **MERCI**:

  Vens, C., Rosso, M.N., & Danchin, E.G.J. (2011). Identifying discriminative classification-based motifs in biological sequences. *Bioinformatics*, 27(9), 1231–1238.  
  DOI: https://doi.org/10.1093/bioinformatics/btr110

- **SeqKit**:

  Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE*, 11(10), e0163962.  
  DOI: https://doi.org/10.1371/journal.pone.0163962

