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

    git clone https://github.com/YOUR_USERNAME/antifp2.git
    cd antifp2

(Optional) Create a virtual environment:

.. code-block:: bash

    python -m venv venv
    source venv/bin/activate

Scripts and Usage
-----------------

1. ESM2-Only Prediction
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python predict_esm2.py --fasta input.fasta --output esm2_results.csv

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

    python predict_esm2.py --help

2. ESM2 + BLAST + MERCI Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python predict_esm2_blast_merci.py \
        --fasta input.fasta \
        --output combo_results.csv \
        --blast blast_results.tsv \
        --merci motif_results.locate \
        --threshold 0.5

**Arguments:**

- ``--fasta``: Input FASTA file containing protein sequences.
- ``--output``: Output CSV file to save predictions.
- ``--blast``: BLAST output in TSV format (qseqid sseqid ...).
- ``--merci``: MERCI motif `.locate` file.
- ``--threshold``: (Optional) Base threshold for classifier (default: 0.5).

**Output:**

A CSV file with the following columns:

- ``ID``: Sequence identifier.
- ``esm2_probability``: Predicted probability from ESM2 model.
- ``adjusted_probability``: Probability adjusted based on BLAST and MERCI hits.
- ``prediction``: Final binary prediction (1 for antifungal, 0 for non-antifungal).

**Scoring Rules:**

- **ESM2 probability** is first calculated using the fine-tuned classifier.
- **BLAST match** to a known positive adds +0.5 to the probability.
- **MERCI motif hit** adds +0.5 to the probability.
- **Final prediction** = 1 if adjusted probability ≥ threshold, else 0.

**Help:**

.. code-block:: bash

    python predict_esm2_blast_merci.py --help

3. MetaPipeline: Prokka + ESM2 Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    python meta_pipeline.py \
        --contigs AM09.contigs.fa \
        --output metapred.csv \
        --prokka_dir tmp_prokka \
        --threshold 0.5

**Arguments:**

- ``--contigs``: Input contigs FASTA file.
- ``--output``: Output CSV file to save predictions.
- ``--prokka_dir``: Directory to store Prokka annotations.
- ``--threshold``: (Optional) Classification threshold (default: 0.5).

**Workflow:**

1. **Prokka Annotation**: Annotates contigs to predict protein-coding sequences.
2. **Sequence Filtering**: Removes sequences with non-standard amino acids and those shorter than 50 amino acids.
3. **ESM2 Prediction**: Classifies the filtered protein sequences using the fine-tuned ESM2 model.

**Output:**

A CSV file with the following columns:

- ``ID``: Sequence identifier.
- ``probability``: Predicted probability of being an antifungal protein.
- ``prediction``: Binary prediction (1 for antifungal, 0 for non-antifungal) based on the threshold.

**Help:**

.. code-block:: bash

    python meta_pipeline.py --help

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


