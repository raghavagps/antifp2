====================================
📘 AntiFP2: Antifungal Protein Prediction Toolkit
====================================

AntiFP2 is a comprehensive toolkit for **antifungal protein prediction** using multiple integrated approaches, combining **deep learning (ESM2)**, **machine learning (PAAC features)**, **sequence similarity (BLAST)**, and **motif analysis (MERCI)**.

It provides both **standalone scripts** and **Dockerized pipelines**, making it fully portable and easy to deploy across Linux and macOS systems.

.. image:: https://img.shields.io/badge/model-huggingface-blue.svg
   :target: https://huggingface.co/raghavagps-group/antifp2
.. image:: https://img.shields.io/badge/license-GPLv3-green.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

---

🔍 Overview
-----------

AntiFP2 predicts antifungal proteins using two complementary frameworks:

* **ESM2 deep learning models** fine-tuned on antifungal protein data.
* **PAAC-based ML models** (XGBoost, Random Forest, etc.) integrated with **BLAST** and **MERCI** for evidence-based refinement.

It supports three primary workflows:

1. **Standalone Predictions** (ESM2-only or ML-only)
2. **Hybrid Predictions** (ESM2 + BLAST + MERCI or ML + BLAST + MERCI)
3. **MetaPipelines** (Prokka + Prediction for contig-level metagenomes)

---

🧠 Core Functionalities
------------------------

* 🧬 **Deep learning** with fine-tuned ESM2 embeddings
* 📈 **Machine learning** using PAAC descriptors
* 🔍 **BLAST integration** for similarity-based adjustments
* 🧩 **MERCI motif scanning** for discriminative motif hits
* 🧫 **Prokka** for genome and metagenome annotation
* 🐳 **Docker support** for fully reproducible environments

---

⚙️ Installation
----------------

You can use **AntiFP2** in three main ways:

🔹 **Option 1 - Install dependencies using conda**  
(Recommended for metagenome or genome pipeline use)

.. code-block:: bash

   conda env create -f environment.yml
   conda activate antifp2

.. code-block:: bash

   git clone https://github.com/patrik-ackerman/antifp2.git
   cd antifp2

🔹 **Option 2 — Install via pip**  
(Recommended for one-by-one prediction use)

.. code-block:: bash

   pip install torch==2.4.1 biopython==1.85 pandas fair-esm==2.0.0 huggingface_hub==0.30.2

.. code-block:: bash

   pip install antifp2

This installs all Python dependencies and makes all scripts available under your Python environment.

🔹 **Option 3 — Run via Docker**  
(Recommended for Linux/macOS)

.. code-block:: bash

   docker pull pratik0297/antifp2

Run the container interactively:

.. code-block:: bash

   docker run -it --rm -v /absolute/path/to/data:/workspace -w /workspace pratik0297/antifp2 bash

You’ll now be inside the AntiFP2 environment, ready to execute any command below.

---

🧭 Usage Overview
------------------

AntiFP2 can be executed using a **master launcher script** or **individual prediction scripts**.

**Master Launcher (Unified Interface):**

.. code-block:: bash

   python3 antifp2.py --method <method-name> [options]

List available methods:

.. code-block:: bash

   python3 antifp2.py --list

**Available methods:**

+------------------+--------------------------------------+
| Method           | Description                          |
+==================+======================================+
| esm2             | ESM2-only prediction                 |
+------------------+--------------------------------------+
| esm2-hybrid      | ESM2 + BLAST + MERCI hybrid          |
+------------------+--------------------------------------+
| esm2-meta        | MetaPipeline (Prokka + ESM2)         |
+------------------+--------------------------------------+
| ml-hybrid        | ML (PAAC) + BLAST + MERCI hybrid     |
+------------------+--------------------------------------+
| ml-hybrid-meta   | MetaPipeline (Prokka + ML hybrid)    |
+------------------+--------------------------------------+

---

🧩 Detailed Script Descriptions
-------------------------------

1️⃣ **ESM2-Only Prediction (esm2)**

Predict antifungal proteins using a fine-tuned **ESM2** model.

.. code-block:: bash

   python3 antifp2.py --method esm2 --input proteins.fasta --output esm2_predictions.csv --threshold 0.5

**Arguments:**

* ``--input``: Input protein FASTA file  
* ``--output``: Output CSV file for results  
* ``--threshold``: Classification cutoff (default: 0.5)  
* ``--no-cleanup``: Keep intermediate files  

**Output CSV Columns:**

+---------------+------------------------------+
| Column        | Description                  |
+===============+==============================+
| ID            | Sequence identifier           |
+---------------+------------------------------+
| probability   | Model-predicted probability   |
+---------------+------------------------------+
| prediction    | Binary label (1/0)            |
+---------------+------------------------------+

---

2️⃣ **ESM2 + BLAST + MERCI Hybrid (esm2-hybrid)**

Combines deep learning with BLAST and MERCI evidence.

.. code-block:: bash

   python3 antifp2.py --method esm2-hybrid --input proteins.fasta --outdir results --threshold 0.6

**Arguments:**

* ``--input``: Input protein FASTA file  
* ``--outdir``: Directory to store results  
* ``--threshold``: Probability cutoff (default: 0.5)  
* ``--no-cleanup``: Retain intermediate BLAST/MERCI files  

**Adjustments Applied:**

* +0.5 → if BLAST hit to known positive  
* −0.5 → if BLAST hit to known negative  
* +0.5 → if MERCI motif hit  

**Outputs:**

* ``*_predictions.csv`` — final probabilities and predictions  
* ``blast_out.csv`` — raw BLAST results  
* ``motifs.locate`` — motif hits  
* ``rejected_log.txt`` — invalid sequence log  

---

3️⃣ **MetaPipeline (Prokka + ESM2) (esm2-meta)**

Annotate contigs using **Prokka**, then predict antifungal proteins via **ESM2**.

.. code-block:: bash

   python3 antifp2.py --method esm2-meta --input contigs.fasta --outdir esm2_meta_out --threads 8 --metagenome

**Arguments:**

* ``--input``: Contigs FASTA file  
* ``--outdir``: Output directory  
* ``--threads``: Threads for Prokka  
* ``--metagenome``: Enable Prokka metagenome mode  

**Workflow:**

1. Run Prokka annotation  
2. Extract predicted proteins  
3. Filter sequences (50–3000 aa)  
4. Predict antifungal proteins using ESM2  

**Outputs:**

* ``*_metapred.csv`` — prediction table  
* ``*_antifp2.fasta`` — antifungal sequences  
* ``prokka/`` — annotation directory  

---

4️⃣ **ML Hybrid (PAAC + BLAST + MERCI) (ml-hybrid)**

Predict using PAAC-based ML classifier (XGBoost), then refine via BLAST and MERCI.

.. code-block:: bash

   python3 antifp2.py --method ml-hybrid --input proteins.fasta --outdir ml_hybrid_out --threshold 0.5

**Arguments:**

* ``--input``: Protein FASTA file  
* ``--outdir``: Output directory  
* ``--threshold``: Cutoff (default: 0.5)  
* ``--envfile``: Path to configuration (DB paths, etc.)  

**Outputs:**

* ``*_predictions.csv`` — final predictions  
* ``blast_out.csv`` and ``motifs.locate`` — evidence data  

---

5️⃣ **ML Hybrid MetaPipeline (Prokka + ML Hybrid) (ml-hybrid-meta)**

Annotate genomes using **Prokka**, then run **ML Hybrid** prediction.

.. code-block:: bash

   python3 antifp2.py --method ml-hybrid-meta --input contigs.fasta --outdir ml_meta_out --threads 8 --threshold 0.5

**Arguments:**

* ``--input``: Contigs FASTA file  
* ``--outdir``: Output directory  
* ``--threads``: CPU threads for Prokka  
* ``--envfile``: Configuration for BLAST and MERCI paths  
* ``--metagenome``: Enable metagenome mode  

**Outputs:**

* ``*_predictions.csv`` — prediction results  
* ``*_antifp2.fasta`` — antifungal protein FASTA  
* ``prokka/`` — annotation files  

---

🧪 Example of Running Docker Full Workflow (Linux)
--------------------------------

Assume your data is in ``/home/user/antifp2_data``.

.. code-block:: bash

   docker run -it --rm -v /home/user/antifp2_data:/workspace -w /workspace pratik0297/antifp2 bash

Then, inside the container:

.. code-block:: bash

   python3 antifp2.py --method ml-hybrid --input proteins.fasta --outdir results

Outputs will be saved in ``results/`` both inside the container and on your host system.

---

🍏 Running Docker on macOS (Apple Silicon)
-----------------------------------

AntiFP2 runs perfectly on macOS via Docker emulation.

1. Install Docker Desktop:  
   https://www.docker.com/products/docker-desktop

2. Enable Rosetta 2 emulation:

   .. code-block:: bash

      softwareupdate --install-rosetta

3. Run the container:

   .. code-block:: bash

      docker run --platform linux/amd64 -it --rm -v ~/Downloads/data:/workspace -w /workspace pratik0297/antifp2 bash

4. Inside container:

   .. code-block:: bash

      python3 antifp2.py --method esm2-meta --input /workspace/contigs.fasta --outdir /workspace/output

---

🧹 Cleanup and Debugging
-------------------------

Use ``--no-cleanup`` to retain intermediate files for inspection (BLAST/MERCI outputs).  
Without this flag, temporary files are removed automatically.

---

📖 Citation
------------

If you use this tool, please cite:

**ESM2 Model:**
  Lin, Z., Akin, H., Rao, R., et al. (2023). *Nature*, 601(7891), 277–284.  
  DOI: https://doi.org/10.1038/s41586-021-03819-2

**Prokka:**
  Seemann, T. (2014). *Bioinformatics*, 30(14), 2068–2069.  
  DOI: https://doi.org/10.1093/bioinformatics/btu153

**BLAST:**
  Camacho, C., et al. (2009). *BMC Bioinformatics*, 10, 421.  
  DOI: https://doi.org/10.1186/1471-2105-10-421

**MERCI:**
  Vens, C., Rosso, M.N., & Danchin, E.G.J. (2011). *Bioinformatics*, 27(9), 1231–1238.  
  DOI: https://doi.org/10.1093/bioinformatics/btr110

**SeqKit:**
  Shen, W., Le, S., Li, Y., & Hu, F. (2016). *PLoS ONE*, 11(10), e0163962.  
  DOI: https://doi.org/10.1371/journal.pone.0163962

---

👨‍💻 Support
-------------

For bug reports or feature requests:

➡️ GitHub: https://github.com/raghavagps/antifp2  
➡️ Email: raghava@iiitd.ac.in
