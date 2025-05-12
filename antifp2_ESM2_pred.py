#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 22:33:19 2025
@author: pratik

Run prediction using fine-tuned ESM model on single or multiple sequences.
Tracks time taken for prediction.
"""

import torch
import esm
import pandas as pd
from pathlib import Path
from torch.utils.data import DataLoader, Dataset
from Bio import SeqIO
import argparse
import time
from huggingface_hub import hf_hub_download

start_time = time.time()
# -------------------- Config --------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"✅ Using device: {device}")


# Download the model files from Hugging Face
repo_id = "raghavagps-group/antifp2"
classifier_path = hf_hub_download(repo_id=repo_id, filename="fine_tuned_classifier.pth")
alphabet_path = hf_hub_download(repo_id=repo_id, filename="esm_alphabet.pth")

# -------------------- Load Model --------------------
torch.serialization.add_safe_globals([esm.data.Alphabet])

alphabet = torch.load(alphabet_path, map_location="cpu", weights_only=False)
batch_converter = alphabet.get_batch_converter()
model, _ = esm.pretrained.esm2_t36_3B_UR50D()

# -------------------- Run --------------------
if __name__ == "__main__":
    

    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Predict with fine-tuned ESM model")
    parser.add_argument("--fasta", type=str,  required=True, help="Path to input FASTA file")
    parser.add_argument("--output", type=str, required=True, help="Path to output CSV file")
    parser.add_argument("--threshold", type=float, default=0.5, help="Prediction threshold (default: 0.5)")
    args = parser.parse_args()

    fasta_input = Path(args.fasta)
    output_csv = Path(args.output)
    threshold = args.threshold

    # --- Define Classes ---
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

    class ProteinDataset(Dataset):
        def __init__(self, fasta_file, batch_converter):
            self.data = []
            self.seq_ids = []
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.data.append((record.id, str(record.seq)))
                self.seq_ids.append(record.id)
            if not self.data:
                raise ValueError(f"No sequences found in {fasta_file}")
            self.batch_converter = batch_converter

        def __len__(self):
            return len(self.data)

        def __getitem__(self, idx):
            return self.data[idx]

        def collate_fn(self, batch):
            _, _, tokens = self.batch_converter(batch)
            return tokens, [seq_id for seq_id, _ in batch]

    # --- Load Fine-tuned Classifier ---
    torch.serialization.add_safe_globals({'ProteinClassifier': ProteinClassifier})
    classifier = torch.load(classifier_path, map_location=device, weights_only=False)
    classifier = classifier.to(device)
    classifier.eval()

    print("✅ Fine-tuned model and alphabet loaded!")

    # --- Predict Function ---
    def predict_with_classifier(fasta_file, output_csv, threshold=0.5, batch_size=4):
        dataset = ProteinDataset(fasta_file, batch_converter)
        dataloader = DataLoader(dataset, batch_size=batch_size, collate_fn=dataset.collate_fn)

        all_probs = []
        all_ids = []

        with torch.no_grad():
            for tokens, seq_ids in dataloader:
                tokens = tokens.to(device)
                logits, _ = classifier(tokens)
                probs = torch.softmax(logits, dim=1)[:, 1].cpu().numpy()
                all_probs.extend(probs)
                all_ids.extend(seq_ids)

        df = pd.DataFrame({"ID": all_ids, "probability": all_probs})
        df["prediction"] = (df["probability"] >= threshold).astype(int)
        df.to_csv(output_csv, index=False)

        elapsed_time = time.time() - start_time
        print(f"✅ Prediction results saved to {output_csv}")
        print(f"⏱️ Time taken for prediction: {elapsed_time:.2f} seconds")

    # --- Run Prediction ---
    predict_with_classifier(fasta_input, output_csv, threshold=threshold)
