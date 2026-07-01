import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

fasta_file = "pep.aligned.mafft.fsa"
nodes_dir = "Models/out"

# Load nodes files
nodes_files = os.listdir(nodes_dir)
nodes_files = [f for f in nodes_files if f.endswith("_ringNodes")]

def find_nodes_file(seq_id):
    c1 = seq_id.lower() + ".cif_ringNodes"
    if c1 in nodes_files: return c1
    c2 = seq_id.replace("_", " ").capitalize() + ".cif_ringNodes"
    if c2 in nodes_files: return c2
    c3 = seq_id.replace("_", " ") + ".cif_ringNodes"
    if c3 in nodes_files: return c3
    seq_id_clean = seq_id.lower().replace("_", "")
    for f in nodes_files:
        f_clean = f.lower().replace("_", "").replace(" ", "")
        if seq_id_clean in f_clean or f_clean.startswith(seq_id_clean):
            return f
    parts = seq_id.split('_')
    if len(parts) >= 2:
        prefix = (parts[0] + parts[1]).lower()
        for f in nodes_files:
            f_clean = f.lower().replace("_", "").replace(" ", "")
            if f_clean.startswith(prefix):
                return f
    return None

def parse_fasta(file_path):
    records = []
    with open(file_path, 'r') as f:
        seq_id = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    records.append((seq_id, "".join(seq)))
                seq_id = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if seq_id:
            records.append((seq_id, "".join(seq)))
    return records

records = parse_fasta(fasta_file)
valid_records = []
species_to_nodes = {}

for seq_id, seq in records:
    match = find_nodes_file(seq_id)
    if match:
        species_to_nodes[seq_id] = os.path.join(nodes_dir, match)
        valid_records.append((seq_id, seq))

species_degrees = {}
for seq_id, fpath in species_to_nodes.items():
    try:
        df = pd.read_csv(fpath, sep='\t')
        df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
        df = df.dropna(subset=['Position'])
        df['Position'] = df['Position'].astype(int)
        pos_to_degree = dict(zip(df['Position'], df['Degree']))
        species_degrees[seq_id] = pos_to_degree
    except Exception as e:
        print(f"Error reading {fpath}: {e}")

# Build feature matrix
alignment_length = len(valid_records[0][1])
feature_matrix = []
labels = []

for seq_id, seq in valid_records:
    if seq_id not in species_degrees:
        continue
        
    protein_features = []
    for i in range(alignment_length):
        char = seq[i]
        if char == '-':
            protein_features.append(0.0) # Gap means 0 degree
        else:
            pos = len(seq[:i+1].replace('-', ''))
            degree = species_degrees[seq_id].get(pos)
            if degree is not None and not np.isnan(degree):
                protein_features.append(degree)
            else:
                protein_features.append(0.0)
                
    feature_matrix.append(protein_features)
    labels.append(seq_id)

X = np.array(feature_matrix)

# We might want to remove columns that are all 0 to avoid zero-variance issues, though standard scaler can handle it if we are careful.
# Standardize features (StandardScaler equivalent)
X_mean = np.mean(X, axis=0)
X_centered = X - X_mean
X_std = np.std(X_centered, axis=0)
X_std[X_std == 0] = 1.0 # prevent division by zero
X_scaled = X_centered / X_std

# --- PCA using numpy ---
cov_matrix = np.cov(X_scaled, rowvar=False)
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

sorted_index = np.argsort(eigenvalues)[::-1]
sorted_eigenvalues = eigenvalues[sorted_index]
sorted_eigenvectors = eigenvectors[:, sorted_index]
explained_variance_ratio = sorted_eigenvalues / np.sum(sorted_eigenvalues)

X_pca = np.dot(X_scaled, sorted_eigenvectors[:, :2])

plt.figure(figsize=(10, 8))
plt.scatter(X_pca[:, 0], X_pca[:, 1], c='blue', edgecolors='k', s=100)

# Add labels
for i, label in enumerate(labels):
    # Simplify label for display
    display_label = label.split('_')[0] + " " + label.split('_')[1] if '_' in label else label
    plt.annotate(display_label, (X_pca[i, 0], X_pca[i, 1]), xytext=(5, 5), textcoords='offset points', fontsize=9)

plt.title(f'PCA of Protein Degrees (Explained Variance: {sum(explained_variance_ratio[:2])*100:.1f}%)')
plt.xlabel(f'PC1 ({explained_variance_ratio[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({explained_variance_ratio[1]*100:.1f}%)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('degree_pca.png', dpi=300)
plt.close()
print("PCA plot saved as degree_pca.png")

# --- HCA (Dendrogram) ---
plt.figure(figsize=(12, 8))
# Using Ward variance minimization algorithm
Z = linkage(X_scaled, 'ward')

# Simplify labels
clean_labels = [l.split('_')[0] + " " + l.split('_')[1] if '_' in l else l for l in labels]

dendrogram(Z, labels=clean_labels, leaf_rotation=90, leaf_font_size=10)
plt.title('HCA Dendrogram of Protein Degrees')
plt.ylabel('Distance (Ward)')
plt.tight_layout()
plt.savefig('degree_hca.png', dpi=300)
plt.close()
print("HCA plot saved as degree_hca.png")
