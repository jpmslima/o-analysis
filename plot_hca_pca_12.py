import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

def main():
    # 1. Read codes
    codes = {}
    with open('12-codes', 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    codes[parts[0]] = " ".join(parts[1:])
    
    # 2. Read distance matrix
    df = pd.read_csv('12-RMSDs.tsv', sep=r'\s+', index_col=0)
    
    # Ensure it's perfectly symmetric and diagonal is zero
    distances = df.values
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)
    
    # 3. HCA
    condensed_dist = squareform(distances)
    Z = linkage(condensed_dist, method='average')
    labels = [codes.get(str(int(idx)), str(idx)) for idx in df.index]
    
    plt.figure(figsize=(10, 7))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=12)
    plt.title('Hierarchical Clustering based on RMSD (12 Models)', fontsize=14)
    plt.ylabel('RMSD Distance', fontsize=12)
    plt.tight_layout()
    plt.savefig('hca_12_dendrogram.png', dpi=300)
    plt.savefig('hca_12_dendrogram.pdf')
    plt.close()
    
    # 4. PCoA / PCA (Manually since sklearn might not be present)
    n = distances.shape[0]
    D_sq = distances ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D_sq).dot(J)
    
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    eigenvalues_2d = np.maximum(eigenvalues[:2], 0)
    eigenvectors_2d = eigenvectors[:, :2]
    coords = eigenvectors_2d * np.sqrt(eigenvalues_2d)
    
    plt.figure(figsize=(10, 7))
    plt.scatter(coords[:, 0], coords[:, 1], c='coral', marker='o', edgecolors='black', s=80, alpha=0.9)
    
    for i, label in enumerate(labels):
        plt.annotate(label, (coords[i, 0], coords[i, 1]), 
                     xytext=(6, 6), textcoords='offset points', 
                     fontsize=11, alpha=0.9)
        
    plt.title('Principal Coordinate Analysis (PCoA) of RMSD Matrix', fontsize=14)
    plt.xlabel('Principal Coordinate 1', fontsize=12)
    plt.ylabel('Principal Coordinate 2', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    plt.savefig('pca_12_plot.png', dpi=300)
    plt.savefig('pca_12_plot.pdf')
    plt.close()
    
    print("Successfully generated HCA and PCA plots for the 12 models dataset.")

if __name__ == '__main__':
    main()
