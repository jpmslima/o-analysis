import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Read the codes mapping Model ID to Name
    codes = {}
    with open('RMSD-structure-codes', 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                # If there are spaces in the name, join them
                codes[parts[0]] = " ".join(parts[1:])

    # Read the distance matrix
    df = pd.read_csv('RMSD-all.csv', sep=r'\s+', index_col=0)
    
    # Ensure it's a square matrix and perfectly symmetric
    distances = df.values
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)
    
    # Perform Principal Coordinate Analysis (PCoA / Classical MDS) manually 
    # since scikit-learn is not available in the environment.
    n = distances.shape[0]
    
    # 1. Squared distance matrix
    D_sq = distances ** 2
    
    # 2. Centering matrix J = I - 1/n * 11^T
    J = np.eye(n) - np.ones((n, n)) / n
    
    # 3. Double-centered distance matrix B
    B = -0.5 * J.dot(D_sq).dot(J)
    
    # 4. Eigendecomposition of B
    # eigh is for symmetric matrices, which B is
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    
    # 5. Sort eigenvalues and eigenvectors in descending order
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # 6. Take top 2 components
    # Handle small negative eigenvalues due to numerical instability by taking max(0, val)
    eigenvalues_2d = np.maximum(eigenvalues[:2], 0)
    eigenvectors_2d = eigenvectors[:, :2]
    
    # Calculate coordinates
    coords = eigenvectors_2d * np.sqrt(eigenvalues_2d)
    
    # Get labels
    labels = [codes[str(int(idx))] for idx in df.index]
    
    # Plot
    plt.figure(figsize=(12, 8))
    plt.scatter(coords[:, 0], coords[:, 1], c='dodgerblue', marker='o', edgecolors='black', s=60, alpha=0.8)
    
    # Add text annotations for each point
    for i, label in enumerate(labels):
        plt.annotate(label, (coords[i, 0], coords[i, 1]), 
                     xytext=(5, 5), textcoords='offset points', 
                     fontsize=9, alpha=0.8)
        
    plt.title('Principal Coordinate Analysis (PCoA / PCA) of RMSD Matrix')
    plt.xlabel('Principal Coordinate 1')
    plt.ylabel('Principal Coordinate 2')
    
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('pca_plot.png', dpi=300)
    plt.savefig('pca_plot.pdf')
    print("Successfully generated and saved PCA/PCoA plot to 'pca_plot.png' and 'pca_plot.pdf'")

if __name__ == '__main__':
    main()
