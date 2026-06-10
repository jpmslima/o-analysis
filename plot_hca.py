import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform

def main():
    # Read the codes mapping Model ID to Name
    codes = {}
    with open('RMSD-structure-codes', 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                # If there are spaces in the name, join them, although the current file uses underscores
                codes[parts[0]] = " ".join(parts[1:])

    # Read the distance matrix
    # Using \s+ separator to handle both spaces and tabs gracefully
    df = pd.read_csv('RMSD-all.csv', sep='\s+', index_col=0)
    
    # Ensure it's a square matrix
    distances = df.values
    
    # Make strictly symmetric to avoid scipy squareform errors
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)
    
    # Convert to condensed distance matrix format required by linkage
    condensed_dist = squareform(distances)
    
    # Perform Hierarchical Clustering Analysis (UPGMA / average linkage)
    Z = linkage(condensed_dist, method='average')
    
    # Get labels for the dendrogram leaves
    labels = [codes[str(int(idx))] for idx in df.index]
    
    # Plot the dendrogram
    plt.figure(figsize=(12, 8))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=10)
    
    plt.title('Hierarchical Clustering of Models based on RMSD (Average Linkage)')
    plt.xlabel('Species / Model')
    plt.ylabel('RMSD')
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('hca_dendrogram.png', dpi=300)
    plt.savefig('hca_dendrogram.pdf')
    print("Successfully generated and saved HCA clustering to 'hca_dendrogram.png' and 'hca_dendrogram.pdf'")

if __name__ == '__main__':
    main()
