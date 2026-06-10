import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import math
import sys

def read_fasta(filename):
    sequences = []
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append("".join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            sequences.append("".join(seq))
    return sequences

def main():
    fasta_file = 'Hs-alignment.fasta'
    
    # Parse the FASTA file
    try:
        sequences = read_fasta(fasta_file)
    except FileNotFoundError:
        print(f"Error: Could not find {fasta_file}")
        sys.exit(1)
        
    if not sequences:
        print("No sequences found in the file.")
        sys.exit(1)
        
    # Ensure all sequences have the same length (alignment)
    max_len = max(len(seq) for seq in sequences)
    padded_seqs = [seq.ljust(max_len, '-') for seq in sequences]
    num_seqs = len(padded_seqs)
    
    # Calculate Shannon Entropy for each column
    entropies = []
    for col in range(max_len):
        # Extract the column
        column_chars = [seq[col].upper() for seq in padded_seqs]
        counts = Counter(column_chars)
        
        # Shannon Entropy formula: H = - sum(p * log2(p))
        entropy = 0.0
        for char, count in counts.items():
            p = count / num_seqs
            entropy -= p * math.log2(p)
            
        entropies.append(entropy)
        
    # Plotting
    plt.figure(figsize=(15, 6))
    
    positions = np.arange(1, max_len + 1)
    
    # We use a bar plot to show the entropy at each position
    plt.bar(positions, entropies, color='cornflowerblue', width=1.0)
    
    plt.title("Shannon's Entropy per Alignment Position", fontsize=16)
    plt.xlabel("Alignment Position", fontsize=14)
    plt.ylabel("Shannon Entropy (bits)", fontsize=14)
    
    # 20 standard amino acids + 1 gap character = 21 possible states
    max_entropy = math.log2(21)
    
    plt.axhline(y=max_entropy, color='crimson', linestyle='--', alpha=0.7, 
                label=f'Max Theoretical Entropy ({max_entropy:.2f} bits)')
    
    # Set axis limits and grid
    plt.xlim(0, max_len + 1)
    
    # Make y-axis dynamic based on the actual max entropy in the data
    actual_max = max(entropies) if entropies else 0
    plt.ylim(0, max(max_entropy * 1.1, actual_max * 1.1))
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save the plot in both PNG and PDF formats
    plt.savefig('shannon_entropy_plot.png', dpi=300)
    plt.savefig('shannon_entropy_plot.pdf')
    
    print(f"Processed {num_seqs} sequences of length {max_len}.")
    print("Successfully generated and saved Shannon Entropy plot to 'shannon_entropy_plot.png' and 'shannon_entropy_plot.pdf'")

if __name__ == '__main__':
    main()
