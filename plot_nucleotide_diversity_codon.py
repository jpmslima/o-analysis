import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
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
    fasta_file = 'dna.aligned.revtrans.fsa'
    
    try:
        sequences = read_fasta(fasta_file)
    except FileNotFoundError:
        print(f"Error: Could not find {fasta_file}")
        sys.exit(1)
        
    if not sequences:
        print("No sequences found in the file.")
        sys.exit(1)
        
    # Ensure all sequences have the same length
    max_len = max(len(seq) for seq in sequences)
    padded_seqs = [seq.ljust(max_len, '-') for seq in sequences]
    num_seqs = len(padded_seqs)
    
    diversities = []
    
    # Calculate Nucleotide Diversity (Expected Heterozygosity: 1 - sum(p_i^2))
    for col in range(max_len):
        column_chars = [seq[col].upper() for seq in padded_seqs]
        counts = Counter(column_chars)
        
        sum_p2 = 0.0
        for char, count in counts.items():
            p = count / num_seqs
            sum_p2 += p * p
            
        diversity = 1.0 - sum_p2
        diversities.append(diversity)
        
    # Plotting
    plt.figure(figsize=(15, 6))
    positions = np.arange(1, max_len + 1)
    
    # Separate the data into 1st, 2nd, and 3rd codon positions
    # 0-indexed column 'col' maps to codon position 'col % 3 + 1'
    colors_map = {
        0: ('#1f77b4', '1st Codon Position'), # Blue
        1: ('#2ca02c', '2nd Codon Position'), # Green
        2: ('#d62728', '3rd Codon Position')  # Red
    }
    
    # Plot each position group separately to assign proper colors and legends
    for remainder in [0, 1, 2]:
        color, label = colors_map[remainder]
        # Find all 0-based indices matching the remainder
        idx = [i for i in range(max_len) if i % 3 == remainder]
        # Extract corresponding 1-based positions and diversity values
        x_vals = positions[idx]
        y_vals = [diversities[i] for i in idx]
        
        plt.bar(x_vals, y_vals, color=color, width=1.0, label=label, alpha=0.9)
    
    plt.title("Nucleotide Diversity Separated by Codon Position", fontsize=16)
    plt.xlabel("Alignment Position (bp)", fontsize=14)
    plt.ylabel(r"Diversity ($1 - \sum p_i^2$)", fontsize=14)
    
    # Theoretical maximum diversity for 5 states
    plt.axhline(y=0.8, color='dimgrey', linestyle='--', alpha=0.7, 
                label='Max Theoretical Diversity (5 states)')
                
    plt.xlim(0, max_len + 1)
    
    actual_max = max(diversities) if diversities else 0
    plt.ylim(0, max(0.85, actual_max * 1.1))
    
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    # Place legend outside or inside
    plt.legend(loc='upper right', framealpha=0.9)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('nucleotide_diversity_codon_plot.png', dpi=300)
    plt.savefig('nucleotide_diversity_codon_plot.pdf')
    
    print(f"Processed {num_seqs} sequences of length {max_len}.")
    print("Successfully generated and saved codon-separated Nucleotide Diversity plot to 'nucleotide_diversity_codon_plot.png' and 'nucleotide_diversity_codon_plot.pdf'")

if __name__ == '__main__':
    main()
