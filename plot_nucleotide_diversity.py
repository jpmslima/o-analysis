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
    
    # Bar plot for nucleotide diversity
    plt.bar(positions, diversities, color='mediumpurple', width=1.0)
    
    plt.title("Nucleotide Diversity per Alignment Position", fontsize=16)
    plt.xlabel("Alignment Position (bp)", fontsize=14)
    plt.ylabel(r"Diversity ($1 - \sum p_i^2$)", fontsize=14)
    
    # Max possible diversity for 5 states (A, C, G, T, -) is 1 - 5*(1/5)^2 = 0.8
    plt.axhline(y=0.8, color='crimson', linestyle='--', alpha=0.7, 
                label='Max Theoretical Diversity (5 states: A,C,G,T,-)')
                
    plt.xlim(0, max_len + 1)
    
    actual_max = max(diversities) if diversities else 0
    # Set y limits nicely
    plt.ylim(0, max(0.85, actual_max * 1.1))
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('nucleotide_diversity_plot.png', dpi=300)
    plt.savefig('nucleotide_diversity_plot.pdf')
    
    print(f"Processed {num_seqs} sequences of length {max_len}.")
    print("Successfully generated and saved Nucleotide Diversity plot to 'nucleotide_diversity_plot.png' and 'nucleotide_diversity_plot.pdf'")

if __name__ == '__main__':
    main()
