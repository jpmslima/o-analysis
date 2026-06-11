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
    fasta_file = 'pep.aligned.mafft.fsa'
    
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
    
    # Calculate Amino Acid Diversity (Expected Heterozygosity: 1 - sum(p_i^2))
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
    
    # Bar plot for amino acid diversity
    plt.bar(positions, diversities, color='darkorange', width=1.0)
    
    plt.title("Amino Acid Diversity per Alignment Position", fontsize=16)
    plt.xlabel("Alignment Position", fontsize=14)
    plt.ylabel(r"Diversity ($1 - \sum p_i^2$)", fontsize=14)
    
    # Max possible diversity for 21 states (20 AAs + 1 gap)
    max_diversity = 1.0 - (1.0 / 21.0)
    plt.axhline(y=max_diversity, color='crimson', linestyle='--', alpha=0.7, 
                label=f'Max Theoretical Diversity (21 states): {max_diversity:.2f}')
                
    plt.xlim(0, max_len + 1)
    
    actual_max = max(diversities) if diversities else 0
    # Set y limits nicely
    plt.ylim(0, max(max_diversity * 1.05, actual_max * 1.1))
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('aa_diversity_plot.png', dpi=300)
    plt.savefig('aa_diversity_plot.pdf')
    
    print(f"Processed {num_seqs} sequences of length {max_len}.")
    print("Successfully generated and saved Amino Acid Diversity plot to 'aa_diversity_plot.png' and 'aa_diversity_plot.pdf'")

if __name__ == '__main__':
    main()
