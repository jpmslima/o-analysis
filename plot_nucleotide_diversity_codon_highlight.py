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
    
    # Calculate Nucleotide Diversity
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
    
    colors_map = {
        0: ('#1f77b4', '1st Codon Position'),
        1: ('#2ca02c', '2nd Codon Position'),
        2: ('#d62728', '3rd Codon Position')
    }
    
    for remainder in [0, 1, 2]:
        color, label = colors_map[remainder]
        idx = [i for i in range(max_len) if i % 3 == remainder]
        x_vals = positions[idx]
        y_vals = [diversities[i] for i in idx]
        
        plt.bar(x_vals, y_vals, color=color, width=1.0, label=label, alpha=0.9)
    
    plt.title("Nucleotide Diversity Separated by Codon Position (Highlighted)", fontsize=16)
    plt.xlabel("Alignment Position (bp)", fontsize=14)
    plt.ylabel(r"Diversity ($1 - \sum p_i^2$)", fontsize=14)
    
    plt.axhline(y=0.8, color='dimgrey', linestyle='--', alpha=0.7, 
                label='Max Theoretical Diversity (5 states)')
                
    # Highlight codon positions 66, 69, 105
    codons_to_highlight = [66, 69, 105]
    for i, codon in enumerate(codons_to_highlight):
        # 1-based start and end bp for the codon
        start_bp = (codon - 1) * 3 + 1
        end_bp = codon * 3
        
        # Add shaded vertical region. We only add the label on the first iteration to avoid legend duplication
        label = 'Highlighted Codons' if i == 0 else None
        plt.axvspan(start_bp - 0.5, end_bp + 0.5, color='gold', alpha=0.4, label=label)
        
        # Add text annotation above the highlighted region
        plt.text((start_bp + end_bp)/2, 0.83, f'Codon {codon}', 
                 ha='center', va='bottom', fontsize=11, fontweight='bold', color='darkgoldenrod')
                
    plt.xlim(0, max_len + 1)
    
    actual_max = max(diversities) if diversities else 0
    # Provide extra headroom for the text annotations
    plt.ylim(0, max(0.9, actual_max * 1.15))
    
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.legend(loc='upper right', framealpha=0.9)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('nucleotide_diversity_codon_plot_highlighted.png', dpi=300)
    plt.savefig('nucleotide_diversity_codon_plot_highlighted.pdf')
    
    print("Successfully generated highlighted plot.")

if __name__ == '__main__':
    main()
