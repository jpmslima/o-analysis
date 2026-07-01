import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from collections import Counter
import math

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

# Read degrees
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
        pass

alignment_length = len(valid_records[0][1])

def calc_entropy(column_chars):
    # Calculate Shannon Entropy for the column (excluding gaps)
    chars = [c for c in column_chars if c != '-']
    if len(chars) == 0:
        return np.nan
    counts = Counter(chars)
    total = len(chars)
    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy

def calc_conservation_identity(column_chars):
    # Percentage of the most frequent amino acid (excluding gaps)
    chars = [c for c in column_chars if c != '-']
    if len(chars) == 0:
        return np.nan
    counts = Counter(chars)
    most_common_count = counts.most_common(1)[0][1]
    return (most_common_count / len(chars)) * 100

mean_degrees = []
entropies = []
identities = []

for i in range(alignment_length):
    col_chars = []
    col_degrees = []
    
    for seq_id, seq in valid_records:
        if seq_id not in species_degrees:
            continue
            
        char = seq[i]
        col_chars.append(char)
        
        if char != '-':
            pos = len(seq[:i+1].replace('-', ''))
            degree = species_degrees[seq_id].get(pos)
            if degree is not None and not np.isnan(degree):
                col_degrees.append(degree)
                
    if len(col_degrees) > 0 and col_chars.count('-') < len(col_chars) * 0.5: # Ignore highly gapped columns
        mean_degrees.append(np.mean(col_degrees))
        entropies.append(calc_entropy(col_chars))
        identities.append(calc_conservation_identity(col_chars))

# Filter NaNs just in case
data = pd.DataFrame({
    'Mean_Degree': mean_degrees,
    'Entropy': entropies,
    'Identity': identities
}).dropna()

# Correlation between Degree and Entropy (Entropy: lower is more conserved)
# If highly conserved residues have higher degree, we expect a NEGATIVE correlation between Degree and Entropy
pearson_corr_ent, p_value_ent = pearsonr(data['Mean_Degree'], data['Entropy'])
spearman_corr_ent, sp_value_ent = spearmanr(data['Mean_Degree'], data['Entropy'])

# Correlation between Degree and Identity (Identity: higher is more conserved)
# We expect a POSITIVE correlation between Degree and Identity
pearson_corr_id, p_value_id = pearsonr(data['Mean_Degree'], data['Identity'])

print(f"Results:")
print(f"Number of valid alignment columns analyzed: {len(data)}")
print(f"Correlation (Degree vs Entropy): Pearson r = {pearson_corr_ent:.3f} (p={p_value_ent:.3e}), Spearman rho = {spearman_corr_ent:.3f} (p={sp_value_ent:.3e})")
print(f"Correlation (Degree vs % Identity): Pearson r = {pearson_corr_id:.3f} (p={p_value_id:.3e})")

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Mean Degree vs Entropy
ax1.scatter(data['Mean_Degree'], data['Entropy'], alpha=0.6, color='coral')
m, b = np.polyfit(data['Mean_Degree'], data['Entropy'], 1)
ax1.plot(data['Mean_Degree'], m*data['Mean_Degree'] + b, color='red', label=f'Trendline (r={pearson_corr_ent:.2f})')
ax1.set_xlabel('Mean Degree')
ax1.set_ylabel('Shannon Entropy (Lower = More Conserved)')
ax1.set_title('Mean Degree vs Conservation (Entropy)')
ax1.legend()
ax1.grid(alpha=0.3)

# Plot 2: Mean Degree vs Identity %
ax2.scatter(data['Mean_Degree'], data['Identity'], alpha=0.6, color='teal')
m2, b2 = np.polyfit(data['Mean_Degree'], data['Identity'], 1)
ax2.plot(data['Mean_Degree'], m2*data['Mean_Degree'] + b2, color='darkblue', label=f'Trendline (r={pearson_corr_id:.2f})')
ax2.set_xlabel('Mean Degree')
ax2.set_ylabel('Conservation (% Identity)')
ax2.set_title('Mean Degree vs Conservation (% Identity)')
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('conservation_vs_degree.png', dpi=300)
print("Plot saved as conservation_vs_degree.png")
