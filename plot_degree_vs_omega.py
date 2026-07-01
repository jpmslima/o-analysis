import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

fasta_file = "pep.aligned.mafft.fsa"
nodes_dir = "Models/out"
fubar_file = "fubar-results.csv"

# Load FUBAR results
fubar = pd.read_csv(fubar_file)
# Sort by Site just in case
fubar = fubar.sort_values(by='Site').reset_index(drop=True)

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
if alignment_length != len(fubar):
    print(f"Warning: Alignment length ({alignment_length}) doesn't match FUBAR rows ({len(fubar)}).")

mean_degrees = []
omegas = []

for i in range(min(alignment_length, len(fubar))):
    col_degrees = []
    
    for seq_id, seq in valid_records:
        if seq_id not in species_degrees:
            continue
            
        char = seq[i]
        
        if char != '-':
            pos = len(seq[:i+1].replace('-', ''))
            degree = species_degrees[seq_id].get(pos)
            if degree is not None and not np.isnan(degree):
                col_degrees.append(degree)
                
    if len(col_degrees) > 0:
        mean_deg = np.mean(col_degrees)
        omega = fubar.loc[i, 'omega']
        # Also let's ignore extreme outlier omegas or nans if any
        if not np.isnan(omega) and not np.isinf(omega):
            mean_degrees.append(mean_deg)
            omegas.append(omega)

data = pd.DataFrame({
    'Mean_Degree': mean_degrees,
    'Omega': omegas
})

# Some omega values might be very high if alpha is near 0. We might want to use log(omega) or limit the axis,
# but let's just plot it as is first.
pearson_corr, p_value_p = pearsonr(data['Mean_Degree'], data['Omega'])
spearman_corr, p_value_s = spearmanr(data['Mean_Degree'], data['Omega'])

print(f"Number of valid alignment columns analyzed: {len(data)}")
print(f"Correlation (Degree vs Omega): Pearson r = {pearson_corr:.3f} (p={p_value_p:.3e}), Spearman rho = {spearman_corr:.3f} (p={p_value_s:.3e})")

plt.figure(figsize=(9, 6))
plt.scatter(data['Mean_Degree'], data['Omega'], alpha=0.6, color='purple')

# Add trendline
m, b = np.polyfit(data['Mean_Degree'], data['Omega'], 1)
plt.plot(data['Mean_Degree'], m*data['Mean_Degree'] + b, color='black', linestyle='--', 
         label=f'Trendline (Spearman rho={spearman_corr:.2f})')

plt.title('Mean Structural Degree vs Evolutionary Rate (Omega = beta/alpha)')
plt.xlabel('Mean Degree')
plt.ylabel('Omega (dN/dS)')
plt.legend()
plt.grid(alpha=0.3)

# If we have extreme omegas, limit y-axis slightly to better see the bulk,
# let's just use the 95th percentile as max if max is way too big.
# if data['Omega'].max() > 10:
#     plt.ylim(-0.1, np.percentile(data['Omega'], 98))

plt.tight_layout()
plt.savefig('degree_vs_omega.png', dpi=300)
print("Plot saved as degree_vs_omega.png")
