import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fasta_file = "pep.aligned.mafft.fsa"
nodes_dir = "Models/out"

# Load nodes files
nodes_files = os.listdir(nodes_dir)
nodes_files = [f for f in nodes_files if f.endswith("_ringNodes")]

def find_nodes_file(seq_id):
    # exact match lower
    c1 = seq_id.lower() + ".cif_ringNodes"
    if c1 in nodes_files: return c1
    # replace _ with space, then Title case
    c2 = seq_id.replace("_", " ").capitalize() + ".cif_ringNodes"
    if c2 in nodes_files: return c2
    # replace _ with space, no cap
    c3 = seq_id.replace("_", " ") + ".cif_ringNodes"
    if c3 in nodes_files: return c3
    
    # partial match
    seq_id_clean = seq_id.lower().replace("_", "")
    for f in nodes_files:
        f_clean = f.lower().replace("_", "").replace(" ", "")
        if seq_id_clean in f_clean or f_clean.startswith(seq_id_clean):
            return f
            
    # Try just the first two words
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
    else:
        print(f"Warning: No nodes file found for {seq_id}")

print(f"Found nodes files for {len(valid_records)}/{len(records)} sequences.")

# Read degrees for each valid record
# The degree is the 7th column (index 6) in the nodes file. 'Position' is the 3rd column.
species_degrees = {}
for seq_id, fpath in species_to_nodes.items():
    try:
        # Some nodes files might have extra spaces or weird formats, but they are tab separated
        df = pd.read_csv(fpath, sep='\t')
        # Filter only RES types, or just take all that have a numeric position
        # Position might be string if it has insertion codes, but usually int
        df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
        df = df.dropna(subset=['Position'])
        df['Position'] = df['Position'].astype(int)
        
        # Create a mapping from position to degree
        pos_to_degree = dict(zip(df['Position'], df['Degree']))
        species_degrees[seq_id] = pos_to_degree
    except Exception as e:
        print(f"Error reading {fpath}: {e}")

if not species_degrees:
    print("No valid degrees loaded. Exiting.")
    exit(1)

# Now, align columns
alignment_length = len(valid_records[0][1])
mean_degrees = []
std_degrees = []

for i in range(alignment_length):
    col_degrees = []
    for seq_id, seq in valid_records:
        if seq_id not in species_degrees:
            continue
            
        char = seq[i]
        if char != '-':
            # Find the position of this residue in the unaligned sequence
            # Position = number of non-gap characters up to this point (1-based)
            pos = len(seq[:i+1].replace('-', ''))
            
            degree = species_degrees[seq_id].get(pos)
            if degree is not None and not np.isnan(degree):
                col_degrees.append(degree)
                
    if len(col_degrees) > 0:
        mean_degrees.append(np.mean(col_degrees))
        std_degrees.append(np.std(col_degrees))
    else:
        mean_degrees.append(np.nan)
        std_degrees.append(np.nan)

# Plot
plt.figure(figsize=(15, 6))
x = np.arange(alignment_length)

# Convert to numpy arrays to handle nans
mean_degrees = np.array(mean_degrees)
std_degrees = np.array(std_degrees)

# Plot valid points
valid_idx = ~np.isnan(mean_degrees)

plt.plot(x[valid_idx], mean_degrees[valid_idx], color='blue', label='Mean Degree')
plt.fill_between(x[valid_idx], 
                 mean_degrees[valid_idx] - std_degrees[valid_idx], 
                 mean_degrees[valid_idx] + std_degrees[valid_idx], 
                 color='blue', alpha=0.2, label='Std Dev')

plt.title('Mean and Standard Deviation of Residue Degrees Based on Alignment')
plt.xlabel('Alignment Position')
plt.ylabel('Degree')
plt.legend()
plt.tight_layout()
plt.savefig('degree_alignment_plot.png', dpi=300)
print("Plot saved as degree_alignment_plot.png")
