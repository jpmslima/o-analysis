import pandas as pd
import numpy as np

# Load FUBAR
fubar = pd.read_csv("fubar-results.csv")
fubar_dict = dict(zip(fubar['Site'], fubar['omega']))

# Align CIF to FUBAR sites
cif_seq = []
with open("Models/Consensus_Model.cif") as f:
    in_seq = False
    for line in f:
        if line.startswith("_entity_poly_seq.hetero"):
            in_seq = True
            continue
        if in_seq:
            if line.startswith("#"):
                break
            parts = line.strip().split()
            if len(parts) >= 3:
                cif_seq.append(parts[2])

aa_map = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'UNK': 'X'}
cif_seq_1L = "".join([aa_map.get(res, 'X') for res in cif_seq])

# Parse fasta to get alignment consensus (with gaps)
records = []
with open("pep.aligned.mafft.fsa", 'r') as f:
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

from collections import Counter
alignment_length = len(records[0][1])
consensus = []
for i in range(alignment_length):
    chars = [rec[1][i] for rec in records if rec[1][i] != '-']
    if chars:
        consensus.append(Counter(chars).most_common(1)[0][0])
    else:
        consensus.append('-')
consensus_seq = "".join(consensus)

def nw_align(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = np.zeros((m+1, n+1))
    for i in range(m+1): score[i][0] = -2 * i
    for j in range(n+1): score[0][j] = -2 * j
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1)
            delete = score[i-1][j] - 2
            insert = score[i][j-1] - 2
            score[i][j] = max(match, delete, insert)
            
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1; j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] - 2:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1
    return align1[::-1], align2[::-1]

a1, a2 = nw_align(consensus_seq, cif_seq_1L)

# Map CIF auth_seq_id (1-based) to FUBAR Site (1-based)
cif_to_omega = {}
site_idx = 1 # FUBAR Site corresponds to consensus_seq index + 1
cif_idx = 1  # CIF auth_seq_id

for char1, char2 in zip(a1, a2):
    if char1 != '-' and char2 != '-':
        if site_idx in fubar_dict:
            cif_to_omega[cif_idx] = fubar_dict[site_idx]
        cif_idx += 1
        site_idx += 1
    elif char1 != '-' and char2 == '-':
        site_idx += 1
    elif char1 == '-' and char2 != '-':
        # CIF has an insertion relative to consensus
        cif_idx += 1
    else:
        # Both gap
        pass

print("Mapped", len(cif_to_omega), "residues.")

# Replace in CIF file
out_lines = []
with open("Models/Consensus_Model.cif") as f:
    for line in f:
        if line.startswith("ATOM "):
            parts = line.split()
            # B_iso_or_equiv is index 17 (0-based) for these files
            auth_seq_id = int(parts[7])
            
            if auth_seq_id in cif_to_omega:
                omega_val = cif_to_omega[auth_seq_id]
                # Format to 3 decimal places
                parts[17] = f"{omega_val:.3f}"
            else:
                parts[17] = "0.000"
                
            # Need to format back. Easiest is just joining with space, 
            # though CIF columns might lose vertical alignment.
            # Let's preserve original formatting roughly by replacing the substring
            # But the original has 45.667. We can replace exactly the 18th whitespace-separated token.
            # Let's do a safe string replace.
            tokens = line.split()
            b_factor_str = tokens[17]
            # Wait, line.replace might replace the wrong thing if B-factor matches coordinate.
            # Let's build the line properly.
            # Actually, `line.split()` and then `.join()` is valid CIF format.
            out_lines.append(" ".join(parts) + "\n")
        else:
            out_lines.append(line)

with open("Models/Consensus_Model-omega.cif", "w") as f:
    f.writelines(out_lines)

print("Saved Consensus_Model-omega.cif")
