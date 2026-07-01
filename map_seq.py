import sys
from collections import Counter

# Parse CIF sequence
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

# Parse fasta
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

alignment_length = len(records[0][1])
consensus = []
for i in range(alignment_length):
    chars = [rec[1][i] for rec in records if rec[1][i] != '-']
    if chars:
        consensus.append(Counter(chars).most_common(1)[0][0])
    else:
        consensus.append('-')

consensus_seq = "".join(consensus)

print(f"CIF length: {len(cif_seq_1L)}")
print(f"Alignment length: {alignment_length}")

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

import numpy as np
cons_nogap = consensus_seq.replace("-", "")
a1, a2 = nw_align(cons_nogap, cif_seq_1L)

print(a1)
print(a2)

# Map CIF index to alignment index
cif_to_aln = {}
aln_idx = 0
cif_idx = 0

for i in range(len(consensus_seq)):
    if consensus_seq[i] != '-':
        # Find which CIF residue this is
        # We can map consensus_seq index directly
        pass

# Actually let's map FUBAR Site to CIF residue
# FUBAR Site = alignment column index + 1
# Let's map CIF auth_seq_id (1-179) to FUBAR Site (1-180)
cif_to_fubar = {}
c_idx = 0
for i in range(len(consensus_seq)):
    if consensus_seq[i] != '-':
        if c_idx < len(cif_seq_1L):
            cif_to_fubar[c_idx + 1] = i + 1  # 1-based indices
            c_idx += 1

print("Mapping (first 10):", list(cif_to_fubar.items())[:10])
print("Mapping (last 10):", list(cif_to_fubar.items())[-10:])
