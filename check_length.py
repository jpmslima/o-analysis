import pandas as pd

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

records = parse_fasta("pep.aligned.mafft.fsa")
alignment_length = len(records[0][1])

fubar = pd.read_csv("fubar-results.csv")

print(f"Alignment length: {alignment_length}")
print(f"FUBAR rows: {len(fubar)}")
