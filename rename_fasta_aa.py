import sys

def main():
    # Load mapping
    mapping = {}
    with open('NCBI-Map-Hs-32.csv', 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                acc = parts[0].strip()
                name = parts[1].strip()
                # Replace spaces with underscores
                name = name.replace(' ', '_')
                mapping[acc] = name
                
                # Map the base accession (without the .version) just in case
                base_acc = acc.split('.')[0]
                if base_acc not in mapping:
                    mapping[base_acc] = name
    
    # Process FASTA
    output = []
    with open('URAD-AA-Hs.fasta', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:] # remove >
                acc = ""
                # Parse the specific efetch format: lcl|XM_1234.1_prot_XP_...
                if header.startswith('lcl|'):
                    parts = header.split('|')
                    if len(parts) > 1:
                        # Extract the part before '_prot'
                        acc = parts[1].split('_prot')[0]
                
                # Fallback if the format is different
                if not acc:
                    acc = header.split()[0]
                
                base_acc = acc.split('.')[0]
                
                new_name = mapping.get(acc)
                if not new_name:
                    new_name = mapping.get(base_acc)
                
                if new_name:
                    output.append(f">{new_name}")
                else:
                    # If not found, just use the extracted accession
                    output.append(f">{acc}")
            else:
                output.append(line)
                
    with open('URAD-AA-Hs.fasta', 'w') as f:
        for line in output:
            f.write(line + '\n')
            
    print("FASTA AA headers updated successfully.")

if __name__ == '__main__':
    main()
