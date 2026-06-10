def merge_and_deduplicate_fasta(input_files, output_file):
    seen_sequences = set()
    unique_records = []

    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                header = ""
                seq = []
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith(">"):
                        # Salva a sequência anterior antes de começar a nova
                        if header:
                            full_seq = "".join(seq)
                            if full_seq not in seen_sequences:
                                seen_sequences.add(full_seq)
                                unique_records.append((header, full_seq))
                        header = line
                        seq = []
                    else:
                        # Remove quebras de linha ou espaços no meio da sequência
                        seq.append(line.replace(" ", ""))
                
                # Salva a última sequência do arquivo
                if header:
                    full_seq = "".join(seq)
                    if full_seq not in seen_sequences:
                        seen_sequences.add(full_seq)
                        unique_records.append((header, full_seq))
        except FileNotFoundError:
            print(f"Erro: O arquivo {file_path} não foi encontrado.")
            return

    # Escreve o novo arquivo FASTA
    with open(output_file, 'w') as out:
        for header, sequence in unique_records:
            out.write(f"{header}\n")
            # Formata a sequência para ter 60 caracteres por linha (padrão FASTA)
            for i in range(0, len(sequence), 60):
                out.write(f"{sequence[i:i+60]}\n")
                
    print(f"Concluído! Arquivo '{output_file}' gerado com {len(unique_records)} sequências únicas.")

# Insira o nome exato dos arquivos que você fez upload
arquivos_entrada = [
    'Dsuniprotkb_uniref_cluster_50_UniRef50_A1_2026_06_09.fasta.txt', 
    'Hs_uniprotkb_uniref_cluster_50_UniRef50_A6_2026_06_09.fasta.txt'
]
arquivo_saida = 'sequencias_unicas_consolidadas.fasta'

merge_and_deduplicate_fasta(arquivos_entrada, arquivo_saida)
