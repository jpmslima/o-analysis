#!/bin/bash

# Cria a pasta de saída "out" caso ela não exista
mkdir -p out/

# Itera sobre todos os arquivos .cif no diretório atual
for file in *.cif; do
    # Verifica se o arquivo realmente existe (útil caso não haja nenhum .cif na pasta)
    if [ -f "$file" ]; then
        echo "Processando: $file"
        ring -i "$file" --all_edges --out_dir out/
    fi
done

echo "Processamento concluído! Resultados salvos na pasta 'out/'."
