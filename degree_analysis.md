# Análise de Degrees dos Resíduos

O script `plot_degree.py` foi criado para realizar as seguintes etapas:
1. **Lê o alinhamento de proteínas** no arquivo `pep.aligned.mafft.fsa`.
2. **Localiza os arquivos `.nodes`** de cada espécie correspondente na pasta `Models/out` e extrai os valores de `Degree` de cada resíduo.
3. **Mapeia as posições alinhadas**, ignorando `gaps` para que os índices correspondam corretamente às posições na estrutura.
4. **Calcula a média e o desvio padrão** para cada coluna do alinhamento.

Foi possível encontrar e associar automaticamente os arquivos `.nodes` para 25 das 29 sequências no alinhamento (alguns nomes diferem entre o FASTA e os arquivos `.cif`, como *Macaca mulatta*, *Lynx rufus*, *Panthera onca* e *Martes martes*, os quais não possuíam arquivos correspondentes).

Abaixo está o gráfico gerado mostrando a média (linha azul) e o desvio padrão (faixa sombreada azul) dos degrees dos resíduos baseando-se nas posições de alinhamento:

![Gráfico da Média e Desvio Padrão do Degree dos Resíduos](/home/jpmslima/.gemini/antigravity/brain/3dbee5f2-ee26-4bd7-b21d-dec0920699b1/artifacts/degree_alignment_plot.png)
