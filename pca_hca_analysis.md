# Análise PCA e HCA de Degree de Proteínas

Com base nos valores de degree de cada resíduo, calculamos matrizes de características para cada proteína, alinhadas de acordo com as posições estruturais (tratando *gaps* do alinhamento como ausência de resíduo, i.e., degree zero). 

Em seguida, aplicamos métodos de redução de dimensionalidade e agrupamento hierárquico para verificar a similaridade entre os perfis topológicos das proteínas.

### Análise de Componentes Principais (PCA)

O gráfico de PCA mostra a distribuição das proteínas em um espaço bidimensional criado para maximizar a variação explicada. Proteínas que aparecem próximas neste gráfico têm perfis de interações inter-resíduos (degrees) muito similares ao longo de toda a cadeia:

![PCA das Proteínas](/home/jpmslima/.gemini/antigravity/brain/3dbee5f2-ee26-4bd7-b21d-dec0920699b1/artifacts/degree_pca.png)

### Agrupamento Hierárquico (HCA)

O dendrograma HCA organiza as proteínas em clusters usando o método de *Ward* (minimização de variância) sobre a distância das características. Quanto menor for a altura em que dois ramos se juntam (no eixo Y), mais similares são as proteínas:

![Dendrograma HCA](/home/jpmslima/.gemini/antigravity/brain/3dbee5f2-ee26-4bd7-b21d-dec0920699b1/artifacts/degree_hca.png)

Essas análises ajudam a revelar grupos evolutivos ou funcionais baseados na topologia estrutural (conectividade dos resíduos) em vez de focar estritamente na sequência de aminoácidos.
