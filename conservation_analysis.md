# Relação entre Conservação de Resíduos e Degree

Para responder à sua pergunta — **"resíduos com maiores degrees são mais conservados?"** —, foi realizado um cálculo relacionando a conservação de cada coluna do alinhamento com a média de *degree* estrutural observada para aquela mesma posição.

Calculamos duas métricas de conservação baseadas nas sequências do arquivo `pep.aligned.mafft.fsa`:
1. **Entropia de Shannon:** Mede a variação de aminoácidos em uma coluna. Valores **menores** indicam **maior conservação** (menos variação).
2. **Identidade (%):** Porcentagem do aminoácido mais frequente na coluna. Valores **maiores** indicam **maior conservação**.

### Resultados Estatísticos
Ao analisarmos 173 posições válidas do alinhamento, obtivemos as seguintes correlações:
- **Degree vs Entropia:** Correlação de Pearson $r = -0.250$ ($p \approx 0.0009$). 
  *(A correlação negativa confirma que maiores degrees estão associados a menores entropias, ou seja, maior conservação).*
- **Degree vs % de Identidade:** Correlação de Pearson $r = 0.230$ ($p \approx 0.0023$). 
  *(A correlação positiva confirma que maiores degrees estão associados a maiores porcentagens de identidade).*

**Conclusão:** 
**Sim.** Existe uma correlação estatisticamente significativa provando que resíduos com um número maior de contatos estruturais (maior *degree*) tendem a ser mais conservados evolutivamente. 

Isso ocorre porque resíduos altamente conectados frequentemente formam o núcleo hidrofóbico da proteína ou sítios funcionais essenciais; mutações nessas posições têm alta probabilidade de desestabilizar a estrutura da proteína, sendo filtradas pela seleção natural.

### Visualização

![Gráficos de Dispersão: Degree vs Conservação](/home/jpmslima/.gemini/antigravity/brain/3dbee5f2-ee26-4bd7-b21d-dec0920699b1/artifacts/conservation_vs_degree.png)
