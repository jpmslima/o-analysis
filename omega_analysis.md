# Relação entre Degree Estrutural e Taxa Evolutiva (Omega)

Para complementar nossa análise anterior, cruzamos o **Degree Médio** de cada posição do alinhamento com os valores de **Omega ($\omega$)** (a razão entre as taxas de mutação não-sinônimas e sinônimas, calculada como $\beta/\alpha$) recém-adicionados ao arquivo do FUBAR.

O valor de $\omega$ reflete o tipo de pressão seletiva atuando no códon:
- **$\omega < 1$**: Seleção purificadora (conservação).
- **$\omega \approx 1$**: Evolução neutra.
- **$\omega > 1$**: Seleção positiva/diversificadora.

### Resultados do Cruzamento
Analisamos todos os 180 sítios mapeados e encontramos as seguintes correlações:
- **Correlação de Spearman ($\rho$):** $-0.276$ ($p \approx 0.00017$)
- **Correlação de Pearson ($r$):** $-0.257$ ($p \approx 0.0005$)

**Conclusão:**
Existe uma **correlação negativa clara e estatisticamente significativa**. Isso significa que à medida que o número de contatos de um resíduo na estrutura 3D (Degree) **aumenta**, a taxa evolutiva $\omega$ **diminui**.

Em termos biológicos: os resíduos centrais e altamente conectados estão sob uma fortíssima pressão seletiva purificadora. Mutações (não-sinônimas) nesses locais são altamente deletérias, mantendo a taxa $\omega$ muito baixa.

### Visualização Gráfica

O gráfico abaixo ilustra a dispersão e a linha de tendência linear dos dados:

![Dispersão Degree vs Omega](/home/jpmslima/.gemini/antigravity/brain/3dbee5f2-ee26-4bd7-b21d-dec0920699b1/artifacts/degree_vs_omega.png)
