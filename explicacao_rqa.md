# Análise de Quantificação de Recorrência (RQA)

Documento dos arquivos `rqa_portable.c`, `rqa_embedded.c` e `rqa_embedded_optimized.c`.

O programa calcula uma medida de organização temporal de um trecho de EEG sintético: o determinismo (DET) da Análise de Quantificação de Recorrência. No experimento proposto, um sinal mais rítmico e repetitivo tende a produzir DET maior que um sinal irregular. O resultado é um demonstrador acadêmico de processamento de sinais e não constitui método clínico de diagnóstico ou predição de epilepsia.

---

## Base teórica

### EEG como série temporal de um sistema complexo

O eletroencefalograma (EEG) é uma série temporal obtida pela medição, no escalpo ou intracranialmente, da atividade elétrica resultante de populações neuronais. Essa atividade apresenta oscilações, ruído, artefatos, acoplamentos e mudanças no tempo. Portanto, métodos espectrais lineares e métodos não lineares devem ser entendidos como abordagens complementares: a análise espectral descreve a distribuição de potência por frequência, enquanto a análise de recorrência investiga quando estados semelhantes reaparecem e por quanto tempo suas evoluções permanecem semelhantes.

A RQA não exige que se suponha previamente um modelo linear para o sinal. Isso não significa, porém, que ela seja imune à não estacionariedade, ao ruído ou à escolha de parâmetros. O tamanho da janela, a normalização, a dimensão de imersão, o atraso, a métrica, o limiar de recorrência e a janela de Theiler afetam diretamente os resultados e precisam ser informados para que o experimento seja reproduzível.

### Reconstrução do espaço de estados

Se apenas uma observação escalar `u(i)` está disponível, uma trajetória pode ser reconstruída por coordenadas atrasadas:

```text
x(i) = [u(i), u(i + τ), ..., u(i + (m − 1)τ)]
```

em que `m` é a dimensão de imersão e `τ` é o atraso temporal. O fundamento dessa reconstrução está no teorema de imersão de Takens: sob hipóteses matemáticas específicas e para uma observação genérica, uma reconstrução de dimensão suficiente preserva propriedades topológicas do atrator original.

Este projeto adota deliberadamente **m = 1** e **τ = 1**. Assim, `x(i) = [u(i)]` e a distância entre estados reduz-se a `|u(i) − u(j)|`. O código mede, portanto, recorrências de amplitude e a persistência temporal dessas recorrências; ele não realiza uma reconstrução multidimensional completa da dinâmica cortical. Essa simplificação reduz memória e processamento e é adequada ao objetivo didático, mas limita interpretações sobre invariantes dinâmicos.

> Takens, F. (1981). Detecting strange attractors in turbulence. In *Dynamical Systems and Turbulence, Warwick 1980*, Lecture Notes in Mathematics, 898, 366–381.

### Matriz de recorrência

Eckmann, Kamphorst e Ruelle (1987) introduziram o *Recurrence Plot* (RP) para representar as recorrências de uma trajetória. Para os vetores de estado `x(i)` e `x(j)`, define-se:

```text
R(i,j) = Θ(ε − ‖x(i) − x(j)‖),    i,j = 0,...,N−1
```

onde `Θ` é a função degrau de Heaviside, `ε` é o limiar de recorrência e `‖·‖` é uma norma. Logo, `R(i,j)=1` quando os dois estados estão dentro de uma vizinhança de raio `ε`, e `R(i,j)=0` caso contrário.

Neste projeto, como `m=1`, a norma euclidiana é simplesmente o valor absoluto:

```text
R(i,j) = 1, se |u(i) − u(j)| ≤ 8
R(i,j) = 0, caso contrário
```

A matriz é simétrica e sua diagonal principal, chamada linha de identidade (LOI), é sempre recorrente porque `R(i,i)=1`. O código exclui essa diagonal do cálculo e usa apenas `j>i`; devido à simetria, isso produz a mesma razão que contar os dois triângulos da matriz.

Linhas diagonais paralelas à LOI têm interpretação temporal: se `R(i+k,j+k)=1` durante vários valores consecutivos de `k`, dois segmentos da trajetória permanecem próximos enquanto evoluem. Pontos isolados indicam recorrências breves. Entretanto, linhas diagonais também podem surgir por periodicidade, forte autocorrelação, suavidade do sinal ou escolha de um `ε` grande; sua presença não prova, isoladamente, caos ou determinismo físico.

> Eckmann, J.-P., Kamphorst, S. O., & Ruelle, D. (1987). Recurrence plots of dynamical systems. *Europhysics Letters*, 4(9), 973–977.

### Quantificação da recorrência e determinismo

Zbilut e Webber (1992) sistematizaram medidas extraídas do RP, dando origem à Análise de Quantificação de Recorrência. Seja `P(l)` o número de linhas diagonais com comprimento exatamente `l`. Para um comprimento mínimo `l_min`, o determinismo é:

```text
             Σ[l=l_min até N] l·P(l)
DET = ─────────────────────────────────
              Σ[l=1 até N] l·P(l)
```

O denominador é o número de pontos recorrentes considerados e o numerador conta somente os pontos pertencentes a linhas com `l ≥ l_min`. No código, `l_min = 2`; portanto, recorrências isoladas entram apenas no denominador. O numerador é implementado pela soma dos comprimentos das linhas aceitas, e o denominador pela contagem de `R(i,j)=1` com `j>i`.

Por construção, `0 ≤ DET ≤ 1`:

- `DET` próximo de 0: a maioria das recorrências é isolada;
- `DET` próximo de 1: a maioria está organizada em segmentos diagonais;
- valores intermediários: combinação das duas estruturas.

O nome “determinismo” é histórico, mas deve ser interpretado com cautela. Um DET alto representa organização diagonal do RP sob os parâmetros escolhidos; ele não demonstra sozinho que o processo gerador é determinístico nem identifica uma patologia. A comparação só é válida quando pré-processamento, `m`, `τ`, `ε`, `l_min` e janela de exclusão são mantidos constantes.

> Zbilut, J. P., & Webber Jr., C. L. (1992). Embeddings and delays as derived from quantification of recurrence plots. *Physics Letters A*, 171(3–4), 199–203.

### Escolha dos parâmetros neste projeto

| Parâmetro | Valor | Função | Consequência |
|---|---:|---|---|
| Frequência de amostragem | 256 Hz | Define a resolução temporal | 3,90625 ms por amostra |
| `N` | 88 ou 240 | Tamanho da janela | ~344 ms ou ~937 ms |
| `m` | 1 | Dimensão de imersão | Recorrência de amplitude, sem reconstrução multidimensional |
| `τ` | 1 | Atraso | Irrelevante quando `m=1` |
| Métrica | Euclidiana em 1D | Distância entre estados | Equivale ao valor absoluto |
| `ε` | 8 unidades de amplitude | Raio da vizinhança | Controla a densidade de recorrências |
| `l_min` | 2 | Menor diagonal aceita | Descarta pontos isolados do numerador |
| Exclusão temporal | Apenas `i=j` | Remove autorrecorrência | Equivale a Theiler = 1 no PyRQA |
| Limiar de classificação | DET > 0,70 | Separa os dois sinais do experimento | Valor empírico, não clínico |

O limiar `ε` é particularmente crítico: valores pequenos tornam a matriz esparsa e valores grandes fazem quase todos os pares recorrerem, elevando artificialmente a taxa de recorrência e podendo alterar o DET. Em estudos comparativos, é comum normalizar os sinais e escolher `ε` por uma fração da dispersão ou fixar uma taxa de recorrência. Aqui, `ε=8` foi calibrado para os sinais sintéticos e deve ser recalibrado antes do uso com dados reais.

A exclusão apenas da LOI é suficiente para reproduzir este experimento e sua validação no PyRQA. Para EEG real, pontos temporalmente próximos podem recorrer apenas por autocorrelação e continuidade do sinal. Uma janela de Theiler maior, escolhida por critério como tempo de autocorrelação, pode ser necessária para excluir essa vizinhança da LOI.

### Aplicação a EEG e epilepsia

Thomasson et al. (2001) aplicaram RQA a EEG epiléptico e relataram que variáveis de recorrência podem distinguir mudanças pré-ictais da atividade de fundo em seus dados. O trabalho é relevante porque emprega RQA em sinais biológicos não estacionários, mas não autoriza concluir que qualquer aumento de DET seja uma crise ou sincronização neuronal.

Lehnertz e Elger (1998) apresentaram evidências de alterações pré-ictais por análise não linear de EEG intracraniano. Trata-se de evidência experimental de que a dinâmica cerebral pode mudar antes da manifestação clínica, e não de prova universal de previsão de crises nem de validação específica do limiar usado neste projeto.

Acharya et al. (2011) extraíram dez medidas de RQA de segmentos de EEG e compararam sete classificadores para separar classes normal, interictal e ictal. O melhor resultado reportado com SVM teve eficiência média de 95,6% no conjunto estudado. A relação com este projeto é metodológica: características de recorrência podem alimentar uma decisão automática. A implementação presente é mais simples, pois utiliza somente DET, sinais sintéticos e um limiar fixo; por isso, sua acurácia clínica não foi avaliada.

Casdagli (1997) e a revisão de Marwan et al. (2007) discutem a interpretação de estruturas recorrentes, aplicações a séries ruidosas e não estacionárias e dificuldades práticas como seleção de limiar, imersão e correções temporais.

### Como relacionar o resultado a uma condição neurológica

O RQA não produz diretamente um diagnóstico. Ele transforma cada janela de EEG em **características quantitativas**. A associação com uma patologia é aprendida ou estimada comparando essas características com um padrão de referência clínico (*ground truth*), obtido por neurologistas a partir do EEG completo, manifestações clínicas e outros exames.

A cadeia de uso é:

```text
EEG do paciente
    → pré-processamento e remoção de artefatos
    → divisão em janelas
    → matriz de recorrência
    → características RQA (DET, RR, Lmax, entropia, LAM etc.)
    → modelo estatístico ou classificador validado
    → probabilidade/alerta
    → interpretação pelo neurologista junto ao contexto clínico
```

#### Exemplo: identificação de atividade epiléptica

Durante determinadas crises, a atividade de populações neuronais pode tornar-se mais rítmica e organizada. Em uma janela de EEG isso pode gerar mais diagonais longas no RP e, consequentemente, aumentar o DET. Um sistema experimental poderia operar assim:

1. registrar EEGs reais e marcar, por revisão de especialistas, janelas **interictais**, **pré-ictais**, **ictais** e com **artefatos**;
2. aplicar o mesmo pré-processamento a todas as janelas;
3. calcular DET e outras medidas de recorrência por canal e por janela;
4. usar parte dos pacientes para ajustar um classificador e seus limiares;
5. testar em pacientes que não participaram do ajuste;
6. apresentar a saída como probabilidade ou alerta de atividade ictal, e não como diagnóstico definitivo.

Exemplo meramente ilustrativo de regra após validação em dados clínicos:

```text
Se DET e outras características permanecerem na região ictal
durante várias janelas consecutivas:
    gerar "suspeita de atividade ictal"
senão:
    manter monitoramento
```

Não se deve transferir para pacientes a regra `DET > 0,70` deste projeto. O limiar clínico teria de ser estimado com uma curva ROC ou procedimento equivalente, equilibrando:

- **sensibilidade**: proporção de crises corretamente detectadas;
- **especificidade**: proporção de períodos não ictais corretamente rejeitados;
- **valor preditivo positivo**: fração dos alertas que correspondem de fato a crises;
- **taxa de falsos alarmes por hora** e atraso de detecção;
- desempenho em pacientes e sessões não utilizados no treinamento.

#### Por que DET isolado não identifica uma patologia

DET alto significa que existem padrões diagonais organizados sob uma configuração específica do RQA. Esse resultado também pode ocorrer em sono rítmico, estimulação periódica, artefato muscular, movimento, interferência elétrica ou outros ritmos fisiológicos. Da mesma forma, nem toda crise produz necessariamente aumento de DET em todos os canais.

Por isso, a associação patológica normalmente utiliza várias informações:

- múltiplas medidas RQA, não apenas DET;
- potência em bandas, morfologia de pontas e outras características do EEG;
- evolução temporal em várias janelas consecutivas;
- distribuição espacial entre canais;
- histórico, sintomas, medicação e avaliação clínica do paciente.

Neste projeto, `anomalia_transicao_estado` deve ser interpretada como **saída de um detector experimental de mudança de regime no sinal sintético**. Para uma aplicação real, a saída mais correta seria algo como `alerta_atividade_ictal`, acompanhado de probabilidade e incerteza. O diagnóstico de epilepsia é mais amplo que detectar uma janela ictal: requer avaliação clínica especializada e não pode ser estabelecido somente por este algoritmo.

### Escopo científico e limitações

Para apresentação acadêmica, as conclusões devem ser delimitadas:

1. o experimento demonstra que a implementação calcula corretamente DET e separa **dois sinais sintéticos previamente construídos**;
2. os valores 0,55 e 0,89 não são estimativas populacionais e não podem ser generalizados para pacientes;
3. “interictal” e “ictal” são rótulos dos modelos sintéticos, não diagnósticos;
4. o limiar 0,70 foi escolhido para o experimento, não derivado de validação clínica;
5. um estudo com EEG real exigiria filtragem e remoção de artefatos, normalização, definição de janelas, seleção de canais, ajuste de parâmetros, separação entre treino e teste e validação em diferentes pacientes;
6. DET isolado não distingue necessariamente crise, sono, artefato muscular ou outra atividade rítmica.

Assim, a contribuição do trabalho é demonstrar, sob restrições severas de memória e processamento, uma implementação reproduzível do núcleo matemático da RQA e sua validação contra uma biblioteca científica independente.

### Ferramentas de referência

- **CRP Toolbox (MATLAB)** — ferramenta para construção e quantificação de recurrence plots, associada ao trabalho de Norbert Marwan e colaboradores.
- **PyRQA (Python)** — implementação paralela baseada em OpenCL para séries longas. Neste projeto, `rqa_validacao.py` usa PyRQA como implementação independente e verifica a equivalência numérica para os mesmos parâmetros; isso valida o cálculo de software, não o desempenho clínico do classificador.

> Rawald, T., Sips, M., & Marwan, N. (2017). PyRQA—Conducting recurrence quantification analysis on very long time series efficiently. *Computers & Geosciences*, 104, 101–108.

---

## Três arquivos, mesmo algoritmo

O projeto possui três versões do mesmo algoritmo, cada uma otimizada para um contexto diferente:

| Arquivo | N | Finalidade | Saída | Representação da matriz | Aritmética do RQA |
|---|---|---|---|---|---|
| `rqa_portable.c` | 88 | Apresentação e validação | `printf` no terminal | `uint8_t [88][88]` (1 byte/elemento) | `float` (fabsf) |
| `rqa_embedded.c` | 88 | Execução em microcontrolador sem FPU | Inspeção por debugger | `uint8_t [88][88]` (1 byte/elemento) | `float` (fabsf) |
| `rqa_embedded_optimized.c` | 240 | Máximo desempenho em microcontrolador | Inspeção por debugger | `uint8_t [240][30]` (1 bit/elemento) | Ponto fixo Q8.8 (`int32_t`) |

O núcleo do algoritmo (`calcular_rqa`) é logicamente idêntico nos três arquivos. Todos respeitam as restrições acadêmicas: sem alocação dinâmica, sem recursão, sem chamada a `malloc`/`calloc`/`free`.

---

## Estruturas de dados

### Versão padrão (`rqa_portable.c` e `rqa_embedded.c`)

As duas estruturas de dados centrais do algoritmo são declaradas como variáveis estáticas globais, alocadas em tempo de compilação no segmento `.bss`:

```c
/* Vetor do sinal EEG — 88 amostras de ponto flutuante */
static float   eeg_sinal[N];                /* 88 × 4 bytes = 352 bytes  */

/* Matriz de recorrência — N×N bytes, cada elemento é 0 ou 1 */
static uint8_t matriz_recorrencia[N][N];     /* 88 × 88 × 1 = 7744 bytes */

/* Total em .bss: 352 + 7744 = 8096 bytes                               */
```

O vetor `eeg_sinal` armazena as N amostras do sinal sintético em ponto flutuante (`float`, 32 bits IEEE 754). As amplitudes são interpretadas em microvolts para manter uma escala plausível, mas não provêm de um sistema de aquisição calibrado.

A `matriz_recorrencia` é uma matriz quadrada N×N onde cada elemento vale 0 ou 1. Ela é preenchida na Etapa 1 do algoritmo e lida na Etapa 2. Os dois laços aninhados (Etapa 1: `for i` × `for j`; Etapa 2: `for i` × `for j` × `for k`) operam diretamente sobre esta matriz.

### Versão otimizada (`rqa_embedded_optimized.c`)

A versão otimizada substitui ambas as estruturas por representações mais compactas:

```c
/* Sinal em ponto fixo Q8.8 — cada amostra é int16_t (2 bytes) */
static int16_t sinal_fixo[N];              /* 240 × 2 bytes =  480 bytes */

/* Matriz bit-packed — 1 bit por elemento, 8 elementos por byte */
static uint8_t matriz_bits[N][BYTES_POR_LINHA];  /* 240 × 30 = 7200 bytes */

/* Total em .bss: 480 + 7200 = 7680 bytes                               */
```

**Ponto fixo Q8.8:** o sinal é convertido de `float` para `int16_t` multiplicando por 256. Os 8 bits superiores representam a parte inteira e os 8 bits inferiores representam a parte fracionária. A comparação de distância na Etapa 1 usa subtração e valor absoluto inteiros nativos de 32 bits, eliminando toda emulação de ponto flutuante do laço O(N²).

**Bit-packing:** cada elemento R(i,j) da matriz ocupa 1 bit em vez de 1 byte. O bit j da linha i é armazenado na posição `(j % 8)` do byte `matriz_bits[i][j / 8]`. O acesso usa operações bitwise:

```c
/* Ligar o bit R(i,j) — operação OR com máscara de 1 bit */
static inline void bit_set(uint_fast16_t i, uint_fast16_t j)
{
    matriz_bits[i][j >> 3] |= (uint8_t)(1u << (j & 7u));
}

/* Ler o bit R(i,j) — shift direito + AND com 1 */
static inline uint_fast8_t bit_get(uint_fast16_t i, uint_fast16_t j)
{
    return (uint_fast8_t)((matriz_bits[i][j >> 3] >> (j & 7u)) & 1u);
}
```

### Variáveis de inspeção

Nas duas versões embarcadas, o resultado do algoritmo é armazenado em variáveis globais `volatile`, acessíveis pelo debugger:

```c
volatile float   determinismo_det;           /* DET calculado (0.0 a 1.0)  */
volatile uint8_t anomalia_transicao_estado;  /* 1 se DET > 0.70, senão 0   */
volatile uint8_t simular_anomalia;           /* 0 = repouso, 1 = crise     */
```

O qualificador `volatile` obriga cada acesso C a ser realizado no objeto em memória, permitindo que o debugger altere `simular_anomalia` enquanto o programa roda e observe os resultados. Ele não torna a operação atômica nem substitui mecanismos de sincronização. Na versão portátil, essas variáveis não são `volatile`, pois são controladas pelo próprio `main()` e impressas no terminal.

### Mapa de memória (8192 bytes de SRAM)

| Estrutura | Versão padrão (N=88) | Versão otimizada (N=240) |
|---|---|---|
| Vetor do sinal | `float [88]` — 352 B | `int16_t [240]` — 480 B |
| Matriz de recorrência | `uint8_t [88][88]` — 7744 B | `uint8_t [240][30]` — 7200 B |
| Variáveis de inspeção | 6 B | 6 B |
| **Total .bss** | **8102 B** | **7686 B** |
| **Sobra para pilha** | **~90 B** | **506 B** |
| Pares analisados (N²) | 7744 | 57600 (7,4× mais) |

---

## Modelo sintético de EEG

O sinal é um modelo didático inspirado em bandas espectrais reconhecidas no EEG de escalpo, amostrado a Fs = 256 Hz. Ele combina senoides e ruído pseudoaleatório e, portanto, não reproduz toda a morfologia, a variabilidade entre pacientes, os artefatos nem as propriedades espaciais de um EEG clínico multicanal.

`simular_anomalia` é o botão do experimento:

### Modo 0 — Modelo interictal sintético

Composição de cinco bandas com amplitudes típicas de escalpo:

| Banda | Frequência | Amplitude | Papel no sinal |
|---|---|---|---|
| Alfa | 10 Hz | 12 | Dominante em repouso com olhos fechados |
| Theta | 6 Hz | 8 | Atividade de fundo |
| Beta | 21 Hz | 5 | Atividade cortical rápida |
| Gama | 35 Hz | 2 | Baixa amplitude no escalpo |
| Delta | 2 Hz | 2 | Oscilação lenta de fundo |
| Ruído | — | 14 | Biológico + instrumental |

A mistura de frequências não harmônicas e o ruído forte produzem um trecho irregular. As amostras adjacentes divergem com frequência, quebrando as linhas diagonais na matriz de recorrência.

### Modo 1 — Modelo ictal sintético

Padrão clássico de **spike-and-wave a 3 Hz**, observado em crises de ausência:

| Componente | Frequência | Amplitude | Papel no sinal |
|---|---|---|---|
| Fundamental | 3 Hz | 22 | Ritmo ictal dominante (onda lenta) |
| 1º harmônico | 6 Hz | 10 | Contribui para a forma pontiaguda (spike) |
| 2º harmônico | 9 Hz | 5 | Acentua o pico do spike |
| Ruído | — | 5 | Menor perturbação aleatória no modelo ictal |

A soma dos harmônicos aproxima uma forma pontiaguda seguida de componente lenta. O ruído menor foi escolhido para representar, de modo simplificado, um regime mais organizado e favorecer a repetição do padrão. Isso produz linhas diagonais mais longas, mas não constitui um modelo fisiológico completo de sincronização cortical.

---

## A matriz de recorrência

Para cada par de instantes `i` e `j`:

```
R(i,j) = 1  se |x(i) - x(j)| ≤ ε
R(i,j) = 0  se |x(i) - x(j)| > ε
```

`EPSILON` vale `8.0`. Dois instantes recorrem quando a amplitude é parecida, mesmo que estejam longe no tempo.

Exemplo feito à mão, série `[0, 1, 0, 1, 0]`, tolerância `0.5`. A matriz vale 1 quando os valores são iguais:

|  | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0 | 1 | 0 | 1 | 0 | 1 |
| 1 | 0 | 1 | 0 | 1 | 0 |
| 2 | 1 | 0 | 1 | 0 | 1 |
| 3 | 0 | 1 | 0 | 1 | 0 |
| 4 | 1 | 0 | 1 | 0 | 1 |

A diagonal principal é sempre 1, porque `|x(i) - x(i)| = 0`. Ela não descreve o sinal. A matriz é simétrica: `R(i,j)` e `R(j,i)` são o mesmo par. O código preenche a matriz inteira. O DET usa só o triângulo superior, `j > i`.

## O determinismo (DET)

Uma linha diagonal é uma sequência de 1s em que os dois índices andam juntos: `R(i,j)`, `R(i+1,j+1)`, `R(i+2,j+2)`, e assim por diante.

Isso significa: a partir de dois instantes parecidos, os instantes seguintes continuam parecidos. O futuro acompanha o passado. É a definição operacional de previsibilidade usada aqui.

No exemplo, a linha que começa em `(0,2)` segue por `(1,3)` e `(2,4)`. Comprimento 3. O ponto `(0,4)` é um 1 isolado, comprimento 1.

`L_MIN` vale 2. Só linhas de comprimento 2 ou mais entram no numerador. O ponto isolado entra só no denominador.

```
DET = 3 / 4 = 0.75
```

No código a mesma conta é:

```
DET = soma_comprimentos / total_pontos_recorrentes
```

Cada linha é medida uma vez, a partir do primeiro 1 cujo antecessor `(i-1, j-1)` vale 0. Cada ponto entra numa linha só, e o DET fica entre 0 e 1. Um DET acima de 1 indicaria a mesma linha contada várias vezes.

Se `determinismo_det > 0.70`, `anomalia_transicao_estado` vai para 1.

## Por que o modelo ictal produz DET alto

Com `simular_anomalia = 1`, o padrão spike-and-wave de 3 Hz se repete a cada ~85 amostras (256/3). Como esse período quase cabe nas 88 amostras, as diagonais de periodicidade são longas. Os harmônicos (6 Hz e 9 Hz) reforçam a repetitividade. O ruído residual de apenas 5 (contra 14 no repouso) preserva as linhas diagonais.

- DET cerca de **0,89** (825 pontos em linhas / 922 recorrências)
- `anomalia_transicao_estado = 1`

No modo interictal (repouso), cinco bandas não harmônicas mais ruído forte produzem divergências frequentes entre amostras adjacentes. A maioria dos pontos recorrentes é isolada.

- DET cerca de **0,55** (597 / 1087)
- `anomalia_transicao_estado = 0`

---

## Entradas de dados e saídas típicas

### Entrada

O programa não recebe dados externos. A entrada é gerada internamente pela função `gerar_sinal_dinamica()`, que produz um vetor de N amostras de EEG sintético. A única variável de controle é `simular_anomalia`:

| Entrada | Tipo | Valores | Significado |
|---|---|---|---|
| `simular_anomalia` | `uint8_t` | `0` | EEG em repouso (interictal) |
| | | `1` | Modelo ictal sintético |

Na versão embarcada, o operador altera `simular_anomalia` pelo debugger em tempo real. Na versão portátil, o `main()` testa os dois valores sequencialmente.

### Saída

O algoritmo produz duas variáveis de saída, estáveis após cada execução de `calcular_rqa()`:

| Saída | Tipo | Faixa | Significado |
|---|---|---|---|
| `determinismo_det` | `float` | 0,0 a 1,0 | Fração de pontos recorrentes organizados em linhas diagonais |
| `anomalia_transicao_estado` | `uint8_t` | `0` ou `1` | Classificação binária: `1` se DET > 0,70 |

### Valores típicos observados

| Modo | Entrada | DET | Recorrências | Pontos em linhas | Classificação |
|---|---|---|---|---|---|
| Repouso (interictal) | `simular_anomalia = 0` | 0,549 | 1087 | 597 | `anomalia = 0` |
| Modelo ictal sintético | `simular_anomalia = 1` | 0,895 | 922 | 825 | `anomalia = 1` |

Nota: os valores acima referem-se à versão com N=88 (`rqa_portable.c` e `rqa_embedded.c`). A versão otimizada (N=240) produz DET interictal ≈ 0,587 e DET ictal ≈ 0,886. A diferença nos valores absolutos deve-se ao maior número de amostras (janela de ~937 ms em vez de ~344 ms), mas a classificação binária é idêntica.

### Variáveis internas de conferência

A versão portátil imprime duas variáveis internas adicionais que permitem conferir a conta:

| Variável | Tipo | Significado |
|---|---|---|
| `total_pontos_recorrentes` | `uint32_t` | Denominador do DET — total de R(i,j) = 1 com j > i |
| `soma_comprimentos` | `uint32_t` | Numerador do DET — soma dos comprimentos de linhas com L ≥ 2 |
| `razao_exata` | `double` (impressão) | Divisão `soma / total` em 64 bits, para conferência |

---

## Procedimentos de teste

### Teste 1 — Validação manual com exemplo 5×5

O código inclui a função `conferir_exemplo_5x5()` que executa o algoritmo RQA completo sobre uma série conhecida `[0, 1, 0, 1, 0]` com ε = 0,5. A resposta correta é calculável à mão:

**Montagem da matriz (Etapa 1):**

```
Série: [0, 1, 0, 1, 0]    ε = 0.5

     j=0  j=1  j=2  j=3  j=4
i=0 [ 1    0    1    0    1  ]
i=1 [ 0    1    0    1    0  ]
i=2 [ 1    0    1    0    1  ]
i=3 [ 0    1    0    1    0  ]
i=4 [ 1    0    1    0    1  ]
```

**Triângulo superior (j > i):** 4 pontos recorrentes — (0,2), (0,4), (1,3), (2,4).

**Detecção de diagonais (Etapa 2):**

- Diagonal iniciando em (0,2): antecessor (−1,1) fora da matriz → início de linha. Segue por (1,3), (2,4). Comprimento = 3 ≥ L_MIN.
- Ponto (0,4): antecessor (−1,3) fora da matriz → início. Comprimento = 1 < L_MIN → descartado.
- Ponto (1,3): antecessor (0,2) = 1 → **não** é início de linha (já contado na diagonal acima).
- Ponto (2,4): antecessor (1,3) = 1 → **não** é início de linha.

**Resultado:** DET = 3 / 4 = **0,75**

**Critério de aprovação:** `total == 4 AND soma == 3`. A função retorna 0 (sucesso) ou 1 (falha).

### Teste 2 — Modo interictal (repouso)

O teste gera 88 amostras de EEG em repouso (cinco bandas + ruído forte) e verifica que o DET fica **abaixo** do limiar de 0,70:

```
simular_anomalia = 0
→ DET = 0.549218 (597 / 1087)
→ anomalia_transicao_estado = 0
```

**Critério de aprovação:** `determinismo_det ≤ 0.70 AND anomalia_transicao_estado == 0`.

### Teste 3 — Modelo ictal sintético

O teste gera 88 amostras de crise epiléptica (spike-and-wave 3 Hz) e verifica que o DET fica **acima** do limiar de 0,70:

```
simular_anomalia = 1
→ DET = 0.894794 (825 / 922)
→ anomalia_transicao_estado = 1
```

**Critério de aprovação:** `determinismo_det > 0.70 AND anomalia_transicao_estado == 1`.

### Teste 4 — Validação cruzada com Python

O script `rqa_validacao.py` compara três fontes independentes:

| Fonte | Descrição | Precisão |
|---|---|---|
| Manual (Python) | Mesma lógica do código C, reimplementada em Python puro | `float64` |
| PyRQA | Biblioteca acadêmica peer-reviewed, com OpenCL | `float64` via OpenCL |
| Código C esperado | Valores compilados de `rqa_portable.exe` | `float32` |

As três fontes devem concordar em contagens inteiras (`soma` e `total`) e classificação binária (`anomalia`). A diferença no DET é tolerada em até 0,01 (arredondamento float32 vs float64).

**Execução:**

```
pip install PyRQA
py -3 rqa_validacao.py
```

**Critério de aprovação:** todas as comparações passam sem divergência.

### Teste 5 — Compilação sem warnings

Os três arquivos C devem compilar sem nenhum warning com flags rigorosas:

```
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Wconversion rqa_portable.c -o rqa_portable.exe -lm
gcc -std=c11 -Os -Wall -Wextra -Wpedantic -Wconversion rqa_embedded.c -o rqa_embedded.exe -lm
gcc -std=c11 -Os -Wall -Wextra -Wpedantic -Wconversion rqa_embedded_optimized.c -o rqa_embedded_optimized.exe -lm
```

**Critério de aprovação:** zero warnings em todos os três arquivos.

### Resumo dos testes

| # | Teste | O que valida | Critério |
|---|---|---|---|
| 1 | Exemplo 5×5 | Lógica do algoritmo (conta manual conferível) | DET = 3/4 = 0,75 |
| 2 | Modo interictal | DET baixo para sinal irregular | DET ≤ 0,70, alarme = 0 |
| 3 | Modo ictal | DET alto para sinal periódico | DET > 0,70, alarme = 1 |
| 4 | Validação cruzada | Consistência entre C, Python e PyRQA | Diferença DET < 0,01 |
| 5 | Compilação limpa | Ausência de bugs detectáveis pelo compilador | 0 warnings |

---

## Como executar a versão de apresentação

Na pasta do projeto, com GCC:

```
gcc -std=c11 -O2 rqa_portable.c -o rqa_portable.exe -lm
.\rqa_portable.exe
```

A saída esperada é:

```
exemplo_5x5 pontos_em_linhas = 3
exemplo_5x5 recorrencias = 4
exemplo_5x5 razao_exata = 0.75

simular_anomalia = 0
pontos_em_linhas = 597
recorrencias = 1087
razao_exata = 0.5492180312787488
determinismo_det = 0.549218059
anomalia_transicao_estado = 0

simular_anomalia = 1
pontos_em_linhas = 825
recorrencias = 922
razao_exata = 0.89479392624728848
determinismo_det = 0.894793928
anomalia_transicao_estado = 1

conferencia: exemplo 5x5 = 3/4 e os dois modos em lados opostos de 0.70
```

`razao_exata` é a divisão dos dois inteiros, impressa com 17 algarismos. `determinismo_det` é o `float` de 32 bits, impresso com 9 algarismos.

---

## Referências

1. Takens, F. (1981). Detecting strange attractors in turbulence. In D. A. Rand & L.-S. Young (Eds.), *Dynamical Systems and Turbulence, Warwick 1980* (Lecture Notes in Mathematics, Vol. 898, pp. 366–381). Springer. https://doi.org/10.1007/BFb0091924

2. Eckmann, J.-P., Kamphorst, S. O., & Ruelle, D. (1987). Recurrence plots of dynamical systems. *Europhysics Letters*, 4(9), 973–977. https://doi.org/10.1209/0295-5075/4/9/004

3. Zbilut, J. P., & Webber Jr., C. L. (1992). Embeddings and delays as derived from quantification of recurrence plots. *Physics Letters A*, 171(3–4), 199–203. https://doi.org/10.1016/0375-9601(92)90426-M

4. Casdagli, M. C. (1997). Recurrence plots revisited. *Physica D: Nonlinear Phenomena*, 108(1–2), 12–44. https://doi.org/10.1016/S0167-2789(97)82003-9

5. Lehnertz, K., & Elger, C. E. (1998). Can epileptic seizures be predicted? Evidence from nonlinear time series analysis of brain electrical activity. *Physical Review Letters*, 80(22), 5019–5022. https://doi.org/10.1103/PhysRevLett.80.5019

6. Thomasson, N., Hoeppner, T. J., Webber Jr., C. L., & Zbilut, J. P. (2001). Recurrence quantification in epileptic EEGs. *Physics Letters A*, 279(1–2), 94–101. https://doi.org/10.1016/S0375-9601(00)00815-X

7. Marwan, N., Romano, M. C., Thiel, M., & Kurths, J. (2007). Recurrence plots for the analysis of complex systems. *Physics Reports*, 438(5–6), 237–329. https://doi.org/10.1016/j.physrep.2006.11.001

8. Acharya, U. R., Sree, S. V., Chattopadhyay, S., Yu, W., & Ang, P. C. A. (2011). Application of recurrence quantification analysis for the automated identification of epileptic EEG signals. *International Journal of Neural Systems*, 21(3), 199–211. https://doi.org/10.1142/S0129065711002808

9. Rawald, T., Sips, M., & Marwan, N. (2017). PyRQA—Conducting recurrence quantification analysis on very long time series efficiently. *Computers & Geosciences*, 104, 101–108. https://doi.org/10.1016/j.cageo.2016.11.016
