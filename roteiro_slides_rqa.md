# Roteiro de slides — RQA em float na NUCLEO-F030R8

Guia para montar a apresentação a partir de `rqa_embedded_serial.c`. Cada bloco é um slide. **No slide** é o texto que a plateia lê. **Falar** desenvolve esse texto.

A tabela é a saída de `rqa_validacao.py`, gravada em `tabela_validacao.txt`. Ela mostra que três implementações da mesma conta concordam. O limiar 0,70 não separa, neste arquivo, janela de crise e janela de fundo.

---

## Slide 1 — Título e objetivo

**No slide**

RQA de EEG na NUCLEO-F030R8. O PC envia uma janela de EEG pela serial. A placa, um STM32F030R8 sem unidade de ponto flutuante e com 8 KiB de RAM, calcula o determinismo em `float` e devolve um número entre 0 e 1.

Esse número responde se amplitudes parecidas reaparecem em sequência. Se ele passa de 0,70, a placa marca alarme. O alarme é o registro dessa comparação. Não é diagnóstico e não prevê crise.

**Falar**

Abrir dizendo o que a placa faz e o que ela não faz. Ela não gera senoide nem lê o arquivo de EEG. Recebe 88 amostras já recortadas, compara cada amostra com as outras e resume a matriz num único valor, o determinismo. A utilidade do trabalho é executar essa conta dentro de 8 KiB e conferir o resultado com a PyRQA. O corte 0,70 só diz se o número ficou alto. A leitura clínica do EEG continua fora deste experimento.

---

## Slide 2 — O que o determinismo mede

**No slide**

A janela é uma sequência `u(0)` … `u(87)`. Com dimensão 1 e atraso 1, cada estado é a própria amplitude: não há vetor de amostras atrasadas.

Dois instantes são recorrentes quando a diferença absoluta das amplitudes não passa de 8:

`R(i,j) = 1` se `|u(i) − u(j)| ≤ 8`, e `0` no caso contrário.

O determinismo é a fração dos pares recorrentes que formam diagonal de pelo menos dois pontos:

`DET = (pontos em diagonais de comprimento ≥ 2) / (pares recorrentes com j > i)`.

O núcleo da primeira comparação, no firmware, é este:

```c
if (fabsf(rqa.eeg_sinal[i] - rqa.eeg_sinal[j]) <= EPSILON) {
    rqa.matriz_recorrencia[i][j] = 1u;
    if (j > i)
        total++;
}
```

`EPSILON` vale 8. `total` é o denominador e só cresce acima da diagonal.

A diagonal principal fica de fora, porque `|u(i) − u(i)|` é sempre zero. Um ponto recorrente isolado entra só no denominador. DET perto de 1 significa que as recorrências vêm em trechos. DET perto de 0 significa que quase todas são isoladas.

**Falar**

Percorrer um exemplo curto. Nas amostras `[0, 1, 20]`, `|0 − 1| = 1` marca recorrência e `|0 − 20| = 20` não marca. O par `(0, 1)` existe, mas não continua em `(1, 2)`, então o comprimento é 1 e o numerador continua 0. Se `(1, 2)` também fosse recorrente, os dois pontos entrariam no numerador uma vez só. `i` e `j` são dois momentos do mesmo sinal, não dois canais.

---

## Slide 3 — Parâmetros da conta

**No slide**

| Parâmetro | Valor | Efeito na conta |
|---|---|---|
| Amostras | 88 | Tamanho da janela e da matriz |
| Duração | 0,34375 s | 88 amostras a 256 Hz |
| Raio | 8 | Até essa diferença as amplitudes contam como parecidas |
| Menor diagonal | 2 | Recorrência de um ponto só não entra no numerador |
| Alarme | DET > 0,70 | Abaixo ou igual a 0,70 o alarme fica 0 |
| Métrica | valor absoluto | Distância entre duas amplitudes |

A placa, o Python e a PyRQA usam esta mesma tabela. O raio 8 só faz sentido depois da normalização, quando o desvio-padrão da janela foi levado a 20.

**Falar**

Dizer que nenhum desses números foi estimado de novo a cada janela. Estão fixos para o teste inteiro poder comparar as três contas. Aumentar o raio enche a matriz de uns. Diminuir a diagonal mínima faz pontos isolados entrarem no determinismo. O 0,70 é um corte empírico deste experimento, aplicado depois que o DET já foi calculado.

---

## Slide 4 — Tipos usados na conta

**No slide**

O tipo acompanha o maior valor da variável e o formato em que esse valor sai pela serial.

Os índices `i`, `j` e `k`, e o comprimento de uma diagonal, são `uint_fast16_t`. Nenhum passa de 88. O tipo pede pelo menos 16 bits e, neste processador, ocupa o registrador de 32 bits, sem cortar o índice a cada soma.

O numerador e o denominador são `uint32_t`. O maior denominador possível é `88 × 87 / 2 = 3828`, que caberia em 16 bits. Eles têm 32 bits porque a resposta envia cada contador em 4 bytes.

O determinismo é `volatile float`: a divisão de dois inteiros, um valor entre 0 e 1, escrito na RAM para os 4 bytes da resposta serem lidos dali. O alarme é `volatile uint8_t`, porque só vale 0 ou 1 e ocupa o último byte. `volatile` obriga essa escrita a chegar na memória antes da transmissão.

A janela é `float[88]`, 352 bytes, o mesmo `float32` recebido. A matriz é `uint8_t[88][88]`, 7744 bytes, porque cada célula guarda apenas 0 ou 1.

**Falar**

Ligar cada tipo a uma pergunta. “Quantas amostras já percorri?” cabe num índice até 88. “Quantos pares recorreram?” é um contador enviado inteiro. “Qual fração formou diagonal?” é o `float`. “Passou de 0,70?” é um byte. “Este par é recorrente?” é um byte na matriz. Um inteiro de 4 bytes por célula não caberia na SRAM.

---

## Slide 5 — Memória

**No slide**

A SRAM do chip tem 8192 bytes. A conta não pede memória durante a execução: o heap reservado é 0.

| O que fica na RAM | Bytes |
|---|---:|
| Determinismo, dois contadores e alarme | 16 |
| 88 amostras em `float` | 352 |
| Matriz 88×88, um byte por par | 7744 |
| Soma da conta | 8112 |
| Pilha reservada | 72 |
| Flash da imagem | 3756 |

A pilha cresce do fim da SRAM para baixo. A subtração em software, usada na comparação das amplitudes, ocupa esses 72 bytes. Por isso os quatro resultados ficam no começo do bloco, longe desse contato. A matriz sozinha toma quase toda a memória. Não há espaço para guardá-la com um tipo maior.

**Falar**

Somar em voz alta: 16 + 352 + 7744 = 8112, mais 72 de pilha, dentro de 8192. O que a plateia precisa reter é que a matriz decide o tamanho do programa. Oito mil bytes de RAM são o motivo de cada recorrência ser um byte, do heap ser zero e de a janela ter 88 amostras nesta versão.

---

## Slide 6 — A RAM acaba no último byte, e o programa ainda assim liga

**No slide**

A objeção é esta: se o programa ocupa os 8192 bytes de SRAM, não deveria compilar nem rodar. Ocupar a memória toda não é o mesmo que passar do fim dela. O ligador recusa a imagem só quando a soma das seções pede 8193 bytes ou mais. Uma imagem de 8192 bytes cabe, e esta cabe.

A conta pede 8112 bytes, já reservados do reset até o fim do programa. O heap é 0: em execução ninguém pede mais RAM. Sobram 72 bytes, e é exatamente essa a pilha reservada. O relógio do CMSIS e o alinhamento fecham os 8 bytes que faltam para 8192. A imagem ligada da configuração Embedded_Serial fica em 8 KiB de RAM, 100%, e em 3756 bytes de Flash, de 64 KiB. O ligador aceita esse mapa.

A matriz não entra na pilha. Se entrasse, 7744 bytes seriam pedidos em cima dos 72 e o enlace falharia, ou a pilha escreveria fora da reserva. Ela é estática:

```c
static struct {
    volatile float determinismo_det;          /* endereço baixo */
    volatile uint32_t pontos_em_linhas;
    volatile uint32_t recorrencias;
    volatile uint8_t anomalia_transicao_estado;
    float eeg_sinal[88];
    uint8_t matriz_recorrencia[88][88];       /* endereço alto, junto da pilha */
} rqa;
```

A pilha cresce do fim da SRAM para baixo, em direção ao final da matriz. Os resultados ficam no começo do bloco. Nos 72 bytes só entram os quadros de chamada. A cadeia mais funda é `main` → `calcular_rqa` → `fabsf`. Não há recursão, não há `malloc` e as variáveis locais são só índices e contadores:

```c
uint_fast16_t i, j, k, comprimento;
uint32_t total = 0u, soma = 0u;
```

No Cortex-M0, sem FPU, a subtração em software de `fabsf` usa esses 72 bytes, e a reserva tem esse tamanho. Por isso a conta roda com a SRAM cheia: tudo o que ela precisa já foi colocado dentro dos 8192 na ligação, e a execução não aloca mais nada.

**Falar**

Responder ao professor com a diferença entre cheio e estourado. Cheio compila. Estourado não. Desenhar a SRAM: 8112 da conta, 72 de pilha, o restante para o relógio e o alinhamento, soma 8192. Dizer que a Flash nem entra nessa conta: sobram dezenas de kibibytes. O programa que rodou na placa é essa imagem, não uma imagem que pediu memória além do chip. Se o bloco de 24 bytes do `crtbegin` entrar no enlace, aí sim a SRAM passa de 8192 e o ligador recusa. Esta configuração deixa esse bloco de fora, de propósito, para a soma fechar.

---

## Slide 7 — Como a conta percorre a matriz

**No slide**

A primeira passada visita os 88×88 pares. Se a diferença absoluta das amplitudes é no máximo 8, a célula vira 1. O denominador só aumenta quando `j > i`, então cada par de instantes diferentes entra uma vez e a diagonal principal não entra.

A segunda passada olha só essa metade de cima. Uma diagonal é a sequência `(i, j)`, `(i+1, j+1)`, `(i+2, j+2)`, enquanto as células continuam em 1. A medida começa no primeiro 1 da linha, aquele cujo antecessor `(i−1, j−1)` é 0. Assim a mesma linha não é somada de novo em cada ponto. Se o comprimento chega a 2, todos esses pontos entram no numerador.

No fim, uma divisão produz o determinismo. Se não houve par recorrente, o determinismo fica 0, para não dividir por zero. O alarme compara esse resultado com 0,70.

```c
if (comprimento >= L_MIN)      /* L_MIN = 2 */
    soma += comprimento;

rqa.determinismo_det = (float)soma / (float)total;
rqa.anomalia_transicao_estado =
    (rqa.determinismo_det > LIMIAR) ? 1u : 0u;   /* LIMIAR = 0,70 */
```

**Falar**

Voltar ao exemplo do slide 2 e apontar as duas passadas. A primeira responde “são parecidos?”. A segunda responde “essa semelhança continua nos instantes seguintes?”. Medir a diagonal só no começo é o que impede um trecho de cinco pontos virar 5+4+3+2+1. A recepção também não acumula quadro: cada byte é lido, e os 88 floats ocupam direto o vetor da janela.

---

## Slide 8 — Dados de teste

**No slide**

Os testes usam um EEG de escalpo da base CHB-MIT, distribuída no PhysioNet, amostrado a 256 Hz.

O registro é o `chb01_03`. O canal lido é o primeiro, a montagem FP1-F7. O arquivo começa às 13:43:04. O sumário do caso anota uma crise de 2996 s a 3036 s, isto é, das 14:33:00 às 14:33:40, quarenta segundos.

Desses sinais saem 20 janelas de 88 amostras. A placa não abre o EDF. O computador lê o registro, escolhe os trechos e manda cada janela pronta.

**Falar**

Situar a base: é EEG pediátrico de escalpo, com as crises marcadas no sumário por início e fim em segundos. Este experimento usa um arquivo e um canal. A anotação diz onde a crise foi marcada; ela não entra no cálculo. Serve depois, na hora de ler o alarme.

---

## Slide 9 — Como as janelas foram preparadas

**No slide**

A preparação acontece no computador, antes de qualquer byte ir para a placa.

1. A crise mais longa do arquivo é o intervalo de 40 s já citado.
2. Dez janelas são colocadas dentro desse intervalo, com as 88 amostras inteiras entre o início e o fim anotados.
3. Dez janelas são colocadas fora, a pelo menos 60 s da crise, também sem cruzar a fronteira.
4. De cada trecho tira-se a média. O desvio-padrão é então levado a 20, em `float32`. Uma janela quase constante vira zeros.
5. Os identificadores 0 a 9 ficam fora da crise. Os identificadores 10 a 19 ficam dentro.
6. O rótulo 0 ou 1 acompanha a janela no texto do teste e não é transmitido.

**Falar**

Explicar por que a janela tem de caber inteira de um lado da fronteira: um trecho que começa antes da crise e termina dentro dela não seria nem “fora” nem “dentro”. A folga de 60 s afasta as janelas de fundo da anotação. A normalização é o passo que permite usar o mesmo raio 8 em trechos de amplitudes diferentes: depois dela, 8 significa oito unidades de uma escala em que o desvio-padrão vale 20.

---

## Slide 10 — Como a janela chega ao microcontrolador

**No slide**

O envio usa a USART2 ligada ao ST-Link, a 115200 baud, 8 bits, sem paridade, um stop, little-endian.

O computador transmite 354 bytes por janela. Os dois primeiros são os caracteres `R` e `Q`. Em seguida vêm as 88 amostras, cada uma em 4 bytes `float32`, na ordem do tempo. A placa descarta tudo o que chega antes de ver `R` imediatamente seguido de `Q`. Um `R` sozinho não abre a janela.

Quando as 88 amostras chegaram, a placa calcula e responde 15 bytes: de novo `R` e `Q`, o determinismo em 4 bytes, o numerador em 4 bytes, o denominador em 4 bytes e o alarme em 1 byte. O alarme é 1 somente se o determinismo é maior que 0,70.

Os bits do `float` viajam como estão. Não há nova escala na serial.

O laço da placa é só isso:

```c
usart2_iniciar();
for (;;) {
    receber_janela();   /* espera 'R''Q' e lê 88 floats */
    calcular_rqa();     /* matriz, DET e alarme */
    enviar_resultado(); /* 'R''Q' + 4 números, 15 bytes */
}
```

**Falar**

Mostrar a viagem de uma janela: texto preparado no PC, bytes `RQ` mais 352 bytes de amostras, espera na placa, resposta de 15 bytes. O `for` não termina: cada volta substitui a janela anterior. O sincronismo existe porque a serial é um fluxo contínuo: sem `RQ`, a placa não saberia qual byte é o começo da primeira amostra. O rótulo da crise fica no computador para a tabela. A placa só vê amplitudes.

---

## Slide 11 — Como o teste decide que a conta confere

**No slide**

As 20 janelas preparadas passam por três cálculos com os mesmos parâmetros:

- o algoritmo em C, executado no computador;
- a mesma sequência de comparações reescrita em Python, em `float32`;
- a PyRQA, com dimensão 1, atraso 1, raio 8, correção de Theiler 1 e diagonal mínima 2.

A PyRQA conta os dois triângulos da matriz. As contagens usadas na comparação são divididas por 2, para ficarem no mesmo triângulo `j > i` da placa.

O trecho só fica OK se os três concordam no numerador, no denominador e no alarme, e se a diferença entre os determinismos é no máximo `1e-5`.

O mesmo texto de janelas pode seguir pela serial. A resposta da placa é então comparada com a PyRQA por esse critério. A tabela dessa conferência, gravada no computador, é o último slide.

**Falar**

Separar duas perguntas. A primeira é “as implementações calculam o mesmo DET?”. A tabela do fim responde essa pergunta. A segunda é “o alarme 1 coincide com a crise anotada?”. Essa pergunta é o slide seguinte. Um trecho pode estar OK e, ainda assim, ter alarme 1 numa janela de fundo.

---

## Slide 12 — O que o alarme fez com o rótulo

**No slide**

O rótulo diz onde a janela foi recortada. O alarme diz se o DET passou de 0,70. São comparações diferentes.

| Onde a janela está | Quantidade | Faixa do DET | Quantas passaram de 0,70 |
|---|---:|---|---:|
| Fora da crise, ids 0 a 9 | 10 | 0,760 a 0,993 | 10 |
| Dentro da crise, ids 10 a 19 | 10 | 0,540 a 0,989 | 8 |

As duas janelas abaixo do corte, ids 13 e 15, estão dentro da crise. Nenhuma janela de fundo ficou abaixo.

Com este raio e esta normalização, o corte 0,70 não separa o trecho anotado como crise do trecho anotado como fundo. O resultado 20/20 mede a concordância entre C, Python e PyRQA. Não é uma taxa de detecção.

**Falar**

Contar com a plateia: das dez janelas de fundo, dez receberam alarme 1. Das dez de dentro, oito também. Um classificador que marcasse “crise” em quase tudo acertaria boa parte das janelas de dentro e erraria todas as de fora. Por isso o slide não apresenta acurácia. O que ficou demonstrado é que as três contas produziram esses valores, inclusive os dois determinismos de 0,54.

---

## Slide 13 — Limites do que foi testado

**No slide**

O teste cobre um canal, um arquivo e vinte janelas. Cada janela dura cerca de um terço de segundo; a crise anotada dura 40 s. Uma janela não representa a crise inteira.

A recorrência é de amplitude, porque a dimensão de imersão é 1. O raio 8 e o corte 0,70 foram fixados para sinais cujo desvio-padrão foi forçado a 20. Outro raio, ou uma janela sem essa normalização, produz outra matriz.

A SRAM está tomada pela matriz de bytes e o heap é zero. Esta versão não aumenta o número de amostras sem mudar a forma de guardar a matriz.

**Falar**

Usar este slide para limitar a frase final. Dá para afirmar que a conta coube na placa e que a referência numérica concordou nas vinte janelas. Não dá para afirmar que o método detecta epilepsia, nem que o mesmo corte valeria para outro canal ou para outra escala de amplitude.

---

## Slide 14 — Fecho

**No slide**

A placa recebe 88 amostras pela serial, monta a matriz de recorrência e devolve o determinismo, os dois contadores e o alarme. A conta ocupa 8112 bytes de um chip com 8192, com heap zero.

Nas vinte janelas do `chb01_03`, o C, o Python e a PyRQA deram o mesmo determinismo e o mesmo alarme. O corte 0,70 marcou as dez janelas de fundo e oito das dez janelas de dentro da crise.

A conta confere. O corte não separa a anotação deste registro.

**Falar**

Repetir as duas frases e, em seguida, mostrar a tabela. A concordância está linha por linha no último slide. A contagem do alarme está no slide 12. Se perguntarem o que mudaria o alarme, a resposta é o raio, a normalização ou o próprio 0,70, não um erro entre as três contas: a partir das mesmas amostras, elas já concordavam.

---

## Slide 15 — Tabela dos 20 trechos

**No slide**

Registro `chb01_03`, canal FP1-F7. Cada janela dura 0,34375 s. As três contas concordam nas 20 linhas.

Ids 0 a 9: fora da crise. Ids 10 a 19: dentro da crise.

| id | DET C | DET Python | DET PyRQA | alarme C | alarme Python | alarme PyRQA |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0,932862 | 0,932862 | 0,932862 | 1 | 1 | 1 |
| 1 | 0,948405 | 0,948405 | 0,948405 | 1 | 1 | 1 |
| 2 | 0,988186 | 0,988186 | 0,988186 | 1 | 1 | 1 |
| 3 | 0,964912 | 0,964912 | 0,964912 | 1 | 1 | 1 |
| 4 | 0,986538 | 0,986538 | 0,986538 | 1 | 1 | 1 |
| 5 | 0,992821 | 0,992821 | 0,992821 | 1 | 1 | 1 |
| 6 | 0,938202 | 0,938202 | 0,938202 | 1 | 1 | 1 |
| 7 | 0,935519 | 0,935519 | 0,935519 | 1 | 1 | 1 |
| 8 | 0,759615 | 0,759615 | 0,759615 | 1 | 1 | 1 |
| 9 | 0,920308 | 0,920308 | 0,920308 | 1 | 1 | 1 |
| 10 | 0,976271 | 0,976271 | 0,976271 | 1 | 1 | 1 |
| 11 | 0,988719 | 0,988719 | 0,988719 | 1 | 1 | 1 |
| 12 | 0,870144 | 0,870144 | 0,870144 | 1 | 1 | 1 |
| 13 | 0,539702 | 0,539702 | 0,539702 | 0 | 0 | 0 |
| 14 | 0,809278 | 0,809278 | 0,809278 | 1 | 1 | 1 |
| 15 | 0,539723 | 0,539723 | 0,539723 | 0 | 0 | 0 |
| 16 | 0,985318 | 0,985318 | 0,985318 | 1 | 1 | 1 |
| 17 | 0,762694 | 0,762694 | 0,762694 | 1 | 1 | 1 |
| 18 | 0,948352 | 0,948352 | 0,948352 | 1 | 1 | 1 |
| 19 | 0,816185 | 0,816185 | 0,816185 | 1 | 1 | 1 |

C, Python e PyRQA: 20/20 OK. O menor DET é 0,539702. O maior é 0,992821.

**Falar**

Este é o slide de encerramento numérico. Ler a tabela como prova de concordância, não como lista de acertos clínicos. Em cada linha, o determinismo do C, o do Python e o da PyRQA são o mesmo número até a sexta casa, e os três alarmes também. O id 8, fora da crise, já está em 0,76, acima do corte. Os ids 13 e 15, dentro da crise, são os únicos abaixo de 0,70. Se a tabela não couber num slide, partir nos ids 0 a 9 e 10 a 19.

---

## Sequência

1. Título e objetivo
2. Definição do determinismo, com o `if` da recorrência
3. Parâmetros, com o efeito de cada um
4. Tipos
5. Memória
6. A RAM acaba no último byte, e o programa ainda assim liga
7. As duas passadas, com a divisão do DET
8. Dados de teste
9. Preparo das janelas
10. Envio à placa, com o `for` de `main`
11. Critério do teste
12. Alarme contra o rótulo
13. Limites
14. Fecho
15. Tabela dos 20 trechos

A conta vem primeiro. Os dados, o preparo e o envio vêm depois dela. A tabela fica por último.
