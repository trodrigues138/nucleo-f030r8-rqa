# RQA em float com janela serial

Documento de `rqa_embedded_serial.c`. A placa recebe uma janela de 88 amostras pela USART2, calcula o determinismo (DET) da Análise de Quantificação de Recorrência e devolve quatro números pela mesma porta.

O programa mede recorrência de amplitude, com dimensão de imersão `m = 1` e atraso `τ = 1`. Cada estado é uma única amostra `u(i)`. O raio `ε = 8` e o corte `DET > 0,70` são empíricos deste experimento. O alarme registra essa comparação; ele não é diagnóstico clínico.

Este arquivo não gera o EEG. A janela já vem pronta em `trechos_eeg.txt`, e o script `enviar_serial.py` a envia. A versão que sintetiza o sinal é `rqa_embedded.c`. A conta do DET é a mesma das duas.

## O que o arquivo contém

O código está na ordem em que os dados andam:

| Trecho | Função |
|---|---|
| Constantes e a struct `rqa` | Janela, matriz e os quatro resultados, todos estáticos |
| `usart2_iniciar`, `ler_byte`, `enviar_byte`, `enviar_u32`, `ler_float` | Um byte de cada vez na USART2 |
| `receber_janela` | Acha `R` `Q` e grava as 88 amostras |
| `calcular_rqa` | Monta a matriz e o DET |
| `enviar_resultado` | Devolve o quadro de 15 bytes |
| `main` | Repete receber, calcular e responder |

Não há `malloc`, recursão nem `printf`. `math.h` entra só por `fabsf`.

## Constantes

```c
#define N        88
#define EPSILON  8.0f
#define L_MIN    2
#define LIMIAR   0.70f

#define BAUD     115200u
#define SYNC_0   0x52u
#define SYNC_1   0x51u
```

`N` é o tamanho da janela. A 256 Hz, 88 amostras cobrem cerca de 344 ms. `EPSILON` é o raio da vizinhança: dois instantes são recorrentes quando a diferença absoluta das amplitudes não passa de 8. `L_MIN` é o menor comprimento de diagonal que entra no numerador do DET. Uma recorrência isolada, de comprimento 1, fica só no denominador. `LIMIAR` é o corte do alarme.

`BAUD` é a velocidade da USART2. O sufixo `u` marca o literal como inteiro sem sinal. `SYNC_0` e `SYNC_1` são os caracteres `'R'` e `'Q'`. Os dois bytes juntos abrem o quadro de entrada e também o de saída. O rótulo da janela, que no arquivo de texto diz se aquele trecho foi marcado como crise, não viaja neste quadro.

`i` e `j`, usados mais adiante, são dois instantes da mesma janela. O programa compara o sinal com ele mesmo.

## Onde os dados ficam

```c
static struct {
    volatile float determinismo_det;
    volatile uint32_t pontos_em_linhas;
    volatile uint32_t recorrencias;
    volatile uint8_t anomalia_transicao_estado;
    float eeg_sinal[N];
    uint8_t matriz_recorrencia[N][N];
} rqa;
```

`static` limita o nome `rqa` a este arquivo e coloca o bloco em memória estática. Nada é criado em tempo de execução.

| Campo | Tamanho | Conteúdo |
|---|---:|---|
| `determinismo_det` | 4 | DET, de 0 a 1 |
| `pontos_em_linhas` | 4 | Numerador: pontos que pertencem a diagonais de comprimento ≥ 2 |
| `recorrencias` | 4 | Denominador: pares recorrentes com `j > i` |
| `anomalia_transicao_estado` | 1 | 1 se o DET passa de 0,70 |
| `eeg_sinal[88]` | 352 | A janela recebida, um `float` por amostra |
| `matriz_recorrencia[88][88]` | 7744 | Um byte por par: 1 se os instantes são recorrentes, 0 caso contrário |

A janela e a matriz somam 8096 bytes. Os quatro resultados vêm antes delas de propósito: a pilha cresce do fim da SRAM para baixo e encosta na matriz. Deixar o DET no começo do bloco afasta esses quatro bytes desse contato.

`volatile` impede o compilador de guardar o resultado só num registrador. O depurador, e a função que transmite a resposta, leem o valor que está na memória depois de cada janela.

Três `_Static_assert` conferem, na compilação, hipóteses de que o resto do arquivo depende: `float` tem 4 bytes, `N` continua 88 e os dois buffers continuam somando 8096. Se alguém alterar um desses números, a compilação para.

A pilha reservada no projeto da IDE, para esta configuração, é `0x30` bytes. A HAL não entra: o BSS dela não cabe junto com a matriz de 7744 bytes.

## Porta serial

A NUCLEO-F030R8 liga a USART2 do microcontrolador ao conversor USB do ST-Link. No PC essa ligação aparece como uma porta COM. PA2 transmite e PA3 recebe. Depois do reset, o CMSIS deste projeto deixa o relógio no HSI de 8 MHz. `SystemInit` não liga o PLL.

### Ligar o periférico

```c
static void usart2_iniciar(void)
{
    uint32_t divisor;

    SystemCoreClockUpdate();
    RCC->AHBENR |= RCC_AHBENR_GPIOAEN;
    RCC->APB1ENR |= RCC_APB1ENR_USART2EN;
    (void)RCC->APB1ENR;
```

`SystemCoreClockUpdate` lê os registradores de relógio e atualiza `SystemCoreClock`. Com o HSI ainda selecionado, esse valor é 8 000 000. Os dois `|=` ligam o clock de GPIOA e o clock da USART2 sem apagar os outros bits do mesmo registrador. A leitura seguinte de `APB1ENR` espera o enable atravessar o barramento antes de usar a USART.

```c
    GPIOA->MODER = (GPIOA->MODER & ~(3u << 4) & ~(3u << 6))
                 | (2u << 4)
                 | (2u << 6);
    GPIOA->AFR[0] = (GPIOA->AFR[0] & ~(0xFu << 8) & ~(0xFu << 12))
                  | (1u << 8)
                  | (1u << 12);
```

Cada pino tem dois bits em `MODER`. O valor 2 escolhe função alternativa. PA2 ocupa os bits 5:4 (`3u << 4` zera esses dois bits) e PA3 ocupa os bits 7:6. `AFR[0]` guarda quatro bits por pino, do PA0 ao PA7. A USART2 neste microcontrolador é a função alternativa 1: o nibble de PA2 são os bits 11:8, e o de PA3 são os bits 15:12. A máscara zera só esses nibbles e escreve 1 no lugar. Os outros pinos de GPIOA permanecem como estavam.

```c
    divisor = (SystemCoreClock + (BAUD / 2u)) / BAUD;
    if (divisor == 0u) {
        divisor = 1u;
    }
    USART2->CR1 = 0u;
    USART2->BRR = divisor;
    USART2->CR1 = USART_CR1_TE | USART_CR1_RE | USART_CR1_UE;
}
```

O divisor inteiro de baud é `clock / 115200`. Somar metade do baud antes da divisão arredonda para o inteiro mais próximo. Com 8 MHz o quociente exato é cerca de 69,44, então `BRR` recebe 69. O teste `divisor == 0` evita programar um divisor nulo se o clock lido fosse menor que o baud. `CR1` é zerado antes de escrever `BRR`, e só então o transmissor (`TE`), o receptor (`RE`) e a USART (`UE`) são ligados. O formato fica 8 bits, sem paridade e com um stop, que é o padrão do registrador.

### Um byte de cada vez

```c
static uint8_t ler_byte(void)
{
    while ((USART2->ISR & USART_ISR_RXNE) == 0u) {
    }
    return (uint8_t)USART2->RDR;
}

static void enviar_byte(uint8_t valor)
{
    while ((USART2->ISR & USART_ISR_TXE) == 0u) {
    }
    USART2->TDR = valor;
}
```

`RXNE` vale 1 quando há um byte novo em `RDR`. `TXE` vale 1 quando `TDR` pode receber outro byte. Os dois `while` são espera ocupada: não há interrupção e não há prazo. Se o PC parar no meio das 88 amostras, a placa permanece em `ler_byte` até o byte chegar. Ler `RDR` baixa `RXNE`. Escrever `TDR` baixa `TXE`.

```c
static void enviar_u32(uint32_t valor)
{
    enviar_byte((uint8_t)(valor & 0xFFu));
    enviar_byte((uint8_t)((valor >> 8) & 0xFFu));
    enviar_byte((uint8_t)((valor >> 16) & 0xFFu));
    enviar_byte((uint8_t)((valor >> 24) & 0xFFu));
}
```

Um `uint32` sai em quatro bytes, o menos significativo primeiro. Essa é a ordem little-endian que `struct.unpack` usa no script. `pontos_em_linhas` e `recorrencias` passam por esta função. O DET não passa, porque ele é `float`: os quatro bytes saem na ordem em que já estão na memória.

```c
static float ler_float(void)
{
    union {
        uint32_t u;
        float f;
    } conv;

    conv.u = (uint32_t)ler_byte()
           | ((uint32_t)ler_byte() << 8)
           | ((uint32_t)ler_byte() << 16)
           | ((uint32_t)ler_byte() << 24);
    return conv.f;
}
```

A `union` ocupa quatro bytes e pode ser lida como inteiro ou como `float`. O programa monta o inteiro com os quatro bytes recebidos e devolve o mesmo padrão de bits visto como `float32`. Não há conversão de texto para número. O primeiro byte lido é o menos significativo, de acordo com o `struct.pack("<f", ...)` do script.

## Recepção da janela

O PC envia, para cada trecho:

```text
'R'  'Q'  u(0)  u(1)  ...  u(87)
```

Cada `u(i)` são quatro bytes. O quadro completo tem 2 + 88 × 4 = 354 bytes.

```c
static void receber_janela(void)
{
    uint8_t anterior = 0u;
    uint_fast16_t i;

    for (;;) {
        uint8_t atual = ler_byte();

        if (anterior == SYNC_0 && atual == SYNC_1) {
            break;
        }
        anterior = atual;
    }

    for (i = 0u; i < N; i++) {
        rqa.eeg_sinal[i] = ler_float();
    }
}
```

O primeiro laço procura a dupla, não um byte isolado. `anterior` começa em 0, que não é `'R'`. A cada byte novo, o programa olha o par `(anterior, atual)`. Só `'R'` seguido imediatamente de `'Q'` interrompe a espera. Um `'R'` sozinho fica guardado em `anterior` e só abre o quadro se o próximo byte for `'Q'`. Bytes que sobraram de uma transmissão anterior são descartados por esse mesmo laço.

Depois do `break`, o segundo laço chama `ler_float` 88 vezes e grava `rqa.eeg_sinal[0]` até `rqa.eeg_sinal[87]`, na ordem do tempo. A função não volta enquanto essas 352 bytes de amostra não tiverem chegado.

## Cálculo do DET

A matriz de recorrência, com `m = 1`, é

```text
R(i,j) = 1, se |u(i) − u(j)| ≤ EPSILON
R(i,j) = 0, caso contrário
```

`u(i)` é `rqa.eeg_sinal[i]`. A matriz é simétrica e a diagonal principal, a linha de identidade, é sempre 1, porque a distância de uma amostra até ela mesma é zero. O DET usa só o triângulo `j > i`.

```text
DET = (soma dos comprimentos das diagonais com L ≥ 2) / (número de pares recorrentes com j > i)
```

```c
static void calcular_rqa(void)
{
    uint_fast16_t i;
    uint_fast16_t j;
    uint_fast16_t k;
    uint_fast16_t comprimento;
    uint32_t total = 0u;
    uint32_t soma = 0u;
```

`total` será o denominador. `soma` será o numerador. Os dois começam em zero a cada janela, então o resultado anterior não vaza para a seguinte.

### Primeira passagem: preencher a matriz

```c
    for (i = 0u; i < N; i++) {
        for (j = 0u; j < N; j++) {
            if (fabsf(rqa.eeg_sinal[i] - rqa.eeg_sinal[j]) <= EPSILON) {
                rqa.matriz_recorrencia[i][j] = 1u;
                if (j > i) {
                    total++;
                }
            } else {
                rqa.matriz_recorrencia[i][j] = 0u;
            }
        }
    }
```

Para cada par de instantes o programa calcula o módulo da diferença com `fabsf` e escreve 0 ou 1 na matriz. O `if (j > i)` deixa de fora a diagonal principal e também o triângulo inferior. Como `R(i,j) = R(j,i)`, contar um lado basta para a razão. `total` é, portanto, o número de recorrências estritamente acima da diagonal.

### Segunda passagem: medir cada diagonal uma vez

Uma diagonal paralela à principal é uma sequência `R(i,j)`, `R(i+1,j+1)`, `R(i+2,j+2)`, … em que os dois instantes avançam juntos. Se essa sequência permanece em 1, dois trechos da janela continuam próximos enquanto o tempo anda. O programa precisa somar o comprimento dessa sequência uma única vez. Se somasse em todo ponto da linha, o numerador contaria a mesma diagonal várias vezes e o DET poderia passar de 1.

```c
    for (i = 0u; i < N; i++) {
        for (j = (uint_fast16_t)(i + 1u); j < N; j++) {
            uint_fast8_t inicio = 0u;

            if (rqa.matriz_recorrencia[i][j] == 1u) {
                if (i == 0u || rqa.matriz_recorrencia[i - 1u][j - 1u] == 0u) {
                    inicio = 1u;
                }
            }
            if (inicio == 0u) {
                continue;
            }
```

`j` já nasce em `i + 1`, então esta passagem nem visita a diagonal principal. `inicio` vale 1 só no primeiro 1 de uma diagonal:

- o ponto `(i, j)` é recorrente;
- e não existe antecessor, porque `i == 0`, ou o antecessor `(i−1, j−1)` é 0.

Se `(i, j)` é 0, ou se ele continua uma linha que já começou na casa anterior, `inicio` permanece 0 e o `continue` pula para o próximo par. O sufixo `u` de `1u` e `0u` é o mesmo marcador de inteiro sem sinal usado nas constantes.

```c
            comprimento = 0u;
            for (k = 0u; (i + k < N) && (j + k < N); k++) {
                if (rqa.matriz_recorrencia[i + k][j + k] == 1u) {
                    comprimento++;
                } else {
                    break;
                }
            }
            if (comprimento >= L_MIN) {
                soma += comprimento;
            }
        }
    }
```

A partir do início achado, `k` caminha pela diagonal enquanto os dois índices couberem na janela e enquanto a matriz continuar em 1. O primeiro 0 encerra a linha. `comprimento` é o número de 1 consecutivos. Como a caminhada começou num 1, esse número é pelo menos 1. A soma só recebe a linha quando o comprimento é pelo menos `L_MIN`, isto é, 2. Uma recorrência isolada aumenta `total` na primeira passagem e não aumenta `soma`.

Um exemplo curto, com uma diagonal de três 1 e depois um 0:

```text
(1,4) = 1    antecessor (0,3) = 0    inicio = 1
(2,5) = 1    antecessor (1,4) = 1    inicio = 0, não mede de novo
(3,6) = 1    antecessor (2,5) = 1    inicio = 0
(4,7) = 0
```

A medição ocorre só em `(1,4)`. O laço de `k` conta 3 e para em `(4,7)`. Como 3 ≥ 2, `soma` recebe 3. Os pontos `(2,5)` e `(3,6)` já entraram nesse 3; contá-los outra vez dobraria o numerador.

### Fechar a razão e o alarme

```c
    rqa.pontos_em_linhas = soma;
    rqa.recorrencias = total;
    if (total == 0u) {
        rqa.determinismo_det = 0.0f;
    } else {
        rqa.determinismo_det = (float)soma / (float)total;
    }
    rqa.anomalia_transicao_estado = (rqa.determinismo_det > LIMIAR) ? 1u : 0u;
}
```

Os dois contadores são gravados antes da divisão, para a resposta serial levar os inteiros exatos e não só a razão. A divisão só acontece com `total > 0`. O cast para `float` faz a divisão em ponto flutuante; uma divisão inteira de `soma / total` perderia a fração. O DET resultante fica entre 0 e 1 porque `soma` reúne um subconjunto dos pontos já contados em `total`.

O alarme é a comparação com 0,70. O operador `>` é estrito: um DET exatamente igual a 0,70 produz alarme 0. O nome do campo registra a passagem desse limiar nesta janela. Ele não classifica o paciente.

## Resposta

```c
static void enviar_resultado(void)
{
    const uint8_t *det = (const uint8_t *)&rqa.determinismo_det;
    uint_fast8_t i;

    enviar_byte(SYNC_0);
    enviar_byte(SYNC_1);
    for (i = 0u; i < 4u; i++) {
        enviar_byte(det[i]);
    }
    enviar_u32(rqa.pontos_em_linhas);
    enviar_u32(rqa.recorrencias);
    enviar_byte(rqa.anomalia_transicao_estado);
    while ((USART2->ISR & USART_ISR_TC) == 0u) {
    }
}
```

O ponteiro `det` olha os quatro bytes já armazenados no `float`, do endereço mais baixo para o mais alto. No Cortex-M0 essa ordem é little-endian, a mesma em que `ler_float` montou as amostras. Em seguida saem os dois `uint32` e o byte do alarme.

O quadro de volta tem 2 + 4 + 4 + 4 + 1 = 15 bytes:

```text
'R'  'Q'  DET  pontos_em_linhas  recorrencias  alarme
      4 bytes float32    4 bytes      4 bytes     1 byte
```

`TXE`, usado dentro de `enviar_byte`, indica que o byte entrou no registrador de transmissão. `TC` indica que esse byte já saiu pelo pino, inclusive o stop bit. O `while` final espera `TC` para a função só retornar quando a resposta inteira tiver deixado a placa. A volta seguinte de `main` pode então bloquear em `ler_byte` sem cortar a transmissão anterior.

## Fluxo de execução

```c
int main(void)
{
    usart2_iniciar();
    for (;;) {
        receber_janela();
        calcular_rqa();
        enviar_resultado();
    }
}
```

Depois do reset, o startup copia `.data`, zera `.bss` (a struct `rqa` nasce zerada), chama `SystemInit` e entra em `main`.

1. `usart2_iniciar` prepara PA2, PA3 e o divisor 69. Ainda não há janela.
2. `receber_janela` bloqueia até o PC mandar `RQ` e as 88 amostras. `rqa.eeg_sinal` passa a ser essa janela.
3. `calcular_rqa` reescreve a matriz inteira e os quatro resultados. A matriz da janela anterior deixa de valer.
4. `enviar_resultado` transmite os 15 bytes e espera o último bit sair.
5. O `for (;;)` volta ao passo 2. `main` não tem `return`: o laço é o modo normal de operação.

Entre uma resposta e a janela seguinte, os quatro campos `volatile` permanecem com o último DET. É esse intervalo que o depurador observa se a execução for pausada depois do envio.

## Teste

A configuração **Embedded_Serial** da STM32CubeIDE compila este arquivo. **Embedded** compila o gerador sintético, que não abre a USART2.

1. Selecione **Embedded_Serial** no menu ao lado do martelo e compile.
2. Grave com a seta verde **Run** e deixe o programa seguir. Se a execução parar em `main`, use **Resume** antes de enviar dados. Com o núcleo parado em `ler_byte`, a placa não consome a porta.
3. Feche o monitor serial da IDE. A COM do ST-Link só aceita um programa por vez.
4. Envie as janelas pelo clique direito em `STM32CubeIDE/Enviar_janelas.launch` e **Run As**, ou, na pasta do repositório, por `py -3 enviar_serial.py`.

O script lê `trechos_eeg.txt`, descobre a porta cujo VID:PID é `0483:374B` e, para cada janela, escreve os mesmos 354 bytes que `receber_janela` espera. Ele então lê exatamente 15 bytes, confere o `RQ` de volta e imprime:

```text
id <id> rotulo <0 ou 1> DET <0 a 1> pontos <numerador> recorrencias <denominador> alarme <0 ou 1>
```

`id` e `rotulo` vêm do arquivo no PC. `DET`, `pontos`, `recorrencias` e `alarme` são `determinismo_det`, `pontos_em_linhas`, `recorrencias` e `anomalia_transicao_estado` desempacotados em little-endian. O valor muda a cada janela. Alarme 1 significa apenas que, com `ε = 8`, `L_MIN = 2` e esta janela, o DET passou de 0,70.

Se a placa não devolver 15 bytes em 5 segundos, o script encerra dizendo quantos bytes chegaram. As causas usuais são a configuração **Embedded** gravada no lugar da **Embedded_Serial**, o núcleo ainda parado no depurador, ou o monitor serial ocupando a COM.

No depurador, com o programa em **Resume**, as mesmas grandezas podem ser lidas em `rqa` depois que uma janela foi respondida.

Pelo PlatformIO, a gravação desta versão é `pio run -e embedded_serial -t upload`, seguida do mesmo `py -3 enviar_serial.py`.
