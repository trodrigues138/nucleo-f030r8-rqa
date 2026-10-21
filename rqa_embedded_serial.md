# RQA em float com janela serial

Documento de `rqa_embedded_serial.c`. A placa recebe uma janela de 88 amostras, calcula o determinismo (DET) da Análise de Quantificação de Recorrência e devolve o resultado pela mesma porta serial.

O programa mede recorrência de amplitude, com dimensão de imersão `m = 1` e atraso `τ = 1`. O raio `ε = 8` e o corte `DET > 0,70` são empíricos deste experimento. O alarme não é diagnóstico clínico.

Este arquivo não gera o EEG. A janela já vem pronta em `trechos_eeg.txt`, e o script `enviar_serial.py` a envia. A versão que sintetiza o sinal é `rqa_embedded.c`.

## Parâmetros

| Parâmetro | Valor | Papel |
|---|---:|---|
| `N` | 88 | Amostras de uma janela, cerca de 344 ms a 256 Hz |
| `m`, `τ` | 1, 1 | Cada estado é uma amplitude |
| `EPSILON` | 8,0 | Raio da vizinhança. `R(i,j) = 1` se `\|u(i) − u(j)\| ≤ 8` |
| `L_MIN` | 2 | Comprimento mínimo de uma diagonal contada no numerador |
| `LIMIAR` | 0,70 | O alarme vale 1 somente se o DET passar deste valor |
| Baud | 115200 8N1 | USART2, clock HSI de 8 MHz após o reset |
| Sincronismo | `R` `Q` | Bytes `0x52` e `0x51` |

`i` e `j` são dois instantes da mesma janela, não dois canais.

## Memória

Tudo é estático. Não há `malloc`, recursão nem `printf`. A struct `rqa` concentra a janela, a matriz e os quatro resultados:

| Campo | Tamanho | Uso |
|---|---:|---|
| `determinismo_det` | 4 bytes | DET, de 0 a 1 |
| `pontos_em_linhas` | 4 bytes | Numerador: pontos em diagonais com comprimento ≥ 2 |
| `recorrencias` | 4 bytes | Denominador: pares recorrentes com `j > i` |
| `anomalia_transicao_estado` | 1 byte | 1 se `DET > 0,70` |
| `eeg_sinal[88]` | 352 bytes | Janela recebida |
| `matriz_recorrencia[88][88]` | 7744 bytes | Um byte por par, 0 ou 1 |

A matriz e a janela somam 8096 bytes. A pilha reservada no projeto da IDE é `0x30`. A HAL não entra: o BSS dela não cabe junto com a matriz.

Os quatro resultados são `volatile` para o depurador lê-los com o programa rodando. A resposta serial usa os mesmos campos.

## Protocolo

O PC envia, em little-endian:

```text
'R'  'Q'  u(0)  u(1)  ...  u(87)
```

Cada `u(i)` é um `float32`. O rótulo da janela não viaja no quadro; ele fica no PC e só aparece na impressão do script.

A placa responde, também em little-endian:

```text
'R'  'Q'  DET  pontos_em_linhas  recorrencias  alarme
      4 B float32    4 B uint32     4 B uint32    1 B
```

São 15 bytes de resposta. O script exige o sincronismo `RQ` de volta e encerra se a placa não devolver o quadro inteiro em 5 segundos.

## Implementação

O arquivo está dividido na ordem em que os dados andam: porta serial, quadro de entrada, conta do DET, quadro de saída e `main`.

### Porta serial

`usart2_iniciar` liga o clock de GPIOA e da USART2, coloca PA2 e PA3 em função alternativa AF1 e programa o divisor de baud. Com o HSI de 8 MHz e 115200 bps, o divisor arredondado é 69. Em seguida habilita transmissor, receptor e a USART.

`ler_byte` e `enviar_byte` esperam a flag do hardware (`RXNE` ou `TXE`) e transferem um byte. Não há interrupção nem timeout: se o PC parar no meio do quadro, a placa permanece nesse laço.

`ler_float` junta quatro bytes na ordem little-endian e os lê como `float32`, sem converter o valor. `enviar_u32` faz o caminho inverso para os contadores.

### Recepção da janela

`receber_janela` guarda o byte anterior. O quadro só abre quando esse byte é `R` e o atual é `Q`. Um `R` isolado não dispara a leitura das amostras. Depois do sincronismo, a função grava as 88 amostras em `rqa.eeg_sinal`, na ordem do tempo.

### Cálculo do DET

`calcular_rqa` faz duas passagens.

A primeira preenche `matriz_recorrencia`. O par `(i, j)` recebe 1 quando a distância absoluta das amplitudes não passa de `EPSILON`. O contador `total` só aumenta quando `j > i`, então a diagonal principal fica fora do denominador. Pela simetria da matriz, o triângulo superior representa a mesma razão que a matriz inteira sem essa diagonal.

A segunda passagem percorre só `j > i` e mede cada diagonal uma vez. Um ponto recorrente inicia uma linha quando não existe antecessor, ou quando o antecessor `(i−1, j−1)` vale 0. A partir daí o laço em `k` conta os 1 consecutivos na direção `(i+k, j+k)`. Se o comprimento é pelo menos `L_MIN`, ele entra em `soma`. Recorrências isoladas ficam só no denominador.

Com `total > 0`:

```text
DET = soma / total
```

`pontos_em_linhas` guarda `soma` e `recorrencias` guarda `total`. Se não houver par recorrente, o DET fica 0. O alarme vale 1 somente quando o DET passa de 0,70.

### Resposta

`enviar_resultado` manda `R`, `Q`, os quatro bytes do DET, os dois contadores e o alarme. Espera a flag `TC` para o último bit sair do pino antes da próxima janela.

## Fluxo de execução

```text
reset
  SystemInit          clock permanece no HSI de 8 MHz
  main
    usart2_iniciar    PA2 = TX, PA3 = RX, 115200 8N1
    repetir
      receber_janela  espera RQ e lê 88 float32
      calcular_rqa    matriz, DET e alarme
      enviar_resultado  devolve RQ e os 13 bytes seguintes
```

`main` não retorna. Cada volta substitui a janela e a matriz anteriores. Os quatro resultados permanecem estáveis até a volta seguinte, tanto na serial quanto para o depurador.

## Teste

A configuração da STM32CubeIDE é **Embedded_Serial**. Ela compila este arquivo, não o gerador sintético.

1. Selecione **Embedded_Serial** no menu ao lado do martelo e compile.
2. Grave com a seta verde **Run** e deixe o programa seguir. Se a execução parar em `main`, use **Resume**. Com o núcleo parado, a placa não responde.
3. Envie as janelas por um destes caminhos, com o monitor serial da IDE fechado:
   - clique direito em `STM32CubeIDE/Enviar_janelas.launch` e escolha **Run As**;
   - ou, na pasta do repositório, `py -3 enviar_serial.py`.

O script localiza a porta do ST-Link, lê `trechos_eeg.txt` e imprime uma linha por janela, neste formato:

```text
id <id> rotulo <0 ou 1> DET <0 a 1> pontos <numerador> recorrencias <denominador> alarme <0 ou 1>
```

`rotulo` vem do arquivo no PC. `DET`, `pontos`, `recorrencias` e `alarme` vêm da placa e mudam a cada janela. O alarme 1 significa apenas que, nesta janela e com estes parâmetros, o DET passou de 0,70.

No depurador, depois do **Resume**, as mesmas grandezas estão em `rqa.determinismo_det`, `rqa.pontos_em_linhas`, `rqa.recorrencias` e `rqa.anomalia_transicao_estado`.

Pelo PlatformIO, a gravação correspondente é `pio run -e embedded_serial -t upload`, seguida do mesmo `py -3 enviar_serial.py`.
