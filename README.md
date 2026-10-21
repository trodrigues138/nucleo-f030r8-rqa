# RQA em EEG para a NUCLEO-F030R8

Trabalho de disciplina. O programa calcula o determinismo (DET) da Análise de Quantificação de Recorrência sobre uma janela de EEG, com imersão `m = 1` e atraso `tau = 1`. O limiar `DET > 0,70` marca uma mudança de regime neste experimento. Não é critério clínico.

- Autor: Tales Rodrigues
- Data: 27/09/2026 (firmware); 07/10/2026 (versão de arquivo e validação com EEG externo)
- Plataforma alvo: placa NUCLEO-F030R8, microcontrolador STM32F030R8, ARM Cortex-M0, sem FPU, 8 KiB de SRAM, sem heap
- Bancada: PC com gcc (C11) e Python 3

A base teórica está em `explicacao_rqa.md`.

## Permissão de uso

Código acadêmico de Tales Rodrigues, 2026. Pode ser lido, compilado e adaptado em atividades de ensino, com citação do autor. Não autoriza uso clínico nem diagnóstico.

Os registros CHB-MIT (`chb01_03.edf`, `chb01-summary.txt` e a cópia em `dados/`) não fazem parte do código e permanecem sob a licença da base no PhysioNet.

## Origem do código

| Classe | O que é |
|---|---|
| Próprio | Os quatro programas em C, o script `rqa_validacao.py` e `explicacao_rqa.md`, escritos para este trabalho. |
| Legado | `calcular_rqa` em `rqa_portable.c` repete a conta de `rqa_embedded.c`. O ruído sintético usa o gerador linear de Numerical Recipes (multiplicador 1664525, incremento 1013904223). A validação consulta a biblioteca PyRQA (Rawald, Sips e Marwan, 2017), que não está neste repositório. |
| Gerado | Não há código de STM32CubeMX, HAL ou de outro gerador. `trechos_eeg.txt`, `resultados_c.txt` e `tabela_validacao.txt` são saídas da execução, não código-fonte. |

## Hierarquia

```text
nucleo-f030r8-rqa/
├── rqa_embedded.c              firmware de referência (N = 88, float)
├── rqa_embedded_serial.c       mesma conta em float, janela pela USART2
├── rqa_embedded_optimized.c    firmware com matriz em bits e Q8.8 (N = 240)
├── rqa_embedded_optimized_serial.c  mesma conta, janela de 88 amostras pela USART2
├── enviar_serial.py            envia trechos_eeg.txt e imprime a resposta
├── rqa_eeg.c                   mesma conta, com conferência no PC
├── rqa_portable.c              lê janelas de um texto e grava o DET
├── rqa_validacao.py            extrai o EDF, chama o C e compara as contas
├── explicacao_rqa.md           teoria, parâmetros e limites do experimento
├── README.md                   este mapa
├── chb01_03.edf                registro externo (dado, não código)
├── chb01-summary.txt           anotações de crise do registro
├── dados/                      cópia opcional do EDF e do sumário
├── trechos_eeg.txt             gerado: entrada do programa em C
├── resultados_c.txt            gerado: saída do programa em C
└── tabela_validacao.txt        gerado: comparação C, Python e PyRQA
```

Cada firmware é um único arquivo, pensado para entrar num projeto STM32 já existente. Nenhum deles inclui startup, HAL ou linker script: esses arquivos de suporte da placa ficam fora deste repositório. No startup do alvo, `Heap_Size = 0`. A pilha fica com o que sobra da SRAM depois da matriz.

## Função de cada fonte

| Arquivo | Papel | Entrada | Saída |
|---|---|---|---|
| `rqa_embedded.c` | Roda no microcontrolador. Gera um EEG sintético e recalcula o DET em laço. | `simular_anomalia` escrito pelo depurador: 0 repouso sintético, 1 ritmo sintético de 3 Hz. | `determinismo_det` (float, 0 a 1) e `anomalia_transicao_estado` (0 ou 1), lidos pelo depurador. |
| `rqa_embedded_serial.c` | A mesma conta em float, sem gerar sinal. A janela chega pela USART2. | 88 float32 little-endian, precedidos dos bytes `R` e `Q`. | O mesmo quadro de `rqa_embedded_optimized_serial.c`. |
| `rqa_embedded_optimized.c` | Mesma experiência, com a matriz em 1 bit por elemento e a comparação em inteiro Q8.8. A janela passa de 88 para 240 amostras. | O mesmo `simular_anomalia`. | As mesmas duas variáveis de inspeção. |
| `rqa_embedded_optimized_serial.c` | A conta otimizada sem gerar sinal. Cada janela do dataset chega pela USART2 e a resposta volta pela mesma porta. | 88 float32 little-endian, precedidos dos bytes `R` e `Q`. | `R` `Q`, DET float32, numerador e denominador uint32, alarme uint8. |
| `rqa_eeg.c` | Variante usada para conferir a conta. Com `-DRQA_LAPTOP`, imprime o DET e um exemplo 5×5. Sem essa macro, é o laço do microcontrolador. | No PC, nenhuma. No alvo, `simular_anomalia`. | No PC, texto no terminal e código de retorno 0 ou 1. No alvo, as variáveis de inspeção. |
| `rqa_portable.c` | Aplica o DET do `rqa_embedded.c` a janelas vindas de arquivo. Não gera sinal. | Texto com a quantidade de trechos, o cabeçalho de cada um e 88 amostras float. | Texto com o cabeçalho repetido, o numerador, o denominador, o DET e o alarme 0/1. |
| `rqa_validacao.py` | Lê um EDF e o sumário de crises, escolhe 10 janelas dentro e 10 fora, grava `trechos_eeg.txt`, compila e executa `rqa_portable.c` e compara com a mesma conta em Python e com a PyRQA. | `.edf` de 16 bits e o sumário CHB (`File Name`, `File Start Time`, `Seizure Start/End Time`). O canal usado é o 0. | `trechos_eeg.txt`, `resultados_c.txt`, `tabela_validacao.txt` e o mesmo texto no terminal. |

`explicacao_rqa.md` não é executado. Ele registra a definição do DET, a escolha de `N`, `EPSILON` e `L_MIN`, e o que o limiar não autoriza concluir.

## Como usar

Versão de arquivo, no PC:

```text
gcc -std=c11 -O2 -Wall -Wextra rqa_portable.c -o rqa_portable.exe -lm
rqa_portable.exe trechos_eeg.txt resultados_c.txt
```

Sem argumentos, `rqa_portable.exe` usa esses dois nomes no diretório corrente.

Validação com o EEG externo:

```text
py -3 rqa_validacao.py
py -3 rqa_validacao.py dados\chb01_03.edf dados\chb01-summary.txt
```

O script procura o `gcc`, compila `rqa_portable.c` e exige NumPy. A coluna PyRQA só aparece se o pacote `pyrqa` estiver instalado.

Conferência da variante com exemplo 5×5:

```text
gcc -std=c11 -DRQA_LAPTOP -O2 rqa_eeg.c -o rqa_eeg.exe -lm
rqa_eeg.exe
```

No microcontrolador, compile `rqa_embedded.c` ou `rqa_embedded_optimized.c` dentro do projeto da placa, com heap zero. A inspeção é pelo depurador (ST-Link / Live Watch): escreva `simular_anomalia` e leia `determinismo_det` e `anomalia_transicao_estado`. Esses dois firmwares não usam `printf`.

A versão serial não gera sinal. Grave `optimized_serial` e envie o dataset pela COM do ST-Link:

```text
pio run -e optimized_serial -t upload
py -3 -m pip install pyserial
py -3 enviar_serial.py COM8 trechos_eeg.txt
```

## Entrada e saída do algoritmo

O núcleo, em todas as versões, recebe uma janela de amplitudes e devolve:

| Dado | Tipo | Significado |
|---|---|---|
| Janela | `float` ou, na versão otimizada, `int16_t` Q8.8 | `N` amostras. 88 no firmware de referência e no programa de arquivo; 240 na versão otimizada. |
| `determinismo_det` | `float` em \[0, 1\] | Fração dos pontos recorrentes do triângulo superior que forma diagonal de comprimento ≥ 2. |
| `anomalia_transicao_estado` | `uint8_t`, 0 ou 1 | 1 quando o DET passa de 0,70. |

Na versão de arquivo, cada linha de `trechos_eeg.txt` ainda carrega `id`, `inicio_amostra`, `inicio_s` e `rotulo`. São inteiro, inteiro longo, número em segundos e o gabarito 0/1. O gabarito é copiado para `resultados_c.txt` e não entra na conta. A linha de resultado acrescenta `pontos_em_linhas` (numerador) e `recorrencias` (denominador).
