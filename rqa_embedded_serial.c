/*
 * =============================================================================
 * RQA em float com janela vinda da serial — sem geracao de sinal
 * =============================================================================
 *
 * Apresentacao
 *   Segunda versao de rqa_embedded.c. A conta e a mesma: matriz de
 *   recorrencia em um byte por par e distancia com fabsf. Este arquivo
 *   nao gera senoide nem ruido. Cada janela chega pela USART2, ja
 *   extraida do dataset (trechos_eeg.txt, 88 amostras float).
 *
 *   A recorrencia usa m = 1 e tau = 1: cada estado e uma unica
 *   amplitude, sem reconstrucao multidimensional. O limiar EPSILON e
 *   o corte DET > 0,70 sao empiricos deste experimento. O alarme nao
 *   e diagnostico clinico.
 *
 * Autor
 *   Tales Rodrigues
 *
 * Data
 *   14/10/2026
 *
 * Contexto
 *   Trabalho de disciplina. A NUCLEO-F030R8 executa o RQA em float
 *   sobre as janelas do dataset, em vez do EEG sintetico.
 *
 * Plataforma alvo
 *   NUCLEO-F030R8, STM32F030R8, ARM Cortex-M0, sem FPU, 8 KiB de SRAM.
 *   USART2 em PA2 (TX) e PA3 (RX), no COM virtual do ST-Link.
 *   Relogio apos o reset: HSI de 8 MHz. Baud 115200, 8N1.
 *   A matriz ocupa 7744 bytes. A pilha reservada fica em 0x30.
 *
 * Permissao de uso (copyright)
 *   Codigo academico de Tales Rodrigues, 2026. Pode ser lido, compilado
 *   e adaptado em atividades de ensino, com citacao do autor.
 *   Nao autoriza uso clinico nem diagnostico.
 *   O EEG CHB-MIT permanece sob a licenca da base no PhysioNet.
 *
 * Origem do codigo
 *   Proprio: USART2 e o sincronismo do quadro.
 *   Legado: calcular_rqa, copiado de rqa_embedded.c. Nao ha seno,
 *   cosseno nem o LCG de Numerical Recipes.
 *   Gerado: nenhum.
 *
 * Como usar
 *   Na STM32CubeIDE, configuracao Embedded_Serial: compilar, Run, e
 *   depois Run As em Enviar_janelas.launch.
 *   No PlatformIO: pio run -e embedded_serial -t upload
 *   e, na pasta do repositorio, py -3 enviar_serial.py
 *
 * Entrada (USART2, little-endian)
 *   Bytes 'R' (0x52) e 'Q' (0x51), depois N = 88 amostras float32.
 *   Sao as mesmas amostras de uma linha de trechos_eeg.txt.
 *   O rotulo da crise nao entra no quadro.
 *
 * Saida (USART2, little-endian)
 *   'R' 'Q'
 *   rqa.determinismo_det          float32, de 0 a 1
 *   rqa.pontos_em_linhas          uint32, numerador do DET
 *   rqa.recorrencias              uint32, denominador do DET
 *   rqa.anomalia_transicao_estado uint8, 1 se o DET passa de 0,70
 *
 * Memoria
 *   rqa.eeg_sinal[88]                 352 bytes
 *   rqa.matriz_recorrencia[88][88]   7744 bytes
 *   Sem malloc, sem recursao e sem printf. math.h entra so por fabsf.
 *
 * Fluxo
 *   usart2_iniciar
 *   repetir: receber_janela, calcular_rqa, enviar_resultado
 * =============================================================================
 */

#include <stdint.h>
#include <math.h>
#include <stm32f0xx.h>

#pragma GCC optimize ("Os")

/* Janela e RQA. m = 1 e tau = 1: cada estado e uma amplitude. */
#define N        88
#define EPSILON  8.0f
#define L_MIN    2
#define LIMIAR   0.70f

/* Quadro da USART2. SYNC_0 e SYNC_1 sao os caracteres 'R' e 'Q'. */
#define BAUD     115200u
#define SYNC_0   0x52u
#define SYNC_1   0x51u

/*
 * A pilha cresce do fim da SRAM para baixo e encosta na matriz.
 * Os resultados ficam no inicio do bloco, longe desse contato.
 * A janela e a matriz ocupam o restante: 352 + 7744 bytes.
 *
 * Os quatro primeiros campos sao volatile porque o depurador os le
 * com o alvo rodando. O envio serial usa os mesmos campos.
 */
static struct {
    volatile float determinismo_det;
    volatile uint32_t pontos_em_linhas;
    volatile uint32_t recorrencias;
    volatile uint8_t anomalia_transicao_estado;
    float eeg_sinal[N];
    uint8_t matriz_recorrencia[N][N];
} rqa;

_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");
_Static_assert(N == 88, "a janela serial acompanha o dataset de 88 amostras");
_Static_assert(sizeof(rqa.eeg_sinal) + sizeof(rqa.matriz_recorrencia) == 8096u,
               "buffers do RQA devem ocupar 8096 bytes");

/* ------------------------------------------------------------------------- */
/* USART2: um byte de cada vez, sem interrupcao e sem printf.                */
/* ------------------------------------------------------------------------- */

/*
 * usart2_iniciar
 *
 * Comportamento
 *   Liga GPIOA e USART2. PA2 transmite e PA3 recebe, funcao alternativa
 *   AF1. O divisor de baud sai de SystemCoreClock. No reset deste CMSIS
 *   o clock e o HSI de 8 MHz, e o divisor de 115200 fica em 69.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Nenhuma do RQA. Escreve RCC, GPIOA e USART2.
 */
static void usart2_iniciar(void)
{
    uint32_t divisor;

    SystemCoreClockUpdate();
    RCC->AHBENR |= RCC_AHBENR_GPIOAEN;
    RCC->APB1ENR |= RCC_APB1ENR_USART2EN;
    (void)RCC->APB1ENR; /* a leitura garante que o enable chegou no barramento */

    /* MODER: dois bits por pino. 2 = funcao alternativa. PA2 usa os bits
     * 5:4 e PA3 usa os bits 7:6. AFR[0]: quatro bits por pino. AF1 da
     * USART2 ocupa os nibbles de PA2 (bits 11:8) e PA3 (bits 15:12). */
    GPIOA->MODER = (GPIOA->MODER & ~(3u << 4) & ~(3u << 6))
                 | (2u << 4)
                 | (2u << 6);
    GPIOA->AFR[0] = (GPIOA->AFR[0] & ~(0xFu << 8) & ~(0xFu << 12))
                  | (1u << 8)
                  | (1u << 12);

    /* Arredonda o divisor para o inteiro mais proximo de clock/baud. */
    divisor = (SystemCoreClock + (BAUD / 2u)) / BAUD;
    if (divisor == 0u) {
        divisor = 1u;
    }
    USART2->CR1 = 0u;
    USART2->BRR = divisor;
    USART2->CR1 = USART_CR1_TE | USART_CR1_RE | USART_CR1_UE;
}

/*
 * ler_byte
 *
 * Comportamento
 *   Espera um byte na USART2. O laco trava ate RXNE: nao ha timeout.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   O byte recebido.
 *
 * Globais afetadas
 *   Nenhuma.
 */
static uint8_t ler_byte(void)
{
    while ((USART2->ISR & USART_ISR_RXNE) == 0u) {
    }
    return (uint8_t)USART2->RDR;
}

/*
 * enviar_byte
 *
 * Comportamento
 *   Espera o registrador de transmissao esvaziar e escreve um byte.
 *
 * Parametros
 *   valor: byte a transmitir.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Nenhuma.
 */
static void enviar_byte(uint8_t valor)
{
    while ((USART2->ISR & USART_ISR_TXE) == 0u) {
    }
    USART2->TDR = valor;
}

/*
 * enviar_u32
 *
 * Comportamento
 *   Transmite um uint32 little-endian, o byte menos significativo primeiro.
 *
 * Parametros
 *   valor: inteiro de 32 bits.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Nenhuma.
 */
static void enviar_u32(uint32_t valor)
{
    enviar_byte((uint8_t)(valor & 0xFFu));
    enviar_byte((uint8_t)((valor >> 8) & 0xFFu));
    enviar_byte((uint8_t)((valor >> 16) & 0xFFu));
    enviar_byte((uint8_t)((valor >> 24) & 0xFFu));
}

/*
 * ler_float
 *
 * Comportamento
 *   Le 4 bytes little-endian e os interpreta como float32, sem conversao
 *   numerica: os bits recebidos sao os bits do float.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   A amostra recebida.
 *
 * Globais afetadas
 *   Nenhuma.
 */
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

/* ------------------------------------------------------------------------- */
/* Quadro de entrada: sincronismo e as 88 amostras.                         */
/* ------------------------------------------------------------------------- */

/*
 * receber_janela
 *
 * Comportamento
 *   Descarta bytes ate ver 'R' seguido de 'Q'. Em seguida grava N
 *   amostras em rqa.eeg_sinal, na ordem do tempo. Um 'R' isolado nao
 *   abre o quadro: ele fica em anterior e so vale se o proximo for 'Q'.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void). Nao retorna enquanto os 88 floats nao chegarem.
 *
 * Globais afetadas
 *   Escreve rqa.eeg_sinal.
 */
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

/* ------------------------------------------------------------------------- */
/* RQA: matriz de recorrencia e determinismo.                               */
/* i e j sao dois instantes da mesma janela, nao dois sinais.               */
/* ------------------------------------------------------------------------- */

/*
 * calcular_rqa
 *
 * Comportamento
 *   R(i,j) = 1 quando |u(i) - u(j)| <= EPSILON. O denominador conta
 *   so o triangulo superior, j > i, e portanto exclui a diagonal
 *   principal. O numerador soma, uma vez, cada diagonal de comprimento
 *   >= L_MIN. DET = soma / total, ou 0 se nao houver par recorrente.
 *   O alarme vale 1 somente quando DET > LIMIAR.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Le rqa.eeg_sinal.
 *   Escreve rqa.matriz_recorrencia, rqa.pontos_em_linhas, rqa.recorrencias,
 *   rqa.determinismo_det e rqa.anomalia_transicao_estado.
 */
static void calcular_rqa(void)
{
    uint_fast16_t i;
    uint_fast16_t j;
    uint_fast16_t k;
    uint_fast16_t comprimento;
    uint32_t total = 0u;
    uint32_t soma = 0u;

    /* Preenche a matriz inteira. total so cresce acima da diagonal. */
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

    /* Cada diagonal paralela a principal e medida a partir do primeiro 1.
     * inicio vale 1 quando (i,j) e recorrente e o antecessor (i-1,j-1)
     * nao e. Assim a linha nao e contada de novo em cada ponto. */
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

    rqa.pontos_em_linhas = soma;
    rqa.recorrencias = total;
    if (total == 0u) {
        rqa.determinismo_det = 0.0f;
    } else {
        rqa.determinismo_det = (float)soma / (float)total;
    }
    rqa.anomalia_transicao_estado = (rqa.determinismo_det > LIMIAR) ? 1u : 0u;
}

/* ------------------------------------------------------------------------- */
/* Quadro de saida: os quatro resultados, na ordem em que o script os le.   */
/* ------------------------------------------------------------------------- */

/*
 * enviar_resultado
 *
 * Comportamento
 *   Devolve 'R' 'Q', o DET em 4 bytes, os dois uint32 e o alarme.
 *   Espera TC para o ultimo bit sair do pino antes da proxima janela.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Le rqa.determinismo_det, rqa.pontos_em_linhas, rqa.recorrencias e
 *   rqa.anomalia_transicao_estado.
 */
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

/*
 * main
 *
 * Comportamento
 *   Inicializa a USART2 e repete tres passos: recebe uma janela,
 *   calcula o DET em float e responde pela mesma porta. Nao ha
 *   geracao de sinal e o laco nao termina.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Declarado como int. O laco nao termina, entao nao devolve valor.
 *
 * Globais afetadas
 *   receber_janela escreve rqa.eeg_sinal. calcular_rqa escreve a matriz
 *   e as globais de resultado.
 */
int main(void)
{
    usart2_iniciar();
    for (;;) {
        receber_janela();
        calcular_rqa();
        enviar_resultado();
    }
}
