/*
 * =============================================================================
 * RQA otimizado com janela vinda da serial — sem geracao de sinal
 * =============================================================================
 *
 * Apresentacao
 *   Segunda versao de rqa_embedded_optimized.c. O calculo e o mesmo:
 *   matriz de recorrencia em 1 bit por elemento e distancia em inteiro
 *   Q8.8. Este arquivo nao gera senoide nem ruido. Cada janela chega
 *   pela USART2, ja extraida do dataset (trechos_eeg.txt, 88 amostras).
 *
 * Autor
 *   Tales Rodrigues
 *
 * Data
 *   14/10/2026
 *
 * Contexto
 *   Trabalho de disciplina. A NUCLEO-F030R8 deixa de simular o EEG e
 *   passa a medir o DET das janelas preparadas no PC.
 *
 * Plataforma alvo
 *   NUCLEO-F030R8, STM32F030R8, ARM Cortex-M0, sem FPU, 8 KiB de SRAM.
 *   USART2 nos pinos PA2 (TX) e PA3 (RX), ligados ao ST-Link. No PC
 *   essa porta aparece como o COM virtual da placa. Relogio apos o
 *   reset: HSI de 8 MHz. Baud 115200, 8 bits, sem paridade, 1 stop.
 *
 * Permissao de uso (copyright)
 *   Codigo academico de Tales Rodrigues, 2026. Pode ser lido, compilado
 *   e adaptado em atividades de ensino, com citacao do autor.
 *   Nao autoriza uso clinico nem diagnostico.
 *   O EEG CHB-MIT permanece sob a licenca da base no PhysioNet.
 *
 * Origem do codigo
 *   Proprio: USART2, sincronismo do quadro e conversao da janela.
 *   Legado: bit_set, bit_get e calcular_rqa, copiados de
 *   rqa_embedded_optimized.c. Nao ha seno, cosseno nem LCG.
 *   Gerado: nenhum. Nao ha HAL nem saida de CubeMX.
 *
 * Como usar
 *   pio run -e optimized_serial -t upload
 *   py -3 enviar_serial.py COM8
 *   Sem argumento, o script usa COM8 e trechos_eeg.txt.
 *
 * Entrada (USART2, little-endian)
 *   Dois bytes de sincronismo: 'R' (0x52) e 'Q' (0x51).
 *   Em seguida, N = 88 amostras float32, as mesmas da linha de
 *   amostras de trechos_eeg.txt. 352 bytes por janela.
 *   O rotulo da crise nao entra no quadro: ele fica no PC.
 *
 * Saida (USART2, little-endian, apos o calculo)
 *   'R' 'Q'
 *   determinismo_det          float32, de 0 a 1
 *   pontos_em_linhas          uint32, numerador do DET
 *   recorrencias              uint32, denominador do DET
 *   anomalia_transicao_estado uint8, 1 se o DET passa de 0,70
 *   O mesmo resultado continua nas globais, para o depurador.
 *
 * Memoria
 *   sinal_fixo[88]              176 bytes
 *   matriz_bits[88][11]         968 bytes
 *   Sem malloc, sem recursao e sem printf.
 * =============================================================================
 */

#include <stdint.h>
#include <stm32f0xx.h>

#pragma GCC optimize ("Os")

#define N               88
#define L_MIN           2

#define ESCALA_FIXO     256        /* float -> Q8.8: valor real = inteiro / 256 */
#define EPSILON_FIXO    2048       /* raio 8.0 na mesma escala do dataset       */
#define LIMIAR_DET_F    0.70f

#define BYTES_POR_LINHA ((N + 7u) >> 3)   /* 11 */

#define BAUD            115200u
#define SYNC_0          0x52u      /* 'R' */
#define SYNC_1          0x51u      /* 'Q' */

/*
 * Janela recebida, em Q8.8. Nao ha segundo buffer: cada float da serial
 * e convertido e guardado direto aqui.
 */
static int16_t sinal_fixo[N];

/* Recurrence plot com 1 bit por par. Bit j da linha i fica em j%8 do byte j/8. */
static uint8_t matriz_bits[N][BYTES_POR_LINHA];

volatile float   determinismo_det;          /* DET da ultima janela, 0 a 1 */
volatile uint8_t anomalia_transicao_estado; /* 1 se DET > 0,70 */
volatile uint32_t pontos_em_linhas;         /* numerador do DET */
volatile uint32_t recorrencias;             /* denominador: pares com j > i */

_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");
_Static_assert(sizeof(int16_t) == 2u, "int16_t deve ter 16 bits");
_Static_assert(N == 88, "a janela serial acompanha o dataset de 88 amostras");
_Static_assert(BYTES_POR_LINHA == 11u, "11 bytes por linha da matriz");
_Static_assert(sizeof(sinal_fixo) == 176u, "vetor EEG Q8.8: 176 bytes");
_Static_assert(sizeof(matriz_bits) == 968u, "matriz R bits: 968 bytes");

/*
 * bit_set
 *
 * Comportamento
 *   Liga o bit R(i,j) sem alterar os outros bits do mesmo byte.
 *   Copiado de rqa_embedded_optimized.c.
 *
 * Parametros
 *   i: linha, de 0 a N-1.
 *   j: coluna, de 0 a N-1.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Escreve um bit de matriz_bits.
 */
static inline void bit_set(uint_fast16_t i, uint_fast16_t j)
{
    matriz_bits[i][j >> 3] |= (uint8_t)(1u << (j & 7u));
}

/*
 * bit_get
 *
 * Comportamento
 *   Le o bit R(i,j).
 *
 * Parametros
 *   i: linha, de 0 a N-1.
 *   j: coluna, de 0 a N-1.
 *
 * Retorno
 *   1 se o par recorre, 0 caso contrario.
 *
 * Globais afetadas
 *   Le matriz_bits. Nao escreve.
 */
static inline uint_fast8_t bit_get(uint_fast16_t i, uint_fast16_t j)
{
    return (uint_fast8_t)((matriz_bits[i][j >> 3] >> (j & 7u)) & 1u);
}

/*
 * usart2_iniciar
 *
 * Comportamento
 *   Liga GPIOA e USART2. PA2 transmite e PA3 recebe, funcao
 *   alternativa AF1. O divisor de baud sai de SystemCoreClock,
 *   que no reset deste CMSIS e o HSI de 8 MHz.
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
    (void)RCC->APB1ENR;

    GPIOA->MODER = (GPIOA->MODER & ~(3u << 4) & ~(3u << 6))
                 | (2u << 4)
                 | (2u << 6);
    GPIOA->AFR[0] = (GPIOA->AFR[0] & ~(0xFu << 8) & ~(0xFu << 12))
                  | (1u << 8)
                  | (1u << 12);

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
 *   Espera um byte na USART2. O polling nao usa interrupcao.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   O byte recebido.
 *
 * Globais afetadas
 *   Nenhuma. Le USART2->RDR.
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
 *   Nenhuma. Escreve USART2->TDR.
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
 *   Transmite um uint32 em little-endian, o mesmo byte order do PC.
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
 *   Le 4 bytes little-endian e os interpreta como float32.
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

/*
 * para_q8_8
 *
 * Comportamento
 *   Multiplica a amplitude por 256 e satura em int16. E a mesma
 *   escala de rqa_embedded_optimized.c, aplicada a cada amostra
 *   que chega, sem guardar o float.
 *
 * Parametros
 *   amostra: amplitude do dataset, em float32.
 *
 * Retorno
 *   Valor Q8.8 saturado.
 *
 * Globais afetadas
 *   Nenhuma.
 */
static int16_t para_q8_8(float amostra)
{
    float escalada = amostra * (float)ESCALA_FIXO;

    if (escalada >= 32767.0f) {
        return 32767;
    }
    if (escalada <= -32768.0f) {
        return (int16_t)(-32768);
    }
    return (int16_t)escalada;
}

/*
 * receber_janela
 *
 * Comportamento
 *   Procura os bytes 'R' e 'Q' e, a partir daí, le N amostras
 *   float32 para sinal_fixo. Um 'R' perdido no meio dos dados
 *   so reabre o quadro se o byte seguinte for 'Q'.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Escreve sinal_fixo.
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
        sinal_fixo[i] = para_q8_8(ler_float());
    }
}

/*
 * calcular_rqa
 *
 * Comportamento
 *   Mesma conta de rqa_embedded_optimized.c. A distancia e inteira.
 *   O triangulo superior (j > i) forma o denominador. Diagonais com
 *   comprimento >= L_MIN formam o numerador. DET = soma / total,
 *   ou 0 se nao houver par recorrente.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Le sinal_fixo.
 *   Escreve matriz_bits, pontos_em_linhas, recorrencias,
 *   determinismo_det e anomalia_transicao_estado.
 */
static void calcular_rqa(void)
{
    uint_fast16_t i;
    uint_fast16_t j;
    uint_fast16_t k;
    uint_fast16_t comprimento;
    uint_fast16_t byte_idx;
    uint32_t total = 0u;
    uint32_t soma = 0u;

    for (i = 0u; i < N; i++) {
        int32_t xi = (int32_t)sinal_fixo[i];

        for (byte_idx = 0u; byte_idx < BYTES_POR_LINHA; byte_idx++) {
            matriz_bits[i][byte_idx] = 0u;
        }

        for (j = 0u; j < N; j++) {
            int32_t diff = xi - (int32_t)sinal_fixo[j];
            int32_t mask = diff >> 31;

            diff = (diff ^ mask) - mask;
            if (diff <= EPSILON_FIXO) {
                bit_set(i, j);
                if (j > i) {
                    total++;
                }
            }
        }
    }

    for (i = 0u; i < N; i++) {
        for (j = (uint_fast16_t)(i + 1u); j < N; j++) {
            if (bit_get(i, j) == 0u) {
                continue;
            }
            if (i > 0u && bit_get((uint_fast16_t)(i - 1u),
                                  (uint_fast16_t)(j - 1u)) != 0u) {
                continue;
            }

            comprimento = 0u;
            for (k = 0u; (i + k < N) && (j + k < N); k++) {
                if (bit_get((uint_fast16_t)(i + k),
                            (uint_fast16_t)(j + k)) != 0u) {
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

    pontos_em_linhas = soma;
    recorrencias = total;
    if (total == 0u) {
        determinismo_det = 0.0f;
    } else {
        determinismo_det = (float)soma / (float)total;
    }
    anomalia_transicao_estado = (determinismo_det > LIMIAR_DET_F) ? 1u : 0u;
}

/*
 * enviar_resultado
 *
 * Comportamento
 *   Devolve o quadro de resposta e espera o ultimo bit sair do pino.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Le determinismo_det, pontos_em_linhas, recorrencias e
 *   anomalia_transicao_estado.
 */
static void enviar_resultado(void)
{
    union {
        float f;
        uint32_t u;
    } conv;

    enviar_byte(SYNC_0);
    enviar_byte(SYNC_1);
    conv.f = determinismo_det;
    enviar_u32(conv.u);
    enviar_u32(pontos_em_linhas);
    enviar_u32(recorrencias);
    enviar_byte(anomalia_transicao_estado);

    while ((USART2->ISR & USART_ISR_TC) == 0u) {
    }
}

/*
 * main
 *
 * Comportamento
 *   Inicializa a USART2 e, para cada quadro completo, calcula o DET
 *   e responde pela mesma porta. Nao ha geracao de sinal.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Declarado como int. O laco nao termina, entao nao devolve valor.
 *
 * Globais afetadas
 *   receber_janela escreve sinal_fixo. calcular_rqa escreve a matriz
 *   e as quatro globais de resultado.
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
