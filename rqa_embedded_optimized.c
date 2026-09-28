/*
 * Analise de Quantificacao de Recorrencia (RQA) para EEG sintetico.
 * Versao ultra-otimizada para ARM Cortex-M0 sem FPU.
 *
 * Usa dimensao de imersao m=1 e atraso tau=1: cada estado e uma unica
 * amplitude, sem reconstrucao multidimensional do espaco de fases.
 * Os modos interictal e ictal sao modelos sinteticos didaticos, e o
 * resultado nao constitui diagnostico clinico.
 *
 * Tres otimizacoes sobre rqa_embedded.c:
 *
 *   1. BIT-PACKING: a matriz de recorrencia armazena 1 bit por elemento
 *      em vez de 1 byte, reduzindo N*N bytes para N*N/8 bytes. Isso
 *      libera SRAM para aumentar N de 88 para 240 (quase 3x).
 *
 *   2. PONTO FIXO Q8.8: o sinal e convertido de float para int16_t
 *      (escala 256). A comparacao |x(i)-x(j)| <= EPSILON na Etapa 1
 *      usa subtracao e valor absoluto em inteiros nativos de 32 bits,
 *      eliminando fabsf e a emulacao de float no laco O(N^2).
 *      Speedup estimado: ~6x na Etapa 1.
 *
 *   3. OPERACOES BITWISE: todo acesso a matriz usa shifts e masks
 *      (>> 3, & 7, |=, >>) em vez de indexacao byte a byte. A geracao
 *      de sinal ainda usa float (sinf/cosf via seno_reduzido) mas a
 *      conversao para ponto fixo acontece so N vezes, enquanto o laco
 *      critico roda N^2 = 57.600 vezes inteiramente em inteiros.
 *
 * Restricoes de compilacao:
 *   - sem alocacao dinamica (malloc, calloc, free): tudo estatico em .bss
 *   - sem recursao: 100% iterativo
 *   - sem stdio, printf ou UART: inspecao exclusivamente por debugger
 *   - math.h apenas com sinf, cosf (precisao simples, via seno_reduzido)
 *
 * Mapa de memoria (8192 bytes de SRAM):
 *   sinal_fixo[240]   (int16_t)     480 bytes
 *   matriz_bits[240][30] (uint8_t)  7200 bytes
 *   determinismo_det      (float)      4 bytes
 *   anomalia_transicao_estado          1 byte
 *   simular_anomalia                   1 byte
 *   total .bss                      7686 bytes
 *   sobra para pilha                 506 bytes
 *
 * Comparativo com rqa_embedded.c:
 *   N:        88  ->  240  (2.7x mais amostras)
 *   Janela:  344 ms -> 937 ms (~1 segundo de EEG a 256 Hz)
 *   Matriz: 7744 B -> 7200 B (bit-packed, 93% do tamanho com 7.4x dados)
 *   Etapa 1: ~310k ciclos (float) -> ~460k ciclos (inteiro, mas 7.4x pares)
 *
 * No startup use Heap_Size = 0 e Stack_Size = 0x1F0 (496 bytes).
 */

#include <stdint.h>
#include <math.h>

#pragma GCC optimize ("Os")

/* ------------------------------------------------------------------ */
/*  Parametros do algoritmo                                            */
/* ------------------------------------------------------------------ */

#define N               240
#define L_MIN           2

#define ESCALA_FIXO     256        /* fator float -> Q8.8             */
#define EPSILON_FIXO    2048       /* epsilon empirico 8.0 em Q8.8   */
#define LIMIAR_DET_F    0.70f      /* limiar empirico, nao clinico    */

#define BYTES_POR_LINHA ((N + 7u) >> 3)   /* ceil(N/8) = 30          */

#define PI_F            3.14159265f
#define DOIS_PI_F       6.28318530f
#define MEIO_PI_F       1.57079632f
#define TRES_MEIO_PI_F  4.71238897f
#define QUARTO_PI_F     0.78539816f

#define ARG_LIBM_MAX    0.78539801f
#define SIN_PI_4_F      0.70710677f

/*
 * Frequencias angulares por amostra (rad/amostra) a Fs = 256 Hz.
 * Bandas clinicas do EEG de escalpo (classificacao IFCN).
 */
#define FS              256.0f

#define W_DELTA         0.04908739f   /*  2 Hz — delta               */
#define W_THETA         0.14726216f   /*  6 Hz — theta               */
#define W_ALPHA         0.24543693f   /* 10 Hz — alfa                */
#define W_BETA          0.51541754f   /* 21 Hz — beta                */
#define W_GAMMA         0.85902924f   /* 35 Hz — gama baixa          */

#define W_CRISE         0.07363108f   /*  3 Hz — ritmo ictal fund.   */
#define W_CRISE_H1      0.14726216f   /*  6 Hz — 1o harmonico        */
#define W_CRISE_H2      0.22089323f   /*  9 Hz — 2o harmonico        */

/* ------------------------------------------------------------------ */
/*  Buffers estaticos e variaveis de inspecao                          */
/* ------------------------------------------------------------------ */

/*
 * sinal_fixo: o EEG em ponto fixo Q8.8 (int16_t).
 *   Valor real = sinal_fixo[i] / 256.0
 *   Range: [-128.0, +127.996], suficiente para o sinal (±35).
 *
 * matriz_bits: matriz de recorrencia empacotada em bits.
 *   Bit j da linha i esta no byte matriz_bits[i][j/8], posicao j%8.
 *   N=240 e divisivel por 8, portanto nenhum bit e desperdicado.
 */
static int16_t sinal_fixo[N];
static uint8_t matriz_bits[N][BYTES_POR_LINHA];

volatile float   determinismo_det;
volatile uint8_t anomalia_transicao_estado;
volatile uint8_t simular_anomalia;

_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");
_Static_assert(sizeof(int16_t) == 2u, "int16_t deve ter 16 bits");
_Static_assert(N == 240, "N fixo em 240 para a versao otimizada");
_Static_assert((N & 7u) == 0u, "N deve ser multiplo de 8 (bit-packing)");
_Static_assert(BYTES_POR_LINHA == 30u, "30 bytes por linha da matriz");
_Static_assert(sizeof(sinal_fixo) == 480u, "vetor EEG Q8.8: 480 bytes");
_Static_assert(sizeof(matriz_bits) == 7200u, "matriz R bits: 7200 bytes");
_Static_assert(sizeof(sinal_fixo) + sizeof(matriz_bits) == 7680u,
               "buffers otimizados: 7680 bytes");
_Static_assert(sizeof(sinal_fixo) + sizeof(matriz_bits)
                   + sizeof(determinismo_det)
                   + sizeof(anomalia_transicao_estado)
                   + sizeof(simular_anomalia)
                   <= 7700u,
               "globais do RQA: 7686 bytes, sobram 506 para pilha");

/* ------------------------------------------------------------------ */
/*  Operacoes bitwise sobre a matriz empacotada                       */
/* ------------------------------------------------------------------ */

/*
 * Cada bit R(i,j) ocupa a posicao (j & 7) do byte matriz_bits[i][j>>3].
 *
 * bit_set:   OR com mascara de 1 bit  -> liga o bit sem tocar os vizinhos
 * bit_get:   shift + AND 1            -> extrai 0 ou 1
 *
 * As operacoes mapeiam diretamente para instrucoes Thumb do M0:
 *   j >> 3   -> LSR r, r, #3    (1 ciclo)
 *   j & 7    -> AND r, r, #7    (via MOV+AND, 2 ciclos)
 *   1 << k   -> MOV+LSL         (2 ciclos)
 *   |= mask  -> LDRB+ORR+STRB   (4 ciclos)
 *   >> + & 1 -> LDRB+LSR+AND    (3 ciclos)
 */
static inline void bit_set(uint_fast16_t i, uint_fast16_t j)
{
    matriz_bits[i][j >> 3] |= (uint8_t)(1u << (j & 7u));
}

static inline uint_fast8_t bit_get(uint_fast16_t i, uint_fast16_t j)
{
    return (uint_fast8_t)((matriz_bits[i][j >> 3] >> (j & 7u)) & 1u);
}

/* ------------------------------------------------------------------ */
/*  Funcoes trigonometricas com reducao de argumento segura            */
/* ------------------------------------------------------------------ */

static float seno_faixa_segura(float angulo)
{
    if (angulo < 0.0f) {
        angulo = -angulo;
        if (angulo <= ARG_LIBM_MAX) {
            return -sinf(angulo);
        }
        return -SIN_PI_4_F;
    }
    if (angulo <= ARG_LIBM_MAX) {
        return sinf(angulo);
    }
    return SIN_PI_4_F;
}

static float cosseno_faixa_segura(float angulo)
{
    if (angulo < 0.0f) {
        angulo = -angulo;
    }
    if (angulo <= ARG_LIBM_MAX) {
        return cosf(angulo);
    }
    return SIN_PI_4_F;
}

static float seno_reduzido(float x)
{
    float voltas;
    float angulo;
    float d;

    if (x < 0.0f) {
        x = -x;
    }

    voltas = x * (1.0f / DOIS_PI_F);
    angulo = (voltas - (float)(int32_t)voltas) * DOIS_PI_F;

    if (angulo < MEIO_PI_F) {
        if (angulo <= QUARTO_PI_F)
            return seno_faixa_segura(angulo);
        return cosseno_faixa_segura(MEIO_PI_F - angulo);
    }
    if (angulo < PI_F) {
        d = PI_F - angulo;
        if (d <= QUARTO_PI_F)
            return seno_faixa_segura(d);
        return cosseno_faixa_segura(angulo - MEIO_PI_F);
    }
    if (angulo < TRES_MEIO_PI_F) {
        d = angulo - PI_F;
        if (d <= QUARTO_PI_F)
            return -seno_faixa_segura(d);
        return -cosseno_faixa_segura(MEIO_PI_F - d);
    }
    d = DOIS_PI_F - angulo;
    if (d <= QUARTO_PI_F)
        return -seno_faixa_segura(d);
    return -cosseno_faixa_segura(angulo - TRES_MEIO_PI_F);
}

static float cosseno_reduzido(float x)
{
    return seno_reduzido(x + MEIO_PI_F);
}

/* ------------------------------------------------------------------ */
/*  Geracao de sinal sintetico (float -> conversao Q8.8)               */
/* ------------------------------------------------------------------ */

/*
 * Modelo didatico: soma de senoides e ruido, nao EEG clinico medido.
 * O modo 0 combina frequencias inspiradas em bandas do EEG e ruido
 * maior. O modo 1 usa 3, 6 e 9 Hz para aproximar uma forma
 * spike-and-wave, com menor ruido para representar um regime mais
 * organizado. Isso nao modela toda a fisiologia de uma crise nem
 * comprova sincronizacao cortical.
 *
 * O sinal e gerado em float (necessario para sinf/cosf) e convertido
 * imediatamente para int16_t em formato Q8.8. Nenhum vetor float e
 * armazenado — o float vive apenas em registradores durante o laco.
 *
 * Conversao: sinal_fixo[i] = (int16_t)(amostra * 256.0f)
 * Range do sinal: ~±35 -> ±8960 em Q8.8 (cabe em int16_t ±32767)
 */
void gerar_sinal_dinamica(void)
{
    uint_fast16_t i;
    uint32_t estado_ruido = 1u;

    if (simular_anomalia == 0u) {
        for (i = 0u; i < N; i++) {
            float t = (float)i;
            float ruido;
            float amostra;

            estado_ruido = (estado_ruido * 1664525u) + 1013904223u;
            ruido = (float)(estado_ruido >> 8) * (1.0f / 16777216.0f);
            ruido = (ruido * 2.0f) - 1.0f;

            amostra  = 12.0f * seno_reduzido(W_ALPHA * t);
            amostra +=  8.0f * cosseno_reduzido(W_THETA * t);
            amostra +=  5.0f * seno_reduzido(W_BETA  * t);
            amostra +=  2.0f * seno_reduzido(W_GAMMA * t);
            amostra +=  2.0f * cosseno_reduzido(W_DELTA * t);
            amostra += 14.0f * ruido;

            sinal_fixo[i] = (int16_t)(amostra * (float)ESCALA_FIXO);
        }
    } else {
        for (i = 0u; i < N; i++) {
            float t = (float)i;
            float ruido;
            float amostra;

            estado_ruido = (estado_ruido * 1664525u) + 1013904223u;
            ruido = (float)(estado_ruido >> 8) * (1.0f / 16777216.0f);
            ruido = (ruido * 2.0f) - 1.0f;

            amostra  = 22.0f * seno_reduzido(W_CRISE    * t);
            amostra += 10.0f * seno_reduzido(W_CRISE_H1 * t);
            amostra +=  5.0f * seno_reduzido(W_CRISE_H2 * t);
            amostra +=  5.0f * ruido;

            sinal_fixo[i] = (int16_t)(amostra * (float)ESCALA_FIXO);
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Calculo do RQA com bit-packing e aritmetica inteira                */
/* ------------------------------------------------------------------ */

/*
 * Etapa 1 — Construcao da matriz bit-packed em aritmetica inteira:
 *
 *   Como m=1, cada estado e somente a amplitude da respectiva amostra.
 *   Para cada par (i,j): calcula diff = sinal_fixo[i] - sinal_fixo[j]
 *   em int32_t nativo do M0. O valor absoluto e obtido por NEG
 *   condicional (sem fabsf). A comparacao com EPSILON_FIXO usa CMP
 *   inteiro. Cada ponto recorrente liga o bit correspondente via OR.
 *
 *   Ciclos por par (Cortex-M0):
 *     LDRSH + LDRSH + SUB + CMP/NEG + CMP + BHI + LDRB + ORR + STRB
 *     ~8-10 ciclos  vs  ~30-50 ciclos com fabsf (emulacao float)
 *
 * Etapa 2 — Deteccao de linhas diagonais com leitura de bits:
 *
 *   Mesma logica do rqa_embedded.c. A leitura de cada R(i,j) usa
 *   bit_get: LDRB + LSR + AND (~3 ciclos).
 *
 * Etapa 3 — DET em float (unica operacao float do RQA):
 *
 *   A divisao soma/total e feita em float para manter a precisao
 *   fracionaria do DET. E uma unica divisao, executada uma vez. O
 *   limiar 0.70 separa apenas os sinais sinteticos deste experimento.
 */
void calcular_rqa(void)
{
    uint_fast16_t i;
    uint_fast16_t j;
    uint_fast16_t k;
    uint_fast16_t comprimento;
    uint_fast16_t byte_idx;
    uint32_t total_pontos_recorrentes = 0u;
    uint32_t soma_comprimentos = 0u;

    /*
     * Etapa 1: construir a matriz bit-packed.
     * Limpa cada linha antes de setar os bits recorrentes.
     */
    for (i = 0u; i < N; i++) {
        int32_t xi = (int32_t)sinal_fixo[i];

        for (byte_idx = 0u; byte_idx < BYTES_POR_LINHA; byte_idx++) {
            matriz_bits[i][byte_idx] = 0u;
        }

        for (j = 0u; j < N; j++) {
            int32_t diff = xi - (int32_t)sinal_fixo[j];

            /* Valor absoluto sem branch (M0: SUB + ASR + EOR + SUB) */
            int32_t mask = diff >> 31;
            diff = (diff ^ mask) - mask;

            if (diff <= EPSILON_FIXO) {
                bit_set(i, j);
                if (j > i) {
                    total_pontos_recorrentes++;
                }
            }
        }
    }

    /* Etapa 2: varrer diagonais do triangulo superior */
    for (i = 0u; i < N; i++) {
        for (j = i + 1u; j < N; j++) {

            if (bit_get(i, j) == 0u) {
                continue;
            }

            if (i > 0u && bit_get((uint_fast16_t)(i - 1u),
                                  (uint_fast16_t)(j - 1u)) != 0u) {
                continue;
            }

            comprimento = 0u;
            for (k = 0u; (i + k < N) && (j + k < N); k++) {
                if (bit_get(i + k, j + k) != 0u) {
                    comprimento++;
                } else {
                    break;
                }
            }

            if (comprimento >= L_MIN) {
                soma_comprimentos += comprimento;
            }
        }
    }

    /* Etapa 3: calcular DET (unica operacao float) e classificar */
    if (total_pontos_recorrentes == 0u) {
        determinismo_det = 0.0f;
    } else {
        determinismo_det = (float)soma_comprimentos
                         / (float)total_pontos_recorrentes;
    }

    if (determinismo_det > LIMIAR_DET_F) {
        anomalia_transicao_estado = 1u;
    } else {
        anomalia_transicao_estado = 0u;
    }
}

/* ------------------------------------------------------------------ */
/*  Ponto de entrada                                                   */
/* ------------------------------------------------------------------ */

int main(void)
{
    for (;;) {
        gerar_sinal_dinamica();
        calcular_rqa();
    }
}
