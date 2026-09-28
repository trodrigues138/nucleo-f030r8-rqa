/*
 * Analise de Quantificacao de Recorrencia (RQA) para um EEG sintetico.
 * Alvo: STM32F030R8 (NUCLEO-F030R8), ARM Cortex-M0, sem FPU, 8 KiB de SRAM.
 *
 * O RQA projeta a serie x[n] num espaco de fase por comparacao direta de
 * amplitudes. A matriz de recorrencia registra os pares de amostras cuja
 * distancia absoluta nao passa de EPSILON. O determinismo (DET) mede a
 * fracao desses pares que se organiza em linhas diagonais longas. Uma
 * crise ritmica (onda quase periodica) eleva o DET; um trecho irregular
 * o mantem baixo.
 *
 * Restricoes da compilacao embarcada:
 *   - sem malloc, calloc ou free (so .bss e escalares de pilha)
 *   - sem recursao
 *   - sem stdio, printf ou UART; a leitura e pelo Live Watch (ST-Link)
 *   - math.h apenas com sinf, cosf e fabsf
 *
 * No laptop, compile com -DRQA_LAPTOP. Essa versao inclui stdio e imprime
 * o DET no terminal. A conta e a mesma.
 *   gcc -std=c11 -DRQA_LAPTOP -O2 rqa_eeg.c -o rqa_eeg.exe -lm
 *
 * Mapa da SRAM (8192 bytes):
 *   eeg_sinal[88]                 352 bytes
 *   matriz_recorrencia[88][88]   7744 bytes
 *   determinismo_det                 4 bytes
 *   anomalia_transicao_estado        1 byte
 *   simular_anomalia                 1 byte
 *   total                         8102 bytes
 *   sobra para a pilha              ~90 bytes
 *
 * No startup_stm32f030x8.s use Heap_Size = 0 e Stack_Size = 0x50.
 * Nao ligue a HAL STM32: o BSS dela nao cabe junto com a matriz.
 * O pragma abaixo compila este arquivo em -Os mesmo que o projeto esteja
 * em -O0. sinf/cosf com argumento grande chamam rem_pio2 e estouram a
 * sobra de pilha; por isso todo angulo e reduzido a cerca de pi/4 antes
 * da chamada, e o codigo gerado precisa ser enxuto.
 *
 * Live Watch (valores estaveis, relidos a cada volta do laco):
 *   simular_anomalia = 0 -> determinismo_det cerca de 0.44, anomalia = 0
 *   simular_anomalia = 1 -> determinismo_det cerca de 0.88, anomalia = 1
 * Escreva simular_anomalia pelo depurador com o alvo rodando.
 */

#include <stdint.h>
#include <math.h>

#ifdef RQA_LAPTOP
#include <stdio.h>
#endif

#pragma GCC optimize ("Os")

#define N 88

#define EPSILON 8.0f
#define L_MIN   2

#define PI_F           3.14159265f
#define DOIS_PI_F      6.28318530f
#define MEIO_PI_F      1.57079632f
#define TRES_MEIO_PI_F 4.71238897f
#define QUARTO_PI_F    0.78539816f
#define OITAVO_PI_F    0.39269908f

/* Maior argumento ainda atendido pelo polinomio de sinf/cosf (newlib),
   sem cair na reducao de Payne-Hanek. pi/4 em float32 passa desse teto. */
#define ARG_LIBM_MAX 0.78539801f
#define SIN_PI_4_F   0.70710677f

/*
 * Buffers grandes em .bss. N = 88 fecha a conta em 8096 bytes:
 * 88 * sizeof(float) + 88 * 88 * sizeof(uint8_t) = 352 + 7744.
 */
static float eeg_sinal[N];
static uint8_t matriz_recorrencia[N][N];

/* Globais de inspecao. volatile: o depurador escreve simular_anomalia
   e le as outras duas enquanto o laco roda. */
volatile float determinismo_det;
volatile uint8_t anomalia_transicao_estado;
volatile uint8_t simular_anomalia;

_Static_assert(sizeof(float) == 4u, "o alvo usa float de 32 bits");
_Static_assert(N == 88, "N fixo em 88 para a SRAM de 8 KiB");
_Static_assert(sizeof(eeg_sinal) == 352u, "vetor EEG: 352 bytes");
_Static_assert(sizeof(matriz_recorrencia) == 7744u, "matriz R: 7744 bytes");
_Static_assert(sizeof(eeg_sinal) + sizeof(matriz_recorrencia) == 8096u,
               "buffers do RQA devem ocupar 8096 bytes");
_Static_assert(sizeof(eeg_sinal) + sizeof(matriz_recorrencia)
                   + sizeof(determinismo_det)
                   + sizeof(anomalia_transicao_estado)
                   + sizeof(simular_anomalia)
                   <= 8128u,
               "globais do RQA precisam deixar pilha na SRAM de 8 KiB");

static float seno_faixa_segura(float angulo);
static float cosseno_faixa_segura(float angulo);
static float seno_reduzido(float x);
static float cosseno_reduzido(float x);

/*
 * sinf apenas com |angulo| no trecho do polinomio. No bordo pi/4,
 * sen e cos valem o mesmo e a constante evita rem_pio2.
 */
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

/*
 * sen(x) para x >= 0, com uma unica chamada a sinf ou cosf.
 * A reducao modulo 2*pi cabe em float32 porque os angulos daqui
 * ficam abaixo de ~40 voltas (indice 0..87).
 */
static float seno_reduzido(float x)
{
    float voltas;
    float angulo;
    float distancia_ao_eixo;

    if (x < 0.0f) {
        x = -x;
    }

    voltas = x * (1.0f / DOIS_PI_F);
    angulo = (voltas - (float)(int32_t)voltas) * DOIS_PI_F;

    if (angulo < MEIO_PI_F) {
        if (angulo <= QUARTO_PI_F) {
            return seno_faixa_segura(angulo);
        }
        return cosseno_faixa_segura(MEIO_PI_F - angulo);
    }

    if (angulo < PI_F) {
        distancia_ao_eixo = PI_F - angulo;
        if (distancia_ao_eixo <= QUARTO_PI_F) {
            return seno_faixa_segura(distancia_ao_eixo);
        }
        return cosseno_faixa_segura(angulo - MEIO_PI_F);
    }

    if (angulo < TRES_MEIO_PI_F) {
        distancia_ao_eixo = angulo - PI_F;
        if (distancia_ao_eixo <= QUARTO_PI_F) {
            return -seno_faixa_segura(distancia_ao_eixo);
        }
        return -cosseno_faixa_segura(MEIO_PI_F - distancia_ao_eixo);
    }

    distancia_ao_eixo = DOIS_PI_F - angulo;
    if (distancia_ao_eixo <= QUARTO_PI_F) {
        return -seno_faixa_segura(distancia_ao_eixo);
    }
    return -cosseno_faixa_segura(angulo - TRES_MEIO_PI_F);
}

/* cos(x) = sen(x + pi/2), com a mesma reducao de argumento. */
static float cosseno_reduzido(float x)
{
    return seno_reduzido(x + MEIO_PI_F);
}

/*
 * simular_anomalia == 0: mistura de senos e cossenos em frequencias
 * diferentes, mais ruido branco de um LCG. O trecho e irregular e o
 * DET fica perto de 0.44.
 * simular_anomalia == 1: seno puro de periodo 16 e amplitude 20.
 * A fase e um contador de 4 bits, entao a onda repete bit a bit e as
 * diagonais de periodicidade atravessam a matriz. DET perto de 0.88.
 * A semente do LCG e rearmada a cada chamada para o Live Watch nao
 * oscilar entre iteracoes.
 */
void gerar_sinal_dinamica(void)
{
    uint16_t i;
    uint32_t estado_ruido = 1u;

    if (simular_anomalia == 0u) {
        for (i = 0u; i < N; i++) {
            float t = (float)i;
            float ruido;
            float amostra;

            estado_ruido = (estado_ruido * 1664525u) + 1013904223u;
            ruido = (float)(estado_ruido >> 8) * (1.0f / 16777216.0f);
            ruido = (ruido * 2.0f) - 1.0f;

            amostra = 16.0f * seno_reduzido(0.35f * t);
            amostra += 11.0f * cosseno_reduzido(1.07f * t);
            amostra += 7.0f * seno_reduzido(2.53f * t);
            amostra += 18.0f * ruido;
            eeg_sinal[i] = amostra;
        }
    } else {
        for (i = 0u; i < N; i++) {
            uint16_t fase = (uint16_t)(i & 15u);
            float unidade;

            switch (fase) {
            case 0:
            case 8:
                unidade = 0.0f;
                break;
            case 4:
                unidade = 1.0f;
                break;
            case 12:
                unidade = -1.0f;
                break;
            case 2:
            case 6:
                unidade = SIN_PI_4_F;
                break;
            case 10:
            case 14:
                unidade = -SIN_PI_4_F;
                break;
            case 1:
            case 7:
                unidade = sinf(OITAVO_PI_F);
                break;
            case 9:
            case 15:
                unidade = -sinf(OITAVO_PI_F);
                break;
            case 3:
            case 5:
                unidade = cosf(OITAVO_PI_F);
                break;
            case 11:
            case 13:
                unidade = -cosf(OITAVO_PI_F);
                break;
            default:
                unidade = 0.0f;
                break;
            }

            eeg_sinal[i] = 20.0f * unidade;
        }
    }
}

/*
 * Etapa 1. Constroi R inteira:
 *   R(i,j) = 1 se |x(i) - x(j)| <= EPSILON, senao 0.
 * O denominador do DET usa so o triangulo superior (j > i). A matriz e
 * simetrica, e a diagonal principal e a linha de identidade (distancia
 * zero para todo i). Contar as duas metades, ou a identidade, duplicaria
 * pares que a etapa 2 varre uma vez so.
 *
 * Etapa 2. Tres lacos percorrem essa metade. Uma linha diagonal comeca
 * num 1 cujo antecessor (i-1, j-1) e 0 ou esta fora da matriz. O terceiro
 * laco mede o comprimento contiguo. Cada ponto recorrente entra em no
 * maximo uma linha, entao DET permanece em [0, 1].
 *   DET = (soma dos comprimentos com L >= L_MIN) / (pontos com R == 1)
 *
 * Etapa 3. DET acima de 0.70 marca a transicao para o regime ritmico.
 */
void calcular_rqa(void)
{
    uint16_t i;
    uint16_t j;
    uint16_t k;
    uint16_t comprimento;
    uint32_t total_pontos_recorrentes = 0u;
    uint32_t soma_comprimentos = 0u;
    float distancia;

    for (i = 0u; i < N; i++) {
        for (j = 0u; j < N; j++) {
            distancia = fabsf(eeg_sinal[i] - eeg_sinal[j]);
            if (distancia <= EPSILON) {
                matriz_recorrencia[i][j] = 1u;
                if (j > i) {
                    total_pontos_recorrentes++;
                }
            } else {
                matriz_recorrencia[i][j] = 0u;
            }
        }
    }

    for (i = 0u; i < N; i++) {
        for (j = (uint16_t)(i + 1u); j < N; j++) {
            uint8_t inicio_de_linha = 0u;

            if (matriz_recorrencia[i][j] == 1u) {
                if (i == 0u) {
                    inicio_de_linha = 1u;
                } else if (matriz_recorrencia[i - 1u][j - 1u] == 0u) {
                    inicio_de_linha = 1u;
                }
            }

            if (inicio_de_linha != 0u) {
                comprimento = 0u;
                for (k = 0u;
                     ((uint16_t)(i + k) < N) && ((uint16_t)(j + k) < N);
                     k++) {
                    if (matriz_recorrencia[i + k][j + k] == 1u) {
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
    }

    if (total_pontos_recorrentes == 0u) {
        determinismo_det = 0.0f;
    } else {
        determinismo_det = (float)soma_comprimentos
                         / (float)total_pontos_recorrentes;
    }

    if (determinismo_det > 0.70f) {
        anomalia_transicao_estado = 1u;
    } else {
        anomalia_transicao_estado = 0u;
    }

#ifdef RQA_LAPTOP
    printf("simular_anomalia = %u\n", (unsigned)simular_anomalia);
    printf("pontos_em_linhas = %lu\n", (unsigned long)soma_comprimentos);
    printf("recorrencias = %lu\n", (unsigned long)total_pontos_recorrentes);
    if (total_pontos_recorrentes == 0u) {
        printf("razao_exata = 0\n");
    } else {
        printf("razao_exata = %.17g\n",
               (double)soma_comprimentos / (double)total_pontos_recorrentes);
    }
    printf("determinismo_det = %.9g\n", (double)determinismo_det);
    printf("anomalia_transicao_estado = %u\n\n",
           (unsigned)anomalia_transicao_estado);
#endif
}

#ifdef RQA_LAPTOP
static int conferir_exemplo_5x5(void)
{
    const float x[5] = {0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
    const float epsilon = 0.5f;
    uint8_t recorrencia[5][5];
    unsigned total = 0u;
    unsigned soma = 0u;
    unsigned i;
    unsigned j;

    for (i = 0u; i < 5u; i++) {
        for (j = 0u; j < 5u; j++) {
            recorrencia[i][j] = (fabsf(x[i] - x[j]) <= epsilon) ? 1u : 0u;
        }
    }

    for (i = 0u; i < 5u; i++) {
        for (j = i + 1u; j < 5u; j++) {
            unsigned comprimento;
            unsigned k;
            uint8_t inicio = 0u;

            if (recorrencia[i][j] == 1u) {
                total++;
                if (i == 0u || recorrencia[i - 1u][j - 1u] == 0u) {
                    inicio = 1u;
                }
            }

            if (inicio == 0u) {
                continue;
            }

            comprimento = 0u;
            for (k = 0u; (i + k) < 5u && (j + k) < 5u; k++) {
                if (recorrencia[i + k][j + k] == 1u) {
                    comprimento++;
                } else {
                    break;
                }
            }
            if (comprimento >= 2u) {
                soma += comprimento;
            }
        }
    }

    printf("exemplo_5x5 pontos_em_linhas = %u\n", soma);
    printf("exemplo_5x5 recorrencias = %u\n", total);
    printf("exemplo_5x5 razao_exata = %.17g\n\n",
           (double)soma / (double)total);

    if (total == 4u && soma == 3u) {
        return 0;
    }
    return 1;
}

int main(void)
{
    int falha = conferir_exemplo_5x5();

    simular_anomalia = 0u;
    gerar_sinal_dinamica();
    calcular_rqa();
    if (!(determinismo_det <= 0.70f) || anomalia_transicao_estado != 0u) {
        falha = 1;
    }

    simular_anomalia = 1u;
    gerar_sinal_dinamica();
    calcular_rqa();
    if (!(determinismo_det > 0.70f) || anomalia_transicao_estado != 1u) {
        falha = 1;
    }

    if (falha == 0) {
        printf("conferencia: exemplo 5x5 = 3/4 e os dois modos em lados opostos de 0.70\n");
    } else {
        printf("conferencia: resultado fora do esperado\n");
    }
    return falha;
}
#else
int main(void)
{
    for (;;) {
        gerar_sinal_dinamica();
        calcular_rqa();
    }
}
#endif
