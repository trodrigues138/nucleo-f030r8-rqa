/*
 * Analise de Quantificacao de Recorrencia (RQA) para EEG sintetico.
 *
 * Esta implementacao usa dimensao de imersao m=1 e atraso tau=1:
 * cada estado e uma unica amplitude, sem reconstrucao multidimensional
 * do espaco de fases. A matriz registra pares cuja distancia absoluta
 * nao ultrapassa EPSILON. O DET mede a fracao desses pontos organizada
 * em linhas diagonais. Os modos interictal e ictal sao modelos sinteticos
 * didaticos; o resultado nao constitui diagnostico clinico.
 *
 * Restricoes academicas:
 *   - sem alocacao dinamica (malloc, calloc, free): tudo estatico em .bss
 *   - sem recursao: 100% iterativo
 *   - math.h apenas com funcoes de precisao simples (sinf, cosf, fabsf)
 *
 * Compilacao:
 *   gcc -std=c11 -O2 rqa_portable.c -o rqa_portable.exe -lm
 *
 * Mapa de memoria dos buffers:
 *   eeg_sinal[88]                 352 bytes
 *   matriz_recorrencia[88][88]   7744 bytes
 *   total                        8096 bytes
 */

#include <stdint.h>
#include <math.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/*  Parametros do algoritmo                                            */
/* ------------------------------------------------------------------ */

#define N       88
#define EPSILON 8.0f  /* valor empirico para os sinais sinteticos */
#define L_MIN   2

#define PI_F      3.14159265f
#define DOIS_PI_F 6.28318530f

/*
 * Taxa de amostragem e frequencias angulares por amostra (rad/amostra).
 * Cada constante W_xxx = 2*pi*f/Fs, onde f e a frequencia em Hz e Fs a
 * taxa de amostragem. Os valores correspondem as bandas clinicas do EEG
 * de escalpo segundo a classificacao da IFCN (International Federation
 * of Clinical Neurophysiology).
 */
#define FS          256.0f

#define W_DELTA     0.04908739f   /*  2 Hz — banda delta             */
#define W_THETA     0.14726216f   /*  6 Hz — banda theta             */
#define W_ALPHA     0.24543693f   /* 10 Hz — banda alfa              */
#define W_BETA      0.51541754f   /* 21 Hz — banda beta              */
#define W_GAMMA     0.85902924f   /* 35 Hz — banda gama baixa        */

#define W_CRISE     0.07363108f   /*  3 Hz — ritmo ictal fundamental */
#define W_CRISE_H1  0.14726216f   /*  6 Hz — 1o harmonico            */
#define W_CRISE_H2  0.22089323f   /*  9 Hz — 2o harmonico            */

/* ------------------------------------------------------------------ */
/*  Buffers estaticos e variaveis de inspecao                          */
/* ------------------------------------------------------------------ */

static float   eeg_sinal[N];
static uint8_t matriz_recorrencia[N][N];

float   determinismo_det;
uint8_t anomalia_transicao_estado;
uint8_t simular_anomalia;

_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");

/* ------------------------------------------------------------------ */
/*  Geracao de sinal sintetico                                         */
/* ------------------------------------------------------------------ */

/*
 * Modelo didatico: soma de senoides e ruido, nao EEG clinico medido.
 *
 * simular_anomalia == 0 (modelo interictal sintetico):
 *   Composicao inspirada em cinco bandas do EEG de escalpo: alfa
 *   dominante, theta e delta de fundo, beta e gama de baixa amplitude,
 *   mais ruido pseudoaleatorio.
 *   A mistura de frequencias nao harmonicas e o ruido produzem um
 *   trecho irregular com DET baixo (~0.55).
 *
 * simular_anomalia == 1 (modelo ictal sintetico):
 *   Ritmo de 3 Hz e harmonicos de 6 Hz e 9 Hz aproximam uma forma
 *   spike-and-wave. O menor ruido representa apenas um regime mais
 *   organizado no modelo e favorece DET alto (~0.89); nao modela toda
 *   a fisiologia de uma crise nem comprova sincronizacao cortical.
 *
 * A semente do LCG e rearmada a cada chamada para garantir
 * reprodutibilidade entre execucoes.
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

            amostra  = 12.0f * sinf(W_ALPHA * t);    /* alfa 10 Hz  */
            amostra +=  8.0f * cosf(W_THETA * t);    /* theta 6 Hz  */
            amostra +=  5.0f * sinf(W_BETA  * t);    /* beta 21 Hz  */
            amostra +=  2.0f * sinf(W_GAMMA * t);    /* gama 35 Hz  */
            amostra +=  2.0f * cosf(W_DELTA * t);    /* delta 2 Hz  */
            amostra += 14.0f * ruido;                 /* ruido bio.  */
            eeg_sinal[i] = amostra;
        }
    } else {
        for (i = 0u; i < N; i++) {
            float t = (float)i;
            float ruido;
            float amostra;

            estado_ruido = (estado_ruido * 1664525u) + 1013904223u;
            ruido = (float)(estado_ruido >> 8) * (1.0f / 16777216.0f);
            ruido = (ruido * 2.0f) - 1.0f;

            amostra  = 22.0f * sinf(W_CRISE    * t); /* 3 Hz fund.  */
            amostra += 10.0f * sinf(W_CRISE_H1 * t); /* 6 Hz harm.  */
            amostra +=  5.0f * sinf(W_CRISE_H2 * t); /* 9 Hz harm.  */
            amostra +=  5.0f * ruido;                 /* ruido resid. */
            eeg_sinal[i] = amostra;
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Calculo do RQA (Recurrence Quantification Analysis)                */
/* ------------------------------------------------------------------ */

/*
 * Etapa 1 — Construcao da matriz de recorrencia:
 *   R(i,j) = 1 se |x(i) - x(j)| <= EPSILON, senao 0.
 *   Como m=1, x(i) e somente a amplitude da amostra i.
 *   O denominador do DET conta apenas o triangulo superior (j > i).
 *   A matriz e simetrica e a diagonal principal (i == j) e trivialmente
 *   1, portanto so o triangulo superior evita contagem dupla e exclui
 *   a linha de identidade.
 *
 * Etapa 2 — Deteccao de linhas diagonais:
 *   Tres lacos percorrem o triangulo superior. Uma linha diagonal comeca
 *   num ponto R(i,j)==1 cujo antecessor (i-1, j-1) vale 0 ou esta fora
 *   da matriz. O terceiro laco mede o comprimento contiguo. Cada ponto
 *   recorrente pertence a no maximo uma linha, mantendo DET em [0, 1].
 *   DET = (soma dos comprimentos com L >= L_MIN) / (total de R(i,j)==1)
 *
 * Etapa 3 — Classificacao:
 *   O limiar empirico 0.70 separa somente os dois sinais deste teste;
 *   nao e um limiar clinico validado.
 */
void calcular_rqa(void)
{
    uint_fast16_t i;
    uint_fast16_t j;
    uint_fast16_t k;
    uint_fast16_t comprimento;
    uint32_t total_pontos_recorrentes = 0u;
    uint32_t soma_comprimentos = 0u;
    float distancia;

    /* Etapa 1: construir a matriz e contar recorrencias no triangulo superior */
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

    /* Etapa 2: varrer diagonais do triangulo superior */
    for (i = 0u; i < N; i++) {
        for (j = i + 1u; j < N; j++) {
            uint_fast8_t inicio_de_linha = 0u;

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
                     (i + k < N) && (j + k < N);
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

    /* Etapa 3: calcular DET e classificar */
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
}

/* ------------------------------------------------------------------ */
/*  Validacao: exemplo 5x5 com conta manual                           */
/* ------------------------------------------------------------------ */

/*
 * Serie [0, 1, 0, 1, 0], epsilon = 0.5.
 * Matriz de recorrencia (1 onde |x(i)-x(j)| <= 0.5):
 *
 *       0  1  2  3  4
 *   0 [ 1  0  1  0  1 ]
 *   1 [ 0  1  0  1  0 ]
 *   2 [ 1  0  1  0  1 ]
 *   3 [ 0  1  0  1  0 ]
 *   4 [ 1  0  1  0  1 ]
 *
 * Triangulo superior: (0,2)=1, (0,4)=1, (1,3)=1, (2,4)=1  -> total = 4
 * Linha diagonal d=2: (0,2),(1,3),(2,4) -> comprimento 3, >= L_MIN
 * Ponto isolado (0,4): comprimento 1, < L_MIN
 *
 * DET = 3 / 4 = 0.75
 */
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

/* ------------------------------------------------------------------ */
/*  Ponto de entrada                                                   */
/* ------------------------------------------------------------------ */

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
        printf("conferencia: exemplo 5x5 = 3/4 e os dois modos em lados"
               " opostos de 0.70\n");
    } else {
        printf("conferencia: resultado fora do esperado\n");
    }
    return falha;
}
