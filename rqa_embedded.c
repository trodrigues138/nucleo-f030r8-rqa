/*
 * Analise de Quantificacao de Recorrencia (RQA) para EEG sintetico.
 * Versao otimizada para microcontrolador ARM Cortex-M0 sem FPU.
 *
 * Esta implementacao usa dimensao de imersao m=1 e atraso tau=1:
 * cada estado e uma unica amplitude, sem reconstrucao multidimensional
 * do espaco de fases. A matriz registra pares cuja distancia absoluta
 * nao ultrapassa EPSILON. O DET mede a fracao desses pontos organizada
 * em linhas diagonais. Os modos interictal e ictal sao modelos sinteticos
 * didaticos; o resultado nao constitui diagnostico clinico.
 *
 * Restricoes de compilacao:
 *   - sem alocacao dinamica (malloc, calloc, free): tudo estatico em .bss
 *   - sem recursao: 100% iterativo
 *   - sem stdio, printf ou UART: inspecao exclusivamente por debugger
 *   - math.h apenas com sinf, cosf e fabsf (precisao simples)
 *   - sinf/cosf com argumento grande invocam rem_pio2 e estouram a
 *     pilha restante; por isso todo angulo e reduzido manualmente a
 *     no maximo pi/4 antes da chamada
 *
 * Mapa de memoria (8192 bytes de SRAM):
 *   eeg_sinal[88]                 352 bytes
 *   matriz_recorrencia[88][88]   7744 bytes
 *   determinismo_det                 4 bytes
 *   anomalia_transicao_estado        1 byte
 *   simular_anomalia                 1 byte
 *   total                         8102 bytes
 *   sobra para a pilha              ~90 bytes
 *
 * No startup use Heap_Size = 0 e Stack_Size = 0x50.
 * Nao ligue HAL: o BSS extra nao cabe junto com a matriz.
 * O pragma abaixo forca -Os neste arquivo mesmo que o projeto use -O0.
 *
 * Variaveis de inspecao por debugger (valores estaveis a cada volta):
 *   simular_anomalia = 0 -> determinismo_det ~ 0.55, anomalia = 0
 *   simular_anomalia = 1 -> determinismo_det ~ 0.89, anomalia = 1
 *   Escreva simular_anomalia pelo debugger com o alvo rodando.
 */

#include <stdint.h>
#include <math.h>

#pragma GCC optimize ("Os")

/* ------------------------------------------------------------------ */
/*  Parametros do algoritmo                                            */
/* ------------------------------------------------------------------ */

#define N       88
#define EPSILON 8.0f  /* valor empirico para os sinais sinteticos */
#define L_MIN   2

#define PI_F           3.14159265f
#define DOIS_PI_F      6.28318530f
#define MEIO_PI_F      1.57079632f

/*
 * Taxa de amostragem e frequencias angulares por amostra (rad/amostra).
 * Cada constante W_xxx = 2*pi*f/Fs, onde f e a frequencia em Hz e Fs a
 * taxa de amostragem. Os valores correspondem as bandas clinicas do EEG
 * de escalpo (classificacao IFCN).
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
#define TRES_MEIO_PI_F 4.71238897f
#define QUARTO_PI_F    0.78539816f
/* Maior argumento atendido pelo polinomio de sinf/cosf (newlib)
   sem cair na reducao de Payne-Hanek. */
#define ARG_LIBM_MAX 0.78539801f
#define SIN_PI_4_F   0.70710677f

/* ------------------------------------------------------------------ */
/*  Buffers estaticos e variaveis de inspecao                          */
/* ------------------------------------------------------------------ */

/*
 * N = 88 fecha a conta em 8096 bytes para os dois buffers:
 * 88 * sizeof(float) + 88 * 88 * sizeof(uint8_t) = 352 + 7744.
 */
static float   eeg_sinal[N];
static uint8_t matriz_recorrencia[N][N];

/* volatile: o debugger escreve simular_anomalia e le as outras duas
   enquanto o laco principal roda. */
volatile float   determinismo_det;
volatile uint8_t anomalia_transicao_estado;
volatile uint8_t simular_anomalia;

_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");
_Static_assert(N == 88, "N fixo em 88 para caber em 8 KiB de SRAM");
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

/* ------------------------------------------------------------------ */
/*  Funcoes trigonometricas com reducao de argumento segura            */
/* ------------------------------------------------------------------ */

/*
 * sinf restrito a |angulo| <= pi/4.  No bordo, sen(pi/4) = cos(pi/4)
 * e a constante SIN_PI_4_F substitui a chamada a rem_pio2.
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
 * Todas as chamadas trigonometricas usam seno_reduzido/cosseno_reduzido
 * para evitar que sinf/cosf invoquem rem_pio2 e estourem a pilha.
 * Maior argumento: W_GAMMA * 87 ~ 74.7 rad ~ 11.9 voltas (seguro).
 *
 * A semente do LCG e rearmada a cada chamada para garantir
 * reprodutibilidade entre iteracoes do laco principal.
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

            amostra  = 12.0f * seno_reduzido(W_ALPHA * t);   /* alfa 10 Hz  */
            amostra +=  8.0f * cosseno_reduzido(W_THETA * t); /* theta 6 Hz  */
            amostra +=  5.0f * seno_reduzido(W_BETA  * t);   /* beta 21 Hz  */
            amostra +=  2.0f * seno_reduzido(W_GAMMA * t);   /* gama 35 Hz  */
            amostra +=  2.0f * cosseno_reduzido(W_DELTA * t); /* delta 2 Hz  */
            amostra += 14.0f * ruido;                          /* ruido bio.  */
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

            amostra  = 22.0f * seno_reduzido(W_CRISE    * t); /* 3 Hz fund.  */
            amostra += 10.0f * seno_reduzido(W_CRISE_H1 * t); /* 6 Hz harm.  */
            amostra +=  5.0f * seno_reduzido(W_CRISE_H2 * t); /* 9 Hz harm.  */
            amostra +=  5.0f * ruido;                           /* ruido resid. */
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
}

/* ------------------------------------------------------------------ */
/*  Ponto de entrada                                                   */
/* ------------------------------------------------------------------ */

/*
 * Laco infinito: gera o sinal e recalcula o RQA continuamente.
 * O debugger le determinismo_det e anomalia_transicao_estado a cada
 * volta, e pode alterar simular_anomalia a qualquer momento.
 */
int main(void)
{
    for (;;) {
        gerar_sinal_dinamica();
        calcular_rqa();
    }
}
