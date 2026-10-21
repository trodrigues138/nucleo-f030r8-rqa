/*
 * =============================================================================
 * RQA portatil — determinismo (DET) de janelas de EEG lidas de arquivo
 * =============================================================================
 *
 * Apresentacao
 *   Calcula o determinismo da Analise de Quantificacao de Recorrencia
 *   com imersao m = 1 e atraso tau = 1. Nessa imersao o estado no
 *   instante i e so a amplitude u(i). Dois instantes recorrem quando
 *   |u(i) - u(j)| <= EPSILON. O DET e a fracao desses pares, contados
 *   uma vez no triangulo superior, que segue proxima ao longo de uma
 *   diagonal de comprimento >= L_MIN:
 *
 *     DET = pontos em diagonais com L >= L_MIN
 *           ------------------------------------
 *           pontos com R(i,j) = 1 e j > i
 *
 *   DET perto de 1: recorrencias alinhadas em trechos parecidos.
 *   DET perto de 0: recorrencias isoladas. O limiar 0,70 separa os
 *   regimes deste experimento. Nao e criterio clinico.
 *
 * Autor
 *  Tales Rodrigues
 *  Ayslana Riene
 *
 * Data
 *   27/09/2026 (nucleo do RQA). 07/10/2026 (esta versao de arquivo).
 *
 * Contexto
 *   Trabalho de disciplina. O algoritmo foi dimensionado para a SRAM
 *   de 8 KiB da placa NUCLEO-F030R8. Este arquivo e a bancada no PC:
 *   confirma o DET em janelas reais antes de levar a mesma conta ao
 *   microcontrolador (rqa_embedded.c).
 *
 * Plataforma alvo
 *   Algoritmo: NUCLEO-F030R8, STM32F030R8, ARM Cortex-M0, sem FPU,
 *   8 KiB de SRAM, sem heap.
 *   Este executavel: PC com gcc em C11. A memoria continua estatica,
 *   como no alvo.
 *
 * Permissao de uso (copyright)
 *   Codigo academico de Tales Rodrigues, 2026. Pode ser lido, compilado
 *   e adaptado em atividades de ensino, com citacao do autor.
 *   Nao autoriza uso clinico nem diagnostico.
 *   Registros CHB-MIT usados pela validacao nao fazem parte deste
 *   fonte e permanecem sob a licenca da base (PhysioNet).
 *
 * Origem do codigo
 *   Proprio: main, processar e o formato de entrada e saida.
 *   Legado: calcular_rqa repete a conta ja escrita em rqa_embedded.c
 *   (mesmo projeto). Nao ha trecho copiado de biblioteca externa.
 *   Gerado: nenhum. Nao ha saida de STM32CubeMX, HAL ou gerador.
 *
 * Como usar
 *   gcc -std=c11 -O2 -Wall -Wextra rqa_portable.c -o rqa_portable.exe -lm
 *   rqa_portable.exe
 *   rqa_portable.exe trechos_eeg.txt resultados_c.txt
 *   rqa_portable.exe <entrada.txt> <saida.txt>
 *
 *   Sem argumentos, le trechos_eeg.txt e grava resultados_c.txt no
 *   diretorio corrente. rqa_validacao.py gera a entrada a partir de
 *   um .edf e chama este programa.
 *
 * Entrada (texto ASCII)
 *   Primeira linha: inteiro, quantidade de trechos (> 0).
 *   Cada trecho ocupa duas linhas:
 *     id inicio_amostra inicio_s rotulo
 *       id              int
 *       inicio_amostra long, indice da primeira amostra no registro
 *       inicio_s        double, instante inicial em segundos
 *       rotulo          int, gabarito da crise: 1 dentro, 0 fora
 *     em seguida, N = 88 amostras float, separadas por espaco.
 *   O rotulo segue para a saida e nao entra no calculo.
 *
 * Saida
 *   Texto ASCII, uma linha por trecho, nesta ordem:
 *     id                    int
 *     inicio_amostra        long
 *     inicio_s              double com 6 casas
 *     rotulo                int (o mesmo da entrada)
 *     pontos_em_linhas      inteiro sem sinal, numerador do DET
 *     recorrencias          inteiro sem sinal, denominador do DET
 *     determinismo_det      float impresso em decimal
 *     anomalia_transicao_estado   0 ou 1 (1 se DET > 0,70)
 *   Sucesso: mensagem em stdout com a quantidade gravada e retorno 0.
 *   Erro de abertura ou de formato: mensagem em stderr e retorno 1.
 *   A saida e aberta em modo escrita e substitui o arquivo anterior.
 *
 * Restricoes herdadas do alvo
 *   Sem malloc, calloc, free ou recursao. math.h entra so por fabsf.
 *   Buffers estaticos: eeg_sinal[88] (352 bytes) e
 *   matriz_recorrencia[88][88] (7744 bytes), 8096 bytes no total.
 * =============================================================================
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#define N       88     /* amostras da janela; a 256 Hz, cerca de 344 ms */
#define EPSILON 8.0f   /* |u(i) - u(j)| dentro deste raio marca recorrencia */
#define L_MIN   2      /* menor diagonal contada no numerador do DET */
#define LIMIAR  0.70f  /* DET acima disto marca mudanca de regime aqui */

/* Janela corrente, lida do arquivo. Unidade: a do texto de entrada. */
static float eeg_sinal[N];

/* Recurrence plot R(i,j): 1 se as amplitudes recorrem, 0 caso contrario. */
static uint8_t matriz_recorrencia[N][N];

/* Ultimo trecho calculado. pontos_em_linhas e o numerador do DET. */
static float determinismo_det;

/* 1 se determinismo_det > LIMIAR, 0 caso contrario. */
static uint8_t anomalia_transicao_estado;

/* Numerador do DET: pontos que pertencem a diagonais com L >= L_MIN. */
static uint32_t pontos_em_linhas;

/* Denominador do DET: pares com R(i,j) = 1 e j > i. */
static uint32_t recorrencias;

/* float de 32 bits e os dois buffers somando 8096 bytes. */
_Static_assert(sizeof(float) == 4u, "float deve ter 32 bits");
_Static_assert(sizeof(eeg_sinal) + sizeof(matriz_recorrencia) == 8096u,
               "buffers do RQA devem ocupar 8096 bytes");

/*
 * calcular_rqa
 *
 * Comportamento
 *   Monta a matriz de recorrencia a partir de eeg_sinal, mede as
 *   diagonais do triangulo superior e grava DET e a classificacao.
 *
 *   A matriz inteira e preenchida: R e simetrica e a diagonal principal
 *   vale sempre 1, porque cada amostra coincide consigo mesma. A conta
 *   usa so j > i, para nao repetir o par (j, i) nem incluir essa diagonal.
 *
 *   Uma linha diagonal comeca em R(i,j) = 1 cujo antecessor R(i-1, j-1)
 *   vale 0 ou cai fora da matriz. Cada ponto recorrente entra em no
 *   maximo uma linha. Um ponto isolado (comprimento 1) fica somente
 *   no denominador. Se nao ha par recorrente, o DET fica 0.
 *
 * Parametros
 *   Nenhum.
 *
 * Retorno
 *   Nenhum (void).
 *
 * Globais afetadas
 *   Le:  eeg_sinal.
 *   Escreve: matriz_recorrencia, pontos_em_linhas, recorrencias,
 *            determinismo_det, anomalia_transicao_estado.
 */
static void calcular_rqa(void)
{
    uint_fast16_t i;           /* linha da matriz */
    uint_fast16_t j;           /* coluna da matriz */
    uint_fast16_t k;           /* passo ao longo da diagonal */
    uint_fast16_t comprimento; /* amostras contiguas da diagonal corrente */
    uint32_t total = 0u;       /* denominador: R(i,j) = 1 com j > i */
    uint32_t soma = 0u;        /* numerador: pontos em diagonais longas */

    /* Etapa 1: R(i,j) = 1 se as amplitudes estao a no maximo EPSILON. */
    for (i = 0u; i < N; i++) {
        for (j = 0u; j < N; j++) {
            if (fabsf(eeg_sinal[i] - eeg_sinal[j]) <= EPSILON) {
                matriz_recorrencia[i][j] = 1u;
                if (j > i) {
                    total++;
                }
            } else {
                matriz_recorrencia[i][j] = 0u;
            }
        }
    }

    /* Etapa 2: so o triangulo superior. O teste de inicio evita contar a mesma diagonal duas vezes. */
    for (i = 0u; i < N; i++) {
        for (j = (uint_fast16_t)(i + 1u); j < N; j++) {
            uint_fast8_t inicio = 0u; /* 1 se (i, j) e o primeiro ponto da diagonal */

            if (matriz_recorrencia[i][j] == 1u) {
                /* Borda da matriz, ou quebra da diagonal no passo anterior. */
                if (i == 0u || matriz_recorrencia[i - 1u][j - 1u] == 0u) {
                    inicio = 1u;
                }
            }

            if (inicio == 0u) {
                continue;
            }

            comprimento = 0u;
            for (k = 0u; (i + k < N) && (j + k < N); k++) {
                if (matriz_recorrencia[i + k][j + k] == 1u) {
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

    /* Etapa 3: DET = soma / total. Sem recorrencias o quociente nao existe. */
    pontos_em_linhas = soma;
    recorrencias = total;
    if (total == 0u) {
        determinismo_det = 0.0f;
    } else {
        determinismo_det = (float)soma / (float)total;
    }
    anomalia_transicao_estado = (determinismo_det > LIMIAR) ? 1u : 0u;
}

/*
 * processar
 *
 * Comportamento
 *   Abre a entrada, le a quantidade de trechos e, para cada um, o
 *   cabecalho e as N amostras. Chama calcular_rqa e acrescenta uma
 *   linha na saida. Para no primeiro trecho incompleto.
 *
 * Parametros
 *   entrada: caminho do texto com as janelas.
 *   saida:   caminho do texto a criar. O modo "w" substitui o arquivo
 *            que ja existir.
 *
 * Retorno
 *   0 se todos os trechos anunciados no cabecalho foram gravados.
 *   1 se um arquivo nao abre ou se o formato esta incompleto. Nesse
 *   caso a saida pode ter ficado vazia ou so com os trechos anteriores.
 *
 * Globais afetadas
 *   Nenhuma diretamente. A cada trecho, calcular_rqa le e escreve
 *   eeg_sinal, matriz_recorrencia e as quatro globais de resultado.
 */
static int processar(const char *entrada, const char *saida)
{
    FILE *in;             /* arquivo de janelas */
    FILE *out;            /* arquivo de resultados */
    int quantidade;       /* numero de trechos anunciado na primeira linha */
    int indice;           /* trecho corrente, de 0 a quantidade-1 */
    int gravados = 0;     /* linhas ja escritas na saida */

    in = fopen(entrada, "r");
    if (in == NULL) {
        fprintf(stderr, "nao foi possivel abrir %s\n", entrada);
        return 1;
    }

    out = fopen(saida, "w");
    if (out == NULL) {
        fprintf(stderr, "nao foi possivel criar %s\n", saida);
        fclose(in);
        return 1;
    }

    if (fscanf(in, "%d", &quantidade) != 1 || quantidade <= 0) {
        fprintf(stderr, "cabecalho invalido em %s\n", entrada);
        fclose(in);
        fclose(out);
        return 1;
    }

    for (indice = 0; indice < quantidade; indice++) {
        long inicio_amostra;      /* indice da primeira amostra no registro */
        double inicio_s;          /* instante inicial, em segundos */
        int id;                   /* identificador do trecho */
        int rotulo;               /* gabarito: 1 dentro da crise, 0 fora */
        uint_fast16_t amostra;    /* posicao de 0 a N-1 dentro da janela */

        if (fscanf(in, "%d %ld %lf %d",
                   &id, &inicio_amostra, &inicio_s, &rotulo) != 4) {
            fprintf(stderr, "cabecalho do trecho %d incompleto\n", indice);
            fclose(in);
            fclose(out);
            return 1;
        }

        for (amostra = 0u; amostra < N; amostra++) {
            if (fscanf(in, "%f", &eeg_sinal[amostra]) != 1) {
                fprintf(stderr, "amostras incompletas no trecho %d\n", id);
                fclose(in);
                fclose(out);
                return 1;
            }
        }

        calcular_rqa();
        fprintf(out,
                "%d %ld %.6f %d %lu %lu %.9g %u\n",
                id,
                inicio_amostra,
                inicio_s,
                rotulo,
                (unsigned long)pontos_em_linhas,
                (unsigned long)recorrencias,
                (double)determinismo_det,
                (unsigned)anomalia_transicao_estado);
        gravados++;
    }

    fclose(in);
    fclose(out);
    printf("gravados %d trechos em %s\n", gravados, saida);
    return 0;
}

/*
 * main
 *
 * Comportamento
 *   Escolhe os caminhos da entrada e da saida e chama processar.
 *
 * Parametros
 *   argc: quantidade de argumentos, incluindo o nome do programa.
 *   argv: argv[1] e a entrada; argv[2] e a saida. Se faltarem, usam
 *         trechos_eeg.txt e resultados_c.txt.
 *
 * Retorno
 *   O codigo de processar: 0 em sucesso, 1 em erro de arquivo ou formato.
 *
 * Globais afetadas
 *   Nenhuma diretamente. processar, por meio de calcular_rqa, altera
 *   as globais do RQA.
 */
int main(int argc, char **argv)
{
    const char *entrada = "trechos_eeg.txt"; /* caminho padrao de leitura */
    const char *saida = "resultados_c.txt";  /* caminho padrao de gravacao */

    if (argc >= 2) {
        entrada = argv[1];
    }
    if (argc >= 3) {
        saida = argv[2];
    }
    return processar(entrada, saida);
}
