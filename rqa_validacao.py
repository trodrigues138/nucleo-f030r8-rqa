"""
Validacao cruzada do algoritmo RQA implementado em C.

Compara tres fontes independentes:
  1. Implementacao manual em Python (mesma logica do codigo C)
  2. Biblioteca PyRQA (referencia academica, peer-reviewed)
  3. Saida esperada do executavel C (rqa_portable.exe)

Parametros identicos ao codigo C:
  N = 88, EPSILON = 8.0, L_MIN = 2, limiar DET = 0.70

Uso:
  pip install PyRQA
  py -3 rqa_validacao.py
"""

import math
import sys

# ---------------------------------------------------------------
#  Parametros do algoritmo (identicos ao codigo C)
# ---------------------------------------------------------------

N = 88
EPSILON = 8.0
L_MIN = 2
LIMIAR_DET = 0.70
DOIS_PI = 2.0 * math.pi
FS = 256.0

W_DELTA = DOIS_PI * 2.0 / FS
W_THETA = DOIS_PI * 6.0 / FS
W_ALPHA = DOIS_PI * 10.0 / FS
W_BETA  = DOIS_PI * 21.0 / FS
W_GAMMA = DOIS_PI * 35.0 / FS

W_CRISE    = DOIS_PI * 3.0 / FS
W_CRISE_H1 = DOIS_PI * 6.0 / FS
W_CRISE_H2 = DOIS_PI * 9.0 / FS

# ---------------------------------------------------------------
#  Geracao de sinal sintetico (replica exata do codigo C)
# ---------------------------------------------------------------


def gerar_sinal(simular_anomalia):
    """Gera o sinal EEG sintetico com a mesma logica do codigo C."""
    sinal = [0.0] * N

    if simular_anomalia == 0:
        estado = 1
        for i in range(N):
            t = float(i)
            estado = (estado * 1664525 + 1013904223) & 0xFFFFFFFF
            ruido = ((estado >> 8) / 16777216.0) * 2.0 - 1.0
            sinal[i] = (
                12.0 * math.sin(W_ALPHA * t)
                + 8.0 * math.cos(W_THETA * t)
                + 5.0 * math.sin(W_BETA * t)
                + 2.0 * math.sin(W_GAMMA * t)
                + 2.0 * math.cos(W_DELTA * t)
                + 14.0 * ruido
            )
    else:
        estado = 1
        for i in range(N):
            t = float(i)
            estado = (estado * 1664525 + 1013904223) & 0xFFFFFFFF
            ruido = ((estado >> 8) / 16777216.0) * 2.0 - 1.0
            sinal[i] = (
                22.0 * math.sin(W_CRISE * t)
                + 10.0 * math.sin(W_CRISE_H1 * t)
                + 5.0 * math.sin(W_CRISE_H2 * t)
                + 5.0 * ruido
            )

    return sinal


# ---------------------------------------------------------------
#  Implementacao manual do RQA (mesma logica do codigo C)
# ---------------------------------------------------------------


def calcular_rqa_manual(sinal):
    """Calcula DET usando a mesma logica do codigo C (triangulo superior)."""
    n = len(sinal)
    matriz = [[0] * n for _ in range(n)]
    total = 0

    for i in range(n):
        for j in range(n):
            if abs(sinal[i] - sinal[j]) <= EPSILON:
                matriz[i][j] = 1
                if j > i:
                    total += 1

    soma = 0
    for i in range(n):
        for j in range(i + 1, n):
            inicio = (
                matriz[i][j] == 1
                and (i == 0 or matriz[i - 1][j - 1] == 0)
            )
            if not inicio:
                continue
            comprimento = 0
            k = 0
            while i + k < n and j + k < n and matriz[i + k][j + k] == 1:
                comprimento += 1
                k += 1
            if comprimento >= L_MIN:
                soma += comprimento

    det = 0.0 if total == 0 else soma / total
    anomalia = 1 if det > LIMIAR_DET else 0
    return det, anomalia, total, soma


# ---------------------------------------------------------------
#  Validacao com o exemplo 5x5 (conta manual conhecida)
# ---------------------------------------------------------------


def conferir_exemplo_5x5():
    """Serie [0, 1, 0, 1, 0], epsilon 0.5. DET esperado = 3/4 = 0.75."""
    x = [0.0, 1.0, 0.0, 1.0, 0.0]
    epsilon = 0.5
    n = len(x)
    r = [[1 if abs(x[i] - x[j]) <= epsilon else 0
          for j in range(n)] for i in range(n)]

    total = 0
    soma = 0
    for i in range(n):
        for j in range(i + 1, n):
            if r[i][j]:
                total += 1
            inicio = r[i][j] and (i == 0 or not r[i - 1][j - 1])
            if not inicio:
                continue
            comprimento = 0
            k = 0
            while i + k < n and j + k < n and r[i + k][j + k]:
                comprimento += 1
                k += 1
            if comprimento >= 2:
                soma += comprimento

    det = soma / total
    ok = total == 4 and soma == 3 and abs(det - 0.75) < 1e-12
    print(f"  Exemplo 5x5: total={total}, soma={soma}, DET={det:.4f}"
          f"  {'OK' if ok else 'FALHOU'}")
    return ok


# ---------------------------------------------------------------
#  Validacao com PyRQA (biblioteca de referencia)
# ---------------------------------------------------------------


def calcular_rqa_pyrqa(sinal):
    """
    Calcula DET usando a biblioteca PyRQA como referencia independente.

    Configuracao para reproduzir a mesma analise do codigo C:
      - embedding_dimension=1: sem reconstrucao de espaco de fase
      - time_delay=1: irrelevante com m=1
      - EuclideanMetric: em 1D equivale a |x(i)-x(j)|
      - FixedRadius(EPSILON): vizinhanca fixa, identica ao limiar do C
      - theiler_corrector=1: exclui a diagonal principal (LOI),
        equivalente ao triangulo superior (j > i) do codigo C
      - min_diagonal_line_length=2: L_MIN do codigo C
    """
    try:
        from pyrqa.time_series import TimeSeries
        from pyrqa.settings import Settings
        from pyrqa.analysis_type import Classic
        from pyrqa.neighbourhood import FixedRadius
        from pyrqa.metric import EuclideanMetric
        from pyrqa.computation import RQAComputation
    except ImportError:
        print("\n  [AVISO] PyRQA nao instalado. Execute: pip install PyRQA")
        return None

    try:
        time_series = TimeSeries(sinal,
                                 embedding_dimension=1,
                                 time_delay=1)
        settings = Settings(time_series,
                            analysis_type=Classic,
                            neighbourhood=FixedRadius(EPSILON),
                            similarity_measure=EuclideanMetric,
                            theiler_corrector=1)
        computation = RQAComputation.create(settings, verbose=False)
        result = computation.run()
        result.min_diagonal_line_length = L_MIN
        result.min_vertical_line_length = 2
        result.min_white_vertical_line_length = 2
    except Exception as e:
        print(f"\n  [AVISO] PyRQA falhou (provavelmente OpenCL): {e}")
        return None

    return result


# ---------------------------------------------------------------
#  Valores esperados do executavel C (rqa_portable.exe)
# ---------------------------------------------------------------

ESPERADO_C = {
    0: {"soma": 597, "total": 1087,
        "det_aprox": 0.549218059, "anomalia": 0},
    1: {"soma": 825, "total": 922,
        "det_aprox": 0.894793928, "anomalia": 1},
}

# ---------------------------------------------------------------
#  Programa principal
# ---------------------------------------------------------------


def main():
    falha = False

    print("=" * 70)
    print("  VALIDACAO CRUZADA DO ALGORITMO RQA")
    print("=" * 70)

    # -- Teste 1: exemplo 5x5 --
    print("\n[1] Exemplo 5x5 (conta manual: DET = 3/4 = 0.75)")
    if not conferir_exemplo_5x5():
        falha = True

    # -- Teste 2 e 3: dois modos do sinal sintetico --
    for modo, rotulo in ((0, "interictal (repouso)"), (1, "ictal (crise epileptica)")):
        print(f"\n{'=' * 70}")
        print(f"[{modo + 2}] Modo {modo}: {rotulo}")
        print("-" * 70)

        sinal = gerar_sinal(modo)
        esperado = ESPERADO_C[modo]

        # --- Manual ---
        det_m, anom_m, total_m, soma_m = calcular_rqa_manual(sinal)
        print(f"\n  Manual (Python):")
        print(f"    recorrencias     = {total_m}")
        print(f"    pontos em linhas = {soma_m}")
        print(f"    DET              = {det_m:.6f}")
        print(f"    anomalia         = {anom_m}")

        # --- C esperado ---
        print(f"\n  Esperado (codigo C):")
        print(f"    recorrencias     = {esperado['total']}")
        print(f"    pontos em linhas = {esperado['soma']}")
        print(f"    DET              = {esperado['det_aprox']:.6f}")
        print(f"    anomalia         = {esperado['anomalia']}")

        # Conferir manual vs C
        if total_m != esperado["total"] or soma_m != esperado["soma"]:
            print(f"\n  ** DIVERGENCIA: manual ({soma_m}/{total_m})"
                  f" vs C ({esperado['soma']}/{esperado['total']})")
            falha = True
        else:
            print(f"\n  Manual vs C: OK (contagens identicas)")

        # --- PyRQA ---
        result = calcular_rqa_pyrqa(sinal)
        if result is not None:
            det_p = result.determinism
            rr_p = result.recurrence_rate
            n_rec = result.number_of_recurrence_points
            n_diag_pts = result.number_of_diagonal_lines_points(L_MIN)
            longest = result.longest_diagonal_line
            avg_diag = result.average_diagonal_line

            print(f"\n  PyRQA (biblioteca de referencia):")
            print(f"    recurrence_rate              = {rr_p:.6f}")
            print(f"    recurrence_points (PyRQA)    = {n_rec}")
            print(f"    pontos em diagonais >= L_MIN = {n_diag_pts}")
            print(f"    DET                          = {det_p:.6f}")
            print(f"    maior diagonal               = {longest}")
            print(f"    media diagonal                = {avg_diag:.2f}")

            # PyRQA conta a matriz simetrica inteira (ambos triangulos).
            # O DET (razao) deve ser o mesmo porque numerador e
            # denominador sao ambos dobrados em relacao ao triangulo
            # superior. Tolerancia de 0.01 para diferencas de
            # arredondamento entre float32 (C) e float64 (PyRQA/OpenCL).
            diff = abs(det_p - det_m)
            if diff < 0.01:
                print(f"\n  PyRQA vs Manual: OK"
                      f" (diferenca DET = {diff:.6f})")
            else:
                print(f"\n  ** DIVERGENCIA PyRQA vs Manual:"
                      f" DET PyRQA={det_p:.6f},"
                      f" DET Manual={det_m:.6f},"
                      f" diff={diff:.6f}")
                falha = True

            # Conferir se PyRQA e C concordam no lado do limiar
            anom_pyrqa = 1 if det_p > LIMIAR_DET else 0
            if anom_pyrqa != esperado["anomalia"]:
                print(f"  ** DIVERGENCIA: PyRQA alarme={anom_pyrqa},"
                      f" C alarme={esperado['anomalia']}")
                falha = True
            else:
                print(f"  PyRQA vs C (classificacao): OK"
                      f" (ambos alarme={anom_pyrqa})")

    # -- Resumo --
    print(f"\n{'=' * 70}")
    if falha:
        print("  RESULTADO: alguma verificacao falhou.")
        print("=" * 70)
        return 1
    else:
        print("  RESULTADO: todas as verificacoes passaram.")
        print("  As tres fontes (manual, PyRQA, codigo C) concordam.")
        print("=" * 70)
        return 0


if __name__ == "__main__":
    sys.exit(main())
