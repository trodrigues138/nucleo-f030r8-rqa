"""RQA no laptop. Mesma conta de rqa_eeg.c, sem placa e sem depurador.

Rode na pasta do projeto:

    py -3 rqa_laptop.py
"""

import math

N = 88
EPSILON = 8.0
L_MIN = 2
LIMIAR_DET = 0.70


def gerar_sinal(simular_anomalia):
    sinal = [0.0] * N
    if simular_anomalia == 0:
        estado = 1
        for i in range(N):
            t = float(i)
            estado = (estado * 1664525 + 1013904223) & 0xFFFFFFFF
            ruido = ((estado >> 8) / 16777216.0) * 2.0 - 1.0
            sinal[i] = (
                16.0 * math.sin(0.35 * t)
                + 11.0 * math.cos(1.07 * t)
                + 7.0 * math.sin(2.53 * t)
                + 18.0 * ruido
            )
    else:
        for i in range(N):
            sinal[i] = 20.0 * math.sin((2.0 * math.pi * i) / 16.0)
    return sinal


def calcular_rqa(sinal):
    matriz = [[0] * N for _ in range(N)]
    total = 0

    for i in range(N):
        for j in range(N):
            if abs(sinal[i] - sinal[j]) <= EPSILON:
                matriz[i][j] = 1
                if j > i:
                    total += 1

    soma = 0
    for i in range(N):
        for j in range(i + 1, N):
            inicio = matriz[i][j] == 1 and (i == 0 or matriz[i - 1][j - 1] == 0)
            if not inicio:
                continue
            comprimento = 0
            k = 0
            while i + k < N and j + k < N and matriz[i + k][j + k] == 1:
                comprimento += 1
                k += 1
            if comprimento >= L_MIN:
                soma += comprimento

    det = 0.0 if total == 0 else soma / total
    anomalia = 1 if det > LIMIAR_DET else 0
    return det, anomalia, total, soma


def conferir_exemplo_5x5():
    """Serie [0, 1, 0, 1, 0], tolerancia 0.5. Conta manual: DET = 3/4."""
    x = [0, 1, 0, 1, 0]
    epsilon = 0.5
    n = len(x)
    r = [[abs(x[i] - x[j]) <= epsilon for j in range(n)] for i in range(n)]
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
    if total != 4 or soma != 3 or abs(det - 0.75) > 1e-12:
        raise SystemExit("Exemplo 5x5 divergiu da conta manual.")
    print("Exemplo 5x5: 3/4 = 0.75. Conta manual conferida.")


def main():
    conferir_exemplo_5x5()
    print()
    print(f"{'modo':<22} {'DET':>8}  {'alarme':>6}  pontos em linhas / recorrencias")

    esperado = {
        0: (0.70, False),
        1: (0.70, True),
    }
    for modo, rotulo in (
        (0, "0 caotico"),
        (1, "1 seno puro"),
    ):
        det, alarme, total, soma = calcular_rqa(gerar_sinal(modo))
        print(f"{rotulo:<22} {det:8.3f}  {alarme:6d}  {soma} / {total}")
        teto, deve_alarmar = esperado[modo]
        alarmou = alarme == 1
        if deve_alarmar and (det <= teto or not alarmou):
            raise SystemExit("O seno puro deveria ultrapassar 0.70 e ligar o alarme.")
        if not deve_alarmar and (det > teto or alarmou):
            raise SystemExit("O sinal caotico deveria ficar em 0.70 ou abaixo, com alarme desligado.")

    print()
    print("Os dois modos ficaram em lados opostos do limiar 0.70.")


if __name__ == "__main__":
    main()
