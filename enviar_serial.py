"""Envia as janelas de trechos_eeg.txt para a USART2 da NUCLEO.

Cada janela sai como os bytes 'R' 'Q' e 88 float32 little-endian,
o quadro que rqa_embedded_optimized_serial.c espera. A resposta e
impressa neste terminal.

Uso:
  py -3 -m pip install pyserial
  py -3 enviar_serial.py
  py -3 enviar_serial.py COM8 trechos_eeg.txt
"""

import struct
import sys
from pathlib import Path

N = 88
SYNC = b"RQ"
RESPOSTA = struct.Struct("<2s f I I B")


def ler_trechos(caminho):
    """Devolve a lista de janelas float32 do texto gerado pela validacao.

    Parametros
        caminho: Path de trechos_eeg.txt.

    Retorno
        Lista de tuplas (id, rotulo, amostras). amostras tem N floats.
        O rotulo nao e enviado; so acompanha a impressao.

    Globais afetadas
        Le N. Nao escreve globais.
    """
    tokens = Path(caminho).read_text(encoding="ascii").split()
    cursor = 1
    quantidade = int(tokens[0])
    trechos = []
    for _ in range(quantidade):
        id_t = int(tokens[cursor])
        rotulo = int(tokens[cursor + 3])
        cursor += 4
        amostras = [float(v) for v in tokens[cursor:cursor + N]]
        cursor += N
        if len(amostras) != N:
            raise SystemExit(f"janela {id_t} nao tem {N} amostras")
        trechos.append((id_t, rotulo, amostras))
    return trechos


def porta_stlink():
    """Acha a COM virtual do ST-Link (VID 0483, PID 374B).

    Retorno
        Nome da porta, por exemplo COM9. Encerra se a placa nao aparece.

    Globais afetadas
        Nenhuma.
    """
    from serial.tools import list_ports

    for porta in list_ports.comports():
        hwid = porta.hwid.upper()
        if "VID:PID=0483:374B" in hwid:
            return porta.device
    raise SystemExit("ST-Link sem porta serial. Portas vistas: "
                     + ", ".join(p.device for p in list_ports.comports()))


def main():
    """Abre a serial, envia cada janela e imprime o DET devolvido.

    Parametros
        Nenhum. sys.argv[1] e a porta. Sem ela, ou se o nome nao
        existir, usa a COM do ST-Link. sys.argv[2] e o texto das
        janelas (padrao trechos_eeg.txt).

    Retorno
        0 se todas as respostas fecharam o sincronismo.

    Globais afetadas
        Nenhuma.
    """
    import serial
    from serial.tools import list_ports

    arquivo = "trechos_eeg.txt"
    if len(sys.argv) > 1 and sys.argv[1].upper().startswith("COM"):
        porta = sys.argv[1]
        if len(sys.argv) > 2:
            arquivo = sys.argv[2]
        existentes = {p.device.upper() for p in list_ports.comports()}
        if porta.upper() not in existentes:
            encontrada = porta_stlink()
            print(f"{porta} nao existe; usando {encontrada}", file=sys.stderr)
            porta = encontrada
    else:
        if len(sys.argv) > 1:
            arquivo = sys.argv[1]
        porta = porta_stlink()
    trechos = ler_trechos(arquivo)

    with serial.Serial(porta, 115200, timeout=5) as canal:
        canal.reset_input_buffer()
        for id_t, rotulo, amostras in trechos:
            canal.write(SYNC + struct.pack("<" + "f" * N, *amostras))
            canal.flush()
            bruto = canal.read(RESPOSTA.size)
            if len(bruto) != RESPOSTA.size:
                raise SystemExit(
                    f"trecho {id_t}: a placa nao respondeu "
                    f"({len(bruto)} de {RESPOSTA.size} bytes)"
                )
            sync, det, pontos, recorrencias, alarme = RESPOSTA.unpack(bruto)
            if sync != SYNC:
                raise SystemExit(f"trecho {id_t}: sincronismo invalido")
            print(
                f"id {id_t} rotulo {rotulo} "
                f"DET {det:.6f} pontos {pontos} "
                f"recorrencias {recorrencias} alarme {alarme}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
