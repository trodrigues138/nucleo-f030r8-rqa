"""Validacao do RQA em 20 trechos de um EEG externo.

Apresentacao
    Le um .edf e o sumario de crises, escolhe 10 janelas dentro da
    crise e 10 fora, grava trechos_eeg.txt e compara o determinismo
    de rqa_portable.c com a mesma conta em Python e com a PyRQA.
    A coluna duracao_s e o tempo das 88 amostras, em segundos.

Autor
    Tales Rodrigues

Data
    07/10/2026

Contexto
    Trabalho de disciplina. Confere, no PC, a conta que o firmware
    da NUCLEO-F030R8 executa na SRAM de 8 KiB.

Plataforma
    PC com Python 3, NumPy e gcc. A PyRQA e opcional. O algoritmo
    comparado e o de rqa_portable.c, dimensionado para o
    STM32F030R8 (Cortex-M0, sem FPU).

Permissao de uso
    Codigo academico de Tales Rodrigues, 2026. Pode ser lido e
    adaptado em atividades de ensino, com citacao do autor.
    Nao autoriza uso clinico. O EDF CHB-MIT permanece sob a licenca
    da base no PhysioNet.

Origem do codigo
    Proprio: leitura do EDF, escolha das janelas, conta manual do
    DET e a montagem da tabela.
    Legado: a conta manual repete calcular_rqa de rqa_portable.c.
    A PyRQA (Rawald, Sips e Marwan, 2017) e biblioteca externa.
    Gerado: nenhum. trechos_eeg.txt, resultados_c.txt e
    tabela_validacao.txt sao saidas da execucao.

Como usar
    py -3 rqa_validacao.py
    py -3 rqa_validacao.py dados\\chb01_03.edf
    py -3 rqa_validacao.py dados\\chb01_03.edf dados\\chb01-summary.txt

Entrada
    argv[1]: caminho do .edf de 16 bits. Sem ele, usa chb01_03.edf
    ao lado deste script.
    argv[2]: sumario CHB com File Name, File Start Time e
    Seizure Start/End Time. Sem ele, usa chb01-summary.txt.
    O canal lido e o de indice 0.

Saida
    trechos_eeg.txt: quantidade de trechos, e para cada um o
    cabecalho (id, inicio_amostra, inicio_s, rotulo) e 88 floats.
    resultados_c.txt: a saida de rqa_portable.c.
    tabela_validacao.txt e o mesmo texto em stdout.
    Retorno do processo: 0 se C, Python e PyRQA concordam; 1 se
    divergem ou se a PyRQA nao puder ser usada.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

N = 88                          # amostras de cada janela
EPSILON = np.float32(8.0)         # raio de recorrencia, igual ao do C
L_MIN = 2                         # menor diagonal do numerador do DET
LIMIAR_DET = 0.70                 # DET acima disto marca o alarme
N_DENTRO = 10                     # janelas com a crise inteira dentro
N_FORA = 10                       # janelas com a crise inteira fora
CANAL = 0                         # primeiro sinal do EDF
AMPLITUDE = 20.0                  # desvio-padrao da janela normalizada

PASTA = Path(__file__).resolve().parent
EDF_PADRAO = PASTA / "chb01_03.edf"
SUMARIO_PADRAO = PASTA / "chb01-summary.txt"
TRECHOS = PASTA / "trechos_eeg.txt"          # entrada gerada para o C
RESULTADOS_C = PASTA / "resultados_c.txt"    # saida gerada pelo C
TABELA = PASTA / "tabela_validacao.txt"      # comparacao das tres contas
FONTE_C = PASTA / "rqa_portable.c"
EXE_C = PASTA / "rqa_portable.exe"


def ler_edf_canal(caminho, canal):
    """Le um canal EDF de 16 bits e devolve a serie em unidade fisica.

    Comportamento
        Interpreta o cabecalho EDF e o corpo de inteiros de 16 bits,
        little-endian. Converte o canal de digital para fisico com
        (pmax - pmin) / (dmax - dmin).

    Parametros
        caminho: Path do arquivo .edf.
        canal: indice do sinal, a partir de 0.

    Retorno
        Tupla (fisico, fs, rotulo). fisico e ndarray float64 na
        unidade do cabecalho; fs e a taxa em Hz; rotulo e o nome
        do canal. Encerra o processo se o arquivo ou o canal falham.

    Globais afetadas
        Nenhuma.
    """
    with caminho.open("rb") as arquivo:
        cabeca = arquivo.read(256)
        if len(cabeca) < 256:
            raise SystemExit(f"EDF invalido: {caminho}")

        def texto(inicio, fim):
            return cabeca[inicio:fim].decode("ascii", errors="replace").strip()

        n_bytes = int(texto(184, 192))
        n_registros = int(texto(236, 244))
        duracao = float(texto(244, 252))
        n_sinais = int(texto(252, 256))
        resto = arquivo.read(n_bytes - 256)

        def lista(offset, tamanho):
            itens = []
            for i in range(n_sinais):
                trecho = resto[offset + i * tamanho:offset + (i + 1) * tamanho]
                itens.append(trecho.decode("latin1").strip())
            return itens

        cursor = 0
        rotulos = lista(cursor, 16)
        cursor += n_sinais * 16
        cursor += n_sinais * 80
        cursor += n_sinais * 8
        pmin = [float(v) for v in lista(cursor, 8)]
        cursor += n_sinais * 8
        pmax = [float(v) for v in lista(cursor, 8)]
        cursor += n_sinais * 8
        dmin = [float(v) for v in lista(cursor, 8)]
        cursor += n_sinais * 8
        dmax = [float(v) for v in lista(cursor, 8)]
        cursor += n_sinais * 8
        cursor += n_sinais * 80
        por_registro = [int(float(v)) for v in lista(cursor, 8)]
        bruto = np.fromfile(arquivo, dtype="<i2")

    if canal < 0 or canal >= n_sinais:
        raise SystemExit(f"canal {canal} fora de 0..{n_sinais - 1}")

    passo = sum(por_registro)
    desloc = sum(por_registro[:canal])
    n_por_canal = por_registro[canal]
    esperado = n_registros * passo
    if bruto.size < esperado:
        raise SystemExit("o EDF terminou antes do numero de amostras do cabecalho")

    indices = (
        np.arange(n_registros)[:, None] * passo
        + desloc
        + np.arange(n_por_canal)
    ).ravel()
    digital = bruto[indices].astype(np.float64)
    escala = (pmax[canal] - pmin[canal]) / (dmax[canal] - dmin[canal])
    fisico = (digital - dmin[canal]) * escala + pmin[canal]
    return fisico, n_por_canal / duracao, rotulos[canal]


def relogio_do_arquivo(caminho, nome_edf):
    """Segundos desde a meia-noite do horario de inicio do EDF no sumario.

    Comportamento
        Procura o bloco File Name igual a nome_edf e le File Start Time.

    Parametros
        caminho: Path do sumario CHB, texto latin1.
        nome_edf: nome do arquivo, como aparece em File Name.

    Retorno
        Inteiro, segundos desde 00:00:00. Encerra o processo se o
        sumario nao traz o inicio desse arquivo.

    Globais afetadas
        Nenhuma.
    """
    atual = None
    for linha in caminho.read_text(encoding="latin1").splitlines():
        linha = linha.strip()
        if linha.startswith("File Name:"):
            atual = linha.split(":", 1)[1].strip()
        elif atual == nome_edf and linha.startswith("File Start Time:"):
            hora, minuto, segundo = linha.split(":", 1)[1].strip().split(":")
            return int(hora) * 3600 + int(minuto) * 60 + int(segundo)
    raise SystemExit(f"o sumario nao informa o inicio de {nome_edf}")


def formatar_horario(segundos):
    """Formata segundos do dia como HH:MM:SS.

    Comportamento
        Arredonda, reduz ao modulo de 24 h e separa hora, minuto
        e segundo.

    Parametros
        segundos: instante em segundos, possivelmente acima de um dia.

    Retorno
        String "HH:MM:SS" com dois digitos em cada campo.

    Globais afetadas
        Nenhuma.
    """
    total = int(round(segundos)) % (24 * 3600)
    hora, resto = divmod(total, 3600)
    minuto, segundo = divmod(resto, 60)
    return f"{hora:02d}:{minuto:02d}:{segundo:02d}"


def crises_do_sumario(caminho, nome_edf):
    """Intervalos [inicio, fim] em segundos para um arquivo do sumario.

    Comportamento
        No bloco File Name de nome_edf, emparelha cada Seizure Start
        Time com o Seizure End Time seguinte.

    Parametros
        caminho: Path do sumario CHB.
        nome_edf: nome do arquivo procurado em File Name.

    Retorno
        Lista de pares (inicio, fim), em segundos a partir do inicio
        do registro. Lista vazia se esse arquivo nao tem crise.

    Globais afetadas
        Nenhuma.
    """
    crises = []
    atual = None
    inicio = None
    for linha in caminho.read_text(encoding="latin1").splitlines():
        linha = linha.strip()
        if linha.startswith("File Name:"):
            atual = linha.split(":", 1)[1].strip()
            inicio = None
        elif atual == nome_edf and linha.startswith("Seizure Start Time:"):
            inicio = float(linha.split(":", 1)[1].replace("seconds", "").strip())
        elif atual == nome_edf and linha.startswith("Seizure End Time:"):
            fim = float(linha.split(":", 1)[1].replace("seconds", "").strip())
            crises.append((inicio, fim))
            inicio = None
    return crises


def escolher_inicios(duracao, fs, crises):
    """Escolhe 10 janelas inteiras fora da crise e 10 inteiras dentro.

    Comportamento
        Usa a crise mais longa. As janelas de dentro cabem no intervalo
        anotado. As de fora ficam pelo menos 60 s afastadas dessa crise.
        Rejeita janela que cruza a fronteira.

    Parametros
        duracao: duracao do registro, em segundos.
        fs: taxa de amostragem, em Hz.
        crises: lista de pares (inicio, fim), em segundos.

    Retorno
        Lista de (t, rotulo). t e o inicio da janela em segundos.
        rotulo e 1 dentro da crise e 0 fora. Encerra o processo se
        nao ha crise ou se nao ha trecho fora longo o bastante.

    Globais afetadas
        Le N, N_DENTRO e N_FORA. Nao escreve globais.
    """
    largura = N / fs
    if not crises:
        raise SystemExit("o sumario nao marca crise neste EDF")

    inicio_crise, fim_crise = max(crises, key=lambda par: par[1] - par[0])
    util = fim_crise - inicio_crise - largura
    if util <= 0.0:
        raise SystemExit("a crise anotada e mais curta que uma janela de 88 amostras")

    dentro = [
        inicio_crise + k * util / (N_DENTRO - 1)
        for k in range(N_DENTRO)
    ]
    limite_fora = inicio_crise - 60.0 - largura
    if limite_fora > 60.0:
        origem = 60.0
        fim = limite_fora
    else:
        origem = fim_crise + 60.0
        fim = duracao - largura - 1.0
    if fim <= origem:
        raise SystemExit("nao ha trecho sem crise longo o bastante fora da anotacao")

    fora = [
        origem + k * (fim - origem) / (N_FORA - 1)
        for k in range(N_FORA)
    ]

    def rotulo(t):
        fim_janela = t + largura
        if t >= inicio_crise and fim_janela <= fim_crise:
            return 1
        if fim_janela <= inicio_crise or t >= fim_crise:
            return 0
        raise SystemExit(f"janela em {t} cruza a fronteira da crise")

    escolhidos = [(t, rotulo(t)) for t in fora + dentro]
    if sum(r for _, r in escolhidos) != N_DENTRO:
        raise SystemExit("a selecao nao ficou com 10 trechos dentro da crise")
    return escolhidos


def normalizar(janela):
    """Tira a media, divide pelo desvio e leva a amplitude 20.

    Comportamento
        Centraliza a janela e escala o desvio-padrao para AMPLITUDE,
        em float32, para o EPSILON do C continuar na mesma faixa.
        Desvio quase nulo vira um vetor de zeros.

    Parametros
        janela: sequencia de N amostras.

    Retorno
        ndarray float32 de tamanho N.

    Globais afetadas
        Le N e AMPLITUDE. Nao escreve globais.
    """
    x = np.asarray(janela, dtype=np.float64)
    desvio = float(x.std())
    if desvio < 1e-8:
        return np.zeros(N, dtype=np.float32)
    return ((x - float(x.mean())) / desvio * AMPLITUDE).astype(np.float32)


def gravar_trechos(caminho, sinal, fs, inicios):
    """Grava as janelas normalizadas no formato lido por rqa_portable.c.

    Comportamento
        A primeira linha e a quantidade. Cada trecho ocupa uma linha
        de cabecalho e uma linha com N amostras em ASCII.

    Parametros
        caminho: Path de trechos_eeg.txt.
        sinal: serie completa do canal.
        fs: taxa de amostragem, em Hz.
        inicios: lista de (t, rotulo) produzida por escolher_inicios.

    Retorno
        Lista de dicionarios com id, inicio_amostra, inicio_s, rotulo
        e o ndarray da janela. Encerra o processo se uma janela passa
        do fim do registro.

    Globais afetadas
        Le N. Escreve o arquivo caminho. Nao altera outras globais.
    """
    linhas = [str(len(inicios))]
    trechos = []
    for indice, (t, rotulo) in enumerate(inicios):
        i0 = int(round(t * fs))
        bloco = normalizar(sinal[i0:i0 + N])
        if bloco.size != N:
            raise SystemExit(f"janela {indice} sai do fim do arquivo")
        cabeca = f"{indice} {i0} {t:.6f} {rotulo}"
        amostras = " ".join(format(float(v), ".16g") for v in bloco)
        linhas.append(cabeca)
        linhas.append(amostras)
        trechos.append({
            "id": indice,
            "inicio_amostra": i0,
            "inicio_s": t,
            "rotulo": rotulo,
            "sinal": bloco,
        })
    caminho.write_text("\n".join(linhas) + "\n", encoding="ascii")
    return trechos


def ler_trechos(caminho):
    """Reler trechos_eeg.txt, para a conta em Python usar o mesmo texto do C.

    Comportamento
        Separa o arquivo em tokens e remonta cada janela em float32.

    Parametros
        caminho: Path do texto gravado por gravar_trechos.

    Retorno
        Lista de dicionarios com id, inicio_amostra, inicio_s, rotulo
        e sinal (ndarray float32 de N amostras).

    Globais afetadas
        Le N. Nao escreve globais.
    """
    tokens = caminho.read_text(encoding="ascii").split()
    cursor = 1
    quantidade = int(tokens[0])
    trechos = []
    for _ in range(quantidade):
        id_t = int(tokens[cursor])
        inicio_amostra = int(tokens[cursor + 1])
        inicio_s = float(tokens[cursor + 2])
        rotulo = int(tokens[cursor + 3])
        cursor += 4
        sinal = np.array([np.float32(float(v)) for v in tokens[cursor:cursor + N]],
                         dtype=np.float32)
        cursor += N
        trechos.append({
            "id": id_t,
            "inicio_amostra": inicio_amostra,
            "inicio_s": inicio_s,
            "rotulo": rotulo,
            "sinal": sinal,
        })
    return trechos


def calcular_rqa_manual(sinal):
    """Mesma conta do C: float32 e triangulo superior, j > i.

    Comportamento
        Monta R com |u(i) - u(j)| <= EPSILON, soma as diagonais de
        comprimento >= L_MIN e divide pelo numero de pares com j > i.
        Codigo proprio que reproduz calcular_rqa de rqa_portable.c.

    Parametros
        sinal: sequencia de N amostras, ja em float32 ou conversivel.

    Retorno
        Tupla (det, alarme, total, soma). det e a razao em float;
        alarme e 1 se det > LIMIAR_DET; total e o denominador;
        soma e o numerador.

    Globais afetadas
        Le N, EPSILON, L_MIN e LIMIAR_DET. Nao escreve globais.
    """
    total = 0
    matriz = [[0] * N for _ in range(N)]
    for i in range(N):
        xi = np.float32(sinal[i])
        for j in range(N):
            if abs(xi - np.float32(sinal[j])) <= EPSILON:
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
    alarme = 1 if det > LIMIAR_DET else 0
    return det, alarme, total, soma


def calcular_rqa_pyrqa(sinal):
    """Pede o DET a biblioteca PyRQA, com os mesmos parametros do C.

    Comportamento
        A PyRQA conta os dois triangulos; a razao do DET e a mesma.
        Este trecho e o unico que chama codigo externo (legado de
        terceiros). O restante da funcao so acerta o Theiler em 1,
        L_MIN e o limiar.

    Parametros
        sinal: sequencia de N amostras.

    Retorno
        Tupla (det, alarme, recorrencias, pontos), no mesmo sentido
        de calcular_rqa_manual. recorrencias e pontos ja estao
        divididos por 2, para comparar com o triangulo superior do C.

    Globais afetadas
        Le EPSILON, L_MIN e LIMIAR_DET. Nao escreve globais.
    """
    from pyrqa.analysis_type import Classic
    from pyrqa.computation import RQAComputation
    from pyrqa.metric import EuclideanMetric
    from pyrqa.neighbourhood import FixedRadius
    from pyrqa.settings import Settings
    from pyrqa.time_series import TimeSeries

    serie = TimeSeries([float(v) for v in sinal],
                       embedding_dimension=1,
                       time_delay=1)
    ajuste = Settings(serie,
                      analysis_type=Classic,
                      neighbourhood=FixedRadius(float(EPSILON)),
                      similarity_measure=EuclideanMetric,
                      theiler_corrector=1)
    resultado = RQAComputation.create(ajuste, verbose=False).run()
    resultado.min_diagonal_line_length = L_MIN
    recorrencias = int(resultado.number_of_diagonal_lines_points(1)) // 2
    pontos = int(resultado.number_of_diagonal_lines_points(L_MIN)) // 2
    det = float(resultado.determinism)
    alarme = 1 if det > LIMIAR_DET else 0
    return det, alarme, recorrencias, pontos


def localizar_gcc():
    """Acha o executavel do gcc.

    Comportamento
        Consulta o PATH e, no Windows, um pacote MinGW instalado
        pelo winget.

    Parametros
        Nenhum.

    Retorno
        Caminho do gcc, como string. Encerra o processo se nao acha.

    Globais afetadas
        Nenhuma.
    """
    encontrado = shutil.which("gcc")
    if encontrado:
        return encontrado
    candidatos = list(Path.home().glob(
        "AppData/Local/Microsoft/WinGet/Packages/**/mingw64/bin/gcc.exe"
    ))
    if candidatos:
        return str(candidatos[0])
    raise SystemExit("gcc nao encontrado no PATH")


def executar_c(trechos, resultados):
    """Compila rqa_portable.c e grava o DET das janelas.

    Comportamento
        Chama gcc -std=c11 -O2 e, em seguida, o executavel com os
        dois caminhos. Propaga o texto de erro do compilador ou do
        programa se o retorno nao e zero.

    Parametros
        trechos: Path da entrada (trechos_eeg.txt).
        resultados: Path da saida (resultados_c.txt).

    Retorno
        A linha que o C imprime em stdout, sem espacos nas pontas.
        Encerra o processo se a compilacao ou a execucao falham.

    Globais afetadas
        Le FONTE_C, EXE_C e PASTA. O subprocesso escreve EXE_C e
        o arquivo resultados.
    """
    gcc = localizar_gcc()
    compilacao = subprocess.run(
        [gcc, "-std=c11", "-O2", "-Wall", "-Wextra",
         str(FONTE_C), "-o", str(EXE_C), "-lm"],
        cwd=PASTA,
        capture_output=True,
        text=True,
    )
    if compilacao.returncode != 0:
        raise SystemExit(compilacao.stderr or compilacao.stdout)
    execucao = subprocess.run(
        [str(EXE_C), str(trechos), str(resultados)],
        cwd=PASTA,
        capture_output=True,
        text=True,
    )
    if execucao.returncode != 0:
        raise SystemExit(execucao.stderr or execucao.stdout)
    return execucao.stdout.strip()


def ler_resultados_c(caminho):
    """Le as linhas gravadas por rqa_portable.c.

    Comportamento
        Cada linha vira um dicionario com os oito campos impressos
        pelo C, ja convertidos para int ou float.

    Parametros
        caminho: Path de resultados_c.txt.

    Retorno
        Lista de dicionarios com id, inicio_amostra, inicio_s, rotulo,
        pontos, recorrencias, det e alarme.

    Globais afetadas
        Nenhuma.
    """
    linhas = []
    for bruta in caminho.read_text(encoding="ascii").splitlines():
        if not bruta.strip():
            continue
        id_t, inicio_amostra, inicio_s, rotulo, pontos, recorrencias, det, alarme = bruta.split()
        linhas.append({
            "id": int(id_t),
            "inicio_amostra": int(inicio_amostra),
            "inicio_s": float(inicio_s),
            "rotulo": int(rotulo),
            "pontos": int(pontos),
            "recorrencias": int(recorrencias),
            "det": float(det),
            "alarme": int(alarme),
        })
    return linhas


def montar_tabela(linhas, fs):
    """Formata a comparacao das tres contas.

    Comportamento
        Uma linha por trecho, com DET e alarme do C, do Python e da
        PyRQA, mais OK ou DIVERGE. A duracao impressa e N/fs e e a
        mesma para todas as janelas.

    Parametros
        linhas: lista de dicionarios montada em main, com det_c,
        det_py, det_rqa, al_c, al_py, al_rqa e impl_ok.
        fs: taxa de amostragem, em Hz.

    Retorno
        Texto da tabela, terminado em quebra de linha.

    Globais afetadas
        Le N. Nao escreve globais.
    """
    duracao = N / fs
    cabeca = (
        f"{'id':>2}  {'duracao_s':>9}  "
        f"{'DET_C':>10}  {'DET_Python':>10}  {'DET_PyRQA':>10}  "
        f"{'al_C':>4}  {'al_Py':>5}  {'al_RQA':>6}  "
        f"{'implementacoes':>15}"
    )
    corpo = [cabeca, "-" * len(cabeca)]
    iguais = 0
    for item in linhas:
        impl = "OK" if item["impl_ok"] else "DIVERGE"
        iguais += item["impl_ok"]
        det_r = f"{item['det_rqa']:.6f}" if item["det_rqa"] is not None else "     n/d"
        al_r = f"{item['al_rqa']}" if item["al_rqa"] is not None else "n/d"
        corpo.append(
            f"{item['id']:2d}  {duracao:9.6f}  "
            f"{item['det_c']:10.6f}  {item['det_py']:10.6f}  {det_r:>10}  "
            f"{item['al_c']:4d}  {item['al_py']:5d}  {al_r:>6}  "
            f"{impl:>15}"
        )
    corpo.append("")
    corpo.append(f"implementacoes concordam: {iguais}/{len(linhas)}")
    corpo.append(f"duracao_s = {N} amostras / {fs:.3f} Hz")
    return "\n".join(corpo) + "\n"


def main():
    """Le o EDF, gera as janelas, roda o C e grava a tabela.

    Comportamento
        Encadeia a leitura, a selecao, a gravacao, a execucao do C
        e as duas contas de conferencia. A PyRQA ausente nao interrompe
        o script: a tabela marca a comparacao como divergente e o
        processo termina com 1.

    Parametros
        Nenhum. Os caminhos vem de sys.argv, descritos no cabecalho
        deste modulo.

    Retorno
        0 se as implementacoes concordam em todos os trechos.
        1 se alguma diverge, se a PyRQA falta, ou o processo e
        encerrado antes por SystemExit em caso de arquivo invalido.

    Globais afetadas
        Le EDF_PADRAO, SUMARIO_PADRAO, TRECHOS, RESULTADOS_C, TABELA,
        CANAL e as constantes do RQA, por meio das funcoes chamadas.
        Escreve os tres arquivos de saida. Nao reatribui essas globais.
    """
    edf = Path(sys.argv[1]) if len(sys.argv) > 1 else EDF_PADRAO
    sumario = Path(sys.argv[2]) if len(sys.argv) > 2 else SUMARIO_PADRAO
    if not edf.is_file():
        raise SystemExit(f"EDF nao encontrado: {edf}")
    if not sumario.is_file():
        raise SystemExit(f"sumario nao encontrado: {sumario}")

    sinal, fs, rotulo_canal = ler_edf_canal(edf, CANAL)
    crises = crises_do_sumario(sumario, edf.name)
    inicio_relogio = relogio_do_arquivo(sumario, edf.name)
    inicios = escolher_inicios(len(sinal) / fs, fs, crises)
    gravar_trechos(TRECHOS, sinal, fs, inicios)
    trechos = ler_trechos(TRECHOS)
    mensagem_c = executar_c(TRECHOS, RESULTADOS_C)
    resultados_c = ler_resultados_c(RESULTADOS_C)
    if len(resultados_c) != len(trechos):
        raise SystemExit("o C gravou uma quantidade de trechos diferente da entrada")

    pyrqa_ok = True
    try:
        calcular_rqa_pyrqa(trechos[0]["sinal"])
    except Exception as erro:
        pyrqa_ok = False
        print(f"PyRQA indisponivel: {erro}", file=sys.stderr)

    linhas = []
    falha_impl = False
    for trecho, c in zip(trechos, resultados_c):
        if trecho["id"] != c["id"] or trecho["rotulo"] != c["rotulo"]:
            raise SystemExit(f"resultado C fora de ordem no trecho {trecho['id']}")
        det_py, al_py, rec_py, pontos_py = calcular_rqa_manual(trecho["sinal"])
        det_rqa = al_rqa = rec_rqa = pontos_rqa = None
        if pyrqa_ok:
            det_rqa, al_rqa, rec_rqa, pontos_rqa = calcular_rqa_pyrqa(trecho["sinal"])

        impl_ok = (
            c["pontos"] == pontos_py
            and c["recorrencias"] == rec_py
            and al_py == c["alarme"]
        )
        if det_rqa is not None:
            impl_ok = impl_ok and (
                pontos_rqa == c["pontos"]
                and rec_rqa == c["recorrencias"]
                and abs(det_rqa - det_py) <= 1e-5
                and al_rqa == c["alarme"]
            )
        else:
            impl_ok = False
        falha_impl = falha_impl or not impl_ok
        linhas.append({
            "id": trecho["id"],
            "inicio_s": trecho["inicio_s"],
            "rotulo": trecho["rotulo"],
            "det_c": c["det"],
            "det_py": det_py,
            "det_rqa": det_rqa,
            "al_c": c["alarme"],
            "al_py": al_py,
            "al_rqa": al_rqa,
            "impl_ok": impl_ok,
        })

    texto = montar_tabela(linhas, fs)
    crise_txt = ", ".join(
        f"{formatar_horario(inicio_relogio + a)}-{formatar_horario(inicio_relogio + b)}"
        for a, b in crises
    )
    prefacio = (
        f"arquivo: {edf.name}\n"
        f"canal: {CANAL} ({rotulo_canal})\n"
        f"inicio do registro: {formatar_horario(inicio_relogio)}\n"
        f"crise anotada: {crise_txt}\n"
        f"c: {mensagem_c}\n\n"
    )
    TABELA.write_text(prefacio + texto, encoding="utf-8")
    print(prefacio + texto, end="")
    return 1 if falha_impl else 0


if __name__ == "__main__":
    sys.exit(main())
