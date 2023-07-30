import numpy as np

# Nome do arquivo que contém os tempos
nome_arquivo = "cuda_sai.txt"

# Lista para armazenar os tempos
tempos = []

# Função para extrair o tempo de uma linha
def extrair_tempo(linha):
    # Procurar a substring que representa o tempo na linha
    inicio = linha.find(":") + 2
    fim = linha.find(" milissegundos")
    tempo = linha[inicio:fim]
    return int(tempo)

# Lendo o arquivo e armazenando os tempos na lista
with open(nome_arquivo, "r") as arquivo:
    for linha in arquivo:
        tempo = extrair_tempo(linha)
        tempos.append(tempo)

# Exibindo a lista de tempos
print(tempos)
print(np.mean(tempos))
print(np.std(tempos))