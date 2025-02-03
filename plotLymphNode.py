import matplotlib.pyplot as plt
import numpy as np
import sys
# sns.set()

# Verifica se os argumentos da linha de comando são passados corretamente
if len(sys.argv) != 3:
    print("Uso: python script.py T_final h_t")
    sys.exit(1)


T_final = int(sys.argv[1])
h_t = float(sys.argv[2])

t = np.linspace(0, T_final, int(T_final/h_t))

# Lê o arquivo e divide em três resultados principais
with open("./result/tCyto.txt", 'r') as f:
    lines = f.readlines()
    

# Divide os dados em três resultados principais
resultados = [lines[i:i+1000] for i in range(0, len(lines), 1000)]

# Verifica se o número de linhas de cada resultado é consistente
for resultado in resultados:
    if len(resultado) != 1000:
        print("Erro: Número de linhas inconsistente nos resultados.")
        sys.exit(1)

# Convertendo os dados para float e plotando
plt.figure(figsize=(10, 6))

# for i, resultado in enumerate(resultados):
#     valores = [float(line.strip()) for line in resultado]
#     plt.plot(t, valores, label=f'Resultado {i+1}')
for i, resultado in enumerate(resultados):
    valores = [float(line.strip()) for line in resultado]
    label = f'{i+1}'
    if i == 0:
        label += " Epslon 0"
    elif i == 1:
        label += " Epslon 0.55"
    elif i == 2:
        label += " Epslon 0.99"
    plt.plot(t, valores, label=label)

# Adicionando pontos experimentais
experimentais = {
    0: [(0, 43.137254901960785), (0, 27.647058823529403)],
    0: [(0, 27.647058823529403)],
    1: [(14, 5.49019607843137), (30, 8.627450980392166)],
    2: [(60, 7.843137254901956)],
     #   (180, 5.098039215686276), (360, 4.705882352941171)]
}

first = True  # Variável para verificar o primeiro conjunto de dados experimentais
for key, values in experimentais.items():
    t_exp, valores_exp = zip(*values)
    if first:
        plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')
        first = False
    else:
        plt.scatter(t_exp, valores_exp, color='red')

plt.title("Lymph node - T $CD8^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cel(Cells/$mm^2$)")
plt.legend()
plt.grid(True)
plt.savefig('result/t_cito_linfonodo.png', dpi = 300)
plt.clf()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import sys

# Parâmetros de entrada
t_start = 0
T_final = int(sys.argv[1])
h_t = float(sys.argv[2])
num_linhas_max = 92  # Número máximo de linhas a ser lido de cada arquivo

# Lista dos nomes dos arquivos de dados
arquivos = ["treatment/dataExecution0.txt", "treatment/dataExecution055.txt", "treatment/dataExecution0935.txt"]

# Cores correspondentes para as curvas
cores = ['blue', 'green', 'orange']

resultados = []

# Lê cada arquivo de dados e limita a 92 linhas
for arquivo in arquivos:
    with open(arquivo, 'r') as f:
        lines = f.readlines()
        # Considera apenas as primeiras 92 linhas (ou menos, se o arquivo tiver menos linhas)
        linhas_lidas = lines[:num_linhas_max]
        resultados.append(linhas_lidas)

# Inicializa listas para armazenar linhas de plotagem e legendas
linhas_plot = []
labels_plot = []

plt.figure(figsize=(10, 6))

# Aumentando o tamanho da fonte
plt.rcParams.update({'font.size': 14})  # Define o tamanho de fonte padrão

# Convertendo os dados para float e plotando as curvas
for i in range(len(arquivos)):
    num_linhas_arquivo = len(resultados[i])
    t = np.linspace(t_start, T_final, num_linhas_arquivo)
    
    valores = [float(line.strip()) for line in resultados[i] if line.strip()]
    
    # Definindo o rótulo e a cor para a curva
    cor = cores[i]
    if i == 0:
        label = "Sem Tratamento, eps = 0"
    elif i == 1:
        label = "Eficácia constante, eps = 0.55"
    elif i == 2:
        label = "Eficácia constante, eps = 0.93"
    
    # Plotando a curva
    linha, = plt.plot(t, valores, label=label, color=cor)
    linhas_plot.append(linha)
    labels_plot.append(label)

# Adicionando pontos experimentais
experimentais = {
    0: [(0, 27.647058823529403)],
    1: [(10, 5.49019607843137)],
    2: [(30, 8.627450980392166)],
    3: [(60, 7.843137254901956)],
    4: [(90, 8.23529411764705)],
    5: [(180, 5.098039215686276)],
    6: [(360, 4.705882352941171)],
}

# Erros de desvio padrão correspondentes aos pontos experimentais
std_err = [12.16427746, 6.031063231, 9.07002932, 7.520979023, 12.28494548, 4.773915361, 4.688042942]

# Extraindo os tempos e valores experimentais
x_exp = [x for (x, _) in sum(experimentais.values(), [])]
y_exp = [y for (_, y) in sum(experimentais.values(), [])]

# Plotando os pontos experimentais com barras de erro
plt.errorbar(x_exp, y_exp, yerr=std_err, fmt='o', color='red', label='Dado Experimental', capsize=5)

# Adiciona a linha de dados experimentais à legenda
plt.legend(fontsize=12)

# Adiciona título e rótulos
plt.title("Total T $CD8^+$", fontsize=18)
plt.xlabel("Tempo (dias)", fontsize=16)
plt.ylabel("Concentração (Células/$mm^2$)", fontsize=16)
plt.ylim(bottom=0)
plt.grid(True)

# Salva o gráfico como uma imagem
plt.savefig('result/Treatment.png', dpi=300)
plt.clf()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import sys

# Parâmetros de entrada
t_start = 0
T_final = int(sys.argv[1])
h_t = float(sys.argv[2])
num_linhas_max = 362  # Número máximo de linhas a ser lido de cada arquivo

# Lista dos nomes dos arquivos de dados
arquivos = ["treatment/dataExecution088.txt"]

# Cores correspondentes para as curvas
cores = ['blue', 'green', 'orange']

resultados = []

# Lê cada arquivo de dados e limita a 362 linhas
for arquivo in arquivos:
    with open(arquivo, 'r') as f:
        lines = f.readlines()
        # Considera apenas as primeiras 362 linhas (ou menos, se o arquivo tiver menos linhas)
        linhas_lidas = lines[:num_linhas_max]
        resultados.append(linhas_lidas)

# Inicializa listas para armazenar linhas de plotagem e legendas
linhas_plot = []
labels_plot = []

plt.figure(figsize=(10, 6))

# Aumentando o tamanho da fonte
plt.rcParams.update({'font.size': 14})  # Define o tamanho de fonte padrão

# Convertendo os dados para float e plotando as curvas
for i in range(len(arquivos)):
    num_linhas_arquivo = len(resultados[i])
    t = np.linspace(t_start, T_final, num_linhas_arquivo)
    
    valores = [float(line.strip()) for line in resultados[i] if line.strip()]
    
    # Definindo o rótulo e a cor para a curva
    cor = cores[i]
    if i == 0:
        label = "Eficácia constante, eps = 0.88"
    
    # Plotando a curva
    linha, = plt.plot(t, valores, label=label, color=cor)
    linhas_plot.append(linha)
    labels_plot.append(label)

# Adicionando pontos experimentais manualmente à legenda
plt.plot([0, 14, 30, 60, 90, 180, 360], [21.39, 7.17, 13.45, 10.42, 5.03, 5.74, 5.59], 'o', color='red', label='Dado Experimental')

# Adiciona a linha de dados experimentais à legenda
plt.legend(fontsize=12)

# Adiciona título e rótulos
plt.title("Total T $CD8^+$ - Paciente 11", fontsize=18)
plt.xlabel("Tempo (dias)", fontsize=16)
plt.ylabel("Concentração (Células/$mm^2$)", fontsize=16)
plt.ylim(bottom=0)
plt.ylim(0, 45)
plt.grid(True)

# Salva o gráfico como uma imagem
plt.savefig('result/TreatmentP11.png', dpi=300)
plt.clf()
plt.show()


# import numpy as np
# import matplotlib.pyplot as plt
# import sys

# # Parâmetros de entrada
# t_start = 0
# T_final = int(sys.argv[1])
# h_t = float(sys.argv[2])
# num_linhas_max = 92  # Número máximo de linhas a ser lido de cada arquivo

# # Lista dos nomes dos arquivos de dados
# arquivos = ["treatment/dataExecution0.txt", "treatment/dataExecution055.txt", "treatment/dataExecution0935.txt"]

# # Cores correspondentes para as curvas
# cores = ['blue', 'green', 'orange']

# resultados = []

# # Lê cada arquivo de dados e limita a 92 linhas
# for arquivo in arquivos:
#     with open(arquivo, 'r') as f:
#         lines = f.readlines()
#         # Considera apenas as primeiras 92 linhas (ou menos, se o arquivo tiver menos linhas)
#         linhas_lidas = lines[:num_linhas_max]
#         resultados.append(linhas_lidas)

# # Inicializa listas para armazenar linhas de plotagem e legendas
# linhas_plot = []
# labels_plot = []

# plt.figure(figsize=(10, 6))

# # Aumentando o tamanho da fonte
# plt.rcParams.update({'font.size': 14})  # Define o tamanho de fonte padrão

# # Convertendo os dados para float e plotando as curvas
# for i in range(len(arquivos)):
#     num_linhas_arquivo = len(resultados[i])
#     t = np.linspace(t_start, T_final, num_linhas_arquivo)
    
#     valores = [float(line.strip()) for line in resultados[i] if line.strip()]
    
#     # Definindo o rótulo e a cor para a curva
#     cor = cores[i]
#     if i == 0:
#         label = "Sem Tratamento, eps = 0"
#     elif i == 1:
#         label = "Eficácia constante, eps = 0.55"
#     elif i == 2:
#         label = "Eficácia constante, eps = 0.93"
    
#     # Plotando a curva
#     linha, = plt.plot(t, valores, label=label, color=cor)
#     linhas_plot.append(linha)
#     labels_plot.append(label)

# # Adicionando pontos experimentais
# experimentais = {
#     0: [(0, 27.647058823529403)],
#     1: [(10, 5.49019607843137)],
#     2: [(30, 8.627450980392166)],
#     3: [(60, 7.843137254901956)],
#     4: [(90, 8.23529411764705)],
#     5: [(180, 5.098039215686276)],
#     6: [(360, 4.705882352941171)],
# }

# # Adiciona pontos experimentais ao gráfico
# scatter = plt.scatter(*zip(*[(x, y) for (x, y) in sum(experimentais.values(), [])]), label='Experimental Data', color='red')

# # Adiciona a linha de dados experimentais à legenda
# linhas_plot.append(scatter)
# labels_plot.append('Dado Experimental')

# # Adiciona a legenda ao gráfico
# plt.legend(linhas_plot, labels_plot, fontsize=12)

# plt.title("Total T $CD8^+$", fontsize=18)
# plt.xlabel("Tempo (dias)", fontsize=16)
# plt.ylabel("Concentração (Células/$mm^2$)", fontsize=16)
# plt.grid(True)
# plt.savefig('result/Treatment.png', dpi=300)
# plt.clf()
# plt.show()


# # Lê cada arquivo de dados
# for arquivo in arquivos:
#     with open(arquivo, 'r') as f:
#         lines = f.readlines()
#         # Verifica se o número de linhas do arquivo é consistente
#         # if len(lines) != 92:
#         #     print(f"Erro: Número de linhas inconsistente no arquivo {arquivo}")
#         #     sys.exit(1)
#         resultados.append(lines)

# # Criando array de tempo baseado no número de linhas nos arquivos de dados
# num_linhas_arquivo = len(resultados[0])  # Considerando que todos os arquivos têm o mesmo número de linhas
# t = np.linspace(0, T_final, num_linhas_arquivo)

# # Convertendo os dados para float e plotando
# plt.figure(figsize=(10, 6))

# for i, resultado in enumerate(resultados):
#     valores = [float(line.strip()) for line in resultado]
#     label = f'{i+1}'
#     if i == 0:
#         label += " Epsilon 0"
#     elif i == 1:
#         label += " Epsilon 0.55"
#     elif i == 2:
#         label += " Epsilon 0.99"
#     elif i == 3:
#         label += " Epsilon variado com o tempo"
#     plt.plot(t, valores, label=label)

# # Adicionando pontos experimentais
# experimentais = {
#     #0: [(0, 43.137254901960785), (0, 27.647058823529403)],
#     1: [(30, 8.627450980392166)],
#     2: [(60, 7.843137254901956)],
#     3: [(90, 8.23529411764705)],
# }

# first = True  # Variável para verificar o primeiro conjunto de dados experimentais
# for key, values in experimentais.items():
#     t_exp, valores_exp = zip(*values)
#     if first:
#         plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')
#         first = False
#     else:
#         plt.scatter(t_exp, valores_exp, color='red')

# plt.title("Total T $CD8^+$")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Cells/$mm^2$)")
# plt.legend()
# plt.grid(True)
# plt.savefig('result/totalTCD8.png', dpi=300)
# plt.clf()
# plt.show()

# Lista para armazenar os valores lidos do arquivo
# mediaCD8_values = []

# # Lendo os valores do arquivo
# with open("dataExecution2.txt", "r") as file:
#     for line in file:
#         mediaCD8_values.append(float(line.strip()))  

# # Pontos experimentais
# experimentais = {
#     0: [(0, 27.647058823529403)],  
#     1: [(14, 5.49019607843137), (30, 8.627450980392166)],
#     2: [(60, 7.843137254901956)],
# }

# # Criando o gráfico
# plt.figure(figsize=(10, 6))
# plt.plot(mediaCD8_values, label='Media T CD8')

# first = True  # Variável para verificar o primeiro conjunto de dados experimentais
# for key, values in experimentais.items():
#     t_exp, valores_exp = zip(*values)
#     if first:
#         plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')
#         first = False
#     else:
#         plt.scatter(t_exp, valores_exp, color='red')
    

# plt.title('Media T CD8 ao longo do tempo')
# plt.xlabel('Time (Days)')
# plt.ylabel("Concentration (Molecules/$mm^2$)")
# plt.legend()
# plt.grid(True)
# plt.savefig('result/mediaTCD8.png', dpi=300)
# plt.show()
# TL_c_vetor = np.zeros(len(t))
# TL_h_vetor = np.zeros(len(t))
# B_vetor = np.zeros(len(t))
# FL_vetor = np.zeros(len(t))
# PL_vetor = np.zeros(len(t))
# DL_vetor = np.zeros(len(t))

# with open("./result/tCyto.txt", 'r') as f:
#     lines = f.readlines()
#     TL_c = [line.rstrip() for line in lines]

# with open("./result/tHelper.txt", 'r') as f:
#     lines = f.readlines()
#     TL_h = [line.rstrip() for line in lines]

# with open("./result/bCell.txt", 'r') as f:
#     lines = f.readlines()
#     BV = [line.rstrip() for line in lines]

# with open("./result/antibody.txt", 'r') as f:
#     lines = f.readlines()
#     FL = [line.rstrip() for line in lines]

# with open("./result/plasmaCell.txt", 'r') as f:
#     lines = f.readlines()
#     PL = [line.rstrip() for line in lines]

# with open("./result/dendritic.txt", 'r') as f:
#     lines = f.readlines()
#     DL = [line.rstrip() for line in lines]


# for i in range(0,len(TL_c)):
#     TL_c_vetor[i] = TL_c[i]
#     TL_h_vetor[i] = TL_h[i]
#     B_vetor[i] = BV[i]
#     FL_vetor[i] = FL[i]
#     PL_vetor[i] = PL[i]
#     DL_vetor[i] = DL[i]
# plt.plot(t,TL_c_vetor, "-r")
# plt.title("Lymph node - T $CD8^+$")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Cells/$mm^2$)")
# plt.savefig('result/t_cito_linfonodo.png', dpi = 300)
# plt.clf()

# plt.plot(t,TL_h_vetor, "-r")
# plt.title("Lymph node - T $CD4^+$")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Cells/$mm^2$)")
# plt.savefig('result/t_helper_linfonodo.png', dpi = 300)
# plt.clf()

# plt.plot(t,B_vetor, "-r")
# plt.title("Lymph node - B Cells")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Cells/$mm^2$)")
# plt.savefig('result/b_cell_linfonodo.png', dpi = 300)
# plt.clf()

# plt.plot(t,FL_vetor, "-r")
# plt.title("Lymph node - Antibodies")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Molecules/$mm^2$)")
# plt.savefig('result/anticorpo_linfonodo.png', dpi = 300)
# plt.clf()

# plt.plot(t,PL_vetor, "-r")
# plt.title("Lymph node - Plasma-cells")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Molecules/$mm^2$)")
# plt.savefig('result/pl_cell_linfonodo.png', dpi = 300)
# plt.clf()

# plt.plot(t,DL_vetor, "-r")
# plt.title("Lymph node - Activated dendritic cells")
# plt.xlabel("Time (days)")
# plt.ylabel("Concentration (Cells/$mm^2$)")
# plt.savefig('result/dc_linfonodo.png', dpi = 300)
# plt.clf()