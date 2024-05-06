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
    
with open("./dataExecution.txt", 'r') as f:
    lines2 = f.readlines()

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
    #2: [(60, 7.843137254901956), (90, 8.23529411764705)]
     #   (180, 5.098039215686276), (360, 4.705882352941171)]
}

for key, values in experimentais.items():
    t_exp, valores_exp = zip(*values)
    plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')   

plt.title("Lymph node - T $CD8^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cel(Cells/$mm^2$)")
plt.legend()
plt.grid(True)
plt.savefig('result/t_cito_linfonodo.png', dpi = 300)
plt.clf()
plt.show()


# Lista para armazenar os valores lidos do arquivo
totalCD8_values = []

# Lendo os valores do arquivo
with open("dataExecution.txt", "r") as file:
    for line in file:
        totalCD8_values.append(float(line.strip()))  # Convertendo para float e adicionando à lista

# Pontos experimentais
experimentais = {
    0: [(0, 27.647058823529403)],  # Adicionando o ponto (0, 27.647058823529403) aqui
    1: [(14, 5.49019607843137), (30, 8.627450980392166)],
}

# Criando o gráfico
plt.figure(figsize=(10, 6))
plt.plot(totalCD8_values, label='Total CD8')

first = True  # Variável para verificar o primeiro conjunto de dados experimentais
for key, values in experimentais.items():
    t_exp, valores_exp = zip(*values)
    if first:
        plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')
        first = False
    else:
        plt.scatter(t_exp, valores_exp, color='red')

plt.title('Total CD8 ao longo do tempo')
plt.xlabel('Dias')
plt.ylabel('Concentração CD8')
plt.legend()
plt.grid(True)
plt.savefig('result/totalCD8.png', dpi=300)
plt.show()


# Lista para armazenar os valores lidos do arquivo
totalCD8_values = []

# Lendo os valores do arquivo
with open("dataExecution2.txt", "r") as file:
    for line in file:
        totalCD8_values.append(float(line.strip()))  # Convertendo para float e adicionando à lista

# Pontos experimentais
experimentais = {
    0: [(0, 43.137254901960785), (0, 27.647058823529403)],
    1: [(14, 5.49019607843137), (30, 8.627450980392166)]
}

# Criando o gráfico
plt.figure(figsize=(10, 6))
plt.plot(totalCD8_values, label='Media CD8')
for key, value in experimentais.items():
    t_exp, valores_exp = zip(*values)
    plt.scatter(t_exp, valores_exp, label=f'Experimental Data', color='red')  
    

plt.title('Media CD8 ao longo do tempo')
plt.xlabel('Dias')
plt.ylabel('Concentração CD8')
plt.legend()
plt.grid(True)
plt.savefig('result/mediaCD8.png', dpi=300)
plt.show()
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