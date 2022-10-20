import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import sys
sns.set()

#Ler L e h_x como parâmetro no terminal e o tempo que é para salvar
L = int(sys.argv[1])
h_x = float(sys.argv[2])
timePrint = float(sys.argv[3])
x = np.linspace(0, L, int(L/h_x))

def createDirectories():
    os.system("mkdir result")
    os.system("mkdir result/odc")
    os.system("mkdir result/mic")
    os.system("mkdir result/dc")
    os.system("mkdir result/tke")
    os.system("mkdir result/da")
    os.system("mkdir result/ant")

if timePrint < 1:
    createDirectories()

populationTitle = {
    "odc": "Destroyed oligodendrocytes",
    "mic": "Microglia",
    "dc": "Conventional dendritic cells",
    "da": "Activated dendritic cells",
    "tke": "$CD8^+$ T",
    "ant": "Antibodies igG"
}

def printMesh(time, population, type):

    x_pts, y_pts = np.meshgrid(x, x)
    max_population = np.max(population)
    if max_population == 0:
        max_population += 1
    levels = np.linspace(0, max_population, 50)

    cp = plt.contourf(x_pts, y_pts,population, levels=levels)
    plt.title(populationTitle[type], fontsize=20)
    matplotlib.rc('xtick', labelsize=15) 
    matplotlib.rc('ytick', labelsize=16)
    plt.rc('axes', labelsize=16)
    plt.rc('font', size=15)
    plt.xlabel("Millimeters")
    plt.ylabel("Millimeters")
    if type == "ant":
        plt.colorbar(cp, label="Concentration (molecules/$mm^2$)")
    else:
        plt.colorbar(cp, label="Concentration (cells/$mm^2$)")
    plt.savefig('result/'+type+'/fig'+'{:.4f}'.format(time)+type+'.png', dpi = 300)
    plt.clf()

mic_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
t_cito_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
olide_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
anticorpo_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_conv_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_ativ_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))

with open("./result/matrix/oligo.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            olide_atual[line][column] = lista[line][column]
    printMesh(timePrint, olide_atual, "odc")

with open("./result/matrix/microglia.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            mic_atual[line][column] = lista[line][column]
    printMesh(timePrint, mic_atual, "mic")

with open("./result/matrix/conventionalDC.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            dendritica_conv_atual[line][column] = lista[line][column]
    printMesh(timePrint, dendritica_conv_atual, "dc")

with open("./result/matrix/activatedDC.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            dendritica_ativ_atual[line][column] = lista[line][column]
    printMesh(timePrint, dendritica_ativ_atual, "da")

with open("./result/matrix/tCyto.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            t_cito_atual[line][column] = lista[line][column]
    printMesh(timePrint, t_cito_atual, "tke")

with open("./result/matrix/antibody.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            anticorpo_atual[line][column] = lista[line][column]
    printMesh(timePrint, anticorpo_atual, "ant")

print("Terminou o plot!!!")