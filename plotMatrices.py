import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
sns.set()

#Ler L e h_x como parâmetro no terminal e o tempo que é para salvar
L = int(sys.argv[1])
h_x = float(sys.argv[2])
finalTime = int(sys.argv[3])
numPlots = int(sys.argv[3])
timesPlots = np.linspace(0, finalTime, numPlots+1)
x = np.linspace(0, L, int(L/h_x))#+1)

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
    plt.savefig('result/'+type+'/fig'+'{:.1f}'.format(time)+type+'.png', dpi = 300)
    plt.clf()

mic_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
t_cito_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
olide_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
anticorpo_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_conv_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_ativ_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))

for t in timesPlots:
    with open("./result/matrix/oligo"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                olide_atual[line][column] = lista[line][column]
        printMesh(t, olide_atual, "odc")

    with open("./result/matrix/microglia"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                mic_atual[line][column] = lista[line][column]
        printMesh(t, mic_atual, "mic")

    with open("./result/matrix/conventionalDC"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                dendritica_conv_atual[line][column] = lista[line][column]
        printMesh(t, dendritica_conv_atual, "dc")

    with open("./result/matrix/activatedDC"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                dendritica_ativ_atual[line][column] = lista[line][column]
        printMesh(t, dendritica_ativ_atual, "da")

    with open("./result/matrix/tCyto"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                t_cito_atual[line][column] = lista[line][column]
        printMesh(t, t_cito_atual, "tke")

    with open("./result/matrix/antibody"+str(t)+".txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                anticorpo_atual[line][column] = lista[line][column]
        printMesh(t, anticorpo_atual, "ant")