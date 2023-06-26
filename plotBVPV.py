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
x = np.linspace(0, L, int(L/h_x))

def printMesh(type, population):

    x_pts, y_pts = np.meshgrid(x, x)
    levels = np.linspace(0, 1, 3)

    cp = plt.contourf(x_pts, y_pts,population, levels=levels)
    matplotlib.rc('xtick', labelsize=15) 
    matplotlib.rc('ytick', labelsize=16)
    plt.rc('axes', labelsize=16)
    plt.rc('font', size=15)
    plt.xlabel("Millimeters")
    plt.ylabel("Millimeters")
    plt.colorbar(cp)
    plt.savefig('result/' + type + 'fig.png', dpi = 300)
    plt.clf()

perivascularSpace = np.zeros(((int(L/h_x)), (int(L/h_x))))
bloodVessel = np.zeros(((int(L/h_x)), (int(L/h_x))))

with open("./result/bv.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            bloodVessel[line][column] = lista[line][column]
    printMesh("bv", bloodVessel)

with open("./result/pv.txt", 'r') as f:
    lista = [line.split(' ')  for line in f]
    for line in range(len(lista[0])-1):
        for column in range(len(lista[0])-1):
            perivascularSpace[line][column] = lista[line][column]
    printMesh("pv", perivascularSpace)