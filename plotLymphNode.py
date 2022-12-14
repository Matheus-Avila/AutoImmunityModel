import matplotlib.pyplot as plt
import numpy as np
import sys
import config
# sns.set()
basePath = config.init()
T_final = int(sys.argv[1])
h_t = float(sys.argv[2])

with open(basePath + "tCyto.txt", 'r') as f:
    lines = f.readlines()
    TL_c = [line.rstrip() for line in lines]

with open(basePath + "tHelper.txt", 'r') as f:
    lines = f.readlines()
    TL_h = [line.rstrip() for line in lines]

with open(basePath + "bCell.txt", 'r') as f:
    lines = f.readlines()
    BV = [line.rstrip() for line in lines]

with open(basePath + "antibody.txt", 'r') as f:
    lines = f.readlines()
    FL = [line.rstrip() for line in lines]

with open(basePath + "plasmaCell.txt", 'r') as f:
    lines = f.readlines()
    PL = [line.rstrip() for line in lines]

with open(basePath + "dendritic.txt", 'r') as f:
    lines = f.readlines()
    DL = [line.rstrip() for line in lines]


t = np.linspace(0, T_final, int(len(TL_c)))

TL_c_vetor = np.zeros(len(TL_c))
TL_h_vetor = np.zeros(len(TL_c))
B_vetor = np.zeros(len(TL_c))
FL_vetor = np.zeros(len(TL_c))
PL_vetor = np.zeros(len(TL_c))
DL_vetor = np.zeros(len(TL_c))

for i in range(0,len(TL_c)):
    TL_c_vetor[i] = TL_c[i]
    TL_h_vetor[i] = TL_h[i]
    B_vetor[i] = BV[i]
    FL_vetor[i] = FL[i]
    PL_vetor[i] = PL[i]
    DL_vetor[i] = DL[i]
plt.plot(t,TL_c_vetor, "-r")
plt.title("Lymph node - T $CD8^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig(basePath + 't_cito_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_h_vetor, "-r")
plt.title("Lymph node - T $CD4^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig(basePath + 't_helper_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,B_vetor, "-r")
plt.title("Lymph node - B Cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig(basePath + 'b_cell_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,FL_vetor, "-r")
plt.title("Lymph node - Antibodies")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig(basePath + 'anticorpo_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,PL_vetor, "-r")
plt.title("Lymph node - Plasma-cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig(basePath + 'pl_cell_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,DL_vetor, "-r")
plt.title("Lymph node - Activated dendritic cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig(basePath + 'dc_linfonodo.png', dpi = 300)
plt.clf()
