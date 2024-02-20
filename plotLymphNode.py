import matplotlib.pyplot as plt
import numpy as np
import sys
# sns.set()

T_final = int(sys.argv[1])
h_t = float(sys.argv[2])

populationTitle = {
    "thp": "Linfonodo - T $CD4^+$",
    "adc": "Linfonodo - Células Dendríticas Ativadas",
    "bcl": "Linfonodo - Células B",
    "pcl": "Linfonodo - Plasmócitos",
    "tke": "Linfonodo - T $CD8^+$",
    "ant": "Linfonodo - Anticorpos IgG"
}

t = np.linspace(0, T_final, int(T_final/h_t))

TL_c_vetor = np.zeros(len(t))
TL_h_vetor = np.zeros(len(t))
B_vetor = np.zeros(len(t))
FL_vetor = np.zeros(len(t))
PL_vetor = np.zeros(len(t))
DL_vetor = np.zeros(len(t))

with open("./result/tCyto.txt", 'r') as f:
    lines = f.readlines()
    TL_c = [line.rstrip() for line in lines]

with open("./result/tHelper.txt", 'r') as f:
    lines = f.readlines()
    TL_h = [line.rstrip() for line in lines]

with open("./result/bCell.txt", 'r') as f:
    lines = f.readlines()
    BV = [line.rstrip() for line in lines]

with open("./result/antibody.txt", 'r') as f:
    lines = f.readlines()
    FL = [line.rstrip() for line in lines]

with open("./result/plasmaCell.txt", 'r') as f:
    lines = f.readlines()
    PL = [line.rstrip() for line in lines]

with open("./result/dendritic.txt", 'r') as f:
    lines = f.readlines()
    DL = [line.rstrip() for line in lines]


for i in range(0,len(TL_c)):
    TL_c_vetor[i] = TL_c[i]
    TL_h_vetor[i] = TL_h[i]
    B_vetor[i] = BV[i]
    FL_vetor[i] = FL[i]
    PL_vetor[i] = PL[i]
    DL_vetor[i] = DL[i]

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,TL_c_vetor, "-r")
plt.title(populationTitle["tke"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/t_cito_linfonodo.png', dpi = 300)
plt.clf()

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,TL_h_vetor, "-r")
plt.title(populationTitle["thp"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/t_helper_linfonodo.png', dpi = 300)
plt.clf()

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,B_vetor, "-r")
plt.title(populationTitle["bcl"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/b_cell_linfonodo.png', dpi = 300)
plt.clf()

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,FL_vetor, "-r")
plt.title(populationTitle["ant"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/anticorpo_linfonodo.png', dpi = 300)
plt.clf()

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,PL_vetor, "-r")
plt.title(populationTitle["pcl"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/pl_cell_linfonodo.png', dpi = 300)
plt.clf()

ax = plt.gca()
ax.set_xticks([0,7,14,21,28])

plt.plot(t,DL_vetor, "-r")
plt.title(populationTitle["adc"])
plt.xlabel("Tempo (dias)")
plt.ylabel("Concentração (Células/$mm^2$)")
plt.savefig('result/dc_linfonodo.png', dpi = 300)
plt.clf()
