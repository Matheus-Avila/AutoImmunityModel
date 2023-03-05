import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from SALib.sample import saltelli
from SALib.analyze import sobol
import math
import pandas as pd
import seaborn as sns


def printResult(indexes, title, fileName):
    
    labels = ["$d_{mic}$",
    "$d_{a}$",
    "$d_{dc}$",
    "$d_{da}$",
    "$d_{t}$",
    "$\chi$",

    "$\mu_{dc}$",
    "$\mu_m$",
    "$r_m$",
    "$r_t$",
    "$\lambda_{fm}$",
    "$b_d$",

    "$\gamma_D$",
    "$\gamma_F$",
    "$\gamma_T$",
    
    "$c_M$",
    "$c_{DC}$",
    "$c_{DA}$",
    "$c_{DL}$",
    "$c_F$",
    
    "$\\alpha_{Th}$",
    "$\\alpha_{Tc}$",
    "$\\alpha_B$",
    "$\\alpha_P$",

    "$b_{Th}$",
    "$b_{Tc}$",
    "$b_{\\rho}$",
    "$b_{\\rho_b}$",
    "$b_{\\rho_p}$",
    "$\\rho_T$",
    "$\\rho_{TC}$",
    "$\\rho_B$",
    "$\\rho_P$",
    "$\\rho_F$"]
    
    # plt.bar(labels, indexes, color ='maroon')
    plt.title(title)
    plt.barh(labels, indexes, color ='maroon')
    plt.xlabel("Sobol's indices")
    # plt.show()
    plt.savefig("./sensitivity_analysis/"+ fileName +".png", dpi = 900)
    plt.clf()


def model(d_mic, d_anti, d_dc, d_da, d_t_cit, chi, mu_dc, mu_m, r_m, r_t, lamb_f_m, b_d, gamma_D, gamma_F, gamma_T,
 cMic, cCDc, cADc, cDl, cF, alpha_T_h, alpha_T_c, alpha_B, alpha_P, b_Th, b_Tc, b_rho, b_rho_b, b_rho_p, rho_Th, rho_Tc, rho_B, rho_P, rho_F):
# def model(mu_m, d_mic):
    #Escrever parametros em um txt
    with open('./sensitivity_analysis/SA_parameters.txt', 'w') as filep:
        filep.write(str(chi)+"\n"+str(d_mic)+"\n"+str(d_dc)+"\n"+str(d_da)+"\n"+str(d_t_cit)+"\n"+str(d_anti)+
        "\n"+str(mu_m)+"\n"+str(r_m)+"\n"+str(lamb_f_m)+"\n"+str(b_d)+"\n"+str(r_t)+
        "\n"+str(mu_dc)+"\n"+str(gamma_D)+"\n"+str(gamma_F)+"\n"+str(gamma_T)+"\n"+str(alpha_T_h)+"\n"+str(alpha_T_c)+
        "\n"+str(alpha_B)+"\n"+str(alpha_P)+"\n"+str(cMic)+"\n"+str(cCDc)+"\n"+str(cADc)+"\n"+str(cDl)+"\n"+str(cF)+
        "\n"+str(b_Th)+"\n"+str(b_Tc)+"\n"+str(b_rho)+"\n"+str(b_rho_b)+"\n"+str(b_rho_p)+"\n"+str(rho_Th)+
        "\n"+str(rho_Tc)+"\n"+str(rho_B)+"\n"+str(rho_P)+"\n"+str(rho_F))
    #Executar codigo em C
    os.system("./mainOMP 2 1")#Executa o modelo C
    #Ler resultado em C
    outPut = 0
    with open("./sensitivity_analysis/SAoutput.txt", 'r') as f:
        outPut = f.readline()
    outPut = float(outPut)
    return outPut


d_mic_mean = 1520*10**-5
d_anti_mean = 1520*10**-4
d_dc_mean = 1520*10**-5
d_da_mean = 1520*10**-5
d_t_cit_mean = 1520*10**-5
chi_mean = 0.298*60*2

mu_dc_mean = 60*24*3*10**-4
mu_m_mean = 60*24*3*10**-6
r_m_mean = 5.702**10**-3
r_t_mean = 0.1
lamb_f_m_mean = 5.702*10**-3
b_d_mean = 0.001

gamma_D_mean = 0.01
gamma_F_mean = 0.3
gamma_T_mean = 2

cMic_mean = 0.1
cCDc_mean = 0.1
cADc_mean = 0.1
cDl_mean  = 0.1
cF_mean   = 0.1

alpha_T_h_mean = 0.1
alpha_T_c_mean = 0.1
alpha_B_mean = 0.1
alpha_P_mean = 1

b_T_mean = 0.17
b_Tc_mean = 0.001
b_rho_mean = 0.6
b_rho_b_mean = 3
b_rho_p_mean = 1.02
rho_T_mean = 2
rho_Tc_mean = 2
rho_B_mean = 11
rho_P_mean = 3
rho_F_mean = 5.1*10**-2

multiplyTerm = .1
upperBound = 1 + multiplyTerm
lowerBound = 1 - multiplyTerm

problem = {
    'num_vars': 34,
    'names': [ 
        'd_mic',
        'd_anti',
        'd_dc',
        'd_da',
        'd_t_cit',
        'chi',

        'mu_dc',
        'mu_m',
        'r_m',
        'r_t',
        'lamb_f_m',
        'b_d',

        'gamma_D',
        'gamma_F',
        'gamma_T',

        'cMic',
        'cCDc',
        'cADc',
        'cDl',
        'cF',
        
        'alpha_T_h',
        'alpha_T_c',
        'alpha_B',
        'alpha_P',

        'b_T',
        'b_Tc',
        'b_rho',
        'b_rho_b',
        'b_rho_p',
        'rho_T',
        'rho_Tc',
        'rho_B',
        'rho_P',
        'rho_F'
    ],
    'bounds': [
        [lowerBound*d_mic_mean, upperBound*d_mic_mean],
        [lowerBound*d_anti_mean, upperBound*d_anti_mean],
        [lowerBound*d_dc_mean, upperBound*d_dc_mean],
        [lowerBound*d_da_mean, upperBound*d_da_mean],
        [lowerBound*d_t_cit_mean, upperBound*d_t_cit_mean],
        [lowerBound*chi_mean, upperBound*chi_mean],

        [lowerBound*mu_dc_mean, upperBound*mu_dc_mean],
        [lowerBound*mu_m_mean, upperBound*mu_m_mean],
        [lowerBound*r_m_mean, upperBound*r_m_mean],
        [lowerBound*r_t_mean, upperBound*r_t_mean],
        [lowerBound*lamb_f_m_mean, upperBound*lamb_f_m_mean],
        [lowerBound*b_d_mean, upperBound*b_d_mean],

        [lowerBound*gamma_D_mean, upperBound*gamma_D_mean],
        [lowerBound*gamma_F_mean, upperBound*gamma_F_mean],
        [lowerBound*gamma_T_mean, upperBound*gamma_T_mean],

        [lowerBound*cMic_mean, upperBound*cMic_mean],
        [lowerBound*cCDc_mean, upperBound*cCDc_mean],
        [lowerBound*cADc_mean, upperBound*cADc_mean],
        [lowerBound*cDl_mean, upperBound*cDl_mean],
        [lowerBound*cF_mean, upperBound*cF_mean],

        [lowerBound*alpha_T_h_mean, upperBound*alpha_T_h_mean],
        [lowerBound*alpha_T_c_mean, upperBound*alpha_T_c_mean],
        [lowerBound*alpha_B_mean, upperBound*alpha_B_mean],
        [lowerBound*alpha_P_mean, upperBound*alpha_P_mean],
        [lowerBound*b_T_mean, upperBound*b_T_mean],
        [lowerBound*b_Tc_mean, upperBound*b_Tc_mean],
        [lowerBound*b_rho_mean, upperBound*b_rho_mean],
        [lowerBound*b_rho_b_mean, upperBound*b_rho_b_mean],
        [lowerBound*b_rho_p_mean, upperBound*b_rho_p_mean],
        [lowerBound*rho_T_mean, upperBound*rho_T_mean],
        [lowerBound*rho_Tc_mean, upperBound*rho_Tc_mean],
        [lowerBound*rho_B_mean, upperBound*rho_B_mean],
        [lowerBound*rho_P_mean, upperBound*rho_P_mean],
        [lowerBound*rho_F_mean, upperBound*rho_F_mean]
    ]
}
'''
problem = {
    'num_vars': 2,
    'names': [
        'mu_m',
        'd_mic'
    ],
    'bounds': [
        [lowerBound*mu_m_mean, upperBound*mu_m_mean],
        [lowerBound*d_mic_mean, upperBound*d_mic_mean]  
    ]
}
'''

def printHeatMap(population, title, fileName):
    df = pd.DataFrame(population)
    matrix = df.corr().round(2)
    label = ["$d_{mic}$", "$d_{a}$", "$d_{dc}$", "$d_{da}$", "$d_{t}$", "$\chi$", "$\mu_{dc}$", "$\mu_m$", "$r_m$", "$r_t$", "$\lambda_{fm}$", "$b_d$", "$\gamma_D$", "$\gamma_F$", "$\gamma_T$", "$c_M$", "$c_{DC}$", "$c_{DA}$", "$c_{DL}$", "$c_F$", "$\\alpha_{Th}$", "$\\alpha_{Tc}$", "$\\alpha_B$", "$\\alpha_P$", "$b_{Th}$", "$b_{Tc}$", "$b_{\\rho}$", "$b_{\\rho_b}$", "$b_{\\rho_p}$", "$\\rho_T$", "$\\rho_{TC}$", "$\\rho_B$", "$\\rho_P$", "$\\rho_F$"]
    sns.heatmap(matrix, xticklabels = label, yticklabels = label, annot=False, vmax=1, vmin=0, center=0.5, linewidths=.2, cmap='vlag')
    plt.title(title)
    plt.savefig("./sensitivity_analysis/"+ fileName +".png", dpi = 900)
    plt.clf()

def printMesh(population, title, fileName):
    x = np.linspace(0, 34, 34)
    x_pts, y_pts = np.meshgrid(x, x)
    for line in population:
        for j in line:
            if math.isnan(j):
                j = 0
    max_population = np.max(population)
    if max_population == 0:
        max_population += 1
    levels = np.linspace(-max_population, max_population, 5)

    cp = plt.contourf(x_pts, y_pts,population, levels=levels)
    plt.title(title)
    label = ["$d_{mic}$",
    "$d_{a}$",
    "$d_{dc}$",
    "$d_{da}$",
    "$d_{t}$",
    "$\chi$",
    "$\mu_{dc}$",
    "$\mu_m$",
    "$r_m$",
    "$r_t$",
    "$\lambda_{fm}$",
    "$b_d$",
    "$\gamma_D$",
    "$\gamma_F$",
    "$\gamma_T$",
    "$c_M$",
    "$c_{DC}$",
    "$c_{DA}$",
    "$c_{DL}$",
    "$c_F$",
    "$\\alpha_{Th}$",
    "$\\alpha_{Tc}$",
    "$\\alpha_B$",
    "$\\alpha_P$",
    "$b_{Th}$",
    "$b_{Tc}$",
    "$b_{\\rho}$",
    "$b_{\\rho_b}$",
    "$b_{\\rho_p}$",
    "$\\rho_T$",
    "$\\rho_{TC}$",
    "$\\rho_B$",
    "$\\rho_P$",
    "$\\rho_F$"]
    # plt.xticks([np.arange(0,33, step = 1)])
    # plt.yticks([np.arange(0,33, step = 1)])
    # plt.colorbar(cp, label="Concentration (cells/$mm^2$)")
    plt.savefig("./sensitivity_analysis/"+ fileName +".png", dpi = 900)
    plt.clf()

def RunSA():
    print("Running Model")
    # open("./sensitivity_analysis/SAalloutput.txt", "w").close()
    numSamples = 2048
    calcSecondOrder=True
    sample = saltelli.sample(problem, numSamples, calc_second_order=calcSecondOrder)
    Y = np.empty([sample.shape[0]])
    if calcSecondOrder:
        numOutPuts = numSamples*(2 * len(sample[0]) + 2)
    else:
        numOutPuts = numSamples*(len(sample[0]) + 2)
    
    # inputFile = open("sample.txt", "w")
    # for samp in sample:
    #     for param in range(len(samp)):
    #         inputFile.write(str(samp[param]) + "\n")
    # inputFile.close()

    numParams = problem["num_vars"]
    with open("sample.txt", "r") as file:
        for i in range(numOutPuts):
            for j in range(numParams):
                sample[i][j] = file.readline()

    # startPos = 0.975
    # startPos = int(numOutPuts * startPos)
    # sample = sample[startPos:]
    # print(len(sample))
    # evaluate the model for eah point in the input sample
    # for i in range(len(Y)):
    #     x = sample[i]
    #     Y[i] = model(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[2], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33])

    output_file = np.zeros(numOutPuts)
    i = 0
    with open("./sensitivity_analysis/SAalloutput.txt", "r") as f:
        # output_file = f.readlines()
        for y in f.read().split("\n"):
            if i == numOutPuts:
                break
            if type(y) == str:
                output_file[i] = float(y)
            i = i + 1

    # print(output_file)

    # estimate the sensitivity indices using Sobol's method
    sensitivity = sobol.analyze(problem, output_file, calc_second_order=calcSecondOrder)

    # firstorder indices
    print("First-order or main effect indices")
    print(sensitivity['S1'])
    print(sensitivity['S1_conf'])

    # higher-order indices
    print("Higher-order or total (interactions) indices")
    print(sensitivity['ST'])
    print(sensitivity['ST_conf'])

    printResult(sensitivity['S1'], 'First Order', 'firstOrder')
    printResult(sensitivity['ST'], 'Total Order', 'totalOrder')
    for i in range(len(sensitivity["S2"][0])):
        for j in range(len(sensitivity["S2"][0])):
            if math.isnan(sensitivity["S2"][i][j]):
                sensitivity["S2"][i][j] = 0

    for i in range(len(sensitivity["S1"])):
        if sensitivity["S1"][i] < 0:
            sensitivity["S1"][i] = 0

    # for line in sensitivity['S2']:
    #     print(line)
    printHeatMap(sensitivity["S2"], "Sobol's Indices Second Order", "heatmap")
    # printMesh(sensitivity["S2"], 'Second Order' , 'secondOrder')

    # count = 0
    # for line in sensitivity['S2']:
    #     print(problem["names"][count])
    #     print(line)
    #     printResult(line, 'Second Order ' +str(problem["names"][count]), 'secondOrder'+str(problem["names"][count]))
    #     count = count + 1

    return sensitivity

sensitivity = RunSA()

