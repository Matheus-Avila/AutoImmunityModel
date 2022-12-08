import numpy as np
import matplotlib.pyplot as plt
import os
from SALib.sample import saltelli
from SALib.analyze import sobol

def modelo(chi, d_mic, mu_m, r_m, d_dc, d_da, d_t_cit, d_anti, lamb_f_m, b_d, r_t, mu_dc, gamma_D, gamma_F, gamma_T, alpha_T_h, alpha_T_c, alpha_B, alpha_P,
cMic, cCDc, cADc, cDl, cF, b_Th, b_Tc, b_rho, b_rho_b, b_rho_p, rho_Th, rho_Tc, rho_B, rho_P, rho_F):
    #Escrever parametros em um txt
    with open('DE_parameters.txt', 'w') as filep:
        filep.write(str(chi)+"\n"+str(d_mic)+"\n"+str(d_dc)+"\n"+str(d_da)+"\n"+str(d_t_cit)+"\n"+str(d_anti)+
        "\n"+str(mu_m)+"\n"+str(r_m)+"\n"+str(lamb_f_m)+"\n"+str(b_d)+"\n"+str(r_t)+
        "\n"+str(mu_dc)+"\n"+str(gamma_D)+"\n"+str(gamma_F)+"\n"+str(gamma_T)+"\n"+str(alpha_T_h)+"\n"+str(alpha_T_c)+
        "\n"+str(alpha_B)+"\n"+str(alpha_P)+"\n"+str(cMic)+"\n"+str(cCDc)+"\n"+str(cADc)+"\n"+str(cDl)+"\n"+str(cF)+
        "\n"+str(b_Th)+"\n"+str(b_Tc)+"\n"+str(b_rho)+"\n"+str(b_rho_b)+"\n"+str(b_rho_p)+"\n"+str(rho_Th)+
        "\n"+str(rho_Tc)+"\n"+str(rho_B)+"\n"+str(rho_P)+"\n"+str(rho_F))
    #Executar codigo em C
    os.system("./mainOMP 4")#Executa o modelo C
    #Ler resultado em C
    outPut = 0
    with open("output.txt", 'r') as f:
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
r_m_mean = 60*24*3.96*10**-6
r_t_mean = 0.1
lamb_f_m_mean = 60*24*3.96*10**-6
b_d_mean = 0.001

gamma_D_mean = 0.01
gamma_F_mean = 0.03
gamma_T_mean = 0.2

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
b_rho_mean = 10**5
b_rho_b_mean = 3
b_rho_p_mean = 1.02
rho_T_mean = 2
rho_Tc_mean = 2
rho_B_mean = 11
rho_P_mean = 3
rho_F_mean = 5.1*10**-2

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
        'rho_F'
    ],
    'bounds': [
        [0.9*d_mic_mean, 1.1*d_mic_mean],
        [0.9*d_anti_mean, 1.1*d_anti_mean],
        [0.9*d_dc_mean, 1.1*d_dc_mean],
        [0.9*d_da_mean, 1.1*d_da_mean],
        [0.9*d_t_cit_mean, 1.1*d_t_cit_mean],
        [0.9*chi_mean, 1.1*chi_mean],

        [0.9*mu_dc_mean, 1.1*mu_dc_mean],
        [0.9*mu_m_mean, 1.1*mu_m_mean],
        [0.9*r_m_mean, 1.1*r_m_mean],
        [0.9*r_t_mean, 1.1*r_t_mean],
        [0.9*lamb_f_m_mean, 1.1*lamb_f_m_mean],
        [0.9*b_d_mean, 1.1*b_d_mean],

        [0.9*gamma_D_mean, 1.1*gamma_D_mean],
        [0.9*gamma_F_mean, 1.1*gamma_F_mean],
        [0.9*gamma_T_mean, 1.1*gamma_T_mean],

        [.9*cMic_mean, 1.1*cMic_mean],
        [.9*cCDc_mean, 1.1*cCDc_mean],
        [.9*cADc_mean, 1.1*cADc_mean],
        [.9*cDl_mean, 1.1*cDl_mean],
        [.9*cF_mean, 1.1*cF_mean],

        [0.9*alpha_T_h_mean, 1.1*alpha_T_h_mean],
        [0.9*alpha_T_c_mean, 1.1*alpha_T_c_mean],
        [0.9*alpha_B_mean, 1.1*alpha_B_mean],
        [0.9*alpha_P_mean, 1.1*alpha_P_mean],
        [0.9*b_T_mean, 1.1*b_T_mean],
        [0.9*b_Tc_mean, 1.1*b_Tc_mean],
        [0.9*b_rho_mean, 1.1*b_rho_mean],
        [0.9*b_rho_b_mean, 1.1*b_rho_b_mean],
        [0.9*b_rho_p_mean, 1.1*b_rho_p_mean],
        [0.9*rho_T_mean, 1.1*rho_T_mean],
        [0.9*rho_Tc_mean, 1.1*rho_Tc_mean],
        [0.9*rho_B_mean, 1.1*rho_B_mean],
        [0.9*rho_P_mean, 1.1*rho_P_mean],
        [0.9*rho_F_mean, 1.1*rho_F_mean]
    ]
}

problem_teste = {
    'num_vars': 3,
    'names': [
        'a',
        'b',
        'c'
    ],
    'bounds': [
        [-1,1],
        [-1,1],
        [-1,1]
    ]
}

print("Running Model")
sample = saltelli.sample(problem, 100, calc_second_order=False)
Y = np.empty([sample.shape[0]])

inputFile = open("sample.txt", "w")
for samp in sample:
    inputFile.write(str(samp) + "\n")
inputFile.close()

# evaluate the model for eah point in the input sample
for i in range(len(Y)):
    x = sample[i]
    Y[i] = modelo(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[2], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33])

output_file = np.zeros(640)
i = 0
with open("returns.txt", "r") as f:
    # output_file = f.readlines()
    for y in f.read().split("\n"):
        print(y)
        if i == 640:
            break
        if type(y) == str:
            output_file[i] = float(y)
        i = i + 1

# print(output_file)

# estimate the sensitivity indices using Sobol's method
sensitivity = sobol.analyze(problem_teste, output_file, calc_second_order=False)

# firstorder indices
print("First-order or main effect indices")
print(sensitivity['S1'])
# interpretation: x1 contributes to half of the 
# total output uncertainty

# higher-order indices
print("Higher-order or total (interactions) indices")
print(sensitivity['ST'])

print("Done!")