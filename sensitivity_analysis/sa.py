import numpy as np
import matplotlib.pyplot as plt
import os
from SALib.sample import saltelli
from SALib.analyze import sobol

def modelo(d_mic, d_anti, d_dc, d_da, d_t_cit, chi, mu_dc, mu_m, r_m, r_t, lamb_f_m, b_d, gamma_D, gamma_F, gamma_T,
cMic, cCDc, cADc, cDl, cF, alpha_T_h, alpha_T_c, alpha_B, alpha_P, b_Th, b_Tc, b_rho, b_rho_b, b_rho_p, rho_Th, rho_Tc, rho_B, rho_P, rho_F):
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

        [.9*cMic_mean, upperBound*cMic_mean],
        [.9*cCDc_mean, upperBound*cCDc_mean],
        [.9*cADc_mean, upperBound*cADc_mean],
        [.9*cDl_mean, upperBound*cDl_mean],
        [.9*cF_mean, upperBound*cF_mean],

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

def RunSA():
    print("Running Model")
    sample = saltelli.sample(problem, 4, calc_second_order=False)
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
    return sensitivity

def ReadSensitivity():
    sensitivity = {
        'S1': [8.10470350e-02, 1.30161888e-06, -1.87588328e-08, -7.98065622e-09,
    1.19157601e-04, 2.15319320e-02, -1.79993643e-07, 8.87319458e-01,
    2.32798953e-05, -1.89251553e-05, -7.30988451e-06, -5.04483338e-07,
    -1.36386018e-08, -4.24069663e-07, -3.79414154e-05, 5.71591451e-03,
    -3.35283511e-08, -5.36640357e-07, -3.59938269e-07 -1.39886999e-08,
    -9.11781785e-09, 6.00908583e-05, -3.75227026e-08, -3.01878232e-06,
    -7.21745691e-08, 1.92672821e-08, -5.38780976e-07, -1.92658619e-08,
    7.56544201e-07, -2.69817677e-07, -5.35171487e-08, 3.66348520e-08,
    -2.89514282e-07, -1.81924591e-06],
        'ST': [8.43824923e-02, 4.91991490e-09, 1.98017551e-09, 2.44353816e-10,
    2.93806199e-06, 2.18838369e-02, 2.86261721e-10, 8.83197561e-01,
    2.24520567e-05, 1.47011920e-05, 3.36275160e-08, 2.20992648e-09,
    2.61491383e-13, 4.33716384e-10, 6.67356909e-06, 3.61700152e-03,
    6.94283658e-13, 2.11023394e-10, 1.09136294e-10, 3.04907957e-13,
    6.27420006e-09, 4.96628909e-05, 1.24207262e-11, 7.78751342e-09,
    1.04881928e-11, 4.06162426e-13, 1.40977732e-08, 4.68120699e-12,
    5.42575796e-09, 3.96713026e-11, 4.99328746e-13, 2.08652193e-09,
    5.66606580e-09, 5.10866780e-09]
    }
    return sensitivity

sensitivity = RunSA()
# sensitivity = ReadSensitivity()

def printResult(indexes):
    labels = ["d_mic",
    "d_anti",
    "d_dc",
    "d_da",
    "d_t_cit",
    "chi",
    "mu_dc",
    "mu_m",
    "r_m",
    "r_t",
    "lamb_f_m",
    "b_d",
    "gamma_D",
    "gamma_F",
    "gamma_T",
    "cMic",
    "cCDc",
    "cADc",
    "cDl",
    "cF",
    "alpha_T_h",
    "alpha_T_c",
    "alpha_B",
    "alpha_P",
    "b_T",
    "b_Tc",
    "b_rho",
    "b_rho_b",
    "b_rho_p",
    "rho_T",
    "rho_Tc",
    "rho_B",
    "rho_F"]
    # plt.bar(labels, indexes, color ='maroon')
    plt.barh(labels, indexes, color ='maroon')
    plt.show()
    # plt.savefig("SAOaT-Results.png", dpi = 900)
    plt.clf()

printResult(sensitivity['S1'])
printResult(sensitivity['ST'])