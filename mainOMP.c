#include "modelOMP.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

structParameters ReadParameters(){
    FILE *file;
    file = fopen("./sensitivity_analysis/SA_parameters.txt","r");
    char lineRead[25];
    structParameters params;
    float fileParameters[34];
    int fileIter = 0;
    while (fgets(lineRead, sizeof(lineRead), file) != NULL){
        // char* valParam = strtok(lineRead, "\n");
        fileParameters[fileIter] = atof(lineRead);
        fileIter++;
    }
    fclose(file);
    //Escrever o vetor de parametros na estrutura parametros
    params.chi = fileParameters[0];
    params.micDiffusion = fileParameters[1];
    params.cDcDiffusion = fileParameters[2];
    params.aDcDiffusion = fileParameters[3];
    params.tCytoDiffusion = fileParameters[4];
    params.antibodyDiffusion = fileParameters[5];
    
    params.muMic = fileParameters[6];
    params.rM = fileParameters[7];
    params.lambAntMic = fileParameters[8];
    params.bD = fileParameters[9];
    params.rT = fileParameters[10];
    params.muCDc = fileParameters[11];
    
    params.gammaD = fileParameters[12];
    params.gammaAntibody = fileParameters[13];
    params.gammaT = fileParameters[14];

    params.alphaTHelper = fileParameters[15];
    params.alphaTCytotoxic = fileParameters[16];
    params.alphaB = fileParameters[17];
    params.alphaP = fileParameters[18];
    params.cMic = fileParameters[19];
    params.cCDc = fileParameters[20];
    params.cADc = fileParameters[21];
    params.cDl = fileParameters[22];
    params.cF = fileParameters[23];
    params.bTHelper = fileParameters[24];
    params.bTCytotoxic = fileParameters[25];
    params.bRho = fileParameters[26];
    params.bRhoB = fileParameters[27];
    params.bRhoP = fileParameters[28];
    params.rhoTHelper = fileParameters[29];
    params.rhoTCytotoxic = fileParameters[30];
    params.rhoB = fileParameters[31];
    params.rhoP = fileParameters[32];
    params.rhoAntibody = fileParameters[33];


    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;
    params.estableTHelper = 84;
    params.estableTCytotoxic = 40;
    params.estableB = 25;
    params.estableP = 2.5;
    params.V_LN = 160;
    params.V_BV = 0;
    params.V_PV = 0;
    return params;
}

structParameters ParametersInitialize(){
    structParameters params;
    params.micDiffusion = 9.6*24*6.6*pow(10,-5);
    params.antibodyDiffusion = 9.6*24*6.6*pow(10,-4);
    params.cDcDiffusion = 9.6*24*6.6*pow(10,-5);
    params.aDcDiffusion = 9.6*24*6.6*pow(10,-5);
    params.tCytoDiffusion = 9.6*24*6.6*pow(10,-5);
    params.chi = 0.298*60*2;
    
    params.muCDc = 60*24*3*pow(10,-4);
    params.muMic = 60*24*3*pow(10,-6);
    params.rM = 60*24*3.96*pow(10,-6);
    params.rT = 0.1;
    params.lambAntMic = 5.702*pow(10,-3);
    params.bD = 0.001;
    
    params.gammaD = 0.01;
    params.gammaAntibody = 0.3;
    params.gammaT = 2;

    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;

    params.cMic = 0.1;
    params.cCDc = 0.1;
    params.cADc = 0.1;
    params.cDl = 0.1;
    params.cF = 0.1;

    params.alphaTHelper = 0.1;
    params.alphaTCytotoxic = 0.1;
    params.alphaB = 0.1;
    params.alphaP = 1;
    
    params.bTHelper = 0.17;
    params.bTCytotoxic = 0.001;
    params.bRho = 0.6;
    params.bRhoB = 3.02;
    params.bRhoP = 1.02;

    params.rhoTHelper = 2;
    params.rhoTCytotoxic = 2;
    params.rhoB = 11;
    params.rhoP = 3;
    params.rhoAntibody = 5.1*pow(10,-2);

    params.estableTHelper = 84;
    params.estableTCytotoxic = 40;
    params.estableB = 25;
    params.estableP = 2.5;
    params.V_LN = 160;
    params.V_BV = 0;
    params.V_PV = 0;

    return params;
}

void clearPhgTxt(){
    system("find ./result/ -name '*.png' -type f -delete");
    system("find ./result/ -name '*.txt' -type f -delete");
}

int main(int argc, char* argv[]){
    clearPhgTxt();
    int tot_thr = strtol(argv[1], NULL, 10);
    int calculateQoI = 0;
    calculateQoI = strtol(argv[2], NULL, 10);
    structParameters parameters;
    if(calculateQoI==0)
        parameters = ParametersInitialize();
    else
        parameters = ReadParameters();
    float ht = 0.000005, hx = 0.2;
    int numFigs = 7, numPointsLN = 1000, time = 14, space = 20;
    structModel model = ModelInitialize(parameters, tot_thr, ht, hx, time, space, numFigs, numPointsLN, calculateQoI);
    float start = omp_get_wtime();
    RunModel(&model);
    float end = omp_get_wtime();
    printf("Work took %f seconds\n", end - start);
    return 0;
}