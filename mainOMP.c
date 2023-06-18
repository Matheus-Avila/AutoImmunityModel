#include "modelOMP.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

structParameters ParametersInitialize(){
    structParameters params;
    params.micDiffusion = 24*6.6*pow(10,-5);
    params.antibodyDiffusion = 24*6.6*pow(10,-4);
    params.cDcDiffusion = 24*6.6*pow(10,-5);
    params.aDcDiffusion = 24*6.6*pow(10,-5);
    params.tCytoDiffusion = 24*6.6*pow(10,-5);
    params.chi = 0.003;
    
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

void WriteTime(float ExecTime){
    FILE *fileTime;
    fileTime = fopen("./ExecsTimes.txt", "a");
    if(fileTime != NULL){
    fprintf(fileTime, "%f\n", ExecTime);
    fclose(fileTime);
    }else{
        printf("Error execution time file\n");
        exit(0);
    }
}

void clearPhgTxt(){
    system("find ./result/ -name '*.png' -type f -delete");
    system("find ./result/ -name '*.txt' -type f -delete");
    system("mkdir result");
    system("mkdir result/matrix");
    system("mkdir result/odc");
    system("mkdir result/mic");
    system("mkdir result/tke");
    system("mkdir result/ant");
    system("mkdir result/da");
    system("mkdir result/dc");
}

int main(int argc, char* argv[]){
    printf("Comecei o main\n");
    clearPhgTxt();
    int tot_thr = strtol(argv[1], NULL, 10);
    float ht = 0.0002, hx = 0.5;
    int numFigs = 28, numPointsLN = 1000, time = 28, space = 20, numStepsLN = 1, saveFigs = 1;
    structParameters parameters = ParametersInitialize();
    structModel model = ModelInitialize(parameters, tot_thr, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);
    printf("Inicializacao feita!!\n\n");
    float start = omp_get_wtime();
    RunModel(&model);
    float end = omp_get_wtime();
    WriteTime(end - start);
    printf("Work took %f seconds\n", end - start);
    return 0;
}