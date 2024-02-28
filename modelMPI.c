#include "modelMPI.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "time.h"

#define DEBUG 1

int xSize = 0;
int tSize = 0;

float* ComposeMatrixToSend(float* oligodendrocyte, float* tCytotoxic, float* microglia, float* antibody, float* conventionalDc, 
float* activatedDc, int numLines, int startLine, int xSize
) {
    float* composed = calloc(xSize * numLines * 6, sizeof(float));
    for(int i = 0; i < numLines; i++) {
        for(int j = 0; j < xSize; j++) {
            composed[(i * xSize) + j] = oligodendrocyte[(startLine + i) * xSize + j];
            composed[(i * xSize) + j + numLines * 1] = tCytotoxic[(startLine + i) * xSize + j];
            composed[(i * xSize) + j + numLines * 2] = microglia[(startLine + i) * xSize + j];
            composed[(i * xSize) + j + numLines * 3] = antibody[(startLine + i) * xSize + j];
            composed[(i * xSize) + j + numLines * 4] = conventionalDc[(startLine + i) * xSize + j];
            composed[(i * xSize) + j + numLines * 5] = activatedDc[(startLine + i) * xSize + j];
        }
    }
    return composed;
}

void DecomposeRecvMatrix(float*composed, float* oligodendrocyte, float* tCytotoxic, float* microglia, float* antibody, float* conventionalDc, 
float* activatedDc, int numLines, int startLine, int xSize
) {
    for(int i = 0; i < numLines; i++) {
        for(int j = 0; j < xSize; j++) {
            oligodendrocyte[(startLine + i) * xSize + j] = composed[(i * xSize) + j];
            tCytotoxic[(startLine + i) * xSize + j] = composed[(i * xSize) + j + numLines * 1];
            microglia[(startLine + i) * xSize + j] = composed[(i * xSize) + j + numLines * 2];
            antibody[(startLine + i) * xSize + j] = composed[(i * xSize) + j + numLines * 3];
            conventionalDc[(startLine + i) * xSize + j] = composed[(i * xSize) + j + numLines * 4];
            activatedDc[(startLine + i) * xSize + j] = composed[(i * xSize) + j + numLines * 5];
        }
    }
}

void InitialConditionTissueMicroglia(structModel* model){
    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        if(pow((i-(int)(model->xSize/2)),2) + pow((j-(int)(model->xSize/2)),2) < 5 / (model->hx * model->hx)){
            model->microglia[0][k] = (float)model->parametersModel.avgMic/3;
        }
    }
}

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN){
    model->dendriticLymphNode[0] = dendriticLN;
    model->tHelperLymphNode[0] = thelperLN;
    model->tCytotoxicLymphNode[0] = tcytotoxicLN;
    model->bCellLymphNode[0] = bcellLN;
    model->plasmaCellLymphNode[0] = plasmacellLN;
    model->antibodyLymphNode[0] = antibodyLN;
}

int VerifyCFL(structParameters parametersModel, float ht, float hx){
    if(parametersModel.micDiffusion*ht/(hx*hx) < 0.25 && parametersModel.cDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.aDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.tCytoDiffusion*ht/(hx*hx) < 0.25 && parametersModel.chi*ht/hx < 1)
        return 1;
    return 0;
}

void WritePopulation(structModel model, float *population, char* fileName, char* bufferTime){
    FILE *file;
    file = fopen(fileName, "w");
    int k = 0;
    while (k < model.xSize*model.xSize){
        int i = k;
        while (i < k + model.xSize){
            fprintf(file, "%f ", population[i]);
            i++;
        }
        fprintf(file,"\n");
        k+=model.xSize;
    }
    fclose(file);
}

void WritePopulationLymphNode(structModel model, float *population, char* fileName){
    FILE *file;
    file = fopen(fileName, "w");
    for(int i=0;i<model.numPointsLN;i++){
        fprintf(file, "%f\n", population[i]);
    }
    fclose(file);
}

void WriteLymphNodeFiles(structModel model, float *dendritic, float *tHelper, float *tCytotoxic, float *bCell, float *plasmaCell, float *antibody){
    clock_t start = clock();
    //printf("Start to write Files\n");
    WritePopulationLymphNode(model, dendritic, "./result/dendritic.txt");
    WritePopulationLymphNode(model, tHelper, "./result/tHelper.txt");
    WritePopulationLymphNode(model, tCytotoxic, "./result/tCyto.txt");
    WritePopulationLymphNode(model, bCell, "./result/bCell.txt");
    WritePopulationLymphNode(model, plasmaCell, "./result/plasmaCell.txt");
    WritePopulationLymphNode(model, antibody, "./result/antibody.txt");
    char buffer[20];
    char command[40] = {};
    strcat(command, "python3 plotLymphNode.py ");
    snprintf(buffer, sizeof(buffer), "%d", model.tFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", (model.tSize/model.numPointsLN)*model.ht);
    strcat(command, buffer);
    system(command);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    //printf("Write files took %lf seconds\n", time_spent);
}

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float  *activatedDC, float time){
    char buffer[10];
    float day = time * model.ht;
    
    snprintf(buffer, sizeof(buffer), "%1.f", day);
    
    char pathOligodendrocytes[200];
    getcwd(pathOligodendrocytes, sizeof(pathOligodendrocytes));
    strcat(pathOligodendrocytes, "/result/matrix/oligo");
    strcat(pathOligodendrocytes, buffer);
    strcat(pathOligodendrocytes, ".txt");
    
    WritePopulation(model, oligodendrocyte, pathOligodendrocytes, buffer);

    char pathMicroglia[200];
    getcwd(pathMicroglia, sizeof(pathMicroglia));
    strcat(pathMicroglia, "/result/matrix/microglia");
    strcat(pathMicroglia, buffer);
    strcat(pathMicroglia, ".txt");
    WritePopulation(model, microglia, pathMicroglia, buffer);

    char pathTCyto[200];
    getcwd(pathTCyto, sizeof(pathTCyto));
    strcat(pathTCyto, "/result/matrix/tCyto");
    strcat(pathTCyto, buffer);
    strcat(pathTCyto, ".txt");
    WritePopulation(model, tCytotoxic, pathTCyto, buffer);

    char pathAntibody[200];
    getcwd(pathAntibody, sizeof(pathAntibody));
    strcat(pathAntibody, "/result/matrix/antibody");
    strcat(pathAntibody, buffer);
    strcat(pathAntibody, ".txt");
    WritePopulation(model, antibody, pathAntibody, buffer);

    char pathConventionalDC[200];
    getcwd(pathConventionalDC, sizeof(pathConventionalDC));
    strcat(pathConventionalDC, "/result/matrix/conventionalDC");
    strcat(pathConventionalDC, buffer);
    strcat(pathConventionalDC, ".txt");
    WritePopulation(model, conventionalDC, pathConventionalDC, buffer);

    char pathActivatedDC[200];
    getcwd(pathActivatedDC, sizeof(pathActivatedDC));
    strcat(pathActivatedDC, "/result/matrix/activatedDC");
    strcat(pathActivatedDC, buffer);
    strcat(pathActivatedDC, ".txt");
    WritePopulation(model, activatedDC, pathActivatedDC, buffer);
}   

void PlotResults(structModel model){
    clock_t start = clock();
    //printf("Ploting results\n");
    char buffer[10];
    char command[200] = {};
    strcat(command, "python3 plotMatrices.py ");
    snprintf(buffer, sizeof(buffer), "%d", model.xFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", model.hx);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", model.tFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", model.intervaloFiguras);
    strcat(command, buffer);
    system(command);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    //printf("Plot results took %lf seconds\n", time_spent);
}

float PreventionOverCrowdingTerm(float populationPoint, float avgValue){
    return populationPoint/(populationPoint + avgValue);
}

float UpDownWind(float frontPoint, float rearPoint, float avgValue){
    return PreventionOverCrowdingTerm(frontPoint, avgValue) - PreventionOverCrowdingTerm(rearPoint, avgValue);
}

float CalculateChemottaxis(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ, float hx){
    float gradientPopulationI, gradientPopulationJ;
    if(gradientOdcI<0)
        gradientPopulationI = UpDownWind(frontIPoint, ijPoint, avgValue)/(float)hx;
    else
        gradientPopulationI = UpDownWind(ijPoint, rearIPoint, avgValue)/(float)hx;
    if(gradientOdcJ<0)
        gradientPopulationJ = UpDownWind(frontJPoint, ijPoint, avgValue)/(float)hx;
    else
        gradientPopulationJ = UpDownWind(ijPoint, rearJPoint, avgValue)/(float)hx;

    return gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

float CalculateDiffusion(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float hx){
    return (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(hx*hx);
}

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
}

void WriteBVPV(structModel *model, float *thetaBV, float *thetaPV){
    FILE *fileBV;
    fileBV = fopen("/result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("/result/pv.txt", "w");
    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        fprintf(fileBV, "%f ", thetaBV[k]);
        fprintf(filePV, "%f ", thetaPV[k]);    
        if(k%model->xSize == 0 && k != 0){
            fprintf(fileBV,"\n");
            fprintf(filePV,"\n");
        }
    }
    fclose(fileBV);
    fclose(filePV);
    char buffer[10];
    char command[70] = {};
    strcat(command, "python3 plotBVPV.py ");
    snprintf(buffer, sizeof(buffer), "%d", model->xFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", model->hx);
    strcat(command, buffer);
    // system(command);
}

void DefineBVPV(structModel *model){
    if(model->my_rank == 0){
        int randomVal;
        for(int i = 0; i < xSize; i++){
            for(int j = 0; j < xSize; j++){
                randomVal = rand() % 100;
                if(randomVal <10){
                    model->parametersModel.V_BV++;
                    model->parametersModel.V_PV++;
                    model->thetaBV[(i * xSize) + j] = 1;
                    if(j != xSize-1)
                        model->thetaPV[(i * xSize) + (j + 1)] = 1;
                    else
                        model->thetaPV[(i * xSize)] = 1;
                }
            }
        }
        model->parametersModel.V_BV = model->parametersModel.V_BV * model->hx * model->hx;
        model->parametersModel.V_PV = model->parametersModel.V_PV * model->hx * model->hx;
        printf("bv = %f, pv = %f \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    }
    // Para cada linha da matriz fazer broadcast
    MPI_Bcast(model->thetaPV, xSize * xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(model->thetaBV, xSize * xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int tempBV = model->parametersModel.V_BV, tempPV = model->parametersModel.V_PV;
    MPI_Bcast(&tempBV, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tempPV, 1, MPI_INT, 0, MPI_COMM_WORLD);
    model->parametersModel.V_BV = tempBV;
    model->parametersModel.V_PV = tempPV;
}

structModel ModelInitialize(structParameters params, float ht, 
float hx, float time, float space, int numFigs, int numPointsLN, int my_rank, int comm_sz, int numStepsLN, int saveFigs){
    structModel model;
    srand(2);
    model.parametersModel = params;
    model.numFigs = numFigs;
    model.numPointsLN = numPointsLN;
    model.ht = ht;
    model.hx = hx;
    model.tFinal = time;
    model.xFinal = space;
    model.tSize = (int)(time/ht);
    model.timeLen = (int)(time/ht);
    model.spaceLen = (int)(space/hx);
    model.xSize = (int)(space/hx);
    model.intervaloFiguras = (int)model.tSize/numFigs;
    xSize = model.xSize;
    tSize = model.tSize;

    model.my_rank = my_rank;
    model.comm_sz = comm_sz;
    
    model.numStepsLN = numStepsLN;
    model.saveFigs = saveFigs;
    
    model.numLines = model.spaceLen/comm_sz;
    model.startLine = my_rank*model.numLines;
    model.endLine = model.startLine + model.numLines-1;
    model.microglia = (float**)malloc(BUFFER * sizeof(float*));
    model.oligodendrocyte = (float**)malloc(BUFFER * sizeof(float*));
    model.tCytotoxic = (float**)malloc(BUFFER * sizeof(float*));
    model.conventionalDc = (float**)malloc(BUFFER * sizeof(float*));
    model.activatedDc = (float**)malloc(BUFFER * sizeof(float*));
    model.antibody = (float**)malloc(BUFFER * sizeof(float*));
    model.startLine = my_rank*model.numLines;
    model.endLine = model.startLine + model.numLines-1;
    //printf("\nProcess %d num line: %d and start line: %d\n", my_rank, model.numLines, model.startLine);
    if(my_rank == comm_sz - 1)
        model.endLine = xSize - 1;
    for (int index=0;index<BUFFER;++index){
        model.microglia[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.oligodendrocyte[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.tCytotoxic[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.conventionalDc[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.activatedDc[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.antibody[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
    }
    //definir BV e PV
    model.thetaBV = (float*)calloc(model.xSize * model.xSize, sizeof(float));
    model.thetaPV = (float*)calloc(model.xSize * model.xSize, sizeof(float));
    DefineBVPV(&model);
    //definir lymph node
    
    model.dendriticLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.tCytotoxicLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.tHelperLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.antibodyLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.bCellLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.plasmaCellLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));

    float dendriticLN = 0.0, thelperLN = params.stableTHelper, tcytotoxicLN = params.stableTCytotoxic, bcellLN = params.stableB, plasmacellLN = 0.0, antibodyLN = 0.0;

    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    InitialConditionTissueMicroglia(&model);
    return model;
}

/*
* Lymphnode
*/
float* EquationsLymphNode(structModel model, float* populationLN, int stepPos){
    float* result =(float *)malloc(sizeof(float)*6);
    
    float dcLN = populationLN[0];
    float tCytoLN = populationLN[1];
    float tHelperLN = populationLN[2];
    float bCellLN = populationLN[3];
    float plasmaCellLN = populationLN[4];
    float antibodyLN = populationLN[5];

    //Describe equations

    //Dendritic cell
    float activatedDcMigration = model.parametersModel.gammaD * (model.activatedDCTissueVessels - dcLN) * (float)(model.parametersModel.V_PV/model.parametersModel.V_LN);
    float activatedDcClearance = model.parametersModel.cDl * dcLN;
    result[0] = activatedDcMigration - activatedDcClearance;

    //T Cytotoxic
    float tCytoActivation = model.parametersModel.bTCytotoxic * (model.parametersModel.rhoTCytotoxic*tCytoLN*dcLN - tCytoLN*dcLN);
    float tCytoHomeostasis = model.parametersModel.alphaTCytotoxic * (model.parametersModel.stableTCytotoxic - tCytoLN);
    float tCytoMigration = model.parametersModel.gammaT * (tCytoLN - model.tCytotoxicTissueVessels) * (float)(model.parametersModel.V_BV/model.parametersModel.V_LN);
    result[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

    //T Helper
    float tHelperActivation = model.parametersModel.bTHelper * (model.parametersModel.rhoTHelper * tHelperLN * dcLN - tHelperLN * dcLN);
    float tHelperHomeostasis = model.parametersModel.alphaTHelper * (model.parametersModel.stableTHelper - tHelperLN);
    float tHelperDispendure = model.parametersModel.bRho * dcLN * tHelperLN * bCellLN;
    result[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

    //B Cell
    float bCellActivation = model.parametersModel.bRhoB * (model.parametersModel.rhoB * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
    float bcellHomeostasis = model.parametersModel.alphaB * (model.parametersModel.stableB - bCellLN);
    result[3] = bcellHomeostasis + bCellActivation;

    //Plasma Cells
    float plasmaActivation = model.parametersModel.bRhoP * (model.parametersModel.rhoP * tHelperLN * dcLN * bCellLN);
    float plasmaHomeostasis = model.parametersModel.alphaP * (model.parametersModel.stableP - plasmaCellLN);
    result[4] = plasmaHomeostasis + plasmaActivation;

    //Antibody
    float antibodyProduction = model.parametersModel.rhoAntibody * plasmaCellLN;
    float antibodyDecayment = model.parametersModel.cF * antibodyLN;
    float antibodyMigration = model.parametersModel.gammaAntibody * (antibodyLN - model.antibodyTissueVessels) * (float)(model.parametersModel.V_BV/model.parametersModel.V_LN);
    result[5] = antibodyProduction - antibodyMigration - antibodyDecayment;

    return result;
}

void SolverLymphNode(structModel *model, int stepPos){
    float populationLN[6];
    int stepKPlus = (stepPos%(2*model->numStepsLN))/model->numStepsLN;
    int stepKMinus = !(stepKPlus && 1);
    populationLN[0] = model->dendriticLymphNode[stepKMinus];
    populationLN[1] = model->tCytotoxicLymphNode[stepKMinus];
    populationLN[2] = model->tHelperLymphNode[stepKMinus];
    populationLN[3] = model->bCellLymphNode[stepKMinus];
    populationLN[4] = model->plasmaCellLymphNode[stepKMinus];
    populationLN[5] = model->antibodyLymphNode[stepKMinus];
    
    float* solutionLN;
    solutionLN = EquationsLymphNode(*model, populationLN, stepPos);

    float htLN = model->ht*model->numStepsLN;
    
    //Execute Euler 
    model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + htLN*solutionLN[0];
    model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + htLN*solutionLN[1];
    model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + htLN*solutionLN[2];
    model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + htLN*solutionLN[3];
    model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + htLN*solutionLN[4];
    model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + htLN*solutionLN[5];
    free(solutionLN);

    int intervalPoints = (int)(model->tSize/model->numPointsLN);
    if(stepPos%intervalPoints){
        int posSave = stepPos/intervalPoints;
        model->dendriticLymphNodeSavedPoints[posSave] = model->dendriticLymphNode[stepKPlus];
        model->tCytotoxicLymphNodeSavedPoints[posSave] = model->tCytotoxicLymphNode[stepKPlus];
        model->tHelperLymphNodeSavedPoints[posSave] = model->tHelperLymphNode[stepKPlus];
        model->bCellLymphNodeSavedPoints[posSave] = model->bCellLymphNode[stepKPlus];
        model->plasmaCellLymphNodeSavedPoints[posSave] = model->plasmaCellLymphNode[stepKPlus];
        model->antibodyLymphNodeSavedPoints[posSave] = model->antibodyLymphNode[stepKPlus];
    }
}

void SavingData(structModel model){
    float totalMic = 0, totalODC = 0, totalCDC = 0, totalADC = 0, totalIGG = 0, totalCD8 = 0;
    for(int kPos = 0; kPos < model.xSize*model.xSize; kPos++){
        totalMic += model.microglia[0][kPos];
        totalODC += model.oligodendrocyte[0][kPos];
        totalCDC += model.conventionalDc[0][kPos];
        totalADC += model.activatedDc[0][kPos];
        totalCD8 += model.tCytotoxic[0][kPos];
        totalIGG += model.antibody[0][kPos];
    }
    FILE *file;
    file = fopen("dataExecution.txt", "w");
    if(file != NULL){
        fprintf(file, "Days = %d - Space = %d - ht = %f, hx = %f, Ht_JumpStep = %d\n", model.tFinal, model.xFinal, model.ht, model.hx, model.numStepsLN);
        fprintf(file, "Lymph node populations\n");
        fprintf(file, "DC = %f, TCD8 = %f, TCD4 = %f, B Cell = %f, Plasma cell = %f, IgG = %f\n", model.dendriticLymphNodeSavedPoints[model.numPointsLN-1],\
         model.tCytotoxicLymphNodeSavedPoints[model.numPointsLN-1], model.tHelperLymphNodeSavedPoints[model.numPointsLN-1], model.bCellLymphNodeSavedPoints[model.numPointsLN-1], model.plasmaCellLymphNodeSavedPoints[model.numPointsLN-1], model.antibodyLymphNodeSavedPoints[model.numPointsLN-1]);
        fprintf(file, "Tissue populations\n");
        fprintf(file, "ODC = %f, Microglia = %f, ConventionalDC = %f, ActivatedDC = %f, TCD8 = %f, IgG = %f\n", totalODC, totalMic, totalCDC, totalADC, totalCD8, totalIGG);    
        fprintf(file, "Parameters\n");
        fprintf(file, "micDiffusion  = %f, antibodyDiffusion = %f, cDcDiffusion = %f, aDcDiffusion = %f, tCytoDiffusion = %f, chi = %f, muCDc = %f, muMic = %f, \
        rM = %f, rT = %f, lambAntMic = %f, bD = %f, gammaD = %f, gammaAntibody = %f, gammaT = %f,  avgT = %f, avgDc = %f, avgMic = %f, avgOdc = %f,  cMic = %f, \
        cCDc = %f, cADc = %f, cDl = %f, cF = %f, alphaTHelper = %f, alphaTCytotoxic = %f, alphaB = %f, alphaP = %f, bTHelper = %f, bTCytotoxic = %f, bRho = %f, \
        bRhoB = %f, bRhoP = %f, rhoTHelper = %f, rhoTCytotoxic = %f, rhoB = %f, rhoP = %f, rhoAntibody = %f, stableTHelper = %f, stableTCytotoxic = %f, \
        stableB = %f, stableP = %f, V_LN = %d, V_BV = %f, V_PV = %f\n",
        model.parametersModel.micDiffusion, model.parametersModel.antibodyDiffusion, model.parametersModel.cDcDiffusion, model.parametersModel.aDcDiffusion, \
        model.parametersModel.tCytoDiffusion, model.parametersModel.chi, model.parametersModel.muCDc, model.parametersModel.muMic, model.parametersModel.rM, \
        model.parametersModel.rT, model.parametersModel.lambAntMic, model.parametersModel.bD, model.parametersModel.gammaD, model.parametersModel.gammaAntibody, \
        model.parametersModel.gammaT,  model.parametersModel.avgT, model.parametersModel.avgDc, model.parametersModel.avgMic, model.parametersModel.avgOdc, \
        model.parametersModel.cMic, model.parametersModel.cCDc, model.parametersModel.cADc, model.parametersModel.cDl, model.parametersModel.cF, \
        model.parametersModel.alphaTHelper, model.parametersModel.alphaTCytotoxic, model.parametersModel.alphaB, model.parametersModel.alphaP, \
        model.parametersModel.bTHelper, model.parametersModel.bTCytotoxic, model.parametersModel.bRho, model.parametersModel.bRhoB, model.parametersModel.bRhoP,\
        model.parametersModel.rhoTHelper, model.parametersModel.rhoTCytotoxic, model.parametersModel.rhoB, model.parametersModel.rhoP,\
        model.parametersModel.rhoAntibody, model.parametersModel.stableTHelper, model.parametersModel.stableTCytotoxic, model.parametersModel.stableB,\
        model.parametersModel.stableP, model.parametersModel.V_LN, model.parametersModel.V_BV, model.parametersModel.V_PV);
        fclose(file);
    }else{
        printf("dataExecution file not found!\n");
        exit(0);
    }
}

void SendBorderThread(structModel *model, int stepKPlus){
    // printf("%d/%d  ::  %d, %d\n",model->my_rank,model->comm_sz, model->startLine, model->endLine);
    if(model->my_rank%2 == 0){
        if(model->my_rank != model->comm_sz-1){
            MPI_Send(&model->oligodendrocyte[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(model->my_rank != 0){
            MPI_Send(&model->oligodendrocyte[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][model->startLine  * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else{
        MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(&model->oligodendrocyte[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
        MPI_Send(&model->tCytotoxic[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
        MPI_Send(&model->microglia[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
        MPI_Send(&model->antibody[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
        MPI_Send(&model->conventionalDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
        MPI_Send(&model->activatedDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
        if(model->my_rank != model->comm_sz-1){                
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(&model->oligodendrocyte[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
            /*MPI_Sendrecv(&model->oligodendrocyte[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, 
            &model->oligodendrocyte[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&model->tCytotoxic[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, 
            &model->tCytotoxic[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&model->microglia[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, 
            &model->microglia[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&model->antibody[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, 
            &model->antibody[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&model->conventionalDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, 
            &model->conventionalDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&model->activatedDc[stepKPlus][(model->endLine+1) * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, 
            &model->activatedDc[stepKPlus][model->endLine * model->xSize], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);*/
        }
        //Bloco de código que sempre será executado
        // if(model->my_rank != 0){                
        //     MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Recv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Recv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Recv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Recv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Recv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //     MPI_Send(&model->oligodendrocyte[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
        //     MPI_Send(&model->tCytotoxic[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
        //     MPI_Send(&model->microglia[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
        //     MPI_Send(&model->antibody[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
        //     MPI_Send(&model->conventionalDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
        //     MPI_Send(&model->activatedDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
        //     MPI_Sendrecv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, 
        //     &model->oligodendrocyte[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, 
        //     &model->tCytotoxic[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, 
        //     &model->microglia[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, 
        //     &model->antibody[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, 
        //     &model->conventionalDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, 
        //     &model->activatedDc[stepKPlus][model->startLine * model->xSize], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // }
    }
}

void RunModel(structModel *model){
    
    if(model->saveFigs && model->my_rank == 0)
        WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
    
    double start = MPI_Wtime();
    //Save IC
    int stepKMinus = 0, stepKPlus, line, column;
    float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float diffusionOdc;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;
    //printf("Process %d spaceLen: %d\n\n", model->my_rank, model->spaceLen);

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0, auxTT = 0.0, auxAT = 0.0, auxDC = 0.0;
    for(int kTime = 1; kTime <= model->timeLen; kTime++){
        if(kTime%model->numStepsLN == 0){
            MPI_Allreduce(&auxTCytotoxicBV, &auxTT, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&auxAntibodyBV, &auxAT, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&auxAdcPV, &auxDC, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

            model->tCytotoxicTissueVessels = auxTT * model->hx * model->hx / model->parametersModel.V_BV;
            model->antibodyTissueVessels = auxAT * model->hx * model->hx / model->parametersModel.V_BV;
            model->activatedDCTissueVessels = auxDC * model->hx * model->hx / model->parametersModel.V_PV;

            // solve lymphnode
            SolverLymphNode(model, kTime);
        }
        auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0, auxTT = 0.0, auxAT = 0.0, auxDC = 0.0;
        
        stepKPlus = kTime%2;
        stepKMinus = !(stepKPlus && 1);
        for(line = model->startLine; line <= model->endLine; line++){
        for(column = 0; column < model->spaceLen; column++){
    
            microgliaKMinus = model->microglia[stepKMinus][line * model->xSize + column];
            conventionalDcKMinus = model->conventionalDc[stepKMinus][line * model->xSize + column];
            activatedDcKMinus = model->activatedDc[stepKMinus][line * model->xSize + column];
            tCytotoxicKMinus = model->tCytotoxic[stepKMinus][line * model->xSize + column];
            antibodyKMinus = model->antibody[stepKMinus][line * model->xSize + column];
            oligodendrocyteKMinus = model->oligodendrocyte[stepKMinus][line * model->xSize + column];
            
            //Define gradient ODCs
            valIPlus = (line != model->xSize-1)? model->oligodendrocyte[stepKMinus][(line+1) * model->xSize + column]: model->oligodendrocyte[stepKMinus][line * model->xSize + column];
            valJPlus = (column != model->xSize-1)? model->oligodendrocyte[stepKMinus][line * model->xSize + column+1]: model->oligodendrocyte[stepKMinus][line * model->xSize + column];
            valIMinus = (line != 0)? model->oligodendrocyte[stepKMinus][(line-1) * model->xSize + column]: model->oligodendrocyte[stepKMinus][line * model->xSize + column];
            valJMinus = (column != 0)? model->oligodendrocyte[stepKMinus][line * model->xSize + column-1]: model->oligodendrocyte[stepKMinus][line * model->xSize + column];
            gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*model->hx);
            gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*model->hx);

	        diffusionOdc = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->oligodendrocyte[stepKMinus][line * model->xSize + column], model->hx);

            //Diffusion and Chemotaxis Mic

            valIPlus  = (line != model->xSize-1)? model->microglia[stepKMinus][(line+1) * model->xSize + column]: model->microglia[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->microglia[stepKMinus][line * model->xSize + column+1]: model->microglia[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->microglia[stepKMinus][(line-1) * model->xSize + column]: model->microglia[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->microglia[stepKMinus][line * model->xSize + column-1]: model->microglia[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);
            
            microgliaDiffusion = model->parametersModel.micDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line * model->xSize + column], model->hx);
            microgliaChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line * model->xSize + column],\
            model->parametersModel.avgMic, gradientOdcI, gradientOdcJ, model->hx) + model->parametersModel.chi*(diffusionOdc * PreventionOverCrowdingTerm(microgliaKMinus, model->parametersModel.avgMic));

            //Diffusion and Chemotaxis CDC

            valIPlus  = (line != model->xSize-1)? model->conventionalDc[stepKMinus][(line+1) * model->xSize + column]: model->conventionalDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->conventionalDc[stepKMinus][line * model->xSize + column+1]: model->conventionalDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->conventionalDc[stepKMinus][(line-1) * model->xSize + column]: model->conventionalDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->conventionalDc[stepKMinus][line * model->xSize + column-1]: model->conventionalDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

            conventionalDcDiffusion = model->parametersModel.cDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line * model->xSize + column], model->hx);
            conventionalDcChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line * model->xSize + column],\
            model->parametersModel.avgDc, gradientOdcI, gradientOdcJ, model->hx) + model->parametersModel.chi * diffusionOdc * PreventionOverCrowdingTerm(conventionalDcKMinus, model->parametersModel.avgDc);;

            //Difussion and Chemotaxis CD8T

            valIPlus  = (line != model->xSize-1)? model->tCytotoxic[stepKMinus][(line+1) * model->xSize + column]: model->tCytotoxic[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->tCytotoxic[stepKMinus][line * model->xSize + column+1]: model->tCytotoxic[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->tCytotoxic[stepKMinus][(line-1) * model->xSize + column]: model->tCytotoxic[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->tCytotoxic[stepKMinus][line * model->xSize + column-1]: model->tCytotoxic[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

            tCytotoxicDiffusion = model->parametersModel.tCytoDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line * model->xSize + column], model->hx);
            tCytotoxicChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line * model->xSize + column],\
            model->parametersModel.avgT, gradientOdcI, gradientOdcJ, model->hx) + model->parametersModel.chi * diffusionOdc * PreventionOverCrowdingTerm(tCytotoxicKMinus, model->parametersModel.avgT);;

            //Difussion ADC

            valIPlus  = (line != model->xSize-1)? model->activatedDc[stepKMinus][(line+1) * model->xSize + column]: model->activatedDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->activatedDc[stepKMinus][line * model->xSize + column+1]: model->activatedDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->activatedDc[stepKMinus][(line-1) * model->xSize + column]: model->activatedDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->activatedDc[stepKMinus][line * model->xSize + column-1]: model->activatedDc[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

            activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][line * model->xSize + column], model->hx);

            //Difussion Antibody

            valIPlus  = (line != model->xSize-1)? model->antibody[stepKMinus][(line+1) * model->xSize + column]: model->antibody[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->antibody[stepKMinus][line * model->xSize + column+1]: model->antibody[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->antibody[stepKMinus][(line-1) * model->xSize + column]: model->antibody[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->antibody[stepKMinus][line * model->xSize + column-1]: model->antibody[stepKMinus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

            antibodyDiffusion = model->parametersModel.antibodyDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKMinus][line * model->xSize + column], model->hx);

            //*******************************************Solving Tissue equations*****************************************************

            //Microglia update
            microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
            microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

            model->microglia[stepKPlus][line * model->xSize + column] = microgliaKMinus + \
            model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

            //Conventional DC update
            conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
            conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
            conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

            model->conventionalDc[stepKPlus][line * model->xSize + column] = conventionalDcKMinus + \
            model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);

            //Activated DC update
            activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
            activatedDcMigration = model->thetaPV[(line * model->xSize) + column]*model->parametersModel.gammaD*(model->dendriticLymphNode[stepKPlus] - activatedDcKMinus);
            
            model->activatedDc[stepKPlus][line * model->xSize + column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);

            //CD8 T update
            tCytotoxicMigration = model->thetaBV[(line * model->xSize) + column]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus);
            
            model->tCytotoxic[stepKPlus][line * model->xSize + column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);

            //Antibody update
            odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
            antibodyMigration = model->thetaBV[(line * model->xSize) + column]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[stepKPlus] - antibodyKMinus);
            
            model->antibody[stepKPlus][line * model->xSize + column] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);

            //Oligodendrocytes update
            odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
            odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

            model->oligodendrocyte[stepKPlus][line * model->xSize + column] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
            if((kTime +1)%model->numStepsLN == 0){
            if(model->thetaBV[(line *model->xSize) + column] == 1){
                auxTCytotoxicBV += model->tCytotoxic[stepKPlus][(line *model->xSize) + column];
                auxAntibodyBV += model->antibody[stepKPlus][(line *model->xSize) + column];
            }
            if(model->thetaPV[(line *model->xSize) + column] == 1){
                auxAdcPV += model->activatedDc[stepKPlus][(line *model->xSize) + column];
            }
            }
        }
        }
        if((kTime%model->intervaloFiguras == 0 && model->saveFigs ==1) || kTime == model->timeLen){
            //Cada myRank manda o resultado para o myrank 0
            if(model->my_rank == 0){
                for(int iterRank = 1; iterRank < model->comm_sz; iterRank++){
                    int startLine = iterRank * model->numLines;
                    int numLines = model->numLines;
                    //MPI_Recv(&numLines, 1, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // MPI_Recv(composed, model->xSize * numLines * 6, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // DecomposeRecvMatrix(composed, model->oligodendrocyte[stepKPlus], model->tCytotoxic[stepKPlus], model->microglia[stepKPlus],
                    //     model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], numLines, 
                    //     startLine, model->xSize
                    // );
                    MPI_Recv(&model->oligodendrocyte[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&model->tCytotoxic[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&model->microglia[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&model->antibody[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&model->conventionalDc[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&model->activatedDc[stepKPlus][startLine * model->xSize], model->xSize * numLines, MPI_FLOAT, iterRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                }
            }else{
                //MPI_Send(&model->numLines, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);
                // float* composed = ComposeMatrixToSend(model->oligodendrocyte[stepKPlus], model->tCytotoxic[stepKPlus], model->microglia[stepKPlus],
                //     model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], model->numLines, model->startLine, model->xSize
                // );
                //MPI_Send(composed, model->xSize * model->numLines * 6, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&model->oligodendrocyte[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&model->tCytotoxic[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
                MPI_Send(&model->microglia[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
                MPI_Send(&model->antibody[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
                MPI_Send(&model->conventionalDc[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
                MPI_Send(&model->activatedDc[stepKPlus][model->startLine * model->xSize], model->xSize * model->numLines, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(model->my_rank == 0 && model->saveFigs == 1){
                WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
            }
        }
       
        if(kTime > 1) {
            SendBorderThread(model, stepKPlus);
        }
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(model->my_rank == 0){
        double end = MPI_Wtime();
        double time_spent = (double)(end - start);
        printf("\nComputation Done in %lf seconds!!\n", time_spent);
        SavingData(*model);
        if(model->saveFigs == 1){
            printf("Saving results...\n\n");
            WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
            PlotResults(*model);
        }
    }
}


