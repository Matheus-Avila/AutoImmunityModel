#include "model.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void InitialConditionTissueMicroglia(structModel* model){
    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        if(pow((i-(int)(model->xSize/2)),2) + pow((j-(int)(model->xSize/2)),2) < 5 / (model->hx * model->hx)){
            model->microglia[0][k] = (float)model->parametersModel.avgMic/3;
        }
    }
}

void InitialConditionLymphNode(structModel *model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN)
{
    model->dendriticLymphNode[0] = dendriticLN;
    model->tHelperLymphNode[0] = thelperLN;
    model->tCytotoxicLymphNode[0] = tcytotoxicLN;
    model->bCellLymphNode[0] = bcellLN;
    model->plasmaCellLymphNode[0] = plasmacellLN;
    model->antibodyLymphNode[0] = antibodyLN;
}

int VerifyCFL(structParameters parametersModel, float ht, float hx){
    if(parametersModel.micDiffusion*ht/(hx*hx) < 0.25 && parametersModel.cDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.aDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.tCytoDiffusion*ht/(hx*hx) < 0.25 && parametersModel.chi*ht/hx < 0.5)
        return 1;
    return 0;
}

void WritePopulation(structModel model, float *population, char *fileName, char *bufferTime)
{
    FILE *file;
    file = fopen(fileName, "w");
    if(file != NULL){
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
    }else{
        printf("Error lymph node file\n");
        exit(0);
    }
}

void WritePopulationLymphNode(structModel model, float *population, char *fileName)
{
    FILE *file;
    file = fopen(fileName, "w");
    if(file != NULL){
        for(int i=0;i<model.numPointsLN;i++){
            fprintf(file, "%f\n", population[i]);
        }
        fclose(file);
    }else{
        printf("Error matrix file\n");
        exit(0);
    }
}

void WriteLymphNodeFiles(structModel model, float *dendritic, float *tHelper, float *tCytotoxic, float *bCell, float *plasmaCell, float *antibody)
{
    WritePopulationLymphNode(model, dendritic, "./result/dendritic.txt");
    WritePopulationLymphNode(model, tHelper, "./result/tHelper.txt");
    WritePopulationLymphNode(model, tCytotoxic, "./result/tCyto.txt");
    WritePopulationLymphNode(model, bCell, "./result/bCell.txt");
    WritePopulationLymphNode(model, plasmaCell, "./result/plasmaCell.txt");
    WritePopulationLymphNode(model, antibody, "./result/antibody.txt");

    char buffer[10];
    char command[40] = {};
    strcat(command, "python3 plotLymphNode.py ");
    snprintf(buffer, sizeof(buffer), "%d", model.tFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", (model.tSize / model.numPointsLN) * model.ht);
    strcat(command, buffer);
    system(command);
}

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float *activatedDC, float time)
{
    char buffer[10];
    float day = time * model.ht;

    snprintf(buffer, sizeof(buffer), "%.1f", day);

    char pathOligodendrocytes[200] = "./result/matrix/oligo";
    strcat(pathOligodendrocytes, buffer);
    strcat(pathOligodendrocytes, ".txt");
    WritePopulation(model, oligodendrocyte, pathOligodendrocytes, buffer);

    char pathMicroglia[200] = "./result/matrix/microglia";
    strcat(pathMicroglia, buffer);
    strcat(pathMicroglia, ".txt");
    WritePopulation(model, microglia, pathMicroglia, buffer);

    char pathTCyto[200] = "./result/matrix/tCyto";
    strcat(pathTCyto, buffer);
    strcat(pathTCyto, ".txt");
    WritePopulation(model, tCytotoxic, pathTCyto, buffer);

    char pathAntibody[200] = "./result/matrix/antibody";
    strcat(pathAntibody, buffer);
    strcat(pathAntibody, ".txt");
    WritePopulation(model, antibody, pathAntibody, buffer);

    char pathConventionalDC[200] = "./result/matrix/conventionalDC";
    strcat(pathConventionalDC, buffer);
    strcat(pathConventionalDC, ".txt");
    WritePopulation(model, conventionalDC, pathConventionalDC, buffer);

    char pathActivatedDC[200] = "./result/matrix/activatedDC";
    strcat(pathActivatedDC, buffer);
    strcat(pathActivatedDC, ".txt");
    WritePopulation(model, activatedDC, pathActivatedDC, buffer);
}

void PlotResults(structModel model)
{
    char buffer[10];
    char command[70] = {};
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
    snprintf(buffer, sizeof(buffer), "%d", model.tSize/model.intervalFigures);
    strcat(command, buffer);
    system(command);
}

__device__ void PreventionOverCrowdingTerm(float populationPoint, float avgValue, float *result)
{
    *result = populationPoint / (populationPoint + avgValue);
}

__device__ void UpDownWind(float frontPoint, float rearPoint, float avgValue, float *result)
{
    float resultF;
    PreventionOverCrowdingTerm(frontPoint, avgValue, result);
    PreventionOverCrowdingTerm(rearPoint, avgValue, &resultF);
    *result = *result - resultF;
}

__device__ void CalculateChemottaxis(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,
                                     float avgValue, float gradientOdcI, float gradientOdcJ, float *result)
{
    float gradientPopulationI, gradientPopulationJ;
    if (gradientOdcI < 0)
    {
        UpDownWind(frontIPoint, ijPoint, avgValue, &gradientPopulationI);
        gradientPopulationI = gradientPopulationI / (float)hx;
    }
    else
    {
        UpDownWind(ijPoint, rearIPoint, avgValue, &gradientPopulationI);
        gradientPopulationI = gradientPopulationI / (float)hx;
    }
    if (gradientOdcJ < 0)
    {
        UpDownWind(frontJPoint, ijPoint, avgValue, &gradientPopulationJ);
        gradientPopulationJ = gradientPopulationJ / (float)hx;
    }
    else
    {
        UpDownWind(ijPoint, rearJPoint, avgValue, &gradientPopulationJ);
        gradientPopulationJ = gradientPopulationJ / (float)hx;
    }

    *result = gradientOdcI * gradientPopulationI + gradientOdcJ * gradientPopulationJ;
}

__device__ void CalculateDiffusion(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float *result)
{
    *result = (float)(frontIPoint + frontJPoint - 4 * ijPoint + rearIPoint + rearJPoint) / (float)(hx * hx);
}

__device__ void fFunc(float valuePopulation, float avgPopulation, float *result)
{
    *result = valuePopulation * valuePopulation / (float)(valuePopulation + avgPopulation);
}

void WriteBVPV(structModel *model, float *thetaBV, float *thetaPV){
    FILE *fileBV;
    fileBV = fopen("./result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./result/pv.txt", "w");
    int k = 0;
    if(fileBV != NULL && filePV != NULL){
        while (k < model->xSize*model->xSize){
            int i = k;
            while (i < k + model->xSize){
                fprintf(fileBV, "%f ", thetaBV[i]);
                fprintf(filePV, "%f ", thetaPV[i]);
                i++;
            }
            fprintf(fileBV,"\n");
            fprintf(filePV,"\n");
            k+=model->xSize;
        }
        fclose(fileBV);
        fclose(filePV);
    }else{
        printf("Error matrix file\n");
        exit(0);
    }
    char buffer[10];
    char command[200] = {};
    strcat(command, "python3 plotBVPV.py ");
    snprintf(buffer, sizeof(buffer), "%d", model->xFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", model->hx);
    strcat(command, buffer);
    // system(command);
}

void DefineBVPV(structModel *model){
    int randomVal;
    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        randomVal = rand() % 100;
        if(randomVal <10){
            model->parametersModel.V_BV++;
            model->parametersModel.V_PV++;
            model->thetaBV[k] = 1;
            if(j != model->xSize-1)
                model->thetaPV[k+1] = 1;
            else
                model->thetaPV[k-model->xSize+1] = 1;
        }
    }
    model->parametersModel.V_BV = model->parametersModel.V_BV * model->hx * model->hx;
    model->parametersModel.V_PV = model->parametersModel.V_PV * model->hx * model->hx;
    printf("bv = %f, pv = %f \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    WriteBVPV(model, model->thetaBV, model->thetaPV);
}

structModel ModelInitialize(structParameters params, float ht, float hx, float time, float space, int numFigs, int numPointsLN, int numStepsLN, int saveFigs){
    structModel model;
    srand(2);
    model.parametersModel = params;
    if(!VerifyCFL(model.parametersModel, ht, hx)){
        printf("Falhou CFL!!\n");
        exit(1);
    }
    model.parametersModel = params;
    model.numFigs = numFigs;
    model.numPointsLN = numPointsLN;
    model.numStepsLN = numStepsLN;

    model.ht = ht;
    model.hx = hx;
    model.tFinal = time;
    model.xFinal = space;
    model.tSize = (int)(time/ht);
    model.xSize = (int)(space/hx);
    model.intervalFigures = (int)model.tSize/numFigs;
    model.saveFigs = saveFigs;

    model.microglia = (float**)malloc(BUFFER * sizeof(float*));
    model.oligodendrocyte = (float**)malloc(BUFFER * sizeof(float*));
    model.tCytotoxic = (float**)malloc(BUFFER * sizeof(float*));
    model.antibody = (float**)malloc(BUFFER * sizeof(float*));
    model.conventionalDc = (float**)malloc(BUFFER * sizeof(float*));
    model.activatedDc = (float**)malloc(BUFFER * sizeof(float*));

    model.activatedDCTissueVessels = 0;
    model.tCytotoxicTissueVessels = 0;
    model.antibodyTissueVessels = 0;

    for (int index=0;index<BUFFER;++index){
        model.microglia[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.oligodendrocyte[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.tCytotoxic[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.conventionalDc[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.activatedDc[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
        model.antibody[index] = (float*)calloc(model.xSize*model.xSize, sizeof(float));
    }
    //definir BV e PV
    model.thetaPV = (float*)calloc(model.xSize*model.xSize, sizeof(float));
    model.thetaBV = (float*)calloc(model.xSize*model.xSize, sizeof(float));
    DefineBVPV(&model);
    //definir lymph node
    model.dendriticLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.tCytotoxicLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.tHelperLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.antibodyLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.bCellLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));
    model.plasmaCellLymphNodeSavedPoints = (float*)calloc(model.numPointsLN, sizeof(float));

    model.dendriticLymphNode = (float*)calloc(2, sizeof(float));
    model.tCytotoxicLymphNode = (float*)calloc(2, sizeof(float));
    model.tHelperLymphNode = (float*)calloc(2, sizeof(float));
    model.antibodyLymphNode = (float*)calloc(2, sizeof(float));
    model.bCellLymphNode = (float*)calloc(2, sizeof(float));
    model.plasmaCellLymphNode = (float*)calloc(2, sizeof(float));    

    float dendriticLN = 0.0, thelperLN = params.stableTHelper, tcytotoxicLN = params.stableTCytotoxic, bcellLN = params.stableB, plasmacellLN = 0.0, antibodyLN = 0.0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    InitialConditionTissueMicroglia(&model);
    return model;
}

void verifyValues(structModel model, float value, int time, char* populationName){
    if(value < 0 ||  isnanf(value)){
        printf("Error: %s = (%f) :: time = %f\n", populationName, value, time*model.ht);
        exit(0);
    }
}

void verifyDerivate(structModel model, float value, int time, char* populationName){
    if(isnanf(value)){
        printf("Error: %s = (%f) :: time = %f\n", populationName, value, time*model.ht);
        exit(0);
    }
}

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
    float tCytoMigration = model.parametersModel.gammaT * (tCytoLN - model.tCytotoxicTissueVessels) * (float)(model.parametersModel.V_BV/model.parametersModel.V_LN) - (1 * model.parametersModel.epslon_x);
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
        fprintf(file, "Execution Time of Kernel = %f secs\n", model.execTimeKernel);
        fprintf(file, "Execution Time of Lymph Node = %f secs\n", model.elapsedTimeLymphNode);
        fprintf(file, "Execution Time of Copies Device to Host = %f secs\n", model.elapsedTimeCopiesDeviceToHost);
        fprintf(file, "Execution Time of Copies Host to Device = %f secs\n", model.elapsedTimeCopiesHostToDevice);
        fprintf(file, "Days = %d - Space = %d - ht = %f, hx = %f, Ht_JumpStep = %d\n", model.tFinal, model.xFinal, model.ht, model.hx, model.numStepsLN);
        fprintf(file, "Lymph node populations\n");
        fprintf(file, "DC = %f, TCD8 = %f, TCD4 = %f, B Cell = %f, Plasma cell = %f, IgG = %f\n", model.dendriticLymphNodeSavedPoints[model.numPointsLN-1], model.tCytotoxicLymphNodeSavedPoints[model.numPointsLN-1], model.tHelperLymphNodeSavedPoints[model.numPointsLN-1], model.bCellLymphNodeSavedPoints[model.numPointsLN-1], model.plasmaCellLymphNodeSavedPoints[model.numPointsLN-1], model.antibodyLymphNodeSavedPoints[model.numPointsLN-1]);
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
    }
}

__device__ __constant__ float upperNeumannBC, lowerNeumannBC, leftNeumannBC, rightNeumannBC, constHx, constHt, consthx2;
__device__ __constant__ int constXSize;
__device__ __constant__ structParameters modelParams;
const int threadsPerBlock = 64;
const int numBlocks = 16;

__global__ void kernelPDE(int kTime, float *tCytoSumVessel, float *activatedDCSumVessel, float *antibodySumVessel, float *devActivatedDCLymphNode, float *devAntibodyLymphNode, float *devTCytotoxicLymphNode, float *devThetaPV, float *devThetaBV, float *devMicrogliaKMinus, float *devMicrogliaKPlus, float *devTCytotoxicKMinus, float *devTCytotoxicKPlus, float *devAntibodyKMinus, float *devAntibodyKPlus, float *devConventionalDCKMinus, float *devConventionalDCKPlus, float *devActivatedDCKMinus, float *devActivatedDCKPlus, float *devOligodendrocyteKMinus, float *devOligodendrocyteKPlus)
{
    int thrIdx = blockIdx.x * blockDim.x + threadIdx.x;
    int vesselIdx = threadIdx.x;
    int line = (int)thrIdx / constXSize;
    int column = thrIdx % constXSize;
    float devOligodendrocyteKMinusThrIdx, devMicrogliaKMinusThrIdx, devConventionalDCKMinusThrIdx, devTCytotoxicKMinusThrIdx, devActivatedDCKMinusThrIdx, devAntibodyKMinusThrIdx;
    float avgOdcMinusODC, diffusionODC;
    __shared__ float tCytoSumVesselBlock[threadsPerBlock];
    __shared__ float activatedDCSumVesselBlock[threadsPerBlock];
    __shared__ float antibodySumVesselBlock[threadsPerBlock];
    for (int i = 0; i < threadsPerBlock; i++)
    {
        tCytoSumVesselBlock[i] = 0;
        activatedDCSumVesselBlock[i] = 0;
        antibodySumVesselBlock[i] = 0;
    }
    while (thrIdx < constXSize * constXSize)
    {
        devOligodendrocyteKMinusThrIdx = devOligodendrocyteKMinus[thrIdx];
        devMicrogliaKMinusThrIdx = devMicrogliaKMinus[thrIdx];
        devConventionalDCKMinusThrIdx = devConventionalDCKMinus[thrIdx];
        devTCytotoxicKMinusThrIdx = devTCytotoxicKMinus[thrIdx];
        devActivatedDCKMinusThrIdx = devActivatedDCKMinus[thrIdx];
        devAntibodyKMinusThrIdx = devAntibodyKMinus[thrIdx];
        line = (int)thrIdx / constXSize;
        column = thrIdx % constXSize;

        // Define gradient ODCs
        float valIPlus = (line != constXSize - 1) ? devOligodendrocyteKMinus[thrIdx + constXSize] : devOligodendrocyteKMinus[thrIdx - constXSize];
        float valJPlus = (column != constXSize - 1) ? devOligodendrocyteKMinus[thrIdx + 1] : devOligodendrocyteKMinus[thrIdx - 1];
        float valIMinus = (line != 0) ? devOligodendrocyteKMinus[thrIdx - constXSize] : devOligodendrocyteKMinus[thrIdx + constXSize];
        float valJMinus = (column != 0) ? devOligodendrocyteKMinus[thrIdx - 1] : devOligodendrocyteKMinus[thrIdx + 1];

        float gradientOdcI = (float)(valIPlus - valIMinus) / (float)(constHx*2);
        float gradientOdcJ = (float)(valJPlus - valJMinus) / (float)(constHx*2);

        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devOligodendrocyteKMinusThrIdx, &diffusionODC);
        // Diffusion and Chemotaxis Mic

        valIPlus = (line != constXSize - 1) ? devMicrogliaKMinus[thrIdx + constXSize] : devMicrogliaKMinus[thrIdx] + ( devMicrogliaKMinus[thrIdx] / (devMicrogliaKMinus[thrIdx] + modelParams.avgMic) ) * constHx * gradientOdcI/modelParams.micDiffusion;
        valJPlus = (column != constXSize - 1) ? devMicrogliaKMinus[thrIdx + 1] : devMicrogliaKMinus[thrIdx] + ( devMicrogliaKMinus[thrIdx] / (devMicrogliaKMinus[thrIdx] + modelParams.avgMic) ) * constHx * gradientOdcJ/modelParams.micDiffusion;
        valIMinus = (line != 0) ? devMicrogliaKMinus[thrIdx - constXSize] : devMicrogliaKMinus[thrIdx] - ( devMicrogliaKMinus[thrIdx] / (devMicrogliaKMinus[thrIdx] + modelParams.avgMic) ) * constHx * gradientOdcI/modelParams.micDiffusion;
        valJMinus = (column != 0) ? devMicrogliaKMinus[thrIdx - 1] : devMicrogliaKMinus[thrIdx] - ( devMicrogliaKMinus[thrIdx] / (devMicrogliaKMinus[thrIdx] + modelParams.avgMic) ) * constHx * gradientOdcJ/modelParams.micDiffusion;

        float microgliaDiffusion = 0;
        float microgliaChemotaxis = 0;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinusThrIdx, &microgliaDiffusion);
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinusThrIdx,
                             modelParams.avgMic, gradientOdcI, gradientOdcJ, &microgliaChemotaxis);
        microgliaChemotaxis += diffusionODC * devMicrogliaKMinusThrIdx / (devMicrogliaKMinusThrIdx + modelParams.avgMic);
        microgliaChemotaxis *= modelParams.chi;
        microgliaDiffusion *= modelParams.micDiffusion;
        // Diffusion and Chemotaxis CDC

        valIPlus = (line != constXSize - 1) ? devConventionalDCKMinus[thrIdx + constXSize] : devConventionalDCKMinus[thrIdx] + ( devConventionalDCKMinus[thrIdx] / (devConventionalDCKMinus[thrIdx] + modelParams.avgDc) ) * constHx * gradientOdcI/modelParams.cDcDiffusion;
        valJPlus = (column != constXSize - 1) ? devConventionalDCKMinus[thrIdx + 1] : devConventionalDCKMinus[thrIdx] + ( devConventionalDCKMinus[thrIdx] / (devConventionalDCKMinus[thrIdx] + modelParams.avgDc) ) * constHx * gradientOdcJ/modelParams.cDcDiffusion;
        valIMinus = (line != 0) ? devConventionalDCKMinus[thrIdx - constXSize] : devConventionalDCKMinus[thrIdx] - ( devConventionalDCKMinus[thrIdx] / (devConventionalDCKMinus[thrIdx] + modelParams.avgDc) ) * constHx * gradientOdcI/modelParams.cDcDiffusion;
        valJMinus = (column != 0) ? devConventionalDCKMinus[thrIdx - 1] : devConventionalDCKMinus[thrIdx] - ( devConventionalDCKMinus[thrIdx] / (devConventionalDCKMinus[thrIdx] + modelParams.avgDc) ) * constHx * gradientOdcJ/modelParams.cDcDiffusion;

        float conventionalDcDiffusion = 0;
        float conventionalDcChemotaxis = 0;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devConventionalDCKMinusThrIdx, &conventionalDcDiffusion);
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devConventionalDCKMinusThrIdx,
                             modelParams.avgDc, gradientOdcI, gradientOdcJ, &conventionalDcChemotaxis);
        conventionalDcChemotaxis += diffusionODC * devConventionalDCKMinusThrIdx / (devConventionalDCKMinusThrIdx + modelParams.avgDc);
        conventionalDcChemotaxis *= modelParams.chi;
        conventionalDcDiffusion *= modelParams.cDcDiffusion;

        // Difussion and Chemotaxis CD8T

        valIPlus = (line != constXSize - 1) ? devTCytotoxicKMinus[thrIdx + constXSize] : devTCytotoxicKMinus[thrIdx] + ( devTCytotoxicKMinus[thrIdx] / (devTCytotoxicKMinus[thrIdx] + modelParams.avgT) ) * constHx * gradientOdcI/modelParams.tCytoDiffusion;
        valJPlus = (column != constXSize - 1) ? devTCytotoxicKMinus[thrIdx + 1] : devTCytotoxicKMinus[thrIdx] + ( devTCytotoxicKMinus[thrIdx] / (devTCytotoxicKMinus[thrIdx] + modelParams.avgT) ) * constHx * gradientOdcJ/modelParams.tCytoDiffusion;
        valIMinus = (line != 0) ? devTCytotoxicKMinus[thrIdx - constXSize] : devTCytotoxicKMinus[thrIdx] - ( devTCytotoxicKMinus[thrIdx] / (devTCytotoxicKMinus[thrIdx] + modelParams.avgT) ) * constHx * gradientOdcI/modelParams.tCytoDiffusion;
        valJMinus = (column != 0) ? devTCytotoxicKMinus[thrIdx - 1] : devTCytotoxicKMinus[thrIdx] - ( devTCytotoxicKMinus[thrIdx] / (devTCytotoxicKMinus[thrIdx] + modelParams.avgT) ) * constHx * gradientOdcJ/modelParams.tCytoDiffusion;

        float tCytotoxicDiffusion = 0;
        float tCytotoxicChemotaxis = 0;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinusThrIdx, &tCytotoxicDiffusion);
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinusThrIdx,
                             modelParams.avgT, gradientOdcI, gradientOdcJ, &tCytotoxicChemotaxis);
        tCytotoxicChemotaxis += diffusionODC * devTCytotoxicKMinusThrIdx / (devTCytotoxicKMinusThrIdx + modelParams.avgT);
        tCytotoxicChemotaxis *= modelParams.chi;
        tCytotoxicDiffusion *= modelParams.tCytoDiffusion;

        // Difussion ADC

        valIPlus = (line != constXSize - 1) ? devActivatedDCKMinus[thrIdx + constXSize] : devActivatedDCKMinus[thrIdx - constXSize] - (float)(constHx*2 * lowerNeumannBC);
        valJPlus = (column != constXSize - 1) ? devActivatedDCKMinus[thrIdx + 1] : devActivatedDCKMinus[thrIdx - 1] - (float)(constHx*2 * rightNeumannBC);
        valIMinus = (line != 0) ? devActivatedDCKMinus[thrIdx - constXSize] : devActivatedDCKMinus[thrIdx + constXSize] - (float)(constHx*2 * upperNeumannBC);
        valJMinus = (column != 0) ? devActivatedDCKMinus[thrIdx - 1] : devActivatedDCKMinus[thrIdx + 1] - (float)(constHx*2 * leftNeumannBC);

        float activatedDCDiffusion = 0;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devActivatedDCKMinusThrIdx, &activatedDCDiffusion);
        activatedDCDiffusion *= modelParams.aDcDiffusion;

        // Difussion Antibody

        valIPlus = (line != constXSize - 1) ? devAntibodyKMinus[thrIdx + constXSize] : devAntibodyKMinus[thrIdx - constXSize] - (float)(constHx*2 * lowerNeumannBC);
        valJPlus = (column != constXSize - 1) ? devAntibodyKMinus[thrIdx + 1] : devAntibodyKMinus[thrIdx - 1] - (float)(constHx*2 * rightNeumannBC);
        valIMinus = (line != 0) ? devAntibodyKMinus[thrIdx - constXSize] : devAntibodyKMinus[thrIdx + constXSize] - (float)(constHx*2 * upperNeumannBC);
        valJMinus = (column != 0) ? devAntibodyKMinus[thrIdx - 1] : devAntibodyKMinus[thrIdx + 1] - (float)(constHx*2 * leftNeumannBC);

        float antibodyDiffusion = 0;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devAntibodyKMinusThrIdx, &antibodyDiffusion);
        antibodyDiffusion *= modelParams.antibodyDiffusion;

        //*******************************************Solving Tissue equations*****************************************************
        
        // Microglia update
        float microgliaReaction = modelParams.muMic * devMicrogliaKMinusThrIdx * (modelParams.avgMic - devMicrogliaKMinusThrIdx);
        float microgliaClearance = modelParams.cMic * devMicrogliaKMinusThrIdx;

        devMicrogliaKPlus[thrIdx] = devMicrogliaKMinusThrIdx +
                                    constHt * (microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

        // Conventional DC update
        float conventionalDcReaction = modelParams.muCDc * devOligodendrocyteKMinusThrIdx * (modelParams.avgDc - devConventionalDCKMinusThrIdx);
        float conventionalDcActivation = modelParams.bD * devConventionalDCKMinusThrIdx * devOligodendrocyteKMinusThrIdx;
        float conventionalDcClearance = modelParams.cCDc * devConventionalDCKMinusThrIdx;

        devConventionalDCKPlus[thrIdx] = devConventionalDCKMinusThrIdx +
                                         constHt * (conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);

        // Activated DC update
        float activatedDcClearance = modelParams.cADc * devActivatedDCKMinusThrIdx;
        float activatedDcMigration = devThetaPV[thrIdx] * modelParams.gammaD * (*devActivatedDCLymphNode - devActivatedDCKMinusThrIdx);

        devActivatedDCKPlus[thrIdx] = devActivatedDCKMinusThrIdx + constHt * (activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);

        // CD8 T update
        float tCytotoxicMigration = devThetaBV[thrIdx] * modelParams.gammaT * (*devTCytotoxicLymphNode - devTCytotoxicKMinusThrIdx);

        devTCytotoxicKPlus[thrIdx] = devTCytotoxicKMinusThrIdx + constHt * (tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);

        // Antibody update
        float resultFFuncMic = 0;
        fFunc(devMicrogliaKMinusThrIdx, modelParams.avgMic, &resultFFuncMic);
        avgOdcMinusODC = modelParams.avgOdc - devOligodendrocyteKMinusThrIdx;
        float odcAntibodyMicrogliaFagocitosis = modelParams.lambAntMic * devAntibodyKMinusThrIdx * avgOdcMinusODC * resultFFuncMic;
        float antibodyMigration = devThetaBV[thrIdx] * modelParams.gammaAntibody * (*devAntibodyLymphNode - devAntibodyKMinusThrIdx);

        devAntibodyKPlus[thrIdx] = devAntibodyKMinusThrIdx + constHt * (antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);

        // Oligodendrocytes update
        float result = 0;
        fFunc(devTCytotoxicKMinusThrIdx, modelParams.avgT, &result);
        float odcMicrogliaFagocitosis = modelParams.rM * resultFFuncMic * avgOdcMinusODC;
        float odcTCytotoxicApoptosis = modelParams.rT * result * avgOdcMinusODC;

        devOligodendrocyteKPlus[thrIdx] = devOligodendrocyteKMinusThrIdx + constHt * (odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);

        if (devThetaBV[thrIdx] == 1)
        {
            tCytoSumVesselBlock[vesselIdx] += devTCytotoxicKPlus[thrIdx];
            antibodySumVesselBlock[vesselIdx] += devAntibodyKPlus[thrIdx];
        }
        if (devThetaPV[thrIdx] == 1)
        {
            activatedDCSumVesselBlock[vesselIdx] += devActivatedDCKPlus[thrIdx];
        }
        thrIdx += gridDim.x * blockDim.x;
    }
    __syncthreads();
    int i = blockDim.x / 2;
    while (i != 0)
    {
        if (vesselIdx < i)
        {
            tCytoSumVesselBlock[vesselIdx] += tCytoSumVesselBlock[vesselIdx + i];
            activatedDCSumVesselBlock[vesselIdx] += activatedDCSumVesselBlock[vesselIdx + i];
            antibodySumVesselBlock[vesselIdx] += antibodySumVesselBlock[vesselIdx + i];
        }
        __syncthreads();
        i /= 2;
    }
    if (vesselIdx == 0)
    {
        tCytoSumVessel[blockIdx.x] = tCytoSumVesselBlock[0];
        activatedDCSumVessel[blockIdx.x] = activatedDCSumVesselBlock[0];
        antibodySumVessel[blockIdx.x] = antibodySumVesselBlock[0];
    }
}

void RunModel(structModel *model)
{
    // Save IC
    if(model->saveFigs)
        WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);

    clock_t start, end;
    float elapsedTimeLymphNode = 0, elapsedTimeCopiesDeviceToHost = 0, elapsedTimeCopiesHostToDevice = 0;
    float elapsedTimeKernel = 0, elapsedTimeKernelAux = 0;
    cudaEvent_t startKernel, stopKernel;
    cudaEventCreate(&startKernel);
    cudaEventCreate(&stopKernel);

    float *activatedDCVessel, *tCytotoxicVessel, *antibodyVessel;

    float *devThetaPV, *devThetaBV, *devActivatedDCVessel, *devTCytotoxicVessel, *devAntibodyVessel, *devActivatedDCLymphNode, *devAntibodyLymphNode, *devTCytotoxicLymphNode, *devMicrogliaKMinus, *devMicrogliaKPlus, *devTCytotoxicKMinus, *devTCytotoxicKPlus, *devAntibodyKMinus, *devAntibodyKPlus, *devConventionalDCKMinus, *devConventionalDCKPlus, *devActivatedDCKMinus, *devActivatedDCKPlus, *devOligodendrocytesDCKMinus, *devOligodendrocytesDCKPlus;

    cudaMalloc((void **)&devThetaPV, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devThetaBV, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devOligodendrocytesDCKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devOligodendrocytesDCKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devMicrogliaKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devMicrogliaKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devTCytotoxicKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devTCytotoxicKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devAntibodyKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devAntibodyKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devConventionalDCKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devConventionalDCKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMalloc((void **)&devActivatedDCKMinus, model->xSize * model->xSize * sizeof(float));
    cudaMalloc((void **)&devActivatedDCKPlus, model->xSize * model->xSize * sizeof(float));

    cudaMemcpy(devThetaBV, model->thetaBV, model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devThetaPV, model->thetaPV, model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(devOligodendrocytesDCKMinus, model->oligodendrocyte[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devMicrogliaKMinus, model->microglia[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devTCytotoxicKMinus, model->tCytotoxic[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devAntibodyKMinus, model->antibody[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devConventionalDCKMinus, model->conventionalDc[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devActivatedDCKMinus, model->activatedDc[0], model->xSize * model->xSize * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&devActivatedDCLymphNode, sizeof(float));
    cudaMalloc((void **)&devAntibodyLymphNode, sizeof(float));
    cudaMalloc((void **)&devTCytotoxicLymphNode, sizeof(float));

    cudaMalloc((void **)&devActivatedDCVessel, numBlocks * sizeof(float));
    cudaMalloc((void **)&devAntibodyVessel, numBlocks * sizeof(float));
    cudaMalloc((void **)&devTCytotoxicVessel, numBlocks * sizeof(float));
    // Inicializar os constant com os valores
    int stepKMinus = 0, stepKPlus;

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;

    float bc = 0.0;
    float hx2 = model->hx * 2;

    cudaMemcpyToSymbol(upperNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(lowerNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(leftNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(rightNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(constHt, &model->ht, sizeof(float));
    cudaMemcpyToSymbol(constHx, &model->hx, sizeof(float));
    cudaMemcpyToSymbol(consthx2, &hx2, sizeof(float));
    
    cudaMemcpyToSymbol(constXSize, &model->xSize, sizeof(float));
    cudaMemcpyToSymbol(modelParams, &model->parametersModel, sizeof(structParameters));

    int devKTime;
    cudaMalloc((void **)&devKTime, sizeof(int));
    activatedDCVessel = (float *)calloc(numBlocks, sizeof(float));
    antibodyVessel = (float *)calloc(numBlocks, sizeof(float));
    tCytotoxicVessel = (float *)calloc(numBlocks, sizeof(float));
    for (int kTime = 1; kTime <= model->tSize; kTime++)
    {	
        // solve lymphnode
        if(kTime%model->numStepsLN == 0){
            // Copia do device para o host as integrais do tecido
            start = clock(); 
            cudaMemcpy(activatedDCVessel, devActivatedDCVessel, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(antibodyVessel, devAntibodyVessel, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(tCytotoxicVessel, devTCytotoxicVessel, numBlocks * sizeof(float), cudaMemcpyDeviceToHost);
            end = clock();

            elapsedTimeCopiesDeviceToHost += ((float) (end - start)) / CLOCKS_PER_SEC;

            auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;
            for (int pos = 0; pos < numBlocks; pos++)
            {
                auxAdcPV += activatedDCVessel[pos];
                auxAntibodyBV += antibodyVessel[pos];
                auxTCytotoxicBV += tCytotoxicVessel[pos];
            }
            model->tCytotoxicTissueVessels = auxTCytotoxicBV * model->hx * model->hx / model->parametersModel.V_BV;
            model->antibodyTissueVessels = auxAntibodyBV * model->hx * model->hx / model->parametersModel.V_BV;
            model->activatedDCTissueVessels = auxAdcPV * model->hx * model->hx / model->parametersModel.V_PV;  

            start = clock(); 
            SolverLymphNode(model, kTime); 
            end = clock();

            elapsedTimeLymphNode += ((float) (end - start)) / CLOCKS_PER_SEC;
        }        
        
        stepKPlus = kTime % 2;
        // copiar LN pra GPU
        if(kTime%model->numStepsLN == 0){
            start = clock(); 
            cudaMemcpy(devActivatedDCLymphNode, &model->dendriticLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(devAntibodyLymphNode, &model->antibodyLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(devTCytotoxicLymphNode, &model->tCytotoxicLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
            end = clock();
            elapsedTimeCopiesHostToDevice += ((float) (end - start)) / CLOCKS_PER_SEC;
        }
        cudaMemcpy(&devKTime, &kTime, sizeof(int), cudaMemcpyHostToDevice);

        cudaEventRecord(startKernel, 0);
        if (stepKPlus % 2 == 1)
            kernelPDE<<<numBlocks, threadsPerBlock>>>(devKTime, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKMinus, devMicrogliaKPlus, devTCytotoxicKMinus, devTCytotoxicKPlus, devAntibodyKMinus, devAntibodyKPlus, devConventionalDCKMinus, devConventionalDCKPlus, devActivatedDCKMinus, devActivatedDCKPlus, devOligodendrocytesDCKMinus, devOligodendrocytesDCKPlus);
        else
            kernelPDE<<<numBlocks, threadsPerBlock>>>(devKTime, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKPlus, devMicrogliaKMinus, devTCytotoxicKPlus, devTCytotoxicKMinus, devAntibodyKPlus, devAntibodyKMinus, devConventionalDCKPlus, devConventionalDCKMinus, devActivatedDCKPlus, devActivatedDCKMinus, devOligodendrocytesDCKPlus, devOligodendrocytesDCKMinus);
        cudaEventRecord(stopKernel, 0);
        cudaEventSynchronize(stopKernel);
        cudaEventElapsedTime(&elapsedTimeKernelAux, startKernel, stopKernel);
        elapsedTimeKernel += elapsedTimeKernelAux;
        if (model->saveFigs && kTime % model->intervalFigures == 0)
        {
            if (stepKPlus % 2 == 1)
            {
                cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->antibody[stepKPlus], devAntibodyKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
            }
            else
            {
                cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->antibody[stepKPlus], devAntibodyKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
            }
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        }else{
            if (kTime == model->tSize){
                if (stepKPlus % 2 == 1)
                {
                    cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->antibody[stepKPlus], devAntibodyKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKPlus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                }
                else
                {
                    cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->antibody[stepKPlus], devAntibodyKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKMinus, model->xSize * model->xSize * sizeof(float), cudaMemcpyDeviceToHost);
                }
                if(model->saveFigs)
                    WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
            }
        }
        stepKMinus += 1;
        stepKMinus = stepKMinus % 2;
    }
    model->elapsedTimeLymphNode = elapsedTimeLymphNode;
    model->elapsedTimeCopiesDeviceToHost = elapsedTimeCopiesDeviceToHost;
    model->elapsedTimeCopiesHostToDevice = elapsedTimeCopiesHostToDevice;
    model->execTimeKernel = elapsedTimeKernel/1000;
    printf("Computation Done!!\n");
    SavingData(*model);
    if(model->saveFigs){
        printf("Saving results...\n\n");
        WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
        PlotResults(*model);
    }
}