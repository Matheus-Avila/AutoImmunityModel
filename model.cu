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
        if(pow((i-(int)(model->xSize/2)),2) + pow((j-(int)(model->xSize/2)),2) < 5){
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
    snprintf(buffer, sizeof(buffer), "%f", (model.tSize/model.numPointsLN)*model.ht);
    strcat(command, buffer);
    system(command);
}

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float  *activatedDC, float time){
    char buffer[10];
    float day = time * model.ht;
    
    snprintf(buffer, sizeof(buffer), "%.1f", day);
    
    char pathOligodendrocytes[50] = "./result/matrix/oligo";
    strcat(pathOligodendrocytes, buffer);
    strcat(pathOligodendrocytes, ".txt");
    WritePopulation(model, oligodendrocyte, pathOligodendrocytes, buffer);

    char pathMicroglia[50] = "./result/matrix/microglia";
    strcat(pathMicroglia, buffer);
    strcat(pathMicroglia, ".txt");
    WritePopulation(model, microglia, pathMicroglia, buffer);

    char pathTCyto[50] = "./result/matrix/tCyto";
    strcat(pathTCyto, buffer);
    strcat(pathTCyto, ".txt");
    WritePopulation(model, tCytotoxic, pathTCyto, buffer);

    char pathAntibody[50] = "./result/matrix/antibody";
    strcat(pathAntibody, buffer);
    strcat(pathAntibody, ".txt");
    WritePopulation(model, antibody, pathAntibody, buffer);

    char pathConventionalDC[50] = "./result/matrix/conventionalDC";
    strcat(pathConventionalDC, buffer);
    strcat(pathConventionalDC, ".txt");
    WritePopulation(model, conventionalDC, pathConventionalDC, buffer);

    char pathActivatedDC[50] = "./result/matrix/activatedDC";
    strcat(pathActivatedDC, buffer);
    strcat(pathActivatedDC, ".txt");
    WritePopulation(model, activatedDC, pathActivatedDC, buffer);
}   

void PlotResults(structModel model){
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
    snprintf(buffer, sizeof(buffer), "%d", model.intervalFigures);
    strcat(command, buffer);
    system(command);
}

float AdvectionTerm(float populationPoint, float avgValue){
    return populationPoint/(populationPoint + avgValue);
}

float UpDownWind(float frontPoint, float rearPoint, float avgValue){
    return AdvectionTerm(frontPoint, avgValue) - AdvectionTerm(rearPoint, avgValue);
}

float CalculateChemottaxis(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ){
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

float CalculateDiffusion(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint){
    return (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(hx*hx);
}

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
}

void WriteBVPV(structModel *model, float *thetaBV, float *thetaPV){
    FILE *fileBV;
    fileBV = fopen("./result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./result/pv.txt", "w");
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
    printf("bv = %d, pv = %d \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    WriteBVPV(model, model->thetaBV, model->thetaPV);
}



structModel ModelInitialize(structParameters params, float ht, float hx, float time, float space, int numFigs, int numPointsLN){
    structModel model;
    srand(2);

    //Pegar os valores pelos parametros
    model.parametersModel = params;
    model.numFigs = numFigs;
    model.numPointsLN = numPointsLN;
    model.ht = ht;
    model.hx = hx;
    model.tFinal = time;
    model.xFinal = space;
    model.tSize = (int)(time/ht);
    model.xSize = (int)(space/hx);
    model.intervalFigures = (int)model.tSize/numFigs;

    //inicializar dinamicamente todos os vetores do tecido
    model.microglia = (float**)malloc(BUFFER * sizeof(float*));
    model.oligodendrocyte = (float**)malloc(BUFFER * sizeof(float*));
    model.tCytotoxic = (float**)malloc(BUFFER * sizeof(float*));
    model.antibody = (float**)malloc(BUFFER * sizeof(float*));
    model.conventionalDc = (float**)malloc(BUFFER * sizeof(float*));
    model.activatedDc = (float**)malloc(BUFFER * sizeof(float*));
    for (int index=0;index<BUFFER;++index){
        model.microglia[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
        model.oligodendrocyte[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
        model.tCytotoxic[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
        model.antibody[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
        model.conventionalDc[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
        model.activatedDc[index] = (float*)malloc(model.xSize*model.xSize * sizeof(float));
    }
    //definir BV e PV
    model.thetaPV = (float*)malloc(model.xSize*model.xSize * sizeof(float));
    model.thetaBV = (float*)malloc(model.xSize*model.xSize * sizeof(float));
    DefineBVPV(&model);
    //definir lymph node
    model.dendriticLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));
    model.tCytotoxicLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));
    model.tHelperLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));
    model.antibodyLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));
    model.bCellLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));
    model.plasmaCellLymphNodeSavedPoints = (float*)malloc(model.numPointsLN * sizeof(float));

    model.dendriticLymphNode = (float*)malloc(2 * sizeof(float));
    model.tCytotoxicLymphNode = (float*)malloc(2 * sizeof(float));
    model.tHelperLymphNode = (float*)malloc(2 * sizeof(float));
    model.antibodyLymphNode = (float*)malloc(2 * sizeof(float));
    model.bCellLymphNode = (float*)malloc(2 * sizeof(float));
    model.plasmaCellLymphNode = (float*)malloc(2 * sizeof(float));    

    float dendriticLN = 0.0, thelperLN = 0.0, tcytotoxicLN = 0.0, bcellLN = 0.0, plasmacellLN = 0.0, antibodyLN = 0.0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    InitialConditionTissueMicroglia(&model);
    return model;
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
    float tCytoHomeostasis = model.parametersModel.alphaTCytotoxic * (model.parametersModel.estableTCytotoxic - tCytoLN);
    float tCytoMigration = model.parametersModel.gammaT * (tCytoLN - model.tCytotoxicTissueVessels) * (float)(model.parametersModel.V_BV/model.parametersModel.V_LN);
    result[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

    //T Helper
    float tHelperActivation = model.parametersModel.bTHelper * (model.parametersModel.rhoTHelper * tHelperLN * dcLN - tHelperLN * dcLN);
    float tHelperHomeostasis = model.parametersModel.alphaTHelper * (model.parametersModel.estableTHelper - tHelperLN);
    float tHelperDispendure = model.parametersModel.bRho * dcLN * tHelperLN * bCellLN;
    result[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

    //B Cell
    float bCellActivation = model.parametersModel.bRhoB * (model.parametersModel.rhoB * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
    float bcellHomeostasis = model.parametersModel.alphaB * (model.parametersModel.estableB - bCellLN);
    result[3] = bcellHomeostasis + bCellActivation;

    //Plasma Cells
    float plasmaActivation = model.parametersModel.bRhoP * (model.parametersModel.rhoP * tHelperLN * dcLN * bCellLN);
    float plasmaHomeostasis = model.parametersModel.alphaP * (model.parametersModel.estableP - plasmaCellLN);
    result[4] = plasmaHomeostasis + plasmaActivation;

    //Antibody
    float antibodyProduction = model.parametersModel.rhoAntibody * plasmaCellLN;
    float antibodyDecayment = model.parametersModel.cF * antibodyLN;
    float antibodyMigration = model.parametersModel.gammaAntibody * (antibodyLN - model.antibodyTissueVessels) * (float)(model.parametersModel.V_BV/model.parametersModel.V_LN);
    result[5] = antibodyProduction - antibodyMigration;

    return result;
}

void verifyValues(structModel model, float value, int time, char* populationName){
    if(value < 0 ||  isnanf(value)){
        printf("Error: %s = (%f) :: time = %f\n", populationName, value, time*model.ht);
        exit(0);
    }
}

void SolverLymphNode(structModel *model, int stepPos){
    float populationLN[6];
    int stepKMinus = stepPos%2;
    int stepKPlus = (stepKMinus+1)%2;
    populationLN[0] = model->dendriticLymphNode[stepKMinus];
    populationLN[1] = model->tCytotoxicLymphNode[stepKMinus];
    populationLN[2] = model->tHelperLymphNode[stepKMinus];
    populationLN[3] = model->bCellLymphNode[stepKMinus];
    populationLN[4] = model->plasmaCellLymphNode[stepKMinus];
    populationLN[5] = model->antibodyLymphNode[stepKMinus];
    
    float* solutionLN;
    solutionLN = EquationsLymphNode(*model, populationLN, stepPos);
    
    //Execute Euler 
    model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + model->ht*solutionLN[0];
    model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + model->ht*solutionLN[1];
    model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + model->ht*solutionLN[2];
    model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + model->ht*solutionLN[3];
    model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + model->ht*solutionLN[4];
    model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + model->ht*solutionLN[5];
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
    verifyValues(*model, model->dendriticLymphNode[stepKPlus], stepPos, "DC lymph node");
    verifyValues(*model, model->tCytotoxicLymphNode[stepKPlus], stepPos, "CD8 T lymph node");
    verifyValues(*model, model->tHelperLymphNode[stepKPlus], stepPos, "CD4 T lymph node");
    verifyValues(*model, model->bCellLymphNode[stepKPlus], stepPos, "B cell lymph node");
    verifyValues(*model, model->plasmaCellLymphNode[stepKPlus], stepPos, "Plasma cell lymph node");
    verifyValues(*model, model->antibodyLymphNode[stepKPlus], stepPos, "Antibody lymph node");
}

__device__ float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0, hx, ht;
__device__ int xSize;
const int threadsPerBlock = 256;
const int numBlocks = 64;

__global__ void kernelPDE(structParameters* devParams, int kTime, int xSize, float* tCytoSumVessel, float* activatedDCSumVessel, float* antibodySumVessel, float *devActivatedDCLymphNode, float *devAntibodyLymphNode, float *devTCytotoxicLymphNode, float *devThetaPV, float *devThetaBV, float *devMicrogliaKMinus, float *devMicrogliaKPlus, float *devTCytotoxicKMinus, float *devTCytotoxicKPlus, float *devAntibodyKMinus, float *devAntibodyKPlus, float *devConventionalDCKMinus, float *devConventionalDCKPlus, float *devActivatedDCKMinus, float *devActivatedDCKPlus, float *devOligodendrocyteKMinus, float *devOligodendrocyteKPlus){
    int thrIdx = blockIdx.x*blockDim.x + threadIdx.x;
    int vesselIdx = threadIdx.x;
    int line = (int)thrIdx/xSize;
    int column = thrIdx%xSize;

    __shared__ float tCytoSumVesselBlock[threadsPerBlock];
    __shared__ float conventionalDCSumVesselBlock[threadsPerBlock];
    __shared__ float antibodySumVesselBlock[threadsPerBlock];
    while(thrIdx < xSize*xSize){
        line = (int)thrIdx/xSize;
        column = thrIdx%xSize;

        //Define gradient ODCs
        float valIPlus = (line != xSize-1)? devOligodendrocyteKMinus[thrIdx + xSize]: devOligodendrocyteKMinus[thrIdx];
        float valJPlus = (column != xSize-1)? devOligodendrocyteKMinus[thrIdx + 1]: devOligodendrocyteKMinus[thrIdx];
        float valIMinus = (line != 0)? devOligodendrocyteKMinus[thrIdx - xSize]: devOligodendrocyteKMinus[thrIdx];
        float valJMinus = (column != 0)? devOligodendrocyteKMinus[thrIdx - 1]: devOligodendrocyteKMinus[thrIdx];
        
        float gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*hx);
        float gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*hx);

        //Diffusion and Chemotaxis Mic

        valIPlus  = (line != xSize-1)? devMicrogliaKMinus[thrIdx + xSize]: devMicrogliaKMinus[thrIdx] - (float)(2*hx*lowerNeumannBC);
        valJPlus  = (column != xSize-1)? devMicrogliaKMinus[thrIdx + 1]: devMicrogliaKMinus[thrIdx] - (float)(2*hx*rightNeumannBC);
        valIMinus = (line != 0)? devMicrogliaKMinus[thrIdx - xSize]: devMicrogliaKMinus[thrIdx] - (float)(2*hx*upperNeumannBC);
        valJMinus = (column != 0)? devMicrogliaKMinus[thrIdx - 1]: devMicrogliaKMinus[thrIdx] - (float)(2*hx*leftNeumannBC);
        
        float microgliaDiffusion = devParams->micDiffusion*CalculateDiffusion(hx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinus[thrIdx]);
        float microgliaChemotaxis = devParams->chi*CalculateChemottaxis(hx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinus[thrIdx],\
        devParams->avgMic, gradientOdcI, gradientOdcJ);

        //Diffusion and Chemotaxis CDC

        valIPlus  = (line != xSize-1)?devConventionalDCKMinus[thrIdx + xSize]:devConventionalDCKMinus[thrIdx] - (float)(2*hx*lowerNeumannBC);
        valJPlus  = (column != xSize-1)?devConventionalDCKMinus[thrIdx + 1]:devConventionalDCKMinus[thrIdx] - (float)(2*hx*rightNeumannBC);
        valIMinus = (line != 0)?devConventionalDCKMinus[thrIdx - xSize]:devConventionalDCKMinus[thrIdx] - (float)(2*hx*upperNeumannBC);
        valJMinus = (column != 0)?devConventionalDCKMinus[thrIdx - 1]:devConventionalDCKMinus[thrIdx] - (float)(2*hx*leftNeumannBC);

        float conventionalDcDiffusion = devParams->cDcDiffusion*CalculateDiffusion(hx, valJPlus, valJMinus, valIPlus, valIMinus,devConventionalDCKMinus[thrIdx]);
        float conventionalDcChemotaxis = devParams->chi*CalculateChemottaxis(hx, valJPlus, valJMinus, valIPlus, valIMinus,devConventionalDCKMinus[thrIdx],\
        devParams->avgDc, gradientOdcI, gradientOdcJ);

        //Difussion and Chemotaxis CD8T

        valIPlus  = (line != xSize-1)? devTCytotoxicKMinus[thrIdx + xSize]: devTCytotoxicKMinus[thrIdx] - (float)(2*hx*lowerNeumannBC);
        valJPlus  = (column != xSize-1)? devTCytotoxicKMinus[thrIdx + 1]: devTCytotoxicKMinus[thrIdx] - (float)(2*hx*rightNeumannBC);
        valIMinus = (line != 0)? devTCytotoxicKMinus[thrIdx - xSize]: devTCytotoxicKMinus[thrIdx] - (float)(2*hx*upperNeumannBC);
        valJMinus = (column != 0)? devTCytotoxicKMinus[thrIdx - 1]: devTCytotoxicKMinus[thrIdx] - (float)(2*hx*leftNeumannBC);

        float tCytotoxicDiffusion = devParams->tCytoDiffusion*CalculateDiffusion(hx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinus[thrIdx]);
        float tCytotoxicChemotaxis = devParams->chi*CalculateChemottaxis(hx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinus[thrIdx],\
        devParams->avgT, gradientOdcI, gradientOdcJ);

        //Difussion ADC

        valIPlus  = (line != xSize-1)? devActivatedDCKMinus[thrIdx + xSize]: devActivatedDCKMinus[thrIdx] - (float)(2*hx*lowerNeumannBC);
        valJPlus  = (column != xSize-1)? devActivatedDCKMinus[thrIdx + 1]: devActivatedDCKMinus[thrIdx] - (float)(2*hx*rightNeumannBC);
        valIMinus = (line != 0)? devActivatedDCKMinus[thrIdx - xSize]: devActivatedDCKMinus[thrIdx] - (float)(2*hx*upperNeumannBC);
        valJMinus = (column != 0)? devActivatedDCKMinus[thrIdx - 1]: devActivatedDCKMinus[thrIdx] - (float)(2*hx*leftNeumannBC);

        float activatedDCDiffusion = devParams->aDcDiffusion*CalculateDiffusion(hx, valJPlus, valJMinus, valIPlus, valIMinus, devActivatedDCKMinus[thrIdx]);

        //Difussion Antibody

        valIPlus  = (line != xSize-1)? devAntibodyKMinus[thrIdx + xSize]: devAntibodyKMinus[thrIdx] - (float)(2*hx*lowerNeumannBC);
        valJPlus  = (column != xSize-1)? devAntibodyKMinus[thrIdx + 1]: devAntibodyKMinus[thrIdx] - (float)(2*hx*rightNeumannBC);
        valIMinus = (line != 0)? devAntibodyKMinus[thrIdx - xSize]: devAntibodyKMinus[thrIdx] - (float)(2*hx*upperNeumannBC);
        valJMinus = (column != 0)? devAntibodyKMinus[thrIdx - 1]: devAntibodyKMinus[thrIdx] - (float)(2*hx*leftNeumannBC);

        float antibodyDiffusion = devParams->antibodyDiffusion*CalculateDiffusion(hx, valJPlus, valJMinus, valIPlus, valIMinus, devAntibodyKMinus[thrIdx]);

        //*******************************************Solving Tissue equations*****************************************************

        //Microglia update
        float microgliaReaction = devParams->muMic*devMicrogliaKMinus[thrIdx]*(devParams->avgMic - devMicrogliaKMinus[thrIdx]);
        float microgliaClearance = devParams->cMic*devMicrogliaKMinus[thrIdx];

        devMicrogliaKPlus[thrIdx] = devMicrogliaKMinus[thrIdx] + \
        ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);
        if((devMicrogliaKPlus[thrIdx]) < 0 || isnanf (devMicrogliaKPlus[thrIdx])){
            printf("Microglia (%f) deu erro no tempo %f\n", devMicrogliaKPlus[thrIdx], kTime*ht);
            exit(0);
        }

        //Conventional DC update
        float conventionalDcReaction = devParams->muCDc*devOligodendrocyteKMinus[thrIdx]*(devParams->avgDc - devConventionalDCKMinus[thrIdx]);
        float conventionalDcActivation = devParams->bD*devConventionalDCKMinus[thrIdx]*devOligodendrocyteKMinus[thrIdx];
        float conventionalDcClearance = devParams->cCDc*devConventionalDCKMinus[thrIdx];

        devConventionalDCKPlus[thrIdx] = devConventionalDCKMinus[thrIdx] + \
        ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);
        if((devConventionalDCKPlus[thrIdx]) < 0 || isnanf (devConventionalDCKPlus[thrIdx])){
            printf("CDC (%f) deu erro no tempo %f\n", devConventionalDCKPlus[thrIdx], kTime*ht);
            exit(0);
        }

        //Activated DC update
        float activatedDcClearance = devParams->cADc*devActivatedDCKMinus[thrIdx];
        float activatedDcMigration = devThetaPV[thrIdx]*devParams->gammaD*(*devActivatedDCLymphNode - devActivatedDCKMinus[thrIdx]);
        
        devActivatedDCKPlus[thrIdx] = devActivatedDCKMinus[thrIdx] + ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
        if((devActivatedDCKPlus[thrIdx]) < 0 || isnanf (devActivatedDCKPlus[thrIdx])){
            printf("ADC (%f) deu erro no tempo %f\n", devActivatedDCKPlus[thrIdx], kTime*ht);
            exit(0);
        }

        //CD8 T update
        float tCytotoxicMigration = devThetaBV[thrIdx]*devParams->gammaT*(*devTCytotoxicLymphNode - devTCytotoxicKMinus[thrIdx]);
        
        devTCytotoxicKPlus[thrIdx] = devTCytotoxicKMinus[thrIdx] + ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
        if((devTCytotoxicKPlus[thrIdx]) < 0 || isnanf (devTCytotoxicKPlus[thrIdx])){
            printf("tCytotoxic (%f) deu erro no tempo %f\n", devTCytotoxicKPlus[thrIdx], kTime*ht);
            exit(0);
        }

        //Antibody update
        float odcAntibodyMicrogliaFagocitosis = devParams->lambAntMic*devAntibodyKMinus[thrIdx]*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx])*fFunc(devMicrogliaKMinus[thrIdx], devParams->avgMic);
        float antibodyMigration = devThetaBV[thrIdx]*devParams->gammaAntibody*(*devAntibodyLymphNode - devAntibodyKMinus[thrIdx]);
        
        devAntibodyKPlus[thrIdx] = devAntibodyKMinus[thrIdx] + ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);
        if((devAntibodyKPlus[thrIdx]) < 0 || isnanf (devAntibodyKPlus[thrIdx])){
            printf("antibody (%.8f) deu erro no tempo %f\n", (devAntibodyKPlus[thrIdx]), kTime*ht);
            exit(0);
        }

        //Oligodendrocytes update
        float odcMicrogliaFagocitosis = devParams->rM*fFunc(devMicrogliaKMinus[thrIdx], devParams->avgMic)*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx]);
        float odcTCytotoxicApoptosis = devParams->rT*fFunc(devTCytotoxicKMinus[thrIdx], devParams->avgT)*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx]);

        devOligodendrocyteKPlus[thrIdx] = devOligodendrocyteKMinus[thrIdx] + ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
        if((devOligodendrocyteKPlus[thrIdx]) < 0 || isnanf (devOligodendrocyteKPlus[thrIdx])){
            printf("oligodendrocyte (%f) deu erro no tempo %f\n", devOligodendrocyteKPlus[thrIdx], kTime*ht);
            exit(0);
        }
        if(devThetaBV[thrIdx] == 1){
            tCytoSumVesselBlock[vesselIdx] += devTCytotoxicKPlus[thrIdx];
            antibodySumVesselBlock[vesselIdx] += devAntibodyKPlus[thrIdx];
        }
        if(devThetaPV[thrIdx] == 1){
            conventionalDCSumVesselBlock[vesselIdx] += devActivatedDCKPlus[thrIdx];
        }
        thrIdx += gridDim.x*blockDim.x;
    }
    __syncthreads();
    int i = blockDim.x / 2;
    while (i != 0) {
        if (vesselIdx < i){
            tCytoSumVesselBlock[vesselIdx] += tCytoSumVesselBlock[vesselIdx + i];
            conventionalDCSumVesselBlock[vesselIdx] += conventionalDCSumVesselBlock[vesselIdx + i];
            antibodySumVesselBlock[vesselIdx] += antibodySumVesselBlock[vesselIdx + i];
        }
        __syncthreads();
        i /= 2;
    }
    if (vesselIdx == 0){
        tCytoSumVessel[blockIdx.x] = tCytoSumVesselBlock[0];
        activatedDCSumVessel[blockIdx.x] = conventionalDCSumVesselBlock[0];
        antibodySumVessel[blockIdx.x] = antibodySumVesselBlock[0];
    }

}

void RunModel(structModel *model){
    //Save IC
    WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
    
    float sumActivatedDCLymphNode, sumAntibodyLymphNode, sumTCytotoxicLymphNode;

    float *devThetaPV, *devThetaBV, *devActivatedDCVessel, *devTCytotoxicVessel, *devAntibodyVessel, *devActivatedDCLymphNode, *devAntibodyLymphNode, *devTCytotoxicLymphNode, *devMicrogliaKMinus, *devMicrogliaKPlus, *devTCytotoxicKMinus, *devTCytotoxicKPlus, *devAntibodyKMinus, *devAntibodyKPlus, *devConventionalDCKMinus, *devConventionalDCKPlus, *devActivatedDCKMinus, *devActivatedDCKPlus, *devOligodendrocytesDCKMinus, *devOligodendrocytesDCKPlus;
    
    cudaMalloc((void**)&devThetaPV, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devThetaBV, model->xSize*model->xSize*sizeof(float));

    cudaMalloc((void**)&devOligodendrocytesDCKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devOligodendrocytesDCKPlus, model->xSize*model->xSize*sizeof(float));

    cudaMalloc((void**)&devMicrogliaKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devMicrogliaKPlus, model->xSize*model->xSize*sizeof(float));
    
    cudaMalloc((void**)&devTCytotoxicKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devTCytotoxicKPlus, model->xSize*model->xSize*sizeof(float));
    
    cudaMalloc((void**)&devAntibodyKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devAntibodyKPlus, model->xSize*model->xSize*sizeof(float));
    
    cudaMalloc((void**)&devConventionalDCKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devConventionalDCKPlus, model->xSize*model->xSize*sizeof(float));
    
    cudaMalloc((void**)&devActivatedDCKMinus, model->xSize*model->xSize*sizeof(float));
    cudaMalloc((void**)&devActivatedDCKPlus, model->xSize*model->xSize*sizeof(float));
    
    cudaMemcpy(devOligodendrocytesDCKMinus, &model->oligodendrocyte[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devMicrogliaKMinus, &model->microglia[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devTCytotoxicKMinus, &model->tCytotoxic[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devAntibodyKMinus, &model->antibody[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devConventionalDCKMinus, &model->conventionalDc[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devActivatedDCKMinus, &model->activatedDc[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    
    structParameters* devParams;
    printf("tamanho parametros = %d", sizeof(devParams));
    //se der errado passar parametro por parametro (tentar com memoria de constantes)
    cudaMalloc((void**)&devParams, sizeof(structParameters));
    cudaMemcpy(devParams, &model->parametersModel, sizeof(structParameters), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&devActivatedDCLymphNode, sizeof(float));
    cudaMalloc((void**)&devAntibodyLymphNode, sizeof(float));
    cudaMalloc((void**)&devTCytotoxicLymphNode, sizeof(float));

    cudaMalloc((void**)&devActivatedDCVessel, numBlocks*sizeof(float));
    cudaMalloc((void**)&devAntibodyVessel, numBlocks*sizeof(float));
    cudaMalloc((void**)&devTCytotoxicVessel, numBlocks*sizeof(float));
    //Inicializar os constant com os valores
    int stepKMinus = 0, stepKPlus;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;

    for(int kTime = 1; kTime <= model->tSize; kTime++){
        auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;
        // solve lymphnode
        SolverLymphNode(model, kTime);

        //copiar LN pra GPU
        cudaMemcpy(devActivatedDCLymphNode, &model->dendriticLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(devAntibodyLymphNode, &model->antibodyLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(devTCytotoxicLymphNode, &model->tCytotoxicLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);        
        stepKPlus = kTime%2;

        if(stepKPlus%2 == 1)
            kernelPDE<<<numBlocks,threadsPerBlock>>>(devParams, kTime, xSize, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKMinus, devMicrogliaKPlus, devTCytotoxicKMinus, devTCytotoxicKPlus, devAntibodyKMinus, devAntibodyKPlus, devConventionalDCKMinus, devConventionalDCKPlus, devActivatedDCKMinus, devActivatedDCKPlus, devOligodendrocytesDCKMinus, devOligodendrocytesDCKPlus);
        else
            kernelPDE<<<numBlocks,threadsPerBlock>>>(devParams, kTime, xSize, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKPlus, devMicrogliaKMinus, devTCytotoxicKPlus, devTCytotoxicKMinus, devAntibodyKPlus, devAntibodyKMinus, devConventionalDCKPlus, devConventionalDCKMinus, devActivatedDCKPlus, devActivatedDCKMinus, devOligodendrocytesDCKPlus, devOligodendrocytesDCKMinus);
        if(kTime%model->intervalFigures == 0 || kTime == model->tSize){
            //Copia tecido para a CPU
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        }
        //Copia do device para o host as integrais do tecido
        model->tCytotoxicTissueVessels = auxTCytotoxicBV/model->parametersModel.V_BV;
        model->antibodyTissueVessels = auxAntibodyBV/model->parametersModel.V_BV;
        model->activatedDCTissueVessels = auxAdcPV/model->parametersModel.V_PV;
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
    }
    printf("Computation Done!!\n");
    // printf("Saving results...\n\n");
    // WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
    // PlotResults(*model);
}