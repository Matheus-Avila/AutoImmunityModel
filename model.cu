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

__device__ void AdvectionTerm(float populationPoint, float avgValue, float* result){
    *result = populationPoint/(populationPoint + avgValue);
}

__device__ void UpDownWind(float frontPoint, float rearPoint, float avgValue, float* result){
    float resultF;
    AdvectionTerm(frontPoint, avgValue, result);
    AdvectionTerm(rearPoint, avgValue, &resultF);
    *result = *result - resultF;
}

__device__ void CalculateChemottaxis(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ, float* result){
    float gradientPopulationI, gradientPopulationJ;
    if(gradientOdcI<0){
        UpDownWind(frontIPoint, ijPoint, avgValue, &gradientPopulationI);
        gradientPopulationI = gradientPopulationI/(float)hx;
    }
    else{
        UpDownWind(ijPoint, rearIPoint, avgValue, &gradientPopulationI);
        gradientPopulationI = gradientPopulationI/(float)hx;
    }
    if(gradientOdcJ<0){
        UpDownWind(frontJPoint, ijPoint, avgValue, &gradientPopulationJ);
        gradientPopulationJ = gradientPopulationJ/(float)hx;
    }
    else{
        UpDownWind(ijPoint, rearJPoint, avgValue, &gradientPopulationJ);
        gradientPopulationJ = gradientPopulationJ/(float)hx;
    }

    *result = gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

__device__ void CalculateDiffusion(float hx, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float* result){
    *result = (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(hx*hx);
}

__device__ void fFunc(float valuePopulation, float avgPopulation, float *result){
    *result = valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
}

void WriteBVPV(structModel *model, float *thetaBV, float *thetaPV){
    FILE *fileBV;
    fileBV = fopen("./result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./result/pv.txt", "w");
    for(int k = 0; k < model->xSize*model->xSize; k++){
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
    model.oligodendrocyte[0][0] = 1;
    model.oligodendrocyte[0][1] = 2;
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
}

__constant__ float upperNeumannBC, lowerNeumannBC, leftNeumannBC, rightNeumannBC, constHx, constHt;
__constant__ int constXSize;
const int threadsPerBlock = 256;
const int numBlocks = 128;

__global__ void kernelPDE(structParameters* devParams, int kTime, float* tCytoSumVessel, float* activatedDCSumVessel, float* antibodySumVessel, float *devActivatedDCLymphNode, float *devAntibodyLymphNode, float *devTCytotoxicLymphNode, float *devThetaPV, float *devThetaBV, float *devMicrogliaKMinus, float *devMicrogliaKPlus, float *devTCytotoxicKMinus, float *devTCytotoxicKPlus, float *devAntibodyKMinus, float *devAntibodyKPlus, float *devConventionalDCKMinus, float *devConventionalDCKPlus, float *devActivatedDCKMinus, float *devActivatedDCKPlus, float *devOligodendrocyteKMinus, float *devOligodendrocyteKPlus){
    int thrIdx = blockIdx.x*blockDim.x + threadIdx.x;
    int vesselIdx = threadIdx.x;
    int line = (int)thrIdx/constXSize;
    int column = thrIdx%constXSize;

    __shared__ float tCytoSumVesselBlock[threadsPerBlock];
    __shared__ float conventionalDCSumVesselBlock[threadsPerBlock];
    __shared__ float antibodySumVesselBlock[threadsPerBlock];
    while(thrIdx < constXSize*constXSize){
        line = (int)thrIdx/constXSize;
        column = thrIdx%constXSize;

        //Define gradient ODCs
        float valIPlus = (line != constXSize-1)? devOligodendrocyteKMinus[thrIdx + constXSize]: devOligodendrocyteKMinus[thrIdx];
        float valJPlus = (column != constXSize-1)? devOligodendrocyteKMinus[thrIdx + 1]: devOligodendrocyteKMinus[thrIdx];
        float valIMinus = (line != 0)? devOligodendrocyteKMinus[thrIdx - constXSize]: devOligodendrocyteKMinus[thrIdx];
        float valJMinus = (column != 0)? devOligodendrocyteKMinus[thrIdx - 1]: devOligodendrocyteKMinus[thrIdx];
        
        float gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*constHx);
        float gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*constHx);

        //Diffusion and Chemotaxis Mic

        valIPlus  = (line != constXSize-1)? devMicrogliaKMinus[thrIdx + constXSize]: devMicrogliaKMinus[thrIdx] - (float)(2*constHx*lowerNeumannBC);
        valJPlus  = (column != constXSize-1)? devMicrogliaKMinus[thrIdx + 1]: devMicrogliaKMinus[thrIdx] - (float)(2*constHx*rightNeumannBC);
        valIMinus = (line != 0)? devMicrogliaKMinus[thrIdx - constXSize]: devMicrogliaKMinus[thrIdx] - (float)(2*constHx*upperNeumannBC);
        valJMinus = (column != 0)? devMicrogliaKMinus[thrIdx - 1]: devMicrogliaKMinus[thrIdx] - (float)(2*constHx*leftNeumannBC);
        
        float microgliaDiffusion;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinus[thrIdx], &microgliaDiffusion);
        float microgliaChemotaxis;
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devMicrogliaKMinus[thrIdx],\
        devParams->avgMic, gradientOdcI, gradientOdcJ, &microgliaChemotaxis);
        microgliaChemotaxis *= devParams->chi;
        microgliaDiffusion *= devParams->micDiffusion;
        //Diffusion and Chemotaxis CDC

        valIPlus  = (line != constXSize-1)?devConventionalDCKMinus[thrIdx + constXSize]:devConventionalDCKMinus[thrIdx] - (float)(2*constHx*lowerNeumannBC);
        valJPlus  = (column != constXSize-1)?devConventionalDCKMinus[thrIdx + 1]:devConventionalDCKMinus[thrIdx] - (float)(2*constHx*rightNeumannBC);
        valIMinus = (line != 0)?devConventionalDCKMinus[thrIdx - constXSize]:devConventionalDCKMinus[thrIdx] - (float)(2*constHx*upperNeumannBC);
        valJMinus = (column != 0)?devConventionalDCKMinus[thrIdx - 1]:devConventionalDCKMinus[thrIdx] - (float)(2*constHx*leftNeumannBC);

        float conventionalDcDiffusion;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus,devConventionalDCKMinus[thrIdx], &conventionalDcDiffusion);
        float conventionalDcChemotaxis;
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus,devConventionalDCKMinus[thrIdx],\
        devParams->avgDc, gradientOdcI, gradientOdcJ, &conventionalDcChemotaxis);
        conventionalDcChemotaxis *= devParams->chi;
        conventionalDcDiffusion *= devParams->cDcDiffusion;

        //Difussion and Chemotaxis CD8T

        valIPlus  = (line != constXSize-1)? devTCytotoxicKMinus[thrIdx + constXSize]: devTCytotoxicKMinus[thrIdx] - (float)(2*constHx*lowerNeumannBC);
        valJPlus  = (column != constXSize-1)? devTCytotoxicKMinus[thrIdx + 1]: devTCytotoxicKMinus[thrIdx] - (float)(2*constHx*rightNeumannBC);
        valIMinus = (line != 0)? devTCytotoxicKMinus[thrIdx - constXSize]: devTCytotoxicKMinus[thrIdx] - (float)(2*constHx*upperNeumannBC);
        valJMinus = (column != 0)? devTCytotoxicKMinus[thrIdx - 1]: devTCytotoxicKMinus[thrIdx] - (float)(2*constHx*leftNeumannBC);

        float tCytotoxicDiffusion;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinus[thrIdx], &tCytotoxicDiffusion);
        float tCytotoxicChemotaxis;
        CalculateChemottaxis(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devTCytotoxicKMinus[thrIdx],\
        devParams->avgT, gradientOdcI, gradientOdcJ, &tCytotoxicChemotaxis);
        tCytotoxicChemotaxis *= devParams->chi;
        tCytotoxicDiffusion *= devParams->tCytoDiffusion;

        //Difussion ADC

        valIPlus  = (line != constXSize-1)? devActivatedDCKMinus[thrIdx + constXSize]: devActivatedDCKMinus[thrIdx] - (float)(2*constHx*lowerNeumannBC);
        valJPlus  = (column != constXSize-1)? devActivatedDCKMinus[thrIdx + 1]: devActivatedDCKMinus[thrIdx] - (float)(2*constHx*rightNeumannBC);
        valIMinus = (line != 0)? devActivatedDCKMinus[thrIdx - constXSize]: devActivatedDCKMinus[thrIdx] - (float)(2*constHx*upperNeumannBC);
        valJMinus = (column != 0)? devActivatedDCKMinus[thrIdx - 1]: devActivatedDCKMinus[thrIdx] - (float)(2*constHx*leftNeumannBC);

        float activatedDCDiffusion;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devActivatedDCKMinus[thrIdx], &activatedDCDiffusion);
        activatedDCDiffusion *= devParams->aDcDiffusion;

        //Difussion Antibody

        valIPlus  = (line != constXSize-1)? devAntibodyKMinus[thrIdx + constXSize]: devAntibodyKMinus[thrIdx] - (float)(2*constHx*lowerNeumannBC);
        valJPlus  = (column != constXSize-1)? devAntibodyKMinus[thrIdx + 1]: devAntibodyKMinus[thrIdx] - (float)(2*constHx*rightNeumannBC);
        valIMinus = (line != 0)? devAntibodyKMinus[thrIdx - constXSize]: devAntibodyKMinus[thrIdx] - (float)(2*constHx*upperNeumannBC);
        valJMinus = (column != 0)? devAntibodyKMinus[thrIdx - 1]: devAntibodyKMinus[thrIdx] - (float)(2*constHx*leftNeumannBC);

        float antibodyDiffusion;
        CalculateDiffusion(constHx, valJPlus, valJMinus, valIPlus, valIMinus, devAntibodyKMinus[thrIdx], &antibodyDiffusion);
        antibodyDiffusion *= devParams->antibodyDiffusion;

        //*******************************************Solving Tissue equations*****************************************************

        //Microglia update
        float microgliaReaction = devParams->muMic*devMicrogliaKMinus[thrIdx]*(devParams->avgMic - devMicrogliaKMinus[thrIdx]);
        float microgliaClearance = devParams->cMic*devMicrogliaKMinus[thrIdx];

        devMicrogliaKPlus[thrIdx] = devMicrogliaKMinus[thrIdx] + \
        constHt*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

        //Conventional DC update
        float conventionalDcReaction = devParams->muCDc*devOligodendrocyteKMinus[thrIdx]*(devParams->avgDc - devConventionalDCKMinus[thrIdx]);
        float conventionalDcActivation = devParams->bD*devConventionalDCKMinus[thrIdx]*devOligodendrocyteKMinus[thrIdx];
        float conventionalDcClearance = devParams->cCDc*devConventionalDCKMinus[thrIdx];

        devConventionalDCKPlus[thrIdx] = devConventionalDCKMinus[thrIdx] + \
        constHt*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);

        //Activated DC update
        float activatedDcClearance = devParams->cADc*devActivatedDCKMinus[thrIdx];
        float activatedDcMigration = devThetaPV[thrIdx]*devParams->gammaD*(*devActivatedDCLymphNode - devActivatedDCKMinus[thrIdx]);
        
        devActivatedDCKPlus[thrIdx] = devActivatedDCKMinus[thrIdx] + constHt*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);

        //CD8 T update
        float tCytotoxicMigration = devThetaBV[thrIdx]*devParams->gammaT*(*devTCytotoxicLymphNode - devTCytotoxicKMinus[thrIdx]);
        
        devTCytotoxicKPlus[thrIdx] = devTCytotoxicKMinus[thrIdx] + constHt*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);

        //Antibody update
        float result;
        fFunc(devMicrogliaKMinus[thrIdx], devParams->avgMic, &result);
        float odcAntibodyMicrogliaFagocitosis = devParams->lambAntMic*devAntibodyKMinus[thrIdx]*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx])*result;
        float antibodyMigration = devThetaBV[thrIdx]*devParams->gammaAntibody*(*devAntibodyLymphNode - devAntibodyKMinus[thrIdx]);
        
        devAntibodyKPlus[thrIdx] = devAntibodyKMinus[thrIdx] + constHt*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);

        //Oligodendrocytes update
        float result1;
        fFunc(devMicrogliaKMinus[thrIdx], devParams->avgMic, &result);
        fFunc(devTCytotoxicKMinus[thrIdx], devParams->avgT, &result1);
        float odcMicrogliaFagocitosis = devParams->rM*result*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx]);
        float odcTCytotoxicApoptosis = devParams->rT*result1*(devParams->avgOdc - devOligodendrocyteKMinus[thrIdx]);

        devOligodendrocyteKPlus[thrIdx] = devOligodendrocyteKMinus[thrIdx] + constHt*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
        // if(devOligodendrocyteKPlus[thrIdx] == 1)
        //     printf("%f: %d - %f\n", kTime*constHt, thrIdx, devOligodendrocyteKPlus[thrIdx]);        
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
    
    float *activatedDCVessel, *tCytotoxicVessel, *antibodyVessel;

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
    
    cudaMemcpy(devOligodendrocytesDCKMinus, model->oligodendrocyte[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devMicrogliaKMinus, model->microglia[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devTCytotoxicKMinus, model->tCytotoxic[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devAntibodyKMinus, model->antibody[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devConventionalDCKMinus, model->conventionalDc[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devActivatedDCKMinus, model->activatedDc[0], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
    
    structParameters* devParams;
    
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
    
    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;

    float bc = 0.0;
    
    cudaMemcpyToSymbol(upperNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(lowerNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(leftNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(rightNeumannBC, &bc, sizeof(float));
    cudaMemcpyToSymbol(constHt, &model->ht, sizeof(float));
    cudaMemcpyToSymbol(constHx, &model->hx, sizeof(float));
    cudaMemcpyToSymbol(constXSize, &model->xSize, sizeof(float));
    for(int kTime = 1; kTime <= model->tSize; kTime++){
        auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;
        // solve lymphnode
        SolverLymphNode(model, kTime);
        //copiar LN pra GPU
        cudaMemcpy(&devActivatedDCLymphNode, &model->dendriticLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(&devAntibodyLymphNode, &model->antibodyLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(&devTCytotoxicLymphNode, &model->tCytotoxicLymphNode[stepKPlus], sizeof(float), cudaMemcpyHostToDevice);        
        stepKPlus = kTime%2;

        if(stepKPlus%2 == 1)
            kernelPDE<<<numBlocks,threadsPerBlock>>>(devParams, kTime, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKMinus, devMicrogliaKPlus, devTCytotoxicKMinus, devTCytotoxicKPlus, devAntibodyKMinus, devAntibodyKPlus, devConventionalDCKMinus, devConventionalDCKPlus, devActivatedDCKMinus, devActivatedDCKPlus, devOligodendrocytesDCKMinus, devOligodendrocytesDCKPlus);
        else
            kernelPDE<<<numBlocks,threadsPerBlock>>>(devParams, kTime, devTCytotoxicVessel, devActivatedDCVessel, devAntibodyVessel, devActivatedDCLymphNode, devAntibodyLymphNode, devTCytotoxicLymphNode, devThetaPV, devThetaBV, devMicrogliaKPlus, devMicrogliaKMinus, devTCytotoxicKPlus, devTCytotoxicKMinus, devAntibodyKPlus, devAntibodyKMinus, devConventionalDCKPlus, devConventionalDCKMinus, devActivatedDCKPlus, devActivatedDCKMinus, devOligodendrocytesDCKPlus, devOligodendrocytesDCKMinus);
        
        //Copia tecido para a CPU
        
        // if(stepKPlus%2 == 1){
        //     cudaMemcpy(&model->oligodendrocyte[stepKPlus], &devOligodendrocytesDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->microglia[stepKPlus], &devMicrogliaKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->tCytotoxic[stepKPlus], &devTCytotoxicKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->antibody[stepKPlus], &devAntibodyKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->conventionalDc[stepKPlus], &devConventionalDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->activatedDc[stepKPlus], &devActivatedDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);

        //     cudaMemcpy(&devOligodendrocytesDCKPlus, &model->oligodendrocyte[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devMicrogliaKPlus, &model->microglia[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devTCytotoxicKPlus, &model->tCytotoxic[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devAntibodyKPlus, &model->antibody[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devConventionalDCKPlus, &model->conventionalDc[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devActivatedDCKPlus, &model->activatedDc[stepKPlus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        // }else{
        //     cudaMemcpy(&model->oligodendrocyte[stepKMinus], &devOligodendrocytesDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->microglia[stepKMinus], &devMicrogliaKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->tCytotoxic[stepKMinus], &devTCytotoxicKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->antibody[stepKMinus], &devAntibodyKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->conventionalDc[stepKMinus], &devConventionalDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
        //     cudaMemcpy(&model->activatedDc[stepKMinus], &devActivatedDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);

        //     cudaMemcpy(&devOligodendrocytesDCKMinus, &model->oligodendrocyte[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devMicrogliaKMinus, &model->microglia[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devTCytotoxicKMinus, &model->tCytotoxic[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devAntibodyKMinus, &model->antibody[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devConventionalDCKMinus, &model->conventionalDc[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        //     cudaMemcpy(&devActivatedDCKMinus, &model->activatedDc[stepKMinus], model->xSize*model->xSize*sizeof(float), cudaMemcpyHostToDevice);
        // }
        
        if(kTime%model->intervalFigures == 0 || kTime == model->tSize){
            if(stepKPlus%2 == 1){
                cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->antibody[stepKPlus], devAntibodyKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKPlus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);

            }else{
                cudaMemcpy(model->oligodendrocyte[stepKPlus], devOligodendrocytesDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->microglia[stepKPlus], devMicrogliaKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->tCytotoxic[stepKPlus], devTCytotoxicKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->antibody[stepKPlus], devAntibodyKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->conventionalDc[stepKPlus], devConventionalDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(model->activatedDc[stepKPlus], devActivatedDCKMinus, model->xSize*model->xSize*sizeof(float), cudaMemcpyDeviceToHost);
            }
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        }
        //Copia do device para o host as integrais do tecido
        activatedDCVessel = (float*)calloc(numBlocks, sizeof(float));
        antibodyVessel = (float*)calloc(numBlocks, sizeof(float));
        tCytotoxicVessel = (float*)calloc(numBlocks, sizeof(float));
        
        cudaMemcpy(&activatedDCVessel, &devActivatedDCVessel, numBlocks*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&antibodyVessel, &devAntibodyVessel, numBlocks*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&tCytotoxicVessel, &devTCytotoxicVessel, numBlocks*sizeof(float), cudaMemcpyDeviceToHost);

        for(int pos = 0; pos < numBlocks; pos++){
            auxAdcPV += activatedDCVessel[pos];
            auxAntibodyBV += antibodyVessel[pos];
            auxTCytotoxicBV += tCytotoxicVessel[pos];
        }
        model->activatedDCTissueVessels = auxAdcPV/model->parametersModel.V_PV;
        model->antibodyTissueVessels = auxAntibodyBV/model->parametersModel.V_BV;
        model->tCytotoxicTissueVessels = auxTCytotoxicBV/model->parametersModel.V_BV;

        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
    }
    printf("Computation Done!!\n");
    printf("Saving results...\n\n");
    WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
    PlotResults(*model);
}