#include "model.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float auxAdcPV = 0.0;
float auxAntibodyBV = 0.0;
float auxTCytotoxicBV = 0.0;

float Min(float a, float b){
    return (a < b) ? a : b;
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
    if(parametersModel.micDiffusion*ht/(hx*hx) < 0.25 && parametersModel.cDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.aDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.tCytoDiffusion*ht/(hx*hx) < 0.25 && parametersModel.chi*ht/hx < 0.5 && parametersModel.chi*ht/(hx*hx) < 0.25)
        return 1;
    return 0;
}

float CalculateMaxGradOdc(structModel model, int stepKMinus){
    float max = 0.0;
    for(int kPos = 0; kPos < model.xSize*model.xSize; kPos++){
        if(model.oligodendrocyte[stepKMinus][kPos] > max)
            max = model.oligodendrocyte[stepKMinus][kPos];
    }    
    return max;
}

float CalculateHt(structModel model, float stepKMinus){
    float maxGradOdc = CalculateMaxGradOdc(model, stepKMinus);
    return Min(model.hx * model.hx / (4 * model.parametersModel.micDiffusion), (4 * model.parametersModel.micDiffusion) / (model.parametersModel.chi * model.parametersModel.chi * maxGradOdc));
}

void WritePopulation(structModel model, float *population, char* fileName, char* bufferTime){
    FILE *file;
    file = fopen(fileName, "w");
    int k = 0;
    if(file != NULL){
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
        printf("Error matrix file\n");
        exit(0);
    }
}

void WritePopulationLymphNode(structModel model, float *population, char* fileName){
    FILE *file;
    file = fopen(fileName, "w");
    if(file != NULL){
        for(int i=0;i<model.numPointsLN;i++){
            fprintf(file, "%f\n", population[i]);
        }
        fclose(file);
    }else{
        printf("Error lymph node file\n");
        exit(0);
    }
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

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float  *activatedDC, int time){
    char buffer[10];
    int day = time;
    
    snprintf(buffer, sizeof(buffer), "%d", day);
    
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
    snprintf(buffer, sizeof(buffer), "%d", model.tSize/model.intervalFigures);
    strcat(command, buffer);
    system(command);
}


float PreventionOverCrowdingTerm(float populationPoint, float avgValue){
    return populationPoint/(populationPoint + avgValue);
}

float UpDownWind(float frontPoint, float rearPoint, float avgValue){
    return PreventionOverCrowdingTerm(frontPoint, avgValue) - PreventionOverCrowdingTerm(rearPoint, avgValue);
}

float CalculateChemotaxis(structModel model, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ){
    float gradientPopulationI, gradientPopulationJ;
    if(gradientOdcI<0)
        gradientPopulationI = UpDownWind(frontIPoint, ijPoint, avgValue)/(float)model.hx;
    else
        gradientPopulationI = UpDownWind(ijPoint, rearIPoint, avgValue)/(float)model.hx;
    if(gradientOdcJ<0)
        gradientPopulationJ = UpDownWind(frontJPoint, ijPoint, avgValue)/(float)model.hx;
    else
        gradientPopulationJ = UpDownWind(ijPoint, rearJPoint, avgValue)/(float)model.hx;
    return gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

float CalculateDiffusion(structModel model, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint){
    return (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(model.hx*model.hx);
}

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
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
    }else{
        printf("dataExecution file not found!\n");
        exit(0);
    }
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

/*
* Lymphnode
*/
float* EquationsLymphNode(structModel* model, int stepPos){
    float* slopes =(float *)malloc(sizeof(float)*6);
    
    int stepKPlus = (stepPos%(2*model->numStepsLN))/model->numStepsLN;
    int stepKMinus = !(stepKPlus && 1);

    float dcLN = model->dendriticLymphNode[stepKMinus];
    float tCytoLN = model->tCytotoxicLymphNode[stepKMinus];
    float tHelperLN = model->tHelperLymphNode[stepKMinus];
    float bCellLN = model->bCellLymphNode[stepKMinus];
    float plasmaCellLN = model->plasmaCellLymphNode[stepKMinus];
    float antibodyLN = model->antibodyLymphNode[stepKMinus];

    //Describe equations

    //Dendritic cell
    float activatedDcMigration = model->parametersModel.gammaD * (model->activatedDCTissueVessels - dcLN) * (float)(model->parametersModel.V_PV/model->parametersModel.V_LN);
    float activatedDcClearance = model->parametersModel.cDl * dcLN;
    slopes[0] = activatedDcMigration - activatedDcClearance;

    //T Cytotoxic
    float tCytoActivation = model->parametersModel.bTCytotoxic * (model->parametersModel.rhoTCytotoxic*tCytoLN*dcLN - tCytoLN*dcLN);
    float tCytoHomeostasis = model->parametersModel.alphaTCytotoxic * (model->parametersModel.stableTCytotoxic - tCytoLN);
    float tCytoMigration = model->parametersModel.gammaT * (tCytoLN - model->tCytotoxicTissueVessels) * (float)(model->parametersModel.V_BV/model->parametersModel.V_LN);
    slopes[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

    //T Helper
    float tHelperActivation = model->parametersModel.bTHelper * (model->parametersModel.rhoTHelper * tHelperLN * dcLN - tHelperLN * dcLN);
    float tHelperHomeostasis = model->parametersModel.alphaTHelper * (model->parametersModel.stableTHelper - tHelperLN);
    float tHelperDispendure = model->parametersModel.bRho * dcLN * tHelperLN * bCellLN;
    slopes[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

    //B Cell
    float bCellActivation = model->parametersModel.bRhoB * (model->parametersModel.rhoB * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
    float bcellHomeostasis = model->parametersModel.alphaB * (model->parametersModel.stableB - bCellLN);
    slopes[3] = bcellHomeostasis + bCellActivation;

    //Plasma Cells
    float plasmaActivation = model->parametersModel.bRhoP * (model->parametersModel.rhoP * tHelperLN * dcLN * bCellLN);
    float plasmaHomeostasis = model->parametersModel.alphaP * (model->parametersModel.stableP - plasmaCellLN);
    slopes[4] = plasmaHomeostasis + plasmaActivation;

    //Antibody
    float antibodyProduction = model->parametersModel.rhoAntibody * plasmaCellLN;
    float antibodyDecayment = model->parametersModel.cF * antibodyLN;
    float antibodyMigration = model->parametersModel.gammaAntibody * (antibodyLN - model->antibodyTissueVessels) * (float)(model->parametersModel.V_BV/model->parametersModel.V_LN);
    slopes[5] = antibodyProduction - antibodyMigration - antibodyDecayment;

    return slopes;
}

/*
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
    
    float* slopeLN;
    slopeLN = EquationsLymphNode(*model, populationLN, stepPos);

    float htLN = model->ht*model->numStepsLN;
    
    //Execute Euler 
    model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + htLN*slopeLN[0];
    model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + htLN*slopeLN[1];
    model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + htLN*slopeLN[2];
    model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + htLN*slopeLN[3];
    model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + htLN*slopeLN[4];
    model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + htLN*slopeLN[5];
    free(slopeLN);

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
*/

derivatives* SlopePDEs(int stepKPlus, float ht, structModel* model){

    derivatives* slopes = calloc(1, sizeof(derivatives));
    slopes->derivativesLymphNode = (float*)calloc(6, sizeof(float));
    slopes->derivativesTissue = (float**)calloc(6, sizeof(float*));
    

    for(int i = 0; i < 6; i++)
        slopes->derivativesTissue[i] = (float*)calloc(model->xSize*model->xSize, sizeof(float));
    /*
     * Solve slope PDEs
    */
    int stepKMinus = !(stepKPlus && 1), line, column;
    
    float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    float diffusionOdc;

    for(int kPos = 0; kPos < model->xSize*model->xSize; kPos++){
        line = (int)kPos/model->xSize;
        column = kPos%model->xSize;
        
        microgliaKMinus = model->microglia[stepKMinus][kPos];
        conventionalDcKMinus = model->conventionalDc[stepKMinus][kPos];
        activatedDcKMinus = model->activatedDc[stepKMinus][kPos];
        tCytotoxicKMinus = model->tCytotoxic[stepKMinus][kPos];
        antibodyKMinus = model->antibody[stepKMinus][kPos];
        oligodendrocyteKMinus = model->oligodendrocyte[stepKMinus][kPos];

        //Define gradient ODCs
            valIPlus = (line != model->xSize-1)? model->oligodendrocyte[stepKMinus][kPos + model->xSize]: model->oligodendrocyte[stepKMinus][kPos - model->xSize];
            valJPlus = (column != model->xSize-1)? model->oligodendrocyte[stepKMinus][kPos + 1]: model->oligodendrocyte[stepKMinus][kPos - 1];
            valIMinus = (line != 0)? model->oligodendrocyte[stepKMinus][kPos - model->xSize]: model->oligodendrocyte[stepKMinus][kPos + model->xSize];
            valJMinus = (column != 0)? model->oligodendrocyte[stepKMinus][kPos - 1]: model->oligodendrocyte[stepKMinus][kPos + 1];

            diffusionOdc = CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->oligodendrocyte[stepKMinus][kPos]);

            gradientOdcI = (float)(valIPlus - valIMinus)/(float)(model->hx);
            gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(model->hx);

            //Diffusion and Chemotaxis Mic

            valIPlus  = (line != model->xSize-1)? model->microglia[stepKMinus][kPos + model->xSize]: model->microglia[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->microglia[stepKMinus][kPos], model->parametersModel.avgMic) * model->hx * gradientOdcI/model->parametersModel.micDiffusion;
            valJPlus  = (column != model->xSize-1)? model->microglia[stepKMinus][kPos + 1]: model->microglia[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->microglia[stepKMinus][kPos], model->parametersModel.avgMic) * model->hx * gradientOdcJ/model->parametersModel.micDiffusion;
            valIMinus = (line != 0)? model->microglia[stepKMinus][kPos - model->xSize]: model->microglia[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->microglia[stepKMinus][kPos], model->parametersModel.avgMic) * model->hx * gradientOdcI/model->parametersModel.micDiffusion;
            valJMinus = (column != 0)? model->microglia[stepKMinus][kPos - 1]: model->microglia[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->microglia[stepKMinus][kPos], model->parametersModel.avgMic) * model->hx * gradientOdcJ/model->parametersModel.micDiffusion;
            
            microgliaDiffusion = model->parametersModel.micDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][kPos]);
            
            microgliaChemotaxis = model->parametersModel.chi * CalculateChemotaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][kPos],\
            model->parametersModel.avgMic, gradientOdcI, gradientOdcJ)\
            + model->parametersModel.chi * diffusionOdc * PreventionOverCrowdingTerm(microgliaKMinus, model->parametersModel.avgMic);

            //Diffusion and Chemotaxis CDC

            valIPlus  = (line != model->xSize-1)? model->conventionalDc[stepKMinus][kPos + model->xSize]: model->conventionalDc[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->conventionalDc[stepKMinus][kPos], model->parametersModel.avgDc) * model->hx * gradientOdcI/model->parametersModel.cDcDiffusion;
            valJPlus  = (column != model->xSize-1)? model->conventionalDc[stepKMinus][kPos + 1]: model->conventionalDc[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->conventionalDc[stepKMinus][kPos], model->parametersModel.avgDc) * model->hx * gradientOdcJ/model->parametersModel.cDcDiffusion;
            valIMinus = (line != 0)? model->conventionalDc[stepKMinus][kPos - model->xSize]: model->conventionalDc[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->conventionalDc[stepKMinus][kPos], model->parametersModel.avgDc) * model->hx * gradientOdcI/model->parametersModel.cDcDiffusion;
            valJMinus = (column != 0)? model->conventionalDc[stepKMinus][kPos - 1]: model->conventionalDc[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->conventionalDc[stepKMinus][kPos], model->parametersModel.avgDc) * model->hx * gradientOdcJ/model->parametersModel.cDcDiffusion;

            conventionalDcDiffusion = model->parametersModel.cDcDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][kPos]);
            
            conventionalDcChemotaxis = model->parametersModel.chi*CalculateChemotaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][kPos],\
            model->parametersModel.avgDc, gradientOdcI, gradientOdcJ)\
            + model->parametersModel.chi * diffusionOdc * PreventionOverCrowdingTerm(conventionalDcKMinus, model->parametersModel.avgDc);

            //Difussion and Chemotaxis CD8T

            valIPlus  = (line != model->xSize-1)? model->tCytotoxic[stepKMinus][kPos + model->xSize]: model->tCytotoxic[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->tCytotoxic[stepKMinus][kPos], model->parametersModel.avgT) * model->hx * gradientOdcI/model->parametersModel.tCytoDiffusion;
            valJPlus  = (column != model->xSize-1)? model->tCytotoxic[stepKMinus][kPos + 1]: model->tCytotoxic[stepKMinus][kPos] + PreventionOverCrowdingTerm(model->tCytotoxic[stepKMinus][kPos], model->parametersModel.avgT) * model->hx * gradientOdcJ/model->parametersModel.tCytoDiffusion;
            valIMinus = (line != 0)? model->tCytotoxic[stepKMinus][kPos - model->xSize]: model->tCytotoxic[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->tCytotoxic[stepKMinus][kPos], model->parametersModel.avgT) * model->hx * gradientOdcI/model->parametersModel.tCytoDiffusion;
            valJMinus = (column != 0)? model->tCytotoxic[stepKMinus][kPos - 1]: model->tCytotoxic[stepKMinus][kPos] - PreventionOverCrowdingTerm(model->tCytotoxic[stepKMinus][kPos], model->parametersModel.avgT) * model->hx * gradientOdcJ/model->parametersModel.tCytoDiffusion;

            tCytotoxicDiffusion = model->parametersModel.tCytoDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][kPos]);
            
            tCytotoxicChemotaxis = model->parametersModel.chi*CalculateChemotaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][kPos],\
            model->parametersModel.avgT, gradientOdcI, gradientOdcJ)\
            + model->parametersModel.chi * diffusionOdc * PreventionOverCrowdingTerm(tCytotoxicKMinus, model->parametersModel.avgT);

            //Difussion ADC

            valIPlus  = (line != model->xSize-1)? model->activatedDc[stepKMinus][kPos + model->xSize]: model->activatedDc[stepKMinus][kPos - model->xSize] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->activatedDc[stepKMinus][kPos + 1]: model->activatedDc[stepKMinus][kPos - 1] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->activatedDc[stepKMinus][kPos - model->xSize]: model->activatedDc[stepKMinus][kPos + model->xSize] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->activatedDc[stepKMinus][kPos - 1]: model->activatedDc[stepKMinus][kPos + 1] - (float)(2*model->hx*leftNeumannBC);

            activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][kPos]);

            //Difussion Antibody

            valIPlus  = (line != model->xSize-1)? model->antibody[stepKMinus][kPos + model->xSize]: model->antibody[stepKMinus][kPos - model->xSize] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->antibody[stepKMinus][kPos + 1]: model->antibody[stepKMinus][kPos - 1] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->antibody[stepKMinus][kPos - model->xSize]: model->antibody[stepKMinus][kPos + model->xSize] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->antibody[stepKMinus][kPos - 1]: model->antibody[stepKMinus][kPos + 1] - (float)(2*model->hx*leftNeumannBC);

            antibodyDiffusion = model->parametersModel.antibodyDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKMinus][kPos]);
            
        //*******************************************Solving Tissue equations*****************************************************

        //Microglia update
        microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
        microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

        slopes->derivativesTissue[0][kPos] = microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance;

        //Conventional DC update
        conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
        conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
        conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

        slopes->derivativesTissue[1][kPos] = conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation;

        //Activated DC update
        activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
        activatedDcMigration = model->thetaPV[kPos]*model->parametersModel.gammaD*(model->dendriticLymphNode[stepKPlus] - activatedDcKMinus);
        
        slopes->derivativesTissue[2][kPos] = activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance;

        //CD8 T update
        tCytotoxicMigration = model->thetaBV[kPos]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus);
        
        slopes->derivativesTissue[3][kPos] = tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration;
        
        //Antibody update
        odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
        antibodyMigration = model->thetaBV[kPos]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[stepKPlus] - antibodyKMinus);
        
        slopes->derivativesTissue[4][kPos] = antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis;
        
        //Oligodendrocytes update
        odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
        odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

        slopes->derivativesTissue[5][kPos] = odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis;
        if(model->thetaBV[kPos] == 1){
            auxTCytotoxicBV += model->tCytotoxic[stepKPlus][kPos];
            auxAntibodyBV += model->antibody[stepKPlus][kPos];
        }
        if(model->thetaPV[kPos] == 1){
            auxAdcPV += model->activatedDc[stepKPlus][kPos];
        }
    }
    return slopes;
}

float Euler(float time, structModel *model, int stepKPlus){
    derivatives* slopeK;
    int stepKMinus = (stepKPlus + 1) % 2;
    float htCalculated = CalculateHt(*model, stepKMinus) * .1; //Multiplicar por uma tolerancia se der errado
    printf("Tempo %f: ht dinamico = %f\n", time, htCalculated);
    
    slopeK = SlopePDEs(stepKPlus, model->ht, model);
    // if(time%model->numStepsLN == 0){
    //     stepKPlus = (time%(2*model->numStepsLN))/model->numStepsLN;
    //     stepKMinus = !(stepKPlus && 1);
    //     // stepKPlus = kTime%2;
    //     // stepKMinus = (stepKPlus + 1) % 2;
    //     model->tCytotoxicTissueVessels = auxTCytotoxicBV * model->hx * model->hx / model->parametersModel.V_BV;
    //     model->antibodyTissueVessels = auxAntibodyBV * model->hx * model->hx / model->parametersModel.V_BV;
    //     model->activatedDCTissueVessels = auxAdcPV * model->hx * model->hx / model->parametersModel.V_PV;
    //     float* slopeLN = EquationsLymphNode(model, time);
    //     slopeK->derivativesLymphNode = slopeLN;

    //     model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[0];
    //     model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[1];
    //     model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[2];
    //     model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[3];
    //     model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[4];
    //     model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + htLymphNode * slopeK->derivativesLymphNode[5];
    // }
    
    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKPlus][spacePoint] = model->microglia[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[0][spacePoint];
        model->conventionalDc[stepKPlus][spacePoint] = model->conventionalDc[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[1][spacePoint];
        model->activatedDc[stepKPlus][spacePoint] = model->activatedDc[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[2][spacePoint];
        model->tCytotoxic[stepKPlus][spacePoint] = model->tCytotoxic[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[3][spacePoint];
        model->antibody[stepKPlus][spacePoint] = model->antibody[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[4][spacePoint];
        model->oligodendrocyte[stepKPlus][spacePoint] = model->oligodendrocyte[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[5][spacePoint];
    }

    free(slopeK);

    // int intervalPoints = (int)(model->tSize/model->numPointsLN);
    // if(time%intervalPoints){
    //     int posSave = time/intervalPoints;
    //     model->dendriticLymphNodeSavedPoints[posSave] = model->dendriticLymphNode[stepKPlus];
    //     model->tCytotoxicLymphNodeSavedPoints[posSave] = model->tCytotoxicLymphNode[stepKPlus];
    //     model->tHelperLymphNodeSavedPoints[posSave] = model->tHelperLymphNode[stepKPlus];
    //     model->bCellLymphNodeSavedPoints[posSave] = model->bCellLymphNode[stepKPlus];
    //     model->plasmaCellLymphNodeSavedPoints[posSave] = model->plasmaCellLymphNode[stepKPlus];
    //     model->antibodyLymphNodeSavedPoints[posSave] = model->antibodyLymphNode[stepKPlus];
    // }

    return htCalculated;

}

void RungeKutta(int kTime, structModel *model){
    derivatives* slopeK1, *slopeK2, *slopeK3 , *slopeK4;
    float htLymphNode = model->ht * model->numStepsLN;

    float** yK;
    yK = (float**) calloc(6, sizeof(float*));
    for(int i = 0; i < 6; i++)
        yK[i] = (float*) calloc(model->xSize * model->xSize, sizeof(float));

    int stepKPlus, stepKMinus;
    if(kTime%model->numStepsLN == 0){
        auxAdcPV = 0.0;
        auxAntibodyBV = 0.0;
        auxTCytotoxicBV = 0.0;
    }
    slopeK1 = SlopePDEs(kTime, model->ht, model);
    if(kTime%model->numStepsLN == 0){
        stepKPlus = (kTime%(2*model->numStepsLN))/model->numStepsLN;
        stepKMinus = !(stepKPlus && 1);
        model->tCytotoxicTissueVessels = auxTCytotoxicBV * model->hx * model->hx / model->parametersModel.V_BV;
        model->antibodyTissueVessels = auxAntibodyBV * model->hx * model->hx / model->parametersModel.V_BV;
        model->activatedDCTissueVessels = auxAdcPV * model->hx * model->hx / model->parametersModel.V_PV;
        float* slopeLN = EquationsLymphNode(model, kTime);
        slopeK1->derivativesLymphNode = slopeLN;

        model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[0];
        model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[1];
        model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[2];
        model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[3];
        model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[4];
        model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + htLymphNode * slopeK1->derivativesLymphNode[5];
    }
    
    stepKPlus = kTime%2;
    stepKMinus = (stepKPlus + 1) % 2;

    //SAlvaro estado setpKminus para usar em cada k

    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        yK[0][spacePoint] = model->microglia[stepKMinus][spacePoint];
        yK[1][spacePoint] = model->conventionalDc[stepKMinus][spacePoint];
        yK[2][spacePoint] = model->activatedDc[stepKMinus][spacePoint];
        yK[3][spacePoint] = model->tCytotoxic[stepKMinus][spacePoint];
        yK[4][spacePoint] = model->antibody[stepKMinus][spacePoint];
        yK[5][spacePoint] = model->oligodendrocyte[stepKMinus][spacePoint];

        model->microglia[stepKMinus][spacePoint] = model->microglia[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[0][spacePoint] / 2;
        model->conventionalDc[stepKMinus][spacePoint] = model->conventionalDc[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[1][spacePoint] / 2;
        model->activatedDc[stepKMinus][spacePoint] = model->activatedDc[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[2][spacePoint] / 2;
        model->tCytotoxic[stepKMinus][spacePoint] = model->tCytotoxic[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[3][spacePoint] / 2;
        model->antibody[stepKMinus][spacePoint] = model->antibody[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[4][spacePoint] / 2;
        model->oligodendrocyte[stepKMinus][spacePoint] = model->oligodendrocyte[stepKMinus][spacePoint] + model->ht * slopeK1->derivativesTissue[5][spacePoint] / 2;
    }

    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKMinus][spacePoint] = yK[0][spacePoint] + model->ht * slopeK1->derivativesTissue[0][spacePoint] / 2;
        model->conventionalDc[stepKMinus][spacePoint] = yK[1][spacePoint] + model->ht * slopeK1->derivativesTissue[1][spacePoint] / 2;
        model->activatedDc[stepKMinus][spacePoint] = yK[2][spacePoint] + model->ht * slopeK1->derivativesTissue[2][spacePoint] / 2;
        model->tCytotoxic[stepKMinus][spacePoint] = yK[3][spacePoint] + model->ht * slopeK1->derivativesTissue[3][spacePoint] / 2;
        model->antibody[stepKMinus][spacePoint] = yK[4][spacePoint] + model->ht * slopeK1->derivativesTissue[4][spacePoint] / 2;
        model->oligodendrocyte[stepKMinus][spacePoint] = yK[5][spacePoint] + model->ht * slopeK1->derivativesTissue[5][spacePoint] / 2;
    }
    
    slopeK2 = SlopePDEs(kTime, model->ht, model);
    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKMinus][spacePoint] = yK[0][spacePoint] + model->ht * slopeK2->derivativesTissue[0][spacePoint] / 2;
        model->conventionalDc[stepKMinus][spacePoint] = yK[1][spacePoint] + model->ht * slopeK2->derivativesTissue[1][spacePoint] / 2;
        model->activatedDc[stepKMinus][spacePoint] = yK[2][spacePoint] + model->ht * slopeK2->derivativesTissue[2][spacePoint] / 2;
        model->tCytotoxic[stepKMinus][spacePoint] = yK[3][spacePoint] + model->ht * slopeK2->derivativesTissue[3][spacePoint] / 2;
        model->antibody[stepKMinus][spacePoint] = yK[4][spacePoint] + model->ht * slopeK2->derivativesTissue[4][spacePoint] / 2;
        model->oligodendrocyte[stepKMinus][spacePoint] = yK[5][spacePoint] + model->ht * slopeK2->derivativesTissue[5][spacePoint] / 2;
    }

    slopeK3 = SlopePDEs(kTime, model->ht, model);
    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKMinus][spacePoint] = yK[0][spacePoint] + model->ht * slopeK3->derivativesTissue[0][spacePoint];
        model->conventionalDc[stepKMinus][spacePoint] = yK[1][spacePoint] + model->ht * slopeK3->derivativesTissue[1][spacePoint];
        model->activatedDc[stepKMinus][spacePoint] = yK[2][spacePoint] + model->ht * slopeK3->derivativesTissue[2][spacePoint];
        model->tCytotoxic[stepKMinus][spacePoint] = yK[3][spacePoint] + model->ht * slopeK3->derivativesTissue[3][spacePoint];
        model->antibody[stepKMinus][spacePoint] = yK[4][spacePoint] + model->ht * slopeK3->derivativesTissue[4][spacePoint];
        model->oligodendrocyte[stepKMinus][spacePoint] = yK[5][spacePoint] + model->ht * slopeK3->derivativesTissue[5][spacePoint];
    }

    slopeK4 = SlopePDEs(kTime, model->ht, model);
    // for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
    //     model->microglia[stepKMinus][spacePoint] = model->microglia[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[0][spacePoint];
    //     model->conventionalDc[stepKMinus][spacePoint] = model->conventionalDc[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[1][spacePoint];
    //     model->activatedDc[stepKMinus][spacePoint] = model->activatedDc[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[2][spacePoint];
    //     model->tCytotoxic[stepKMinus][spacePoint] = model->tCytotoxic[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[3][spacePoint];
    //     model->antibody[stepKMinus][spacePoint] = model->antibody[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[4][spacePoint];
    //     model->oligodendrocyte[stepKMinus][spacePoint] = model->oligodendrocyte[stepKMinus][spacePoint] + model->ht * slopeK4->derivativesTissue[5][spacePoint];
    // }
    
    // for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
    //     model->microglia[stepKMinus][spacePoint] = model->microglia[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[0][spacePoint] / 2;
    //     model->conventionalDc[stepKMinus][spacePoint] = model->conventionalDc[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[1][spacePoint] / 2;
    //     model->activatedDc[stepKMinus][spacePoint] = model->activatedDc[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[2][spacePoint] / 2;
    //     model->tCytotoxic[stepKMinus][spacePoint] = model->tCytotoxic[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[3][spacePoint] / 2;
    //     model->antibody[stepKMinus][spacePoint] = model->antibody[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[4][spacePoint] / 2;
    //     model->oligodendrocyte[stepKMinus][spacePoint] = model->oligodendrocyte[stepKPlus][spacePoint] + model->ht * slopeK1->derivativesTissue[5][spacePoint] / 2;
    // }

    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKPlus][spacePoint] = yK[0][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[0][spacePoint] + 2 * slopeK2->derivativesTissue[0][spacePoint] + 2 * slopeK3->derivativesTissue[0][spacePoint] + slopeK4->derivativesTissue[0][spacePoint]);
        model->conventionalDc[stepKPlus][spacePoint] = yK[1][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[1][spacePoint] + 2 * slopeK2->derivativesTissue[1][spacePoint] + 2 * slopeK3->derivativesTissue[1][spacePoint] + slopeK4->derivativesTissue[1][spacePoint]);
        model->activatedDc[stepKPlus][spacePoint] = yK[2][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[2][spacePoint] + 2 * slopeK2->derivativesTissue[2][spacePoint] + 2 * slopeK3->derivativesTissue[2][spacePoint] + slopeK4->derivativesTissue[2][spacePoint]);
        model->tCytotoxic[stepKPlus][spacePoint] = yK[3][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[3][spacePoint] + 2 * slopeK2->derivativesTissue[3][spacePoint] + 2 * slopeK3->derivativesTissue[3][spacePoint] + slopeK4->derivativesTissue[3][spacePoint]);
        model->antibody[stepKPlus][spacePoint] = yK[4][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[4][spacePoint] + 2 * slopeK2->derivativesTissue[4][spacePoint] + 2 * slopeK3->derivativesTissue[4][spacePoint] + slopeK4->derivativesTissue[4][spacePoint]);
        model->oligodendrocyte[stepKPlus][spacePoint] = yK[5][spacePoint] + (model->ht / 6) * ( slopeK1->derivativesTissue[5][spacePoint] + 2 * slopeK2->derivativesTissue[5][spacePoint] + 2 * slopeK3->derivativesTissue[5][spacePoint] + slopeK4->derivativesTissue[5][spacePoint]);
    }
    
    
    free(slopeK1->derivativesLymphNode);
    free(slopeK2->derivativesLymphNode);
    free(slopeK3->derivativesLymphNode);
    free(slopeK4->derivativesLymphNode);

    free(slopeK1->derivativesTissue[0]);
    free(slopeK1->derivativesTissue[1]);
    free(slopeK1->derivativesTissue[2]);
    free(slopeK1->derivativesTissue[3]);
    free(slopeK1->derivativesTissue[4]);
    free(slopeK1->derivativesTissue[5]);
    
    free(slopeK2->derivativesTissue[0]);
    free(slopeK2->derivativesTissue[1]);
    free(slopeK2->derivativesTissue[2]);
    free(slopeK2->derivativesTissue[3]);
    free(slopeK2->derivativesTissue[4]);
    free(slopeK2->derivativesTissue[5]);

    free(slopeK3->derivativesTissue[0]);
    free(slopeK3->derivativesTissue[1]);
    free(slopeK3->derivativesTissue[2]);
    free(slopeK3->derivativesTissue[3]);
    free(slopeK3->derivativesTissue[4]);
    free(slopeK3->derivativesTissue[5]);

    free(slopeK4->derivativesTissue[0]);
    free(slopeK4->derivativesTissue[1]);
    free(slopeK4->derivativesTissue[2]);
    free(slopeK4->derivativesTissue[3]);
    free(slopeK4->derivativesTissue[4]);
    free(slopeK4->derivativesTissue[5]);
    
    for(int i = 0; i < 6; i++)
        free(yK[i]);
    
    free(yK);
    free(slopeK1);
    free(slopeK2);
    free(slopeK3);
    free(slopeK4);


    int intervalPoints = (int)(model->tSize/model->numPointsLN);
    if(kTime%intervalPoints){
        int posSave = kTime/intervalPoints;
        model->dendriticLymphNodeSavedPoints[posSave] = model->dendriticLymphNode[stepKPlus];
        model->tCytotoxicLymphNodeSavedPoints[posSave] = model->tCytotoxicLymphNode[stepKPlus];
        model->tHelperLymphNodeSavedPoints[posSave] = model->tHelperLymphNode[stepKPlus];
        model->bCellLymphNodeSavedPoints[posSave] = model->bCellLymphNode[stepKPlus];
        model->plasmaCellLymphNodeSavedPoints[posSave] = model->plasmaCellLymphNode[stepKPlus];
        model->antibodyLymphNodeSavedPoints[posSave] = model->antibodyLymphNode[stepKPlus];
    }
}

/*
 * Central Nervous System - PDEs - Finite Differences
 */
void RunModel(structModel *model){

    int stepKPlus = 1, stepKMinus = 0;

    //Save IC
    if(model->saveFigs)
        WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
    
    int kTime = 1; 
    float time = 0.0;
    float htDynamic = 0.0;

    while(time < model->tFinal){
        htDynamic = Euler(time, model, stepKPlus);
        time += htDynamic;
        
        if(model->saveFigs && ( (int)time - kTime)){
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
            printf("%d!!\n", (int)(kTime));
        }
        kTime = time;
        stepKMinus = stepKPlus;
        stepKPlus = !stepKMinus;
        // printf("%f!!\n", time);
    }    
        
    WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
    
    printf("Computation Done!!\n");
    SavingData(*model);
    for(int index=0;index<BUFFER;++index){
        free(model->microglia[index]);
        free(model->oligodendrocyte[index]);
        free(model->tCytotoxic[index]);
        free(model->antibody[index]);
        free(model->conventionalDc[index]);
        free(model->activatedDc[index]);
    }

    free(model->microglia);
    free(model->oligodendrocyte);
    free(model->tCytotoxic);
    free(model->antibody);
    free(model->conventionalDc);
    free(model->activatedDc);
    free(model->thetaPV);
    free(model->thetaBV);

    free(model->dendriticLymphNode);
    free(model->tCytotoxicLymphNode);
    free(model->tHelperLymphNode);
    free(model->antibodyLymphNode);
    free(model->bCellLymphNode);
    free(model->plasmaCellLymphNode);

    if(model->saveFigs){
        printf("Saving results...\n\n");
        WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
        PlotResults(*model);
    }
    free(model->dendriticLymphNodeSavedPoints);
    free(model->tCytotoxicLymphNodeSavedPoints);
    free(model->tHelperLymphNodeSavedPoints);
    free(model->antibodyLymphNodeSavedPoints);
    free(model->bCellLymphNodeSavedPoints);
    free(model->plasmaCellLymphNodeSavedPoints);
}