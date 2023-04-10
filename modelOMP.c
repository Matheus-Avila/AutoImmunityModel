#include "modelOMP.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>

void InitialConditionTissueMicroglia(structModel* model){
    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        if(pow((i-(int)(0.2*model->xSize)),2) + pow((j-(int)(0.2*model->xSize)),2) < 5 / (model->hx * model->hx)){
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
    if(parametersModel.micDiffusion*ht/(hx*hx) < 0.25 && parametersModel.cDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.aDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.tCytoDiffusion*ht/(hx*hx) < 0.25 && parametersModel.chi*ht/hx < 0.5)
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
    
    char pathOligodendrocytes[70];
    getcwd(pathOligodendrocytes, sizeof(pathOligodendrocytes));
    strcat(pathOligodendrocytes, "/result/matrix/oligo");
    strcat(pathOligodendrocytes, buffer);
    strcat(pathOligodendrocytes, ".txt");
    WritePopulation(model, oligodendrocyte, pathOligodendrocytes, buffer);

    char pathMicroglia[70];
    getcwd(pathMicroglia, sizeof(pathMicroglia));
    strcat(pathMicroglia, "/result/matrix/microglia");
    strcat(pathMicroglia, buffer);
    strcat(pathMicroglia, ".txt");
    WritePopulation(model, microglia, pathMicroglia, buffer);

    char pathTCyto[70];
    getcwd(pathTCyto, sizeof(pathTCyto));
    strcat(pathTCyto, "/result/matrix/tCyto");
    strcat(pathTCyto, buffer);
    strcat(pathTCyto, ".txt");
    WritePopulation(model, tCytotoxic, pathTCyto, buffer);

    char pathAntibody[70];
    getcwd(pathAntibody, sizeof(pathAntibody));
    strcat(pathAntibody, "/result/matrix/antibody");
    strcat(pathAntibody, buffer);
    strcat(pathAntibody, ".txt");
    WritePopulation(model, antibody, pathAntibody, buffer);

    char pathConventionalDC[70];
    getcwd(pathConventionalDC, sizeof(pathConventionalDC));
    strcat(pathConventionalDC, "/result/matrix/conventionalDC");
    strcat(pathConventionalDC, buffer);
    strcat(pathConventionalDC, ".txt");
    WritePopulation(model, conventionalDC, pathConventionalDC, buffer);

    char pathActivatedDC[70];
    getcwd(pathActivatedDC, sizeof(pathActivatedDC));
    strcat(pathActivatedDC, "/result/matrix/activatedDC");
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

float CalculateChemottaxis(structModel model, float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
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
    model->parametersModel.V_BV = 160;//model->parametersModel.V_BV * model->hx * model->hx;
    model->parametersModel.V_PV = 160;//model->parametersModel.V_PV * model->hx * model->hx;
    printf("bv = %f, pv = %f \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    WriteBVPV(model, model->thetaBV, model->thetaPV);
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

structModel ModelInitialize(structParameters params, int totThr, float ht, float hx, float time, float space, int numFigs, int numPointsLN){
    structModel model;
    srand(2);
    model.parametersModel = params;
    if(!VerifyCFL(model.parametersModel, ht, hx)){
        printf("Falhou CFL!!\n");
        exit(0);
    }        
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
    
    model.totalThreads = totThr;
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
    result[5] = antibodyProduction - antibodyMigration - antibodyDecayment;

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

void RunModel(structModel *model){
    //Save IC
    WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
    
    int stepKMinus = 0, stepKPlus, line, column;

    float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;

    int kTime;
    int tid;
    #pragma omp parallel num_threads(model->totalThreads) default(none) shared(model, auxTCytotoxicBV, auxAntibodyBV, auxAdcPV, lowerNeumannBC, rightNeumannBC, upperNeumannBC, leftNeumannBC)\
            private(tid, kTime, stepKMinus, stepKPlus, line, column, microgliaKMinus, conventionalDcKMinus,\
            activatedDcKMinus, tCytotoxicKMinus, antibodyKMinus, oligodendrocyteKMinus, valIPlus, valJPlus, valIMinus, valJMinus, gradientOdcI,\
            gradientOdcJ, microgliaDiffusion, microgliaChemotaxis, conventionalDcDiffusion, conventionalDcChemotaxis, tCytotoxicDiffusion, tCytotoxicChemotaxis,\
            activatedDCDiffusion, antibodyDiffusion, microgliaReaction, microgliaClearance, conventionalDcReaction, conventionalDcActivation, conventionalDcClearance,\
            activatedDcClearance, activatedDcMigration, tCytotoxicMigration, odcAntibodyMicrogliaFagocitosis, antibodyMigration, odcMicrogliaFagocitosis, odcTCytotoxicApoptosis)
    for(kTime = 1; kTime <= model->tSize; kTime++){
        tid = omp_get_thread_num();
        
        model->tCytotoxicTissueVessels = auxTCytotoxicBV * model->hx * model->hx / model->parametersModel.V_BV;
        model->antibodyTissueVessels = auxAntibodyBV * model->hx * model->hx / model->parametersModel.V_BV;
        model->activatedDCTissueVessels = auxAdcPV * model->hx * model->hx / model->parametersModel.V_PV;
   
        if(tid == 0){
            auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;
            SolverLymphNode(model, kTime);
        }
        stepKPlus = kTime%2;
        #pragma omp for reduction(+:auxTCytotoxicBV, auxAntibodyBV, auxAdcPV)
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
            valIPlus = (line != model->xSize-1)? model->oligodendrocyte[stepKMinus][kPos + model->xSize]: model->oligodendrocyte[stepKMinus][kPos];
            valJPlus = (column != model->xSize-1)? model->oligodendrocyte[stepKMinus][kPos + 1]: model->oligodendrocyte[stepKMinus][kPos];
            valIMinus = (line != 0)? model->oligodendrocyte[stepKMinus][kPos - model->xSize]: model->oligodendrocyte[stepKMinus][kPos];
            valJMinus = (column != 0)? model->oligodendrocyte[stepKMinus][kPos - 1]: model->oligodendrocyte[stepKMinus][kPos];
            
            gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*model->hx);
            gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*model->hx);

            //Diffusion and Chemotaxis Mic

            valIPlus  = (line != model->xSize-1)? model->microglia[stepKMinus][kPos + model->xSize]: model->microglia[stepKMinus][kPos] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->microglia[stepKMinus][kPos + 1]: model->microglia[stepKMinus][kPos] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->microglia[stepKMinus][kPos - model->xSize]: model->microglia[stepKMinus][kPos] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->microglia[stepKMinus][kPos - 1]: model->microglia[stepKMinus][kPos] - (float)(2*model->hx*leftNeumannBC);
            
            microgliaDiffusion = model->parametersModel.micDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][kPos]);
            microgliaChemotaxis = model->parametersModel.chi*CalculateChemottaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][kPos],\
            model->parametersModel.avgMic, gradientOdcI, gradientOdcJ);

            //Diffusion and Chemotaxis CDC

            valIPlus  = (line != model->xSize-1)? model->conventionalDc[stepKMinus][kPos + model->xSize]: model->conventionalDc[stepKMinus][kPos] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->conventionalDc[stepKMinus][kPos + 1]: model->conventionalDc[stepKMinus][kPos] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->conventionalDc[stepKMinus][kPos - model->xSize]: model->conventionalDc[stepKMinus][kPos] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->conventionalDc[stepKMinus][kPos - 1]: model->conventionalDc[stepKMinus][kPos] - (float)(2*model->hx*leftNeumannBC);

            conventionalDcDiffusion = model->parametersModel.cDcDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][kPos]);
            conventionalDcChemotaxis = model->parametersModel.chi*CalculateChemottaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][kPos],\
            model->parametersModel.avgDc, gradientOdcI, gradientOdcJ);

            //Difussion and Chemotaxis CD8T

            valIPlus  = (line != model->xSize-1)? model->tCytotoxic[stepKMinus][kPos + model->xSize]: model->tCytotoxic[stepKMinus][kPos] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->tCytotoxic[stepKMinus][kPos + 1]: model->tCytotoxic[stepKMinus][kPos] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->tCytotoxic[stepKMinus][kPos - model->xSize]: model->tCytotoxic[stepKMinus][kPos] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->tCytotoxic[stepKMinus][kPos - 1]: model->tCytotoxic[stepKMinus][kPos] - (float)(2*model->hx*leftNeumannBC);

            tCytotoxicDiffusion = model->parametersModel.tCytoDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][kPos]);
            tCytotoxicChemotaxis = model->parametersModel.chi*CalculateChemottaxis(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][kPos],\
            model->parametersModel.avgT, gradientOdcI, gradientOdcJ);

            //Difussion ADC

            valIPlus  = (line != model->xSize-1)? model->activatedDc[stepKMinus][kPos + model->xSize]: model->activatedDc[stepKMinus][kPos] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->activatedDc[stepKMinus][kPos + 1]: model->activatedDc[stepKMinus][kPos] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->activatedDc[stepKMinus][kPos - model->xSize]: model->activatedDc[stepKMinus][kPos] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->activatedDc[stepKMinus][kPos - 1]: model->activatedDc[stepKMinus][kPos] - (float)(2*model->hx*leftNeumannBC);

            activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][kPos]);

            //Difussion Antibody

            valIPlus  = (line != model->xSize-1)? model->antibody[stepKMinus][kPos + model->xSize]: model->antibody[stepKMinus][kPos] - (float)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->antibody[stepKMinus][kPos + 1]: model->antibody[stepKMinus][kPos] - (float)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->antibody[stepKMinus][kPos - model->xSize]: model->antibody[stepKMinus][kPos] - (float)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->antibody[stepKMinus][kPos - 1]: model->antibody[stepKMinus][kPos] - (float)(2*model->hx*leftNeumannBC);

            antibodyDiffusion = model->parametersModel.antibodyDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKMinus][kPos]);

            //*******************************************Solving Tissue equations*****************************************************

            //Microglia update
            microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
            microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

            model->microglia[stepKPlus][kPos] = microgliaKMinus + \
            model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

            //Conventional DC update
            conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
            conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
            conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

            model->conventionalDc[stepKPlus][kPos] = conventionalDcKMinus + \
            model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);

            //Activated DC update
            activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
            activatedDcMigration = model->thetaPV[kPos]*model->parametersModel.gammaD*(model->dendriticLymphNode[stepKPlus] - activatedDcKMinus);
            
            model->activatedDc[stepKPlus][kPos] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);

            //CD8 T update
            tCytotoxicMigration = model->thetaBV[kPos]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus);
            
            model->tCytotoxic[stepKPlus][kPos] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
            
            //Antibody update
            odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
            antibodyMigration = model->thetaBV[kPos]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[stepKPlus] - antibodyKMinus);
            
            model->antibody[stepKPlus][kPos] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);
            
            //Oligodendrocytes update
            odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
            odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

            model->oligodendrocyte[stepKPlus][kPos] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
            if(model->thetaBV[kPos] == 1){
                auxTCytotoxicBV += model->tCytotoxic[stepKPlus][kPos];
                auxAntibodyBV += model->antibody[stepKPlus][kPos];
            }
            if(model->thetaPV[kPos] == 1){
                auxAdcPV += model->activatedDc[stepKPlus][kPos];
            }
        }
        if(kTime%model->intervalFigures == 0 || kTime == model->tSize)
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;

    }
        
    printf("Computation Done!!\n");

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

    printf("Saving results...\n\n");
    WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
    
    free(model->dendriticLymphNodeSavedPoints);
    free(model->tCytotoxicLymphNodeSavedPoints);
    free(model->tHelperLymphNodeSavedPoints);
    free(model->antibodyLymphNodeSavedPoints);
    free(model->bCellLymphNodeSavedPoints);
    free(model->plasmaCellLymphNodeSavedPoints);

    PlotResults(*model);
}