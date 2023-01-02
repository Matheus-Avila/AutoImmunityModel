#include "modelMPI.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "log.h"

#define DEBUG 1

void InitialConditionTissueMicroglia(structModel* model){
    for(int i = 0; i < model->xSize; i++){
        for(int j = 0; j < model->xSize; j++){
            if(pow((i-(int)(model->xSize/2)),2) + pow((j-(int)(model->xSize/2)),2) < 5){
                model->microglia[0][i * model->xSize + j] = (float)model->parametersModel.avgMic/3;
            }
        }
    }
}

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN){
    printf("Initialize lymph node\n");
    model->dendriticLymphNode[0] = dendriticLN;
    model->tHelperLymphNode[0] = thelperLN;
    model->tCytotoxicLymphNode[0] = tcytotoxicLN;
    model->bCellLymphNode[0] = bcellLN;
    model->plasmaCellLymphNode[0] = plasmacellLN;
    model->antibodyLymphNode[0] = antibodyLN;
    printf("Finish initial lymph Node\n");
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

float CalculateChemottaxis(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ, float HX){
    float gradientPopulationI, gradientPopulationJ;
    if(gradientOdcI<0)
        gradientPopulationI = UpDownWind(frontIPoint, ijPoint, avgValue)/(float)HX;
    else
        gradientPopulationI = UpDownWind(ijPoint, rearIPoint, avgValue)/(float)HX;
    if(gradientOdcJ<0)
        gradientPopulationJ = UpDownWind(frontJPoint, ijPoint, avgValue)/(float)HX;
    else
        gradientPopulationJ = UpDownWind(ijPoint, rearJPoint, avgValue)/(float)HX;

    return gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

float CalculateDiffusion(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float HX){
    return (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(HX*HX);
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
    if(model->my_rank == 0){
        int randomVal;
        for(int k = 0; k < model->xSize; k++){
            int i = k / model->xSize;
            int j = k % model->xSize;
            randomVal = rand() % 100;
            if(randomVal <10){
                model->parametersModel.V_BV++;
                model->parametersModel.V_PV++;
                model->thetaBV[i * model->xSize + j] = 1;
                if(j != model->xSize-1)
                    model->thetaPV[i * model->xSize + (j + 1)] = 1;
                else
                    model->thetaPV[i * model->xSize] = 1;
            }
        }
        printf("bv = %d, pv = %d \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
        // WriteBVPV(model->thetaBV, model->thetaPV);
    }
    // Para cada linha da matriz fazer broadcast
    MPI_Bcast(model->thetaPV, model->xSize * model->xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(model->thetaBV, model->xSize * model->xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&model->parametersModel.V_BV, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&model->parametersModel.V_PV, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

structModel ModelInitialize(structParameters params, float ht, float hx, float time, float space, int numFigs, int numPointsLN, int my_rank, int comm_sz
){
    structModel model;
    srand(2);
    model.parametersModel = params;
    
    model.my_rank = my_rank;
    model.comm_sz = comm_sz;
    model.ht = ht;
    model.hx = hx;
    model.tFinal = time;
    model.xFinal = space;
    model.timeLen = (int)(time/ht);
    model.spaceLen = (int)(space/hx);
    model.tSize = (int)(time/ht);
    model.xSize = (int)(space/hx);
    model.intervalFigures = (int)model.tSize/numFigs;
    model.numLines = model.spaceLen/comm_sz;
    model.startLine = my_rank*model.numLines;
    model.endLine = model.startLine + model.numLines-1;
    model.microglia = (float**)malloc(BUFFER * sizeof(float*));
    model.oligodendrocyte = (float**)malloc(BUFFER * sizeof(float*));
    model.tCytotoxic = (float**)malloc(BUFFER * sizeof(float*));
    model.antibody = (float**)malloc(BUFFER * sizeof(float*));
    model.conventionalDc = (float**)malloc(BUFFER * sizeof(float*));
    model.activatedDc = (float**)malloc(BUFFER * sizeof(float*));
    if(my_rank == comm_sz - 1)
        model.endLine = model.xSize - 1;
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
    //definir lymph node
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
    result[5] = antibodyProduction - antibodyDecayment - antibodyMigration;

    return result;
}

void verifyValues(float value, int time, char* populationName){
    if(value < 0 ||  isnanf(value)){
        //printf("Error: %s = (%f) :: time = %f\n", populationName, value, time*HT);
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
    verifyValues(model->dendriticLymphNode[stepKPlus], stepPos, "DC lymph node");
    verifyValues(model->tCytotoxicLymphNode[stepKPlus], stepPos, "CD8 T lymph node");
    verifyValues(model->tHelperLymphNode[stepKPlus], stepPos, "CD4 T lymph node");
    verifyValues(model->bCellLymphNode[stepKPlus], stepPos, "B cell lymph node");
    verifyValues(model->plasmaCellLymphNode[stepKPlus], stepPos, "Plasma cell lymph node");
    verifyValues(model->antibodyLymphNode[stepKPlus], stepPos, "Antibody lymph node");
    
}

void SendBorderThread(structModel *model, int stepKPlus){
    // printf("%d/%d  ::  %d, %d\n",model->my_rank,model->comm_sz, model->startLine, model->endLine);
    if(model->my_rank%2 == 0){
        if(model->my_rank != model->comm_sz-1){
            MPI_Send(&model->oligodendrocyte[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][model->endLine * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(model->my_rank != 0){
            MPI_Send(&model->oligodendrocyte[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
    }else{
        if(model->my_rank != model->comm_sz-1){                
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->endLine+1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(&model->oligodendrocyte[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][(model->endLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
        }
        if(model->my_rank != 0){                
            MPI_Recv(&model->oligodendrocyte[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->tCytotoxic[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->microglia[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->antibody[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->conventionalDc[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&model->activatedDc[stepKPlus][(model->startLine-1) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(&model->oligodendrocyte[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
            MPI_Send(&model->tCytotoxic[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
            MPI_Send(&model->microglia[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
            MPI_Send(&model->antibody[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
            MPI_Send(&model->conventionalDc[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
            MPI_Send(&model->activatedDc[stepKPlus][(model->startLine) * model->xSize], model->xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
        }
    }
}

void RunModel(structModel *model){
    //Save IC
    double initTime, endTime;
    if(model->my_rank == 0) {
        initTime = MPI_Wtime();
    }
    printf("Initialize files");
    WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
    printf("Files writed\n");
    int stepKMinus = 0, stepKPlus, line, column;

    float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0, auxTT = 0.0, auxAT = 0.0, auxDC = 0.0;
    
    for(int kTime = 1; kTime <= model->timeLen; kTime++){
        auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0, auxTT = 0.0, auxAT = 0.0, auxDC = 0.0;
        // solve lymphnode
        SolverLymphNode(model, kTime);
        stepKPlus = kTime%2;
        for(line = model->startLine; line <= model->endLine; line++){
            for(column = 0; column < model->spaceLen; column++){
                
                microgliaKMinus = model->microglia[stepKPlus][line * model->xSize + column];
                conventionalDcKMinus = model->conventionalDc[stepKPlus][line * model->xSize + column];
                activatedDcKMinus = model->activatedDc[stepKPlus][line * model->xSize + column];
                tCytotoxicKMinus = model->tCytotoxic[stepKPlus][line * model->xSize + column];
                antibodyKMinus = model->antibody[stepKPlus][line * model->xSize + column];
                oligodendrocyteKMinus = model->oligodendrocyte[stepKPlus][line * model->xSize + column];

                //Define gradient ODCs
                valIPlus = (line != model->spaceLen-1) ? model->oligodendrocyte[stepKMinus][(line + 1) * model->xSize + column]: model->oligodendrocyte[stepKPlus][line * model->xSize + column];
                valJPlus = (column != model->spaceLen-1) ? model->oligodendrocyte[stepKMinus][(line * model->xSize + (column + 1))]: model->oligodendrocyte[stepKPlus][line * model->xSize + column];
                valIMinus = (line != 0) ? model->oligodendrocyte[stepKMinus][(line - 1) * model->xSize + column]: model->oligodendrocyte[stepKPlus][line * model->xSize + column];
                valJMinus = (column != 0) ? model->oligodendrocyte[stepKMinus][line * model->xSize + (column - 1)]: model->oligodendrocyte[stepKPlus][line * model->xSize + column];
                
                gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*model->hx);
                gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*model->hx);

                //Diffusion and Chemotaxis Mic

                valIPlus  = (line != model->spaceLen-1) ? model->microglia[stepKMinus][(line + 1) * model->xSize + column]: model->microglia[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1) ? model->microglia[stepKMinus][(line * model->xSize + (column + 1))]: model->microglia[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
                valIMinus = (line != 0) ? model->microglia[stepKMinus][(line - 1) * model->xSize + column]: model->microglia[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
                valJMinus = (column != 0) ? model->microglia[stepKMinus][line * model->xSize + (column - 1)]: model->microglia[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);
                
                microgliaDiffusion = model->parametersModel.micDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKPlus][line * model->xSize + column], model->hx);
                microgliaChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKPlus][line * model->xSize + column],\
                model->parametersModel.avgMic, gradientOdcI, gradientOdcJ, model->hx);

                //Diffusion and Chemotaxis CDC

                valIPlus  = (line != model->spaceLen-1) ? model->conventionalDc[stepKMinus][(line + 1) * model->xSize + column]: model->conventionalDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1) ? model->conventionalDc[stepKMinus][(line * model->xSize + (column + 1))]: model->conventionalDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
                valIMinus = (line != 0) ? model->conventionalDc[stepKMinus][(line - 1) * model->xSize + column]: model->conventionalDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
                valJMinus = (column != 0) ? model->conventionalDc[stepKMinus][line * model->xSize + (column - 1)]: model->conventionalDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

                conventionalDcDiffusion = model->parametersModel.cDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKPlus][line * model->xSize + column], model->hx);
                conventionalDcChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKPlus][line * model->xSize + column],\
                model->parametersModel.avgDc, gradientOdcI, gradientOdcJ, model->hx);

                //Difussion and Chemotaxis CD8T

                valIPlus  = (line != model->spaceLen-1) ? model->tCytotoxic[stepKMinus][(line + 1) * model->xSize + column]: model->tCytotoxic[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1) ? model->tCytotoxic[stepKMinus][(line * model->xSize + (column + 1))]: model->tCytotoxic[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
                valIMinus = (line != 0) ? model->tCytotoxic[stepKMinus][(line - 1) * model->xSize + column]: model->tCytotoxic[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
                valJMinus = (column != 0) ? model->tCytotoxic[stepKMinus][line * model->xSize + (column - 1)]: model->tCytotoxic[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

                tCytotoxicDiffusion = model->parametersModel.tCytoDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKPlus][line * model->xSize + column], model->hx);
                tCytotoxicChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKPlus][line * model->xSize + column],\
                model->parametersModel.avgT, gradientOdcI, gradientOdcJ, model->hx);

                //Difussion ADC

                valIPlus  = (line != model->spaceLen-1) ? model->activatedDc[stepKMinus][(line + 1) * model->xSize + column]: model->activatedDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1) ? model->activatedDc[stepKMinus][(line * model->xSize + (column + 1))]: model->activatedDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
                valIMinus = (line != 0) ? model->activatedDc[stepKMinus][(line - 1) * model->xSize + column]: model->activatedDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
                valJMinus = (column != 0) ? model->activatedDc[stepKMinus][line * model->xSize + (column - 1)]: model->activatedDc[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

                activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKPlus][line * model->xSize + column], model->hx);

                //Difussion Antibody

                valIPlus  = (line != model->spaceLen-1) ? model->antibody[stepKMinus][(line + 1) * model->xSize + column]: model->antibody[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1) ? model->antibody[stepKMinus][(line * model->xSize + (column + 1))]: model->antibody[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*rightNeumannBC);
                valIMinus = (line != 0) ? model->antibody[stepKMinus][(line - 1) * model->xSize + column]: model->antibody[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*upperNeumannBC);
                valJMinus = (column != 0) ? model->antibody[stepKMinus][line * model->xSize + (column - 1)]: model->antibody[stepKPlus][line * model->xSize + column] - (float)(2*model->hx*leftNeumannBC);

                antibodyDiffusion = model->parametersModel.antibodyDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKPlus][line * model->xSize + column], model->hx);

                //*******************************************Solving Tissue equations*****************************************************

                //Microglia update
                microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
                microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

                model->microglia[stepKPlus][line * model->xSize + column] = microgliaKMinus + \
                model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);
                if((model->microglia[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->microglia[stepKPlus][line * model->xSize + column])){
                    //printf("Microglia (%f) deu erro no tempo %f\n", model->microglia[stepKPlus][line * model->xSize + column], kTime*HT);
                    exit(0);
                }

                //Conventional DC update
                conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
                conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
                conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

                model->conventionalDc[stepKPlus][line * model->xSize + column] = conventionalDcKMinus + \
                model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);
                if((model->conventionalDc[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->conventionalDc[stepKPlus][line * model->xSize + column])){
                    //printf("CDC (%f) deu erro no tempo %f\n", model->conventionalDc[stepKPlus][line * model->xSize + column], kTime*HT);
                    exit(0);
                }

                //Activated DC update
                activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
                activatedDcMigration = model->thetaPV[line * model->xSize + column]*model->parametersModel.gammaD*(model->dendriticLymphNode[stepKPlus] - activatedDcKMinus);
                
                model->activatedDc[stepKPlus][line * model->xSize + column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
                if((model->activatedDc[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->activatedDc[stepKPlus][line * model->xSize + column])){
                    //printf("ADC (%f) deu erro no tempo %f\n", model->activatedDc[stepKPlus][line * model->xSize + column], kTime*HT);
                    exit(0);
                }

                //CD8 T update
                tCytotoxicMigration = model->thetaBV[line * model->xSize + column]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus);
                
                model->tCytotoxic[stepKPlus][line * model->xSize + column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
                if((model->tCytotoxic[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->tCytotoxic[stepKPlus][line * model->xSize + column])){
                    //printf("tCytotoxic (%f) deu erro no tempo %f\n", model->tCytotoxic[stepKPlus][line * model->xSize + column], kTime*HT);
                    exit(0);
                }

                //Antibody update
                odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
                antibodyMigration = model->thetaBV[line * model->xSize + column]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[stepKPlus] - antibodyKMinus);
                
                model->antibody[stepKPlus][line * model->xSize + column] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);
                if((model->antibody[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->antibody[stepKPlus][line * model->xSize + column])){
                    //printf("antibody (%.8f) deu erro no tempo %f\n", (model->antibody[stepKPlus][line * model->xSize + column]), kTime*HT);
                    exit(0);
                }

                //Oligodendrocytes update
                odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
                odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

                model->oligodendrocyte[stepKPlus][line * model->xSize + column] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
                if((model->oligodendrocyte[stepKPlus][line * model->xSize + column]) < 0 || isnanf (model->oligodendrocyte[stepKPlus][line * model->xSize + column])){
                    //printf("oligodendrocyte (%f) deu erro no tempo %f\n", model->oligodendrocyte[stepKPlus][line * model->xSize + column], kTime*HT);
                    exit(0);
                }
                if(model->thetaBV[line * model->xSize + column] == 1){
                    auxTCytotoxicBV += model->tCytotoxic[stepKPlus][line * model->xSize + column];
                    auxAntibodyBV += model->antibody[stepKPlus][line * model->xSize + column];
                }
                if(model->thetaPV[line * model->xSize + column] == 1){
                    auxAdcPV += model->activatedDc[stepKPlus][line * model->xSize + column];
                }
            }
            
        }
        if(kTime%model->intervalFigures == 0 || kTime == model->timeLen){
            //Cada myRank manda o resultado para o myrank 0
            if(model->my_rank == 0){
                for(int iterRank = 1; iterRank < model->comm_sz; iterRank++){
                    for(int line = model->numLines*iterRank; line < model->numLines*(iterRank + 1) - 1; line++){
                        MPI_Recv(&model->oligodendrocyte[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->tCytotoxic[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->microglia[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->antibody[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->conventionalDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->activatedDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, iterRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                //Se tiver mais alguma linha pra ser enviada pelo ultimo processo
                int lastLineRead = model->numLines*model->comm_sz - 1;
                if(lastLineRead < model->xSize - 1){
                    int recvRank = model->comm_sz - 1;
                    for(int line = lastLineRead; line < model->xSize; line++){
                        MPI_Recv(&model->oligodendrocyte[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->tCytotoxic[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->microglia[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->antibody[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->conventionalDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&model->activatedDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, recvRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }else{
                for(int line = model->startLine; line < model->endLine; line ++){
                    MPI_Send(&model->oligodendrocyte[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&model->tCytotoxic[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&model->microglia[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
                    MPI_Send(&model->antibody[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
                    MPI_Send(&model->conventionalDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
                    MPI_Send(&model->activatedDc[stepKPlus][line * model->xSize], model->xSize, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
                }
            }
            if(model->my_rank == 0)
                WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        }
        MPI_Allreduce(&auxTCytotoxicBV, &auxTT, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&auxAntibodyBV, &auxAT, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&auxAdcPV, &auxDC, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        model->tCytotoxicTissueVessels = (float)(auxTT/model->parametersModel.V_BV);
        model->antibodyTissueVessels = (float)(auxAT/model->parametersModel.V_BV);
        model->activatedDCTissueVessels = (float)(auxDC/model->parametersModel.V_PV);

        SendBorderThread(model, stepKPlus);
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(model->my_rank == 0){
        endTime = MPI_Wtime();
        double totalTime = endTime - initTime;
        printf("Time to launch code in MPI is: %.2lf\n", totalTime);
        printf("Computation Done!!\n");
        WriteLymphNodeFiles(*model, model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
        PlotResults(*model);
    }
}


