#include "modelMPI.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "log.h"

#define DEBUG 1

const int xSize = (int)(LENGTH/HX);
const int tSize = (int)(TIME/HT);

void InitialConditionTissueMicroglia(structModel* model){
    for(int i = 0; i < xSize; i++){
        for(int j = 0; j < xSize; j++){
            if(pow((i-(int)(xSize/2)),2) + pow((j-(int)(xSize/2)),2) < 5){
                model->microglia[0][i][j] = (float)model->parametersModel.avgMic/3;
            }
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

int VerifyCFL(structParameters parametersModel){
    if(parametersModel.micDiffusion*HT/(HX*HX) < 0.25 && parametersModel.cDcDiffusion*HT/(HX*HX) < 0.25 && parametersModel.aDcDiffusion*HT/(HX*HX) < 0.25 && parametersModel.tCytoDiffusion*HT/(HX*HX) < 0.25 && parametersModel.chi*HT/HX < 1)
        return 1;
    return 0;
}

void WritePopulation(float population[xSize][xSize], char* fileName, char* bufferTime){
    FILE *file;
    file = fopen(fileName, "w");
    for(int i=0;i<xSize;i++) {
        for(int j=0;j<xSize;j++){
            fprintf(file, "%f ", population[i][j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

void WritePopulationLymphNode(float population[NUMPOINTSLYMPHNODE], char* fileName){
    char absolutePath[70];
    getcwd(absolutePath, sizeof(absolutePath));
    strcat(absolutePath, fileName);
    FILE *file;
    file = fopen(absolutePath, "w");
    for(int i=0;i<NUMPOINTSLYMPHNODE;i++) {
        fprintf(file, "%f\n", population[i]);
    }
    fclose(file);
}


void WriteLymphNodeFiles(float dendritic[NUMPOINTSLYMPHNODE], float tHelper[NUMPOINTSLYMPHNODE], float tCytotoxic[NUMPOINTSLYMPHNODE], float bCell[NUMPOINTSLYMPHNODE], float plasmaCell[NUMPOINTSLYMPHNODE], float antibody[NUMPOINTSLYMPHNODE]){
    WritePopulationLymphNode(dendritic, "/modelmpi/result/dendritic.txt");
    WritePopulationLymphNode(tHelper, "/modelmpi/result/tHelper.txt");
    WritePopulationLymphNode(tCytotoxic, "/modelmpi/result/tCyto.txt");
    WritePopulationLymphNode(bCell, "/modelmpi/result/bCell.txt");
    WritePopulationLymphNode(plasmaCell, "/modelmpi/result/plasmaCell.txt");
    WritePopulationLymphNode(antibody, "/modelmpi/result/antibody.txt");

    char buffer[10];
    char command[40] = {};
    strcat(command, "python3 plotLymphNode.py ");
    snprintf(buffer, sizeof(buffer), "%d", TIME);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", (tSize/NUMPOINTSLYMPHNODE)*HT);
    strcat(command, buffer);
    system(command);
}

void WriteFiles(structModel model, float oligodendrocyte[xSize][xSize], float microglia[xSize][xSize], float tCytotoxic[xSize][xSize], float antibody[xSize][xSize], float conventionalDC[xSize][xSize], float  activatedDC[xSize][xSize], float time){
    char buffer[10];
    float day = time*HT;
    snprintf(buffer, sizeof(buffer), "%.1f", day);
    char pathOligodendrocytes[70];
    getcwd(pathOligodendrocytes, sizeof(pathOligodendrocytes));
    strcat(pathOligodendrocytes, "/modelmpi/result/matrix/oligo");
    strcat(pathOligodendrocytes, buffer);
    strcat(pathOligodendrocytes, ".txt");
    printf("%s\n", pathOligodendrocytes);
    WritePopulation(oligodendrocyte, pathOligodendrocytes, buffer);
    char pathMicroglia[70];
    getcwd(pathMicroglia, sizeof(pathMicroglia));
    strcat(pathMicroglia, "/modelmpi/result/matrix/microglia");
    strcat(pathMicroglia, buffer);
    strcat(pathMicroglia, ".txt");
    WritePopulation(microglia, pathMicroglia, buffer);
    char pathTCyto[70];
    getcwd(pathTCyto, sizeof(pathTCyto));
    strcat(pathTCyto, "/modelmpi/result/matrix/tCyto");
    strcat(pathTCyto, buffer);
    strcat(pathTCyto, ".txt");
    WritePopulation(tCytotoxic, pathTCyto, buffer);
    char pathAntibody[70];
    getcwd(pathAntibody, sizeof(pathAntibody));
    strcat(pathAntibody, "/modelmpi/result/matrix/antibody");
    strcat(pathAntibody, buffer);
    strcat(pathAntibody, ".txt");
    WritePopulation(antibody, pathAntibody, buffer);
    char pathConventionalDC[70];
    getcwd(pathConventionalDC, sizeof(pathConventionalDC));
    strcat(pathConventionalDC, "/modelmpi/result/matrix/conventionalDC");
    strcat(pathConventionalDC, buffer);
    strcat(pathConventionalDC, ".txt");
    WritePopulation(conventionalDC, pathConventionalDC, buffer);
    char pathActivatedDC[70];
    getcwd(pathActivatedDC, sizeof(pathActivatedDC));
    strcat(pathActivatedDC, "/modelmpi/result/matrix/activatedDC");
    strcat(pathActivatedDC, buffer);
    strcat(pathActivatedDC, ".txt");
    WritePopulation(activatedDC, pathActivatedDC, buffer);
}      

void PlotResults(){
    printf("Saving results...\n\n");
    char buffer[10];
    char command[70] = {};
    strcat(command, "python3 plotMatrices.py ");
    snprintf(buffer, sizeof(buffer), "%d", LENGTH);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", HX);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", TIME);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", NUMFIGS);
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
 float avgValue, float gradientOdcI, float gradientOdcJ){
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

float CalculateDiffusion(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint){
    return (float)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(float)(HX*HX);
}

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
}

void WriteBVPV(float thetaBV[xSize*xSize], float thetaPV[xSize*xSize]){
    FILE *fileBV;
    fileBV = fopen("./modelmpi/result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./modelmpi/result/pv.txt", "w");
    for(int i = 0; i < xSize; i++){        
    for(int j = 0; j < xSize; j++){
        fprintf(fileBV, "%f ", thetaBV[i * xSize + j]);
        fprintf(filePV, "%f ", thetaPV[i * xSize + j]);    
    }
    fprintf(fileBV,"\n");
    fprintf(filePV,"\n");
    }
    fclose(fileBV);
    fclose(filePV);
    char buffer[10];
    char command[70] = {};
    strcat(command, "python3 plotBVPV.py ");
    snprintf(buffer, sizeof(buffer), "%d", LENGTH);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", HX);
    strcat(command, buffer);
    // system(command);
}

void DefineBVPV(structModel *model){
    if(model->my_rank == 0){
        int randomVal;
        for(int k = 0; k < xSize; k++){
            int i = k / xSize;
            int j = k % xSize;
            randomVal = rand() % 100;
            if(randomVal <10){
                model->parametersModel.V_BV++;
                model->parametersModel.V_PV++;
                model->thetaBV[i * xSize + j] = 1;
                if(j != xSize-1)
                    model->thetaPV[i * xSize + (j + 1)] = 1;
                else
                    model->thetaPV[i * xSize] = 1;
            }
        }
        printf("bv = %d, pv = %d \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
        // WriteBVPV(model->thetaBV, model->thetaPV);
    }
    // Para cada linha da matriz fazer broadcast
    MPI_Bcast(model->thetaPV, xSize*xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(model->thetaBV, xSize*xSize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&model->parametersModel.V_BV, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&model->parametersModel.V_PV, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

structModel ModelInitialize(structParameters params, int my_rank, int comm_sz){
    structModel model;
    srand(2);
    model.parametersModel = params;
    model.intervaloFiguras = (int)tSize/NUMFIGS;
    
    model.my_rank = my_rank;
    model.comm_sz = comm_sz;
    
    model.ht = HT;
    model.hx = HX;
    model.tFinal = TIME;
    model.xFinal = LENGTH;
    model.timeLen = (int)(TIME/HT);
    model.spaceLen = (int)(LENGTH/HX);

    model.numLines = model.spaceLen/comm_sz;
    model.startLine = my_rank*model.numLines;
    model.endLine = model.startLine + model.numLines-1;
    if(my_rank == comm_sz - 1)
        model.endLine = xSize - 1;
    //definir BV e PV
    DefineBVPV(&model);
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
        printf("Error: %s = (%f) :: time = %f\n", populationName, value, time*HT);
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

    int intervalPoints = (int)(tSize/NUMPOINTSLYMPHNODE);
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
            MPI_Send(model->oligodendrocyte[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(model->tCytotoxic[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(model->microglia[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(model->antibody[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(model->conventionalDc[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(model->activatedDc[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(model->oligodendrocyte[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->tCytotoxic[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->microglia[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->antibody[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->conventionalDc[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->activatedDc[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(model->my_rank != 0){
            MPI_Send(model->oligodendrocyte[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
            MPI_Send(model->tCytotoxic[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
            MPI_Send(model->microglia[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
            MPI_Send(model->antibody[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
            MPI_Send(model->conventionalDc[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
            MPI_Send(model->activatedDc[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
            
            MPI_Recv(model->oligodendrocyte[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->tCytotoxic[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->microglia[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->antibody[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->conventionalDc[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->activatedDc[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
    }else{
        if(model->my_rank != model->comm_sz-1){                
            MPI_Recv(model->oligodendrocyte[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->tCytotoxic[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->microglia[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->antibody[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->conventionalDc[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->activatedDc[stepKPlus][model->endLine+1], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(model->oligodendrocyte[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(model->tCytotoxic[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 1, MPI_COMM_WORLD);
            MPI_Send(model->microglia[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 2, MPI_COMM_WORLD);
            MPI_Send(model->antibody[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 3, MPI_COMM_WORLD);
            MPI_Send(model->conventionalDc[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 4, MPI_COMM_WORLD);
            MPI_Send(model->activatedDc[stepKPlus][model->endLine], xSize, MPI_FLOAT, model->my_rank+1, 5, MPI_COMM_WORLD);
        }
        if(model->my_rank != 0){                
            MPI_Recv(model->oligodendrocyte[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->tCytotoxic[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->microglia[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->antibody[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->conventionalDc[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(model->activatedDc[stepKPlus][model->startLine-1], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(model->oligodendrocyte[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 0, MPI_COMM_WORLD);
            MPI_Send(model->tCytotoxic[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 1, MPI_COMM_WORLD);
            MPI_Send(model->microglia[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 2, MPI_COMM_WORLD);
            MPI_Send(model->antibody[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 3, MPI_COMM_WORLD);
            MPI_Send(model->conventionalDc[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 4, MPI_COMM_WORLD);
            MPI_Send(model->activatedDc[stepKPlus][model->startLine], xSize, MPI_FLOAT, model->my_rank-1, 5, MPI_COMM_WORLD);
        }
    }
}

void RunModel(structModel *model){
    //Save IC
    double initTime, endTime;
    if(model->my_rank == 0) {
        initTime = MPI_Wtime();
    }
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
                
                microgliaKMinus = model->microglia[stepKMinus][line][column];
                conventionalDcKMinus = model->conventionalDc[stepKMinus][line][column];
                activatedDcKMinus = model->activatedDc[stepKMinus][line][column];
                tCytotoxicKMinus = model->tCytotoxic[stepKMinus][line][column];
                antibodyKMinus = model->antibody[stepKMinus][line][column];
                oligodendrocyteKMinus = model->oligodendrocyte[stepKMinus][line][column];

                //Define gradient ODCs
                valIPlus = (line != model->spaceLen-1)? model->oligodendrocyte[stepKMinus][line+1][column]: model->oligodendrocyte[stepKMinus][line][column];
                valJPlus = (column != model->spaceLen-1)? model->oligodendrocyte[stepKMinus][line][column+1]: model->oligodendrocyte[stepKMinus][line][column];
                valIMinus = (line != 0)? model->oligodendrocyte[stepKMinus][line-1][column]: model->oligodendrocyte[stepKMinus][line][column];
                valJMinus = (column != 0)? model->oligodendrocyte[stepKMinus][line][column-1]: model->oligodendrocyte[stepKMinus][line][column];
                
                gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2*HX);
                gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2*HX);

                //Diffusion and Chemotaxis Mic

                valIPlus  = (line != model->spaceLen-1)? model->microglia[stepKMinus][line+1][column]: model->microglia[stepKMinus][line][column] - (float)(2*HX*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1)? model->microglia[stepKMinus][line][column+1]: model->microglia[stepKMinus][line][column] - (float)(2*HX*rightNeumannBC);
                valIMinus = (line != 0)? model->microglia[stepKMinus][line-1][column]: model->microglia[stepKMinus][line][column] - (float)(2*HX*upperNeumannBC);
                valJMinus = (column != 0)? model->microglia[stepKMinus][line][column-1]: model->microglia[stepKMinus][line][column] - (float)(2*HX*leftNeumannBC);
                
                microgliaDiffusion = model->parametersModel.micDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line][column]);
                microgliaChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line][column],\
                model->parametersModel.avgMic, gradientOdcI, gradientOdcJ);

                //Diffusion and Chemotaxis CDC

                valIPlus  = (line != model->spaceLen-1)? model->conventionalDc[stepKMinus][line+1][column]: model->conventionalDc[stepKMinus][line][column] - (float)(2*HX*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1)? model->conventionalDc[stepKMinus][line][column+1]: model->conventionalDc[stepKMinus][line][column] - (float)(2*HX*rightNeumannBC);
                valIMinus = (line != 0)? model->conventionalDc[stepKMinus][line-1][column]: model->conventionalDc[stepKMinus][line][column] - (float)(2*HX*upperNeumannBC);
                valJMinus = (column != 0)? model->conventionalDc[stepKMinus][line][column-1]: model->conventionalDc[stepKMinus][line][column] - (float)(2*HX*leftNeumannBC);

                conventionalDcDiffusion = model->parametersModel.cDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line][column]);
                conventionalDcChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line][column],\
                model->parametersModel.avgDc, gradientOdcI, gradientOdcJ);

                //Difussion and Chemotaxis CD8T

                valIPlus  = (line != model->spaceLen-1)? model->tCytotoxic[stepKMinus][line+1][column]: model->tCytotoxic[stepKMinus][line][column] - (float)(2*HX*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1)? model->tCytotoxic[stepKMinus][line][column+1]: model->tCytotoxic[stepKMinus][line][column] - (float)(2*HX*rightNeumannBC);
                valIMinus = (line != 0)? model->tCytotoxic[stepKMinus][line-1][column]: model->tCytotoxic[stepKMinus][line][column] - (float)(2*HX*upperNeumannBC);
                valJMinus = (column != 0)? model->tCytotoxic[stepKMinus][line][column-1]: model->tCytotoxic[stepKMinus][line][column] - (float)(2*HX*leftNeumannBC);

                tCytotoxicDiffusion = model->parametersModel.tCytoDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line][column]);
                tCytotoxicChemotaxis = model->parametersModel.chi*CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line][column],\
                model->parametersModel.avgT, gradientOdcI, gradientOdcJ);

                //Difussion ADC

                valIPlus  = (line != model->spaceLen-1)? model->activatedDc[stepKMinus][line+1][column]: model->activatedDc[stepKMinus][line][column] - (float)(2*HX*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1)? model->activatedDc[stepKMinus][line][column+1]: model->activatedDc[stepKMinus][line][column] - (float)(2*HX*rightNeumannBC);
                valIMinus = (line != 0)? model->activatedDc[stepKMinus][line-1][column]: model->activatedDc[stepKMinus][line][column] - (float)(2*HX*upperNeumannBC);
                valJMinus = (column != 0)? model->activatedDc[stepKMinus][line][column-1]: model->activatedDc[stepKMinus][line][column] - (float)(2*HX*leftNeumannBC);

                activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][line][column]);

                //Difussion Antibody

                valIPlus  = (line != model->spaceLen-1)? model->antibody[stepKMinus][line+1][column]: model->antibody[stepKMinus][line][column] - (float)(2*HX*lowerNeumannBC);
                valJPlus  = (column != model->spaceLen-1)? model->antibody[stepKMinus][line][column+1]: model->antibody[stepKMinus][line][column] - (float)(2*HX*rightNeumannBC);
                valIMinus = (line != 0)? model->antibody[stepKMinus][line-1][column]: model->antibody[stepKMinus][line][column] - (float)(2*HX*upperNeumannBC);
                valJMinus = (column != 0)? model->antibody[stepKMinus][line][column-1]: model->antibody[stepKMinus][line][column] - (float)(2*HX*leftNeumannBC);

                antibodyDiffusion = model->parametersModel.antibodyDiffusion*CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKMinus][line][column]);

                //*******************************************Solving Tissue equations*****************************************************

                //Microglia update
                microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
                microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

                model->microglia[stepKPlus][line][column] = microgliaKMinus + \
                model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);
                if((model->microglia[stepKPlus][line][column]) < 0 || isnanf (model->microglia[stepKPlus][line][column])){
                    printf("Microglia (%f) deu erro no tempo %f\n", model->microglia[stepKPlus][line][column], kTime*HT);
                    exit(0);
                }

                //Conventional DC update
                conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
                conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
                conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

                model->conventionalDc[stepKPlus][line][column] = conventionalDcKMinus + \
                model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);
                if((model->conventionalDc[stepKPlus][line][column]) < 0 || isnanf (model->conventionalDc[stepKPlus][line][column])){
                    printf("CDC (%f) deu erro no tempo %f\n", model->conventionalDc[stepKPlus][line][column], kTime*HT);
                    exit(0);
                }

                //Activated DC update
                activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
                activatedDcMigration = model->thetaPV[line * xSize + column]*model->parametersModel.gammaD*(model->dendriticLymphNode[stepKPlus] - activatedDcKMinus);
                
                model->activatedDc[stepKPlus][line][column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
                if((model->activatedDc[stepKPlus][line][column]) < 0 || isnanf (model->activatedDc[stepKPlus][line][column])){
                    printf("ADC (%f) deu erro no tempo %f\n", model->activatedDc[stepKPlus][line][column], kTime*HT);
                    exit(0);
                }

                //CD8 T update
                tCytotoxicMigration = model->thetaBV[line * xSize + column]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus);
                
                model->tCytotoxic[stepKPlus][line][column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
                if((model->tCytotoxic[stepKPlus][line][column]) < 0 || isnanf (model->tCytotoxic[stepKPlus][line][column])){
                    printf("tCytotoxic (%f) deu erro no tempo %f\n", model->tCytotoxic[stepKPlus][line][column], kTime*HT);
                    exit(0);
                }

                //Antibody update
                odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
                antibodyMigration = model->thetaBV[line * xSize + column]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[stepKPlus] - antibodyKMinus);
                
                model->antibody[stepKPlus][line][column] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);
                if((model->antibody[stepKPlus][line][column]) < 0 || isnanf (model->antibody[stepKPlus][line][column])){
                    printf("antibody (%.8f) deu erro no tempo %f\n", (model->antibody[stepKPlus][line][column]), kTime*HT);
                    exit(0);
                }

                //Oligodendrocytes update
                odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
                odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

                model->oligodendrocyte[stepKPlus][line][column] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
                if((model->oligodendrocyte[stepKPlus][line][column]) < 0 || isnanf (model->oligodendrocyte[stepKPlus][line][column])){
                    printf("oligodendrocyte (%f) deu erro no tempo %f\n", model->oligodendrocyte[stepKPlus][line][column], kTime*HT);
                    exit(0);
                }
                if(model->thetaBV[line * xSize + column] == 1){
                    auxTCytotoxicBV += model->tCytotoxic[stepKPlus][line][column];
                    auxAntibodyBV += model->antibody[stepKPlus][line][column];
                }
                if(model->thetaPV[line * xSize + column] == 1){
                    auxAdcPV += model->activatedDc[stepKPlus][line][column];
                }
            }
            
        }
        if(kTime%model->intervaloFiguras == 0 || kTime == model->timeLen){
            //Cada myRank manda o resultado para o myrank 0
            if(model->my_rank == 0){
                for(int iterRank = 1; iterRank < model->comm_sz; iterRank++){
                    for(int line = model->numLines*iterRank; line < model->numLines*(iterRank + 1) - 1; line++){
                        MPI_Recv(model->oligodendrocyte[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->tCytotoxic[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->microglia[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->antibody[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->conventionalDc[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->activatedDc[stepKPlus][line], xSize, MPI_FLOAT, iterRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                //Se tiver mais alguma linha pra ser enviada pelo ultimo processo
                int lastLineRead = model->numLines*model->comm_sz - 1;
                if(lastLineRead < xSize - 1){
                    int recvRank = model->comm_sz - 1;
                    for(int line = lastLineRead; line < xSize; line++){
                        MPI_Recv(model->oligodendrocyte[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->tCytotoxic[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->microglia[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->antibody[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->conventionalDc[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(model->activatedDc[stepKPlus][line], xSize, MPI_FLOAT, recvRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }else{
                for(int line = model->startLine; line < model->endLine; line ++){
                    MPI_Send(model->oligodendrocyte[stepKPlus][line], xSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(model->tCytotoxic[stepKPlus][line], xSize, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(model->microglia[stepKPlus][line], xSize, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
                    MPI_Send(model->antibody[stepKPlus][line], xSize, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
                    MPI_Send(model->conventionalDc[stepKPlus][line], xSize, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
                    MPI_Send(model->activatedDc[stepKPlus][line], xSize, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
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
        WriteLymphNodeFiles(model->dendriticLymphNodeSavedPoints, model->tHelperLymphNodeSavedPoints, model->tCytotoxicLymphNodeSavedPoints, model->bCellLymphNodeSavedPoints, model->plasmaCellLymphNodeSavedPoints, model->antibodyLymphNodeSavedPoints);
        PlotResults();
    }
}


