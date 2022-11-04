#include "model.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int VerifyCFL(structParameters parametersModel, float ht, float hx){

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

void WritePopulationLymphNode(float population[tSize], char* fileName){
    FILE *file;
    file = fopen(fileName, "w");
    for(int i=0;i<tSize;i++) {
        fprintf(file, "%f\n", population[i]);
    }
    fclose(file);
}

void WriteLymphNodeFiles(float dendritic[tSize], float tHelper[tSize], float tCytotoxic[tSize], float bCell[tSize], float plasmaCell[tSize], float antibody[tSize]){
    WritePopulationLymphNode(dendritic, "./result/dendritic.txt");
    WritePopulationLymphNode(tHelper, "./result/tHelper.txt");
    WritePopulationLymphNode(tCytotoxic, "./result/tCyto.txt");
    WritePopulationLymphNode(bCell, "./result/bCell.txt");
    WritePopulationLymphNode(plasmaCell, "./result/plasmaCell.txt");
    WritePopulationLymphNode(antibody, "./result/antibody.txt");

    char buffer[10];
    char command[40] = {};
    strcat(command, "python3 plotLymphNode.py ");
    snprintf(buffer, sizeof(buffer), "%d", TIME);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", HT);
    strcat(command, buffer);
    system(command);
}

void WriteFiles(structModel model, float oligodendrocyte[xSize][xSize], float microglia[xSize][xSize], float tCytotoxic[xSize][xSize], float antibody[xSize][xSize], float conventionalDC[xSize][xSize], float  activatedDC[xSize][xSize], float time){
    char buffer[10];
    float day = time*HT;
    snprintf(buffer, sizeof(buffer), "%f", day);
    WritePopulation(oligodendrocyte, "./result/matrix/oligo.txt", buffer);
    WritePopulation(microglia, "./result/matrix/microglia.txt", buffer);
    WritePopulation(tCytotoxic, "./result/matrix/tCyto.txt", buffer);
    WritePopulation(antibody, "./result/matrix/antibody.txt", buffer);
    WritePopulation(conventionalDC, "./result/matrix/conventionalDC.txt", buffer);
    WritePopulation(activatedDC, "./result/matrix/activatedDC.txt", buffer);
    //Fazer chamada de sistema passando o tempo como argumento para o python e no python salvar a matriz das populações
    char command[70] = {};
    strcat(command, "python3 plotMatrices.py ");
    snprintf(buffer, sizeof(buffer), "%d", LENGTH);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", HX);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", day);
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

void DefineBVPV(structModel *model){
    int randomVal;
    for(int i = 0; i < xSize; i++){
        for(int j = 0; j < xSize; j++){
            randomVal = rand() % 100;
            if(randomVal <10){
                model->parametersModel.V_BV++;
                model->parametersModel.V_PV++;
                model->thetaBV[i][j] = 1;
                if(j != xSize-1)
                    model->thetaPV[i][j+1] = 1;
                else
                    model->thetaPV[i][0] = 1;
            }
        }
    }
    printf("bv = %d, pv = %d \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    WriteBVPV(model->thetaBV, model->thetaPV);
}

void WriteBVPV(float thetaBV[xSize][xSize], float thetaPV[xSize][xSize]){
    FILE *fileBV;
    fileBV = fopen("./result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./result/pv.txt", "w");
    for(int i = 0; i < xSize; i++){        
    for(int j = 0; j < xSize; j++){
        fprintf(fileBV, "%f ", thetaBV[i][j]);
        fprintf(filePV, "%f ", thetaPV[i][j]);    
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

structModel ModelInitialize(structParameters params){
    structModel model;
    srand(2);
    
    model.parametersModel = params;
    model.intervaloFiguras = tSize/NUMFIGURAS;
    
    model.ht = HT;
    model.hx = HX;
    model.tFinal = TIME;
    model.xFinal = LENGTH;
    model.timeLen = (int)(TIME/HT);
    model.spaceLen = (int)(LENGTH/HX);
    
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
    float activatedDcClearance = 0.0;//model.parametersModel.cDl * dcLN;
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

    populationLN[0] = model->dendriticLymphNode[stepPos-1];
    populationLN[1] = model->tCytotoxicLymphNode[stepPos-1];
    populationLN[2] = model->tHelperLymphNode[stepPos-1];
    populationLN[3] = model->bCellLymphNode[stepPos-1];
    populationLN[4] = model->plasmaCellLymphNode[stepPos-1];
    populationLN[5] = model->antibodyLymphNode[stepPos-1];
    
    float* solutionLN;
    solutionLN = EquationsLymphNode(*model, populationLN, stepPos);
    
    //Execute Euler (or RungeKutta4ThOrder)
    model->dendriticLymphNode[stepPos] = model->dendriticLymphNode[stepPos-1] + model->ht*solutionLN[0];
    model->tCytotoxicLymphNode[stepPos] = model->tCytotoxicLymphNode[stepPos-1] + model->ht*solutionLN[1];
    model->tHelperLymphNode[stepPos] = model->tHelperLymphNode[stepPos-1] + model->ht*solutionLN[2];
    model->bCellLymphNode[stepPos] = model->bCellLymphNode[stepPos-1] + model->ht*solutionLN[3];
    model->plasmaCellLymphNode[stepPos] = model->plasmaCellLymphNode[stepPos-1] + model->ht*solutionLN[4];
    model->antibodyLymphNode[stepPos] = model->antibodyLymphNode[stepPos-1] + model->ht*solutionLN[5];
    free(solutionLN); 

    if(model->dendriticLymphNode[stepPos] < 0 ||  isnanf (model->dendriticLymphNode[stepPos])){
        printf("DC lymph node (%f) deu erro no tempo %f\n", model->dendriticLymphNode[stepPos], stepPos*HT);
        exit(0);
    }

    if(model->tCytotoxicLymphNode[stepPos] < 0 ||  isnanf (model->tCytotoxicLymphNode[stepPos])){
        printf("CD8 T lymph node (%f) deu erro no tempo %f\n", model->tCytotoxicLymphNode[stepPos], stepPos*HT);
        exit(0);
    }

    if(model->tHelperLymphNode[stepPos] < 0 ||  isnanf (model->tHelperLymphNode[stepPos])){
        printf("CD4 T lymph node (%f) deu erro no tempo %f\n", model->tHelperLymphNode[stepPos], stepPos*HT);
        exit(0);
    }

    if(model->bCellLymphNode[stepPos] < 0 || isnanf (model->bCellLymphNode[stepPos])){
        printf("B cell lymph node (%f) deu erro no tempo %f\n", model->bCellLymphNode[stepPos], stepPos*HT);
        exit(0);
    }

    if(model->plasmaCellLymphNode[stepPos] < 0 || isnanf (model->plasmaCellLymphNode[stepPos])){
        printf("Plasma cell lymph node (%f) deu erro no tempo %f\n", model->plasmaCellLymphNode[stepPos], stepPos*HT);
        exit(0);
    }

    if(model->antibodyLymphNode[stepPos] < 0 || isnanf (model->antibodyLymphNode[stepPos])){
        printf("Antibody lymph node (%f) deu erro no tempo %f\n", model->antibodyLymphNode[stepPos], stepPos*HT);
        exit(0);
    }
}

void RunModel(structModel *model){
    //Save IC
    
    WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);

    int stepKMinus = 0, stepKPlus;

    float upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    float valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    float microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    float microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    float microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    float auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;

    for(int kTime = 1; kTime <= model->timeLen; kTime++){
        auxAdcPV = 0.0, auxAntibodyBV = 0.0, auxTCytotoxicBV = 0.0;
        // solve lymphnode
        // SolverLymphNode(model, kTime);
        stepKPlus = kTime%2;
        for(int line = 0; line < model->spaceLen; line++){//Iterando todas as colunas de uma linha antes de ir pra proxima linha
        for(int column = 0; column < model->spaceLen; column++){            
            
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
            activatedDcMigration = model->thetaPV[line][column]*model->parametersModel.gammaD*(model->dendriticLymphNode[kTime] - activatedDcKMinus);
            
            model->activatedDc[stepKPlus][line][column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
            if((model->activatedDc[stepKPlus][line][column]) < 0 || isnanf (model->activatedDc[stepKPlus][line][column])){
                printf("ADC (%f) deu erro no tempo %f\n", model->activatedDc[stepKPlus][line][column], kTime*HT);
                exit(0);
            }

            //CD8 T update
            tCytotoxicMigration = model->thetaBV[line][column]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[kTime] - tCytotoxicKMinus);
            
            model->tCytotoxic[stepKPlus][line][column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
            if((model->tCytotoxic[stepKPlus][line][column]) < 0 || isnanf (model->tCytotoxic[stepKPlus][line][column])){
                printf("tCytotoxic (%f) deu erro no tempo %f\n", model->tCytotoxic[stepKPlus][line][column], kTime*HT);
                exit(0);
            }

            //Antibody update
            odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
            antibodyMigration = model->thetaBV[line][column]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[kTime] - antibodyKMinus);
            
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
            if(model->thetaBV[line][column] == 1){
                auxTCytotoxicBV += model->tCytotoxic[stepKPlus][line][column];
                auxAntibodyBV += model->antibody[stepKPlus][line][column];
            }
            if(model->thetaPV[line][column] == 1){
                auxAdcPV += model->activatedDc[stepKPlus][line][column];
            }
        }
        }
        if(kTime%model->intervaloFiguras == 0 || kTime == model->timeLen)
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
        model->tCytotoxicTissueVessels = auxTCytotoxicBV;
        model->antibodyTissueVessels = auxAntibodyBV;
        model->activatedDCTissueVessels = auxAdcPV;
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
    }

    WriteLymphNodeFiles(model->dendriticLymphNode, model->tHelperLymphNode, model->tCytotoxicLymphNode, model->bCellLymphNode, model->plasmaCellLymphNode, model->antibodyLymphNode);
}