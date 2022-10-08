#include "model.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void InitialConditionTissueMicroglia(structModel* model){
    for(int i = 0; i < model->xFinal; i++){
        for(int j = 0; j < model->xFinal; j++){
            if(pow((i-(int)(model->xFinal/model->hx)),2) + pow((j-(int)(model->xFinal/model->hx)),2) < 5)
                model->microglia[0][i][j] = model->parametersModel.avgMic/3;
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

void WritePopulation(float** population, char* fileName, char* bufferTime, int spaceLen){
    printf("Populaçao: %f\n", population[0][5]);
    FILE *file;
    file = fopen(fileName, "w");
    for(int i=0;i<spaceLen;i++) {
        for(int j=0;j<spaceLen;j++) {
            fprintf(file,"%f ",population[0][5]);
            printf("comecando\n");
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

void WriteFiles(structModel model, float** oligodendrocyte, float** microglia, float** tCytotoxic, float** antibody, float** conventionalDC, float ** activatedDC, float time){
    char buffer[6];
    sprintf(buffer, "%f", time);
    // gctv(time, 6,buffer);
    printf("Salvando...\n");
    WritePopulation(oligodendrocyte, "/result/matrix/oligo.txt", buffer, model.spaceLen);
    WritePopulation(microglia, "/result/matrix/microglia.txt", buffer, model.spaceLen);
    WritePopulation(tCytotoxic, "/result/matrix/tCyto.txt", buffer, model.spaceLen);
    WritePopulation(antibody, "/result/matrix/antibody.txt", buffer, model.spaceLen);
    WritePopulation(conventionalDC, "/result/matrix/conventionalDC.txt", buffer, model.spaceLen);
    WritePopulation(activatedDC, "/result/matrix/activatedDC.txt", buffer, model.spaceLen);    

    //Fazer chamada de sistema passando o tempo como argumento para o python e no python salvar a matriz das populações

}

float AdvectionTerm(float populationPoint, float avgValue){
    return populationPoint/(populationPoint + avgValue);
}

float UpDownWind(float frontPoint, float rearPoint, float avgValue){
    return AdvectionTerm(frontPoint, avgValue) - AdvectionTerm(rearPoint, avgValue);
}

float CalculateChemottaxis(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ, float hx){
    float gradientPopulationI = 0, gradientPopulationJ = 0;
    if(gradientOdcI<0)
        gradientPopulationI = UpDownWind(frontIPoint, ijPoint, avgValue)/hx;
    else
        gradientPopulationI = UpDownWind(ijPoint, rearIPoint, avgValue)/hx;
    if(gradientOdcJ<0)
        gradientPopulationJ = UpDownWind(frontJPoint, ijPoint, avgValue)/hx;
    else
        gradientPopulationJ = UpDownWind(ijPoint, rearJPoint, avgValue)/hx;

    return gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

float CalculateDiffusion(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float hx){
    return (frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(hx*hx);
}

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(valuePopulation + avgPopulation);
}

void DefineBVPV(structModel *model){
    int randomVal;
    for(int i = 0; i < model->xFinal; i++){
        for(int j = 0; j < model->xFinal; j++){
            randomVal = rand() % 100;
            if(randomVal <10){
                model->thetaBV[i][j] = 1;
                if(j != model->xFinal-1)
                    model->thetaPV[i][j+1] = 1;
                else
                    model->thetaPV[i][0] = 1;
            }
        }
    }
}

structModel ModelInitialize(structParameters params, float dt, float dx, float tFinal, float xFinal, int numFiguras){
    structModel model;
    srand(time(NULL));
    
    model.parametersModel = params;
    model.numFiguras = numFiguras;
    
    model.ht = dt;
    model.hx = dx;
    model.tFinal = tFinal;
    model.xFinal = xFinal;
    model.timeLen = (int)(tFinal/dt);
    model.spaceLen = (int)(xFinal/dx);
    
    float mesh[2][model.spaceLen][model.spaceLen]; // verificar se começa com zero!
    model.microglia = (float***)malloc(2*sizeof(float**));
    model.microglia[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.microglia[1] = (float**)malloc(model.spaceLen*sizeof(float*));
    
    model.oligodendrocyte = (float***)malloc(2*sizeof(float**));
    model.oligodendrocyte[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.oligodendrocyte[1] = (float**)malloc(model.spaceLen*sizeof(float*));
    
    model.conventionalDc = (float***)malloc(2*sizeof(float**));
    model.conventionalDc[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.conventionalDc[1] = (float**)malloc(model.spaceLen*sizeof(float*));

    model.activatedDc = (float***)malloc(2*sizeof(float**));
    model.activatedDc[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.activatedDc[1] = (float**)malloc(model.spaceLen*sizeof(float*));

    model.tCytotoxic = (float***)malloc(2*sizeof(float**));
    model.tCytotoxic[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.tCytotoxic[1] = (float**)malloc(model.spaceLen*sizeof(float*));

    model.antibody = (float***)malloc(2*sizeof(float**));
    model.antibody[0] = (float**)malloc(model.spaceLen*sizeof(float*));
    model.antibody[1] = (float**)malloc(model.spaceLen*sizeof(float*));
    
    for (int i = 0; i < model.spaceLen; i++){
        model.microglia[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.microglia[1][i] = (float*)malloc(model.spaceLen * sizeof(float));

        model.oligodendrocyte[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.oligodendrocyte[1][i] = (float*)malloc(model.spaceLen * sizeof(float));

        model.conventionalDc[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.conventionalDc[1][i] = (float*)malloc(model.spaceLen * sizeof(float));

        model.activatedDc[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.activatedDc[1][i] = (float*)malloc(model.spaceLen * sizeof(float));

        model.tCytotoxic[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.tCytotoxic[1][i] = (float*)malloc(model.spaceLen * sizeof(float));

        model.antibody[0][i] = (float*)malloc(model.spaceLen * sizeof(float));
        model.antibody[1][i] = (float*)malloc(model.spaceLen * sizeof(float));
    }
    
    model.thetaBV = (float**)malloc(model.spaceLen*sizeof(float*));
    for (int i = 0; i < model.spaceLen; i++)
        model.thetaBV[i] = (float*)malloc(model.spaceLen * sizeof(float));
    
    model.thetaPV = (float**)malloc(model.spaceLen*sizeof(float*));
    for (int i = 0; i < model.spaceLen; i++)
        model.thetaPV[i] = (float*)malloc(model.spaceLen * sizeof(float));
    InitialConditionTissueMicroglia(&model);

    model.dendriticLymphNode  = (float*)malloc(model.spaceLen * sizeof(float));
    model.tHelperLymphNode    = (float*)malloc(model.spaceLen * sizeof(float));
    model.tCytotoxicLymphNode = (float*)malloc(model.spaceLen * sizeof(float));
    model.bCellLymphNode      = (float*)malloc(model.spaceLen * sizeof(float));
    model.plasmaCellLymphNode = (float*)malloc(model.spaceLen * sizeof(float));
    model.antibodyLymphNode   = (float*)malloc(model.spaceLen * sizeof(float));

    //definir BV e PV
    DefineBVPV(&model);
    printf("BV e PV feitos!!\n");
    //definir lymph node
    float dendriticLN = 0, thelperLN = 0, tcytotoxicLN = 0, bcellLN = 0, plasmacellLN = 0, antibodyLN = 0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    
    return model;

}

float* EquationsLymphNode(structModel model, float* populationLN, int stepPos){
    static float result[6];
    
    float dcLN = populationLN[0];
    float tCytoLN = populationLN[1];
    float tHelperLN = populationLN[2];
    float bCellLN = populationLN[3];
    float plasmaCellLN = populationLN[4];
    float antibodyLN = populationLN[5];

    //Describe equations

    //Dendritic cell
    float activatedDcMigration = model.parametersModel.gammaD * (model.activatedDCTissueVessels - dcLN) * (model.parametersModel.V_LV/model.parametersModel.V_LN);
    float activatedDcClearance = model.parametersModel.cDl * dcLN;
    result[0] = activatedDcMigration - activatedDcClearance;

    //T Cytotoxic
    float tCytoActivation = model.parametersModel.bTCytotoxic * (model.parametersModel.rhoTCytotoxic*tCytoLN*dcLN - tCytoLN*tCytoLN*dcLN/model.parametersModel.estableTCytotoxic);
    float tCytoHomeostasis = model.parametersModel.alphaTCytotoxic * (model.parametersModel.estableTCytotoxic - tCytoLN);
    float tCytoMigration = model.parametersModel.gammaT * (tCytoLN - model.tCytotoxicTissueVessels) * model.parametersModel.V_BV/model.parametersModel.V_LN;
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
    float antibodyMigration = model.parametersModel.gammaAntibody * (antibodyLN - model.antibodyTissueVessels) * (model.parametersModel.V_BV/model.parametersModel.V_LN);
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
    
    float *solutionLN = (float*)malloc(6*sizeof(float));
    solutionLN = EquationsLymphNode(*model, populationLN, stepPos);
    
    //Execute Euler (or RungeKutta4ThOrder)
    model->dendriticLymphNode[stepPos] = model->dendriticLymphNode[stepPos-1] + model->ht*solutionLN[0];
    model->tCytotoxicLymphNode[stepPos] = model->tCytotoxicLymphNode[stepPos-1] + model->ht*solutionLN[1];
    model->tHelperLymphNode[stepPos] = model->tHelperLymphNode[stepPos-1] + model->ht*solutionLN[2];
    model->bCellLymphNode[stepPos] = model->bCellLymphNode[stepPos-1] + model->ht*solutionLN[3];
    model->plasmaCellLymphNode[stepPos] = model->plasmaCellLymphNode[stepPos-1] + model->ht*solutionLN[4];
    model->antibodyLymphNode[stepPos] = model->antibodyLymphNode[stepPos-1] + model->ht*solutionLN[5];
    // free(solutionLN);
}

void RunModel(structModel *model){
    //Save IC
    
    int stepKMinus = 0, stepKPlus;

    float upperNeumannBC = 0, lowerNeumannBC = 0, leftNeumannBC = 0, rightNeumannBC = 0;
    
    float valIPlus, valIMinus, valJPlus, valJMinus, gradientOdcI, gradientOdcJ;

    float microgliaChemotaxis, tCytotoxicChemotaxis, conventionalDcChemotaxis,\
     microgliaDiffusion, tCytotoxicDiffusion, conventionalDcDiffusion, activatedDCDiffusion, antibodyDiffusion ;

    float microgliaReaction, microgliaClearance, tCytotoxicMigration, odcAntibodyMicrogliaFagocitosis, \
    odcMicrogliaFagocitosis, odcTCytotoxicApoptosis, conventionalDcReaction, conventionalDcClearance, conventionalDcActivation, \
    activatedDcClearance, activatedDcMigration, antibodyMigration;

    float microgliaKMinus, conventionalDcKMinus, activatedDcKMinus, tCytotoxicKMinus, antibodyKMinus, oligodendrocyteKMinus;

    float auxAdcPV, auxAntibodyBV, auxTCytotoxicBV;

    for(int kTime = 1; kTime <= model->timeLen; kTime++){
        auxAdcPV = 0, auxAntibodyBV = 0, auxTCytotoxicBV = 0;
        // solve lymphnode
        SolverLymphNode(model, kTime);
        
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
            valIPlus = (line != 0)? model->oligodendrocyte[stepKMinus][line-1][column]: model->oligodendrocyte[stepKMinus][line][column];
            valJPlus = (column != 0)? model->oligodendrocyte[stepKMinus][line][column-1]: model->oligodendrocyte[stepKMinus][line][column];
            
            gradientOdcI = (valIPlus - valIMinus)/(2*model->hx);
            gradientOdcJ = (valJPlus - valJMinus)/(2*model->hx);

            //Diffusion and Chemotaxis Mic

            valIPlus = (line != model->spaceLen-1)? model->microglia[stepKMinus][line+1][column]: model->microglia[stepKMinus][line][column];
            valJPlus = (column != model->spaceLen-1)? model->microglia[stepKMinus][line][column+1]: model->microglia[stepKMinus][line][column];
            valIPlus = (line != 0)? model->microglia[stepKMinus][line-1][column]: model->microglia[stepKMinus][line][column];
            valJPlus = (column != 0)? model->microglia[stepKMinus][line][column-1]: model->microglia[stepKMinus][line][column];

            microgliaDiffusion = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line][column], model->hx);
            microgliaChemotaxis = CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[stepKMinus][line][column],\
            model->parametersModel.avgMic, gradientOdcI, gradientOdcJ, model->hx);

            //Diffusion and Chemotaxis CDC

            valIPlus = (line != model->spaceLen-1)? model->conventionalDc[stepKMinus][line+1][column]: model->conventionalDc[stepKMinus][line][column];
            valJPlus = (column != model->spaceLen-1)? model->conventionalDc[stepKMinus][line][column+1]: model->conventionalDc[stepKMinus][line][column];
            valIPlus = (line != 0)? model->conventionalDc[stepKMinus][line-1][column]: model->conventionalDc[stepKMinus][line][column];
            valJPlus = (column != 0)? model->conventionalDc[stepKMinus][line][column-1]: model->conventionalDc[stepKMinus][line][column];

            conventionalDcDiffusion = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line][column], model->hx);
            conventionalDcChemotaxis = CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->conventionalDc[stepKMinus][line][column],\
            model->parametersModel.avgDc, gradientOdcI, gradientOdcJ, model->hx);

            //Difussion and Chemotaxis CD8T

            valIPlus = (line != model->spaceLen-1)? model->tCytotoxic[stepKMinus][line+1][column]: model->tCytotoxic[stepKMinus][line][column];
            valJPlus = (column != model->spaceLen-1)? model->tCytotoxic[stepKMinus][line][column+1]: model->tCytotoxic[stepKMinus][line][column];
            valIPlus = (line != 0)? model->tCytotoxic[stepKMinus][line-1][column]: model->tCytotoxic[stepKMinus][line][column];
            valJPlus = (column != 0)? model->tCytotoxic[stepKMinus][line][column-1]: model->tCytotoxic[stepKMinus][line][column];

            tCytotoxicDiffusion = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line][column], model->hx);
            tCytotoxicChemotaxis = CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->tCytotoxic[stepKMinus][line][column],\
            model->parametersModel.avgT, gradientOdcI, gradientOdcJ, model->hx);

            //Difussion ADC

            valIPlus = (line != model->spaceLen-1)? model->activatedDc[stepKMinus][line+1][column]: model->activatedDc[stepKMinus][line][column];
            valJPlus = (column != model->spaceLen-1)? model->activatedDc[stepKMinus][line][column+1]: model->activatedDc[stepKMinus][line][column];
            valIPlus = (line != 0)? model->activatedDc[stepKMinus][line-1][column]: model->activatedDc[stepKMinus][line][column];
            valJPlus = (column != 0)? model->activatedDc[stepKMinus][line][column-1]: model->activatedDc[stepKMinus][line][column];

            activatedDCDiffusion = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][line][column], model->hx);

            //Difussion Antibody

            valIPlus = (line != model->spaceLen-1)? model->antibody[stepKMinus][line+1][column]: model->antibody[stepKMinus][line][column];
            valJPlus = (column != model->spaceLen-1)? model->antibody[stepKMinus][line][column+1]: model->antibody[stepKMinus][line][column];
            valIPlus = (line != 0)? model->antibody[stepKMinus][line-1][column]: model->antibody[stepKMinus][line][column];
            valJPlus = (column != 0)? model->antibody[stepKMinus][line][column-1]: model->antibody[stepKMinus][line][column];

            antibodyDiffusion = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->antibody[stepKMinus][line][column], model->hx);

            //Solve Tissue equations

            //Microglia update
            microgliaReaction = model->parametersModel.muMic*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
            microgliaClearance = model->parametersModel.cMic*microgliaKMinus;

            model->microglia[stepKPlus][line][column] = microgliaKMinus + \
            model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

            //Conventional DC update
            conventionalDcReaction = model->parametersModel.muCDc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
            conventionalDcActivation = model->parametersModel.bD*conventionalDcKMinus*oligodendrocyteKMinus;
            conventionalDcClearance = model->parametersModel.cCDc*conventionalDcKMinus;

            model->conventionalDc[stepKPlus][line][column] = conventionalDcKMinus + \
            model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);
            
            //Activated DC update
            activatedDcClearance = model->parametersModel.cADc*activatedDcKMinus;
            activatedDcMigration = model->thetaPV[line][column]*model->parametersModel.gammaD*(model->dendriticLymphNode[kTime] - activatedDcKMinus);
            
            model->activatedDc[stepKPlus][line][column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
            
            //CD8 T update
            tCytotoxicMigration = model->thetaBV[line][column]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[kTime] - tCytotoxicKMinus);
            
            model->tCytotoxic[stepKPlus][line][column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
            

            //Antibody update
            odcAntibodyMicrogliaFagocitosis = model->parametersModel.lambAntMic*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
            antibodyMigration = model->thetaBV[line][column]*model->parametersModel.gammaAntibody*(model->antibodyLymphNode[kTime] - antibodyKMinus);
            
            model->antibody[stepKPlus][line][column] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);

            //Oligodendrocytes update
            odcMicrogliaFagocitosis = model->parametersModel.rM*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
            odcTCytotoxicApoptosis = model->parametersModel.rT*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

            model->oligodendrocyte[stepKPlus][line][column] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
            
            if(model->thetaBV[line][column] == 1){
                auxTCytotoxicBV += model->tCytotoxic[stepKPlus][line][column];
                auxAntibodyBV += model->antibody[stepKPlus][line][column];
            }
            if(model->thetaPV[line][column] == 1){
                auxAdcPV += model->activatedDc[stepKPlus][line][column];
            }
            
            WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);


        }   
        }
        model->tCytotoxicTissueVessels = auxTCytotoxicBV;
        model->antibodyTissueVessels = auxAntibodyBV;
        model->activatedDCTissueVessels = auxAdcPV;
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
    }
}