#include "model.h"
#include <time.h>
#include <math.h>
// definir estrutura com os parametros do modelo


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

int VerifyCFL(structParameters parametersModel){

    return 0;
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

structModel ModelInitialize(structParameters params, int dt, int dx, int tFinal, int xFinal, int numFiguras){
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
    model.microglia = mesh;
    model.oligodendrocyte = mesh;
    model.convecionalDC = mesh;
    model.activatedDC = mesh;
    model.antibody = mesh;
    model.tCytotoxic = mesh;
    model.thetaBV = mesh[0];
    model.thetaPV = mesh[0];
    InitialConditionTissueMicroglia(&model);

    //definir BV e PV
    DefineBVPV(&model);
    
    //definir lymph node
    float dendriticLN = 0, thelperLN = 0, tcytotoxicLN = 0, bcellLN = 0, plasmacellLN = 0, antibodyLN = 0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    
    return model;

}

float* EquationsLymphNode(float* populationLN, float step){
    //Describe equations
}

void SolverLymphNode(structModel *model, float step){
    float *solutionLN;

    //Execute Euler (or RungeKutta4ThOrder)

}

void RunModel(structModel *model){
    //Save IC


    float upperNeumannBC = 0, lowerNeumannBC = 0, leftNeumannBC = 0, rightNeumannBC = 0;
    
    float valIPlus, valIMinus, valJPlus, valJMinus, gradientOdcI, gradientOdcJ;
    float chemotaxisMic, chemotaxisTCytotoxic, chemotaxisCDC,\
     diffusionMic, diffusionTCytotoxic, diffusionCDC, diffusionADC, diffusionAntibody ;

    for(int kTime = 0; kTime < model->timeLen; kTime++){
        // solve lymphnode
        SolverLymphNode(model, (float)kTime/model->ht);
        
        //Iterando todas as colunas de uma linha antes de ir pra proxima
        for(int line = 0; line < model->spaceLen; line++){
        for(int column = 0; column < model->spaceLen; column++){
            
            //Define gradient ODCs
            valIPlus = (line != model->spaceLen-1)? model->oligodendrocyte[0][line+1][column]: model->oligodendrocyte[0][line][column];
            valJPlus = (column != model->spaceLen-1)? model->oligodendrocyte[0][line][column+1]: model->oligodendrocyte[0][line][column];
            valIPlus = (line != 0)? model->oligodendrocyte[0][line-1][column]: model->oligodendrocyte[0][line][column];
            valJPlus = (column != 0)? model->oligodendrocyte[0][line][column-1]: model->oligodendrocyte[0][line][column];
            
            gradientOdcI = (valIPlus - valIMinus)/(2*model->hx);
            gradientOdcJ = (valJPlus - valJMinus)/(2*model->hx);

            //Diffusion and Chemotaxis Mic

            valIPlus = (line != model->spaceLen-1)? model->microglia[0][line+1][column]: model->microglia[0][line][column];
            valJPlus = (column != model->spaceLen-1)? model->microglia[0][line][column+1]: model->microglia[0][line][column];
            valIPlus = (line != 0)? model->microglia[0][line-1][column]: model->microglia[0][line][column];
            valJPlus = (column != 0)? model->microglia[0][line][column-1]: model->microglia[0][line][column];

            diffusionMic = CalculateDiffusion(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[0][line][column], model->hx);
            chemotaxisMic = CalculateChemottaxis(valJPlus, valJMinus, valIPlus, valIMinus, model->microglia[0][line][column],\
            model->parametersModel.avgMic, gradientOdcI, gradientOdcJ, model->hx);

            
            


            //Solve diffusions


            //Solve Chemotaxis


            //Solve Tissue equations


            //Save Results


        }   
        }
    }
}