#include "model.h"
#include <time.h>
#include <math.h>
// definir estrutura com os parametros do modelo


void InitialConditionTissueMicroglia(structModel* model){
    for(int i = 0; i < model->xFinal; i++){
        for(int j = 0; j < model->xFinal; j++){
            if(pow((i-(int)(model->xFinal/model->hx)),2) + pow((j-(int)(model->xFinal/model->hx)),2) < 5)
                model->microglia[0][i][j] = model->parametersModel.mic_media/3;
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

float CalculaQuimiotaxia(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual,\
 float valor_medio, float gradiente_odc_i, float gradiente_odc_j, float hx){
    return 0;
}

float CalculaDifusao(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual, float hx){
    return 0;
}

void DefineBVPV(structModel *model, int randomvessels){
    int randomVal;
    for(int i = 0; i < model->xFinal; i++){
        for(int j = 0; j < model->xFinal; j++){
            randomVal = rand() % 100;
            if(randomVal <10){
                model->matrixBV[i][j] = 1;
                model->matrixBV[i][j+1] = (j != model->xFinal-1) ?  1: 0;
            }
        }
    }
}

structModel ModelInitialize(structParameters params, int dt, int dx, int tFinal, int xFinal, int randomvessels){
    structModel model;
    if (randomvessels)
        srand(time(NULL));
    model.parametersModel = params;

    model.ht = dt;
    model.hx = dx;
    model.tFinal = tFinal;
    model.xFinal = xFinal;
    model.timeLen = (int)(tFinal/dt);
    model.spaceLen = (int)(xFinal/dx);
    
    float mesh[2][model.spaceLen][model.spaceLen]; // verificar se começa com zero!
    model.oligodendrocyte = mesh;
    model.convecionalDC = mesh;
    model.activatedDC = mesh;
    model.antibody = mesh;
    model.tCytotoxic = mesh;
    model.matrixBV = mesh[0];
    model.matrixPV = mesh[0];
    InitialConditionTissueMicroglia(&model);

    //definir BV e PV
    DefineBVPV(&model, randomvessels);
    
    //definir lymph node
    float dendriticLN = 0, thelperLN = 0, tcytotoxicLN = 0, bcellLN = 0, plasmacellLN = 0, antibodyLN = 0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    
    return model;

}



// estrutura com os principais campos(ht, hx, matrizes tecido, vetores)