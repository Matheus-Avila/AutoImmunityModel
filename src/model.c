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
    model.conventionalDc = mesh;
    model.activatedDc = mesh;
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

float* EquationsLymphNode(structModel model, float* populationLN, float step){
    float result[6];
    
    float dcLN = populationLN[0];
    float tCytoLN = populationLN[1];
    float tHelperLN = populationLN[2];
    float bCellLN = populationLN[3];
    float plasmaCellLN = populationLN[4];
    float iggGLN = populationLN[5];

    //Describe equations

    //Dendritic cell
    float activatedDcMigration = model.parametersModel.gamma_D * (model.activatedDCTissueVessels - dcLN) * (model.parametersModel.V_LV/model.parametersModel.V_LN);
    float activatedDcClearance = model.parametersModel.c_dl * dcLN;
    result[0] = activatedDcMigration - activatedDcClearance;

    //T Cytotoxic
    float tCytoActivation = model.parametersModel.b_Tc * (model.parametersModel.rho_Tc*tCytoLN*dcLN - tCytoLN*tCytoLN*dcLN/model.parametersModel.estable_T_c);
    float tCytoHomeostasis = model.parametersModel.alpha_T_c * (model.parametersModel.estable_T_c - tCytoLN);
    float tCytoMigration = model.parametersModel.gamma_T * (tCytoLN - model.tCytotoxicTissueVessels) * model.parametersModel.V_BV/model.parametersModel.V_LN;
    result[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

    //T Helper
    float tHelperActivation = model.parametersModel.b_T * (model.parametersModel.rho_T * tHelperLN * dcLN - tHelperLN * dcLN);
    float tHelperHomeostasis = model.parametersModel.alpha_T_h * (model.parametersModel.estable_T_h - tHelperLN);
    float tHelperDispendure = model.parametersModel.b_rho * dcLN * tHelperLN * bCellLN;
    result[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

    //B Cell
    float bCellActivation = model.parametersModel.b_rho_b * (model.parametersModel.rho_B * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
    float bcellHomeostasis = model.parametersModel.alpha_B * (model.parametersModel.estable_B - bCellLN);
    result[3] = bcellHomeostasis + bCellActivation;

    //Plasma Cells
    float plasmaActivation = model.parametersModel.b_rho_p * (model.parametersModel.rho_P * tHelperLN * dcLN * bCellLN);
    float plasmaHomeostasis = model.parametersModel.alpha_P * (model.parametersModel.estable_P - plasmaCellLN);
    result[4] = plasmaHomeostasis + plasmaActivation;

    //Antibody
    float antibodyProduction = model.parametersModel.rho_F * plasmaCellLN;
    float antibodyMigration = model.parametersModel.gamma_F * (iggGLN - model.antibodyTissueVessels) * (model.parametersModel.V_BV/model.parametersModel.V_LN);
    result[5] = antibodyProduction - antibodyMigration;

    return result;
}

void SolverLymphNode(structModel *model, float step){
    float solutionLN[6];


    //Execute Euler (or RungeKutta4ThOrder)

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
    conventionalDcClearance, activatedDcClearance, activatedDcMigration, antibodyMigration;

    float microgliaKMinus, conventionalDcKMinus, activatedDcKMinus, tCytotoxicKMinus, antibodyKMinus, oligodendrocyteKMinus;

    for(int kTime = 1; kTime <= model->timeLen; kTime++){
        // solve lymphnode
        SolverLymphNode(model, (float)kTime*model->ht);
        
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
            microgliaReaction = model->parametersModel.mu_m*microgliaKMinus*(model->parametersModel.avgMic - microgliaKMinus);
            microgliaClearance = model->parametersModel.microgliaClearance*microgliaKMinus;

            model->microglia[stepKPlus][line][column] = microgliaKMinus + \
            model->ht*(microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance);

            //Conventional DC update
            conventionalDcReaction = model->parametersModel.mu_dc*oligodendrocyteKMinus*(model->parametersModel.avgDc - conventionalDcKMinus);
            conventionalDcActivation = model->parametersModel.b_d*conventionalDcKMinus;
            conventionalDcClearance = model->parametersModel.c_dc*conventionalDcKMinus;

            model->conventionalDc[stepKPlus][line][column] = conventionalDcKMinus + \
            model->ht*(conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation);
            
            //Activated DC update
            activatedDcClearance = model->parametersModel.c_da*activatedDcKMinus;
            activatedDcMigration = model->thetaPV[line][column]*model->parametersModel.gamma_D*(model->dendriticLymphNode[kTime] - activatedDcKMinus);
            
            model->activatedDc[stepKPlus][line][column] = activatedDcKMinus + model->ht*(activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance);
            
            //CD8 T update
            tCytotoxicMigration = model->thetaBV[line][column]*model->parametersModel.gamma_T*(model->tCytotoxicLymphNode[kTime] - tCytotoxicKMinus);
            
            model->tCytotoxic[stepKPlus][line][column] = tCytotoxicKMinus + model->ht*(tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration);
            

            //Antibody update
            odcAntibodyMicrogliaFagocitosis = model->parametersModel.lamb_f_m*antibodyKMinus*(model->parametersModel.avgOdc - oligodendrocyteKMinus)*fFunc(microgliaKMinus, model->parametersModel.avgMic);
            antibodyMigration = model->thetaBV[line][column]*model->parametersModel.gamma_F*(model->antibodyLymphNode[kTime] - antibodyKMinus);
            
            model->antibody[stepKPlus][line][column] = antibodyKMinus + model->ht*(antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis);

            //Oligodendrocytes update
            odcMicrogliaFagocitosis = model->parametersModel.r_m*fFunc(microgliaKMinus, model->parametersModel.avgMic)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);
            odcTCytotoxicApoptosis = model->parametersModel.r_t*fFunc(tCytotoxicKMinus, model->parametersModel.avgT)*(model->parametersModel.avgOdc - oligodendrocyteKMinus);

            model->oligodendrocyte[stepKPlus][line][column] = oligodendrocyteKMinus + model->ht*(odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis);
            //Save Results



        }   
        }
        stepKMinus += 1;
        stepKMinus = stepKMinus%2;
    }
}