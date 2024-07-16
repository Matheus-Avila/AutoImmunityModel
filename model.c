#include "model.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double auxAdcPV = 0.0;
double auxAntibodyBV = 0.0;
double auxTCytotoxicBV = 0.0;

double Min(double a, double b){
    return (a < b) ? a : b;
}

double Max(double a, double b){
    return (a > b) ? a : b;
}

void InitialConditionTissueMicroglia(structModel* model){

    for(int k = 0; k < model->xSize*model->xSize; k++){
        int i = (int)k/model->xSize;
        int j = k%model->xSize;
        if(pow((i-(int)(model->xSize/2)),2) + pow((j-(int)(model->xSize/2)),2) < 5 / (model->hx * model->hx)){
            model->microglia[0][k] = (double)model->parametersModel.avgMic/3;
        }
        
        //inicializar apenas com vaso sanguineo
        //model->tCytotoxic[0][k] = model->thetaBV[i] * 28.4;   
        
       
    }
    
}

void InitialConditionLymphNode(structModel* model, double dendriticLN, double thelperLN, double tcytotoxicLN, double bcellLN, double plasmacellLN, double antibodyLN){
    model->dendriticLymphNodeSavedPoints[0] = model->dendriticLymphNode[0] = dendriticLN;
    model->tHelperLymphNodeSavedPoints[0] = model->tHelperLymphNode[0] = thelperLN;
    model->tCytotoxicLymphNodeSavedPoints[0] = model->tCytotoxicLymphNode[0] = tcytotoxicLN;
    model->bCellLymphNodeSavedPoints[0] = model->bCellLymphNode[0] = bcellLN;
    model->plasmaCellLymphNodeSavedPoints[0] = model->plasmaCellLymphNode[0] = plasmacellLN;
    model->antibodyLymphNodeSavedPoints[0] = model->antibodyLymphNode[0] = antibodyLN;
}

int VerifyCFL(structParameters parametersModel, double ht, double hx){
    if(parametersModel.micDiffusion*ht/(hx*hx) < 0.25 && parametersModel.cDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.aDcDiffusion*ht/(hx*hx) < 0.25 && parametersModel.tCytoDiffusion*ht/(hx*hx) < 0.25 && parametersModel.chi*ht/hx < 0.5 && parametersModel.chi*ht/(hx*hx) < 0.25)
        return 1;
    return 0;
}

double CalculateMaxGradOdc(structModel model, int stepKMinus){
    double max = 0.0;
    double valIPlus;
    double valJPlus;
    double valIMinus;
    double valJMinus;
    double auxMax;
    int line, column;
    for(int kPos = 0; kPos < model.xSize*model.xSize; kPos++){
        line = (int)kPos/model.xSize;
        column = kPos%model.xSize;
        valIPlus = (line != model.xSize-1)? model.oligodendrocyte[stepKMinus][kPos + model.xSize]: model.oligodendrocyte[stepKMinus][kPos - model.xSize];
        valJPlus = (column != model.xSize-1)? model.oligodendrocyte[stepKMinus][kPos + 1]: model.oligodendrocyte[stepKMinus][kPos - 1];
        valIMinus = (line != 0)? model.oligodendrocyte[stepKMinus][kPos - model.xSize]: model.oligodendrocyte[stepKMinus][kPos + model.xSize];
        valJMinus = (column != 0)? model.oligodendrocyte[stepKMinus][kPos - 1]: model.oligodendrocyte[stepKMinus][kPos + 1];
        auxMax = Max((double)(valIPlus - valIMinus)/(double)(2 * model.hx), (double)(valJPlus - valJMinus)/(double)(2 * model.hx));
        if(auxMax > max)
            max = auxMax;
    }    
    return max;
}

double CalculateHt(structModel model, double stepKMinus){
    double maxGradOdc = CalculateMaxGradOdc(model, stepKMinus);
    return Min(model.hx * model.hx / (4 * model.parametersModel.antibodyDiffusion), (4 * model.parametersModel.antibodyDiffusion) / (model.parametersModel.chi * model.parametersModel.chi * maxGradOdc));
}

void WritePopulation(structModel model, double *population, char* fileName, char* bufferTime){
    FILE *file;
    file = fopen(fileName, "w");
    int k = 0;
    if(file != NULL){
        while (k < model.xSize*model.xSize){
            int i = k;
            while (i < k + model.xSize){
                fprintf(file, "%lf ", population[i]);
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

void WritePopulationLymphNode(structModel model, double *population, char* fileName){
    FILE *file;
    file = fopen(fileName, "a");
    if(file != NULL){
        for(int i=0;i<model.numPointsLN;i++){
            fprintf(file, "%lf\n", population[i]);
        }
        fclose(file);
    }else{
        printf("Error lymph node file\n");
        exit(0);
    }
    // FILE *dendritic;
    // dendritic = fopen("./result/dendritic.txt", "w");
    // FILE *tHelper;
    // tHelper = fopen("./result/tHelper.txt", "w");
    // FILE *tCyto;
    // tCyto = fopen("./result/tCyto.txt", "w");
    // FILE *bCell;
    // bCell = fopen("./result/bCell.txt", "w");
    // FILE *plasmaCell;
    // plasmaCell = fopen("./result/plasmaCell.txt", "w");
    // FILE *antibody;
    // antibody = fopen("./result/antibody.txt", "w");
    // int k = 0;
    // if(dendritic != NULL && tHelper != NULL && tCyto != NULL && bCell != NULL && plasmaCell != NULL && antibody != NULL){
    //     while (k < model->xSize*model->xSize){
    //         int i = k;
    //         while (i < k + model->xSize){
    //             fprintf(dendritic, "%lf\n ", population[i]);
    //             fprintf(tHelper, "%lf\n ", population[i]);
    //             fprintf(tCyto, "%lf\n ", population[i]);
    //             fprintf(bCell, "%lf\n ", population[i]);
    //             fprintf(plasmaCell, "%lf\n ", population[i]);
    //             fprintf(antibody, "%lf\n ", population[i]);
    //             i++;
    //         }
    //         // fprintf(dendritic,"\n");
    //         // fprintf(tHelper,"\n");
    //         // fprintf(tCyto,"\n");
    //         // fprintf(bCell,"\n");
    //         // fprintf(plasmaCell,"\n");
    //         // fprintf(antibody,"\n");
    //         // k+=model->xSize;
    //     }
    //     fclose(dendritic);
    //     fclose(tHelper);
    //     fclose(tCyto);
    //     fclose(bCell);
    //     fclose(plasmaCell);
    //     fclose(antibody);
    // }else{
    //     printf("Error matrix file\n");
    //     exit(0);
    // }
}

void WriteLymphNodeFiles(structModel model, double *dendritic, double *tHelper, double *tCytotoxic, double *bCell, double *plasmaCell, double *antibody){
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
    snprintf(buffer, sizeof(buffer), "%lf", (model.tSize/model.numPointsLN)*model.ht);
    strcat(command, buffer);
    system(command);
}

void WriteFiles(structModel model, double *oligodendrocyte, double *microglia, double *tCytotoxic, double *antibody, double *conventionalDC, double  *activatedDC, int Ktime){
    char buffer[10];
    int day = Ktime;
    
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

    //SavingData(model, Ktime);
}   

void PlotResults(structModel model){
    char buffer[10];
    char command[200] = {};
    strcat(command, "python3 plotMatrices.py ");
    snprintf(buffer, sizeof(buffer), "%d", model.xFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%lf", model.hx);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", model.tFinal);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%d", model.tSize/model.intervalFigures);
    strcat(command, buffer);
    system(command);
}


double PreventionOverCrowdingTerm(double populationPoint, double avgValue){
    return populationPoint/(populationPoint + avgValue);
}

double UpDownWind(double frontPoint, double rearPoint, double avgValue){
    return PreventionOverCrowdingTerm(frontPoint, avgValue) - PreventionOverCrowdingTerm(rearPoint, avgValue);
}

double CalculateChemotaxis(structModel model, double frontJPoint, double rearJPoint, double frontIPoint, double rearIPoint, double ijPoint,\
 double avgValue, double gradientOdcI, double gradientOdcJ){
    double gradientPopulationI, gradientPopulationJ;
    if(gradientOdcI<0)
        gradientPopulationI = UpDownWind(frontIPoint, ijPoint, avgValue)/(double)model.hx;
    else
        gradientPopulationI = UpDownWind(ijPoint, rearIPoint, avgValue)/(double)model.hx;
    if(gradientOdcJ<0)
        gradientPopulationJ = UpDownWind(frontJPoint, ijPoint, avgValue)/(double)model.hx;
    else
        gradientPopulationJ = UpDownWind(ijPoint, rearJPoint, avgValue)/(double)model.hx;
    return gradientOdcI*gradientPopulationI + gradientOdcJ*gradientPopulationJ;
}

double CalculateDiffusion(structModel model, double frontJPoint, double rearJPoint, double frontIPoint, double rearIPoint, double ijPoint){
    return (double)(frontIPoint + frontJPoint - 4*ijPoint + rearIPoint + rearJPoint)/(double)(model.hx*model.hx);
}

double fFunc(double valuePopulation, double avgPopulation){
    return valuePopulation*valuePopulation/(double)(valuePopulation + avgPopulation);
}
void WriteBVPV(structModel *model, double *thetaBV, double *thetaPV){
    FILE *fileBV;
    fileBV = fopen("./result/bv.txt", "w");
    FILE *filePV;
    filePV = fopen("./result/pv.txt", "w");
    int k = 0;
    if(fileBV != NULL && filePV != NULL){
        while (k < model->xSize*model->xSize){
            int i = k;
            while (i < k + model->xSize){
                fprintf(fileBV, "%lf ", thetaBV[i]);
                fprintf(filePV, "%lf ", thetaPV[i]);
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
    snprintf(buffer, sizeof(buffer), "%lf", model->hx);
    strcat(command, buffer);
    // system(command);
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
    printf("bv = %lf, pv = %lf \n", model->parametersModel.V_BV, model->parametersModel.V_PV);
    WriteBVPV(model, model->thetaBV, model->thetaPV);
}

structParameters ParametersInitialize();
void SavingData(structModel model, int Ktime){
    double totalCD8m = 0, totalMic = 0, totalODC = 0, totalCDC = 0, totalADC = 0, totalIGG = 0, totalCD8 = 0;
    FILE *file;
    structParameters parameters = ParametersInitialize();
    if(parameters.epslon_x  == 0){
        file = fopen("dataExecution0.txt", "a");
    }else if(parameters.epslon_x  == 0.55){
        file = fopen("dataExecution055.txt", "a");
    }else if(parameters.epslon_x  == 0.99){
        file = fopen("dataExecution099.txt", "a");
    }
    
    
    if(file == NULL){
            printf("dataExecution file not found!\n");
            exit(0);
    }
    for(int kPos = 0; kPos < model.xSize*model.xSize; kPos++){
        totalMic += model.microglia[0][kPos];
        totalODC += model.oligodendrocyte[0][kPos];
        totalCDC += model.conventionalDc[0][kPos];
        totalADC += model.activatedDc[0][kPos];
        totalCD8 += model.tCytotoxic[0][kPos];
        totalCD8m += (model.tCytotoxic[0][kPos] * model.thetaBV[kPos]); /*model.thetaBV[KPos]*/
        totalIGG += model.antibody[0][kPos];
    }
        //if(file != NULL){
        fprintf(file, "%lf\n", totalCD8/(model.xFinal*model.xFinal));
        
        fclose(file);

        FILE *file2;
        file2 = fopen("dataExecution2.txt", "a");
        
        if(file2 == NULL){
            printf("dataExecution file not found!\n");
            exit(0);
        }
            float media = totalCD8m/model.parametersModel.V_BV;
            fprintf(file2, "%lf\n", media);
            //fprintf(file, "%d\n", kPos);
            //fprintf(file, "Days = %d - Space = %d - ht = %lf, hx = %lf, Ht_JumpStep = %d\n", model.tFinal, model.xFinal, model.ht, model.hx, model.numStepsLN);
            //fprintf(file, "Lymph node populations\n");
            //fprintf(file, "DC = %lf, TCD8 = %lf, TCD4 = %lf, B Cell = %lf, Plasma cell = %lf, IgG = %lf\n", model.dendriticLymphNodeSavedPoints[model.numPointsLN-1], model.tCytotoxicLymphNodeSavedPoints[model.numPointsLN-1], model.tHelperLymphNodeSavedPoints[model.numPointsLN-1], model.bCellLymphNodeSavedPoints[model.numPointsLN-1], model.plasmaCellLymphNodeSavedPoints[model.numPointsLN-1], model.antibodyLymphNodeSavedPoints[model.numPointsLN-1]);
            //fprintf(file, "Tissue populations\n");
            //fprintf(file, "ODC = %lf, Microglia = %lf, ConventionalDC = %lf, ActivatedDC = %lf, TCD8 = %lf, IgG = %lf\n", totalODC, totalMic, totalCDC, totalADC, totalCD8, totalIGG);    
            //fprintf(file, "Parameters\n");
            /*fprintf(file, "micDiffusion  = %lf, antibodyDiffusion = %lf, cDcDiffusion = %lf, aDcDiffusion = %lf, tCytoDiffusion = %lf, chi = %lf, muCDc = %lf, muMic = %lf, \
            rM = %lf, rT = %lf, lambAntMic = %lf, bD = %lf, gammaD = %lf, gammaAntibody = %lf, gammaT = %lf,  avgT = %lf, avgDc = %lf, avgMic = %lf, avgOdc = %lf,  cMic = %lf, \
            cCDc = %lf, cADc = %lf, cDl = %lf, cF = %lf, alphaTHelper = %lf, alphaTCytotoxic = %lf, alphaB = %lf, alphaP = %lf, bTHelper = %lf, bTCytotoxic = %lf, bRho = %lf, \
            bRhoB = %lf, bRhoP = %lf, rhoTHelper = %lf, rhoTCytotoxic = %lf, rhoB = %lf, rhoP = %lf, rhoAntibody = %lf, stableTHelper = %lf, stableTCytotoxic = %lf, \
            stableB = %lf, stableP = %lf, V_LN = %d, V_BV = %lf, V_PV = %lf\n",
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
            */
        
        fclose(file2);
    
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

    model.microglia = (double**)malloc(BUFFER * sizeof(double*));
    model.oligodendrocyte = (double**)malloc(BUFFER * sizeof(double*));
    model.tCytotoxic = (double**)malloc(BUFFER * sizeof(double*));
    model.antibody = (double**)malloc(BUFFER * sizeof(double*));
    model.conventionalDc = (double**)malloc(BUFFER * sizeof(double*));
    model.activatedDc = (double**)malloc(BUFFER * sizeof(double*));

    model.activatedDCTissueVessels = 0;
    model.tCytotoxicTissueVessels = 0;
    model.antibodyTissueVessels = 0;

    for (int index=0;index<BUFFER;++index){
        model.microglia[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
        model.oligodendrocyte[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
        model.tCytotoxic[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
        model.conventionalDc[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
        model.activatedDc[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
        model.antibody[index] = (double*)calloc(model.xSize*model.xSize, sizeof(double));
    }
    //definir BV e PV
    model.thetaPV = (double*)calloc(model.xSize*model.xSize, sizeof(double));
    model.thetaBV = (double*)calloc(model.xSize*model.xSize, sizeof(double));
    DefineBVPV(&model);
    //definir lymph node
    model.dendriticLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));
    model.tCytotoxicLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));
    model.tHelperLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));
    model.antibodyLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));
    model.bCellLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));
    model.plasmaCellLymphNodeSavedPoints = (double*)calloc(model.numPointsLN, sizeof(double));

    model.dendriticLymphNode = (double*)calloc(2, sizeof(double));
    model.tCytotoxicLymphNode = (double*)calloc(2, sizeof(double));
    model.tHelperLymphNode = (double*)calloc(2, sizeof(double));
    model.antibodyLymphNode = (double*)calloc(2, sizeof(double));
    model.bCellLymphNode = (double*)calloc(2, sizeof(double));
    model.plasmaCellLymphNode = (double*)calloc(2, sizeof(double));

    double dendriticLN = 0.0, thelperLN = params.stableTHelper, tcytotoxicLN = params.stableTCytotoxic, bcellLN = params.stableB, plasmacellLN = 0.0, antibodyLN = 0.0;
    InitialConditionLymphNode(&model, dendriticLN, thelperLN, tcytotoxicLN, bcellLN, plasmacellLN, antibodyLN);
    InitialConditionTissueMicroglia(&model);
    return model;
}

/*
* Lymphnode
*/
double* EquationsLymphNode(structModel* model, int stepPos){
    double* slopes =(double *)malloc(sizeof(double)*6);
    
    int stepKPlus = (stepPos%(2*model->numStepsLN))/model->numStepsLN;
    int stepKMinus = !(stepKPlus && 1);

    double dcLN = model->dendriticLymphNode[stepKMinus];
    double tCytoLN = model->tCytotoxicLymphNode[stepKMinus];
    double tHelperLN = model->tHelperLymphNode[stepKMinus];
    double bCellLN = model->bCellLymphNode[stepKMinus];
    double plasmaCellLN = model->plasmaCellLymphNode[stepKMinus];
    double antibodyLN = model->antibodyLymphNode[stepKMinus];

    //Describe equations

    //Dendritic cell
    double activatedDcMigration = model->parametersModel.gammaD * (model->activatedDCTissueVessels - dcLN) * (double)(model->parametersModel.V_PV/model->parametersModel.V_LN);
    double activatedDcClearance = model->parametersModel.cDl * dcLN;
    slopes[0] = activatedDcMigration - activatedDcClearance;

    //T Cytotoxic
    //printf("%.2lf\n",model->parametersModel.eps_new);
    double tCytoActivation = model->parametersModel.bTCytotoxic * (model->parametersModel.rhoTCytotoxic*tCytoLN*dcLN - tCytoLN*dcLN);
    double tCytoHomeostasis = model->parametersModel.alphaTCytotoxic * (model->parametersModel.stableTCytotoxic - tCytoLN);
    double tCytoMigration = model->parametersModel.gammaT * (tCytoLN - model->tCytotoxicTissueVessels) * (double)(model->parametersModel.V_BV/model->parametersModel.V_LN)* (1 - model->parametersModel.eps_new);
    slopes[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

    //T Helper
    double tHelperActivation = model->parametersModel.bTHelper * (model->parametersModel.rhoTHelper * tHelperLN * dcLN - tHelperLN * dcLN);
    double tHelperHomeostasis = model->parametersModel.alphaTHelper * (model->parametersModel.stableTHelper - tHelperLN);
    double tHelperDispendure = model->parametersModel.bRho * dcLN * tHelperLN * bCellLN;
    slopes[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

    //B Cell
    double bCellActivation = model->parametersModel.bRhoB * (model->parametersModel.rhoB * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
    double bcellHomeostasis = model->parametersModel.alphaB * (model->parametersModel.stableB - bCellLN);
    slopes[3] = bcellHomeostasis + bCellActivation;

    //Plasma Cells
    double plasmaActivation = model->parametersModel.bRhoP * (model->parametersModel.rhoP * tHelperLN * dcLN * bCellLN);
    double plasmaHomeostasis = model->parametersModel.alphaP * (model->parametersModel.stableP - plasmaCellLN);
    slopes[4] = plasmaHomeostasis + plasmaActivation;

    //Antibody
    double antibodyProduction = model->parametersModel.rhoAntibody * plasmaCellLN;
    double antibodyDecayment = model->parametersModel.cF * antibodyLN;
    double antibodyMigration = model->parametersModel.gammaAntibody * (antibodyLN - model->antibodyTissueVessels) * (double)(model->parametersModel.V_BV/model->parametersModel.V_LN);
    slopes[5] = antibodyProduction - antibodyMigration - antibodyDecayment;

    return slopes;
}

/*
void SolverLymphNode(structModel *model, int stepPos){
    double populationLN[6];
    int stepKPlus = (stepPos%(2*model->numStepsLN))/model->numStepsLN;
    int stepKMinus = !(stepKPlus && 1);
    populationLN[0] = model->dendriticLymphNode[stepKMinus];
    populationLN[1] = model->tCytotoxicLymphNode[stepKMinus];
    populationLN[2] = model->tHelperLymphNode[stepKMinus];
    populationLN[3] = model->bCellLymphNode[stepKMinus];
    populationLN[4] = model->plasmaCellLymphNode[stepKMinus];
    populationLN[5] = model->antibodyLymphNode[stepKMinus];
    
    double* slopeLN;
    slopeLN = EquationsLymphNode(*model, populationLN, stepPos);

    double htLN = model->ht*model->numStepsLN;
    
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
derivatives* SlopePDEs(int stepKPlus, double ht, structModel* model){

    derivatives* slopes = (derivatives*)calloc(1, sizeof(derivatives));
    slopes->derivativesLymphNode = (double*)calloc(6, sizeof(double));
    slopes->derivativesTissue = (double**)calloc(6, sizeof(double*));
    
    
    for(int i = 0; i < 6; i++)
        slopes->derivativesTissue[i] = (double*)calloc(model->xSize*model->xSize, sizeof(double));
    /*
     * Solve slope PDEs
    */
    int stepKMinus = !(stepKPlus && 1), line, column;
    
    double upperNeumannBC = 0.0, lowerNeumannBC = 0.0, leftNeumannBC = 0.0, rightNeumannBC = 0.0;
    
    double valIPlus = 0.0, valIMinus = 0.0, valJPlus = 0.0, valJMinus = 0.0, gradientOdcI = 0.0, gradientOdcJ = 0.0;

    double microgliaChemotaxis = 0.0, tCytotoxicChemotaxis = 0.0, conventionalDcChemotaxis = 0.0,\
     microgliaDiffusion = 0.0, tCytotoxicDiffusion = 0.0, conventionalDcDiffusion = 0.0, activatedDCDiffusion = 0.0, antibodyDiffusion = 0.0;

    double microgliaReaction = 0.0, microgliaClearance = 0.0, tCytotoxicMigration = 0.0, odcAntibodyMicrogliaFagocitosis = 0.0, \
    odcMicrogliaFagocitosis = 0.0, odcTCytotoxicApoptosis = 0.0, conventionalDcReaction = 0.0, conventionalDcClearance = 0.0, conventionalDcActivation = 0.0, \
    activatedDcClearance = 0.0, activatedDcMigration = 0.0, antibodyMigration = 0.0;

    double microgliaKMinus = 0.0, conventionalDcKMinus = 0.0, activatedDcKMinus = 0.0, tCytotoxicKMinus = 0.0, antibodyKMinus = 0.0, oligodendrocyteKMinus = 0.0;

    double diffusionOdc;
    //float eps_new;
    
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

            gradientOdcI = (float)(valIPlus - valIMinus)/(float)(2 * model->hx);
            gradientOdcJ = (float)(valJPlus - valJMinus)/(float)(2 * model->hx);

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

            valIPlus  = (line != model->xSize-1)? model->activatedDc[stepKMinus][kPos + model->xSize]: model->activatedDc[stepKMinus][kPos - model->xSize] - (double)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->activatedDc[stepKMinus][kPos + 1]: model->activatedDc[stepKMinus][kPos - 1] - (double)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->activatedDc[stepKMinus][kPos - model->xSize]: model->activatedDc[stepKMinus][kPos + model->xSize] - (double)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->activatedDc[stepKMinus][kPos - 1]: model->activatedDc[stepKMinus][kPos + 1] - (double)(2*model->hx*leftNeumannBC);

            activatedDCDiffusion = model->parametersModel.aDcDiffusion*CalculateDiffusion(*model, valJPlus, valJMinus, valIPlus, valIMinus, model->activatedDc[stepKMinus][kPos]);

            //Difussion Antibody

            valIPlus  = (line != model->xSize-1)? model->antibody[stepKMinus][kPos + model->xSize]: model->antibody[stepKMinus][kPos - model->xSize] - (double)(2*model->hx*lowerNeumannBC);
            valJPlus  = (column != model->xSize-1)? model->antibody[stepKMinus][kPos + 1]: model->antibody[stepKMinus][kPos - 1] - (double)(2*model->hx*rightNeumannBC);
            valIMinus = (line != 0)? model->antibody[stepKMinus][kPos - model->xSize]: model->antibody[stepKMinus][kPos + model->xSize] - (double)(2*model->hx*upperNeumannBC);
            valJMinus = (column != 0)? model->antibody[stepKMinus][kPos - 1]: model->antibody[stepKMinus][kPos + 1] - (double)(2*model->hx*leftNeumannBC);

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
        
        tCytotoxicMigration = model->thetaBV[kPos]*model->parametersModel.gammaT*(model->tCytotoxicLymphNode[stepKPlus] - tCytotoxicKMinus) * (1 - model->parametersModel.eps_new) ;
        
        slopes->derivativesTissue[3][kPos] = tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration - (model->parametersModel.ct * tCytotoxicKMinus);
        
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

double Euler(double time, structModel *model, int stepKPlus, int *posSave){
    derivatives* slopeK;
    int stepKMinus = (stepKPlus + 1) % 2;
    double htCalculated = Min(CalculateHt(*model, stepKMinus), 0.25 * (model->hx*model->hx) / model->parametersModel.chi) * .01;
    // printf("Tempo %lf: ht dinamico = %.16lf\n", time, htCalculated);
    
    slopeK = SlopePDEs(stepKPlus, model->ht, model);
    model->tCytotoxicTissueVessels = auxTCytotoxicBV * model->hx * model->hx / model->parametersModel.V_BV;
    model->antibodyTissueVessels = auxAntibodyBV * model->hx * model->hx / model->parametersModel.V_BV;
    model->activatedDCTissueVessels = auxAdcPV * model->hx * model->hx / model->parametersModel.V_PV;
    double* slopeLN = EquationsLymphNode(model, time);
    slopeK->derivativesLymphNode = slopeLN;

    model->dendriticLymphNode[stepKPlus] = model->dendriticLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[0];
    model->tCytotoxicLymphNode[stepKPlus] = model->tCytotoxicLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[1];
    model->tHelperLymphNode[stepKPlus] = model->tHelperLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[2];
    model->bCellLymphNode[stepKPlus] = model->bCellLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[3];
    model->plasmaCellLymphNode[stepKPlus] = model->plasmaCellLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[4];
    model->antibodyLymphNode[stepKPlus] = model->antibodyLymphNode[stepKMinus] + htCalculated * slopeK->derivativesLymphNode[5];
    
    for(int spacePoint = 0; spacePoint < model->xSize * model->xSize; spacePoint++){
        model->microglia[stepKPlus][spacePoint] = model->microglia[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[0][spacePoint];
        model->conventionalDc[stepKPlus][spacePoint] = model->conventionalDc[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[1][spacePoint];
        model->activatedDc[stepKPlus][spacePoint] = model->activatedDc[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[2][spacePoint];
        model->tCytotoxic[stepKPlus][spacePoint] = model->tCytotoxic[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[3][spacePoint];
        model->antibody[stepKPlus][spacePoint] = model->antibody[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[4][spacePoint];
        model->oligodendrocyte[stepKPlus][spacePoint] = model->oligodendrocyte[stepKMinus][spacePoint] + htCalculated * slopeK->derivativesTissue[5][spacePoint];
    }

    for(int i = 0; i < 6; i++)
        free(slopeK->derivativesTissue[i]);
    free(slopeK->derivativesLymphNode);
    free(slopeK->derivativesTissue);
    free(slopeK);

    auxTCytotoxicBV = 0.0;
    auxAntibodyBV = 0.0;
    auxAdcPV = 0.0;
    double intervalPoints = (double) ( (double) model->tFinal / (double) model->numPointsLN);
    if( time > intervalPoints * *posSave){
        model->dendriticLymphNodeSavedPoints[*posSave] = model->dendriticLymphNode[stepKPlus];
        model->tCytotoxicLymphNodeSavedPoints[*posSave] = model->tCytotoxicLymphNode[stepKPlus];
        model->tHelperLymphNodeSavedPoints[*posSave] = model->tHelperLymphNode[stepKPlus];
        model->bCellLymphNodeSavedPoints[*posSave] = model->bCellLymphNode[stepKPlus];
        model->plasmaCellLymphNodeSavedPoints[*posSave] = model->plasmaCellLymphNode[stepKPlus];
        model->antibodyLymphNodeSavedPoints[*posSave] = model->antibodyLymphNode[stepKPlus];
        *posSave = *posSave + 1;
        // printf("%d\n", *posSave);
    }

    return htCalculated;

}

void RungeKutta(int kTime, structModel *model){
    derivatives* slopeK1, *slopeK2, *slopeK3 , *slopeK4;
    double htLymphNode = model->ht * model->numStepsLN;

    double** yK;
    yK = (double**) calloc(6, sizeof(double*));
    for(int i = 0; i < 6; i++)
        yK[i] = (double*) calloc(model->xSize * model->xSize, sizeof(double));

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
        double* slopeLN = EquationsLymphNode(model, kTime);
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

int isIn(int ktime, int vec[], int size) {
    for(int i = 0; i < size; i++) {
        if(ktime - 1 == vec[i]) {
            //printf("entrou \n");
            return 1;
        } 
    }
    return 0;
}
/*
 * Central Nervous System - PDEs - Finite Differences
 */

float RunModel(structModel *model, int* save_times, int size, float* points_values, int* targetDays,int targetSize){

    int stepKPlus = 1, stepKMinus = 0;
    int posSave = 1;

    float sum = 0.f;
    int currentIndex = 0;
    //Save IC
    if(model->saveFigs)
        WriteFiles(*model, model->oligodendrocyte[0], model->microglia[0], model->tCytotoxic[0], model->antibody[0], model->conventionalDc[0], model->activatedDc[0], 0);
        SavingData(*model, 0);
    int kTime = 1; 
    double time = 0.0;
    double htDynamic = 0.0;
    int days = 0;
    //double eps_new;
   

    while(time < model->tFinal){

        // if(time <= 30){
        //     model->parametersModel.eps_new = 0;
        // }else if(time > 30){
        //     model->parametersModel.eps_new = model->parametersModel.epslon_x;
        // }

        htDynamic = Euler(time, model, stepKPlus, &posSave);
        time += htDynamic;
        if(!((int)time - kTime)){
            if(model->saveFigs) {
                WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);
                SavingData(*model, kTime);
            }

            days++;
            // printf("%d!!%lf\n", kTime, time);
            kTime++;
        }

        // Verifique se `days` é igual a qualquer um dos dias em `targetDays`
        for (int i = 0; i < targetSize; ++i) {
            if (days == targetDays[i]) {
                sum += (model->tCytotoxicLymphNode[stepKPlus] - points_values[i]);
            }
        }

        stepKMinus = stepKPlus;
        stepKPlus = !stepKMinus;
        // printf("%lf!!\n", time);
        // printf("%d \n", kTime);   
    }    
        
    WriteFiles(*model, model->oligodendrocyte[stepKPlus], model->microglia[stepKPlus], model->tCytotoxic[stepKPlus], model->antibody[stepKPlus], model->conventionalDc[stepKPlus], model->activatedDc[stepKPlus], kTime);

    printf("Computation Done!!\n");
    SavingData(*model, kTime);
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

    //sum = (model->tCytotoxicLymphNode[stepKPlus] - points_values[0]);

    free(model->dendriticLymphNodeSavedPoints);
    free(model->tCytotoxicLymphNodeSavedPoints);
    free(model->tHelperLymphNodeSavedPoints);
    free(model->antibodyLymphNodeSavedPoints);
    free(model->bCellLymphNodeSavedPoints);
    free(model->plasmaCellLymphNodeSavedPoints);

    
    return sum;
}