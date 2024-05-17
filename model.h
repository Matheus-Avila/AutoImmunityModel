#ifndef _MODEL_H_
#define _MODEL_H_

#define BUFFER 2
#include <unistd.h> 

typedef struct derivatives
{
    double** derivativesTissue;
    double* derivativesLymphNode;
}derivatives;

typedef struct structParameters
{
    double muMic;
    double rM;

    double chi;
    double micDiffusion;
    double cDcDiffusion;
    double aDcDiffusion;
    double tCytoDiffusion;
    double antibodyDiffusion;
    
    double lambAntMic;
    double bD;
    double rT;
    double cMic;
    double muCDc;

    double gammaD;
    double gammaAntibody;
    double gammaT;

    double avgT;
    double avgDc;
    double avgMic;
    double avgOdc;
    
    double cCDc;
    double cADc;
    double cDl;
    double cF;

    double alphaTHelper;
    double alphaTCytotoxic;
    double alphaB;
    double alphaP;
    
    double bTHelper;
    double bTCytotoxic;
    double bRho;
    double bRhoB;
    double bRhoP;
    
    double rhoTHelper;
    double rhoTCytotoxic;
    double rhoB;
    double rhoP;
    double rhoAntibody;
    double stableTHelper;
    double stableTCytotoxic;
    double stableB;
    double stableP;
    double V_PV;
    double V_BV;
    int V_LN;
    double epslon_x;

}structParameters;

typedef struct structModel
{
    double ht;
    double hx;

    int tFinal;
    int xFinal;
    
    int numStepsLN;
    int intervalFigures;
    int numPointsLN;
    int numFigs;
    int saveFigs;

    int tSize;
    int xSize;

    double *thetaBV;
    double *thetaPV;

    double activatedDCTissueVessels;
    double antibodyTissueVessels;
    double tCytotoxicTissueVessels;
    
    double **microglia;
    double **oligodendrocyte;
    double **conventionalDc;
    double **activatedDc;
    double **antibody;
    double **tCytotoxic;

    double *dendriticLymphNode;
    double *tHelperLymphNode;
    double *tCytotoxicLymphNode;
    double *bCellLymphNode;
    double *plasmaCellLymphNode;
    double *antibodyLymphNode;

    double *dendriticLymphNodeSavedPoints;
    double *tHelperLymphNodeSavedPoints;
    double *tCytotoxicLymphNodeSavedPoints;
    double *bCellLymphNodeSavedPoints;
    double *plasmaCellLymphNodeSavedPoints;
    double *antibodyLymphNodeSavedPoints;

    double eps_new;

    structParameters parametersModel;

}structModel;

int VerifyCFL(structParameters parametersModel, double ht, double hx);

double PreventionOverCrowdingTerm(double populationPoint, double avgValue);

double UpDownWind(double frontIPoint, double ijPoint, double avgValue);

double CalculateChemotaxis(structModel model, double ponto_posterior_j, double ponto_anterior_j, double ponto_posterior_i, double ponto_anterior_i, double ponto_atual,\
 double valor_medio, double gradiente_odc_i, double gradiente_odc_j);

double CalculateDiffusion(structModel model, double ponto_posterior_j, double ponto_anterior_j, double ponto_posterior_i, double ponto_anterior_i, double ponto_atual);

double fFunc(double valuePopulation, double avgPopulation);

double* EquationsLymphNode(structModel* model, int stepPos);

double CalculateHt(structModel model, double stepKMinus);

void SolverLymphNode(structModel *model, int stepPos);

structModel ModelInitialize(structParameters params, float ht, float hx, float time, float space, int numFigs, int numPointsLN, int numStepsLN, int saveFigs);

void DefineBVPV(structModel *model);

void InitialConditionTissueMicroglia(structModel* model);

void InitialConditionLymphNode(structModel* model, double dendriticLN, double thelperLN, double tcytotoxicLN, double bcellLN, double plasmacellLN, double antibodyLN);

float RunModel(structModel *model,int* save_times, int size, float* points_values);

void WritePopulation(structModel model, double *population, char* fileName, char* bufferTime);

void WriteFiles(structModel model, double *oligodendrocyte, double *microglia, double *tCytotoxic, double *antibody, double *conventionalDC, double  *activatedDC, int time);

#endif