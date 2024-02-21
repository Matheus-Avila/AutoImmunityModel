#ifndef _MODELOMP_H_
#define _MODELOMP_H_

#define BUFFER 2

typedef struct structParameters
{
    float muMic;
    float rM;

    float chi;
    float micDiffusion;
    float cDcDiffusion;
    float aDcDiffusion;
    float tCytoDiffusion;
    float antibodyDiffusion;
    
    float lambAntMic;
    float bD;
    float rT;
    float cMic;
    float muCDc;

    float gammaD;
    float gammaAntibody;
    float gammaT;
    
    float cCDc;
    float cADc;
    float cDl;
    float cF;

    float alphaTHelper;
    float alphaTCytotoxic;
    float alphaB;
    float alphaP;
    
    float bTHelper;
    float bTCytotoxic;
    float bRho;
    float bRhoB;
    float bRhoP;
    
    float rhoTHelper;
    float rhoTCytotoxic;
    float rhoB;
    float rhoP;
    float rhoAntibody;

    float avgT;
    float avgDc;
    float avgMic;
    float avgOdc;

    float stableTHelper;
    float stableTCytotoxic;
    float stableB;
    float stableP;
    float V_PV;
    float V_BV;
    int V_LN;
}structParameters;

typedef struct structModel
{
    int totalThreads;

    float ht;
    float hx;

    int tFinal;
    int xFinal;
    
    int numStepsLN;
    int intervalFigures;
    int numPointsLN;
    int numFigs;
    int saveFigs;

    int tSize;
    int xSize;

    float *thetaBV;
    float *thetaPV;

    float activatedDCTissueVessels;
    float antibodyTissueVessels;
    float tCytotoxicTissueVessels;
    
    float **microglia;
    float **oligodendrocyte;
    float **conventionalDc;
    float **activatedDc;
    float **antibody;
    float **tCytotoxic;

    float *dendriticLymphNode;
    float *tHelperLymphNode;
    float *tCytotoxicLymphNode;
    float *bCellLymphNode;
    float *plasmaCellLymphNode;
    float *antibodyLymphNode;

    float *dendriticLymphNodeSavedPoints;
    float *tHelperLymphNodeSavedPoints;
    float *tCytotoxicLymphNodeSavedPoints;
    float *bCellLymphNodeSavedPoints;
    float *plasmaCellLymphNodeSavedPoints;
    float *antibodyLymphNodeSavedPoints;

    structParameters parametersModel;

}structModel;

int VerifyCFL(structParameters parametersModel, float ht, float hx);

float PreventionOverCrowdingTerm(float populationPoint, float avgValue);

float UpDownWind(float frontIPoint, float ijPoint, float avgValue);

float CalculateChemottaxis(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint,\
 float avgValue, float gradientOdcI, float gradientOdcJ, float hx);

float CalculateDiffusion(float frontJPoint, float rearJPoint, float frontIPoint, float rearIPoint, float ijPoint, float hx);

float fFunc(float valuePopulation, float avgPopulation);

float* EquationsLymphNode(structModel model, float* populationLN, int stepPos);

void SolverLymphNode(structModel *model, int stepPos);

structModel ModelInitialize(structParameters params, int totThr, float ht, float hx, float time, float space, int numFigs, int numPointsLN, int numStepsLN, int saveFigs);

void DefineBVPV(structModel *model);

void InitialConditionTissueMicroglia(structModel* model);

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN);

void RunModel(structModel *model);

void WritePopulation(structModel model, float *population, char* fileName, char* bufferTime);

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float  *activatedDC, float time);
#endif