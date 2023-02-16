#ifndef _MODEL_H_
#define _MODEL_H_

#define BUFFER 2
#define TIME 28
#define HT 0.0002
#define LENGTH 20
#define HX 0.5
#define NUMFIGS 28
#define NUMPOINTSLYMPHNODE 1000

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

    float avgT;
    float avgDc;
    float avgMic;
    float avgOdc;
    
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
    float estableTHelper;
    float estableTCytotoxic;
    float estableB;
    float estableP;
    int V_PV;
    int V_BV;
    int V_LN;
}structParameters;

typedef struct structModel
{
    int my_rank;
    int comm_sz;
    
    int numLines;
    int numPointsLN;
    int startLine;
    int endLine;
    

    float ht;
    float hx;

    int tFinal;
    int tSize;
    int xFinal;
    int xSize;
    int intervaloFiguras;
    int numFigs;

    // int *timeVec;
    // int *spaceVec;

    int timeLen;
    int spaceLen;

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

    float dendriticLymphNode[BUFFER];
    float tHelperLymphNode[BUFFER];
    float tCytotoxicLymphNode[BUFFER];
    float bCellLymphNode[BUFFER];
    float plasmaCellLymphNode[BUFFER];
    float antibodyLymphNode[BUFFER];

    float dendriticLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];
    float tHelperLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];
    float tCytotoxicLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];
    float bCellLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];
    float plasmaCellLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];
    float antibodyLymphNodeSavedPoints[NUMPOINTSLYMPHNODE];

    structParameters parametersModel;

}structModel;

int VerifyCFL(structParameters parametersModel, float ht, float hx);

float AdvectionTerm(float populationPoint, float avgValue);

float UpDownWind(float frontIPoint, float ijPoint, float avgValue);

float CalculateChemottaxis(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual,\
 float valor_medio, float gradiente_odc_i, float gradiente_odc_j, float hx);

float CalculateDiffusion(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual, float hx);

float fFunc(float valuePopulation, float avgPopulation);

float* EquationsLymphNode(structModel model, float* populationLN, int stepPos);

void SolverLymphNode(structModel *model, int stepPos);

structModel ModelInitialize(structParameters params, float ht, float hx, float time, float space, int numFigs, int numPointsLN, int my_rank, int comm_sz);

void DefineBVPV(structModel *model);

void InitialConditionTissueMicroglia(structModel* model);

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN);

void RunModel(structModel *model);

void WritePopulation(structModel model, float *population, char* fileName, char* bufferTime);

void WriteFiles(structModel model, float *oligodendrocyte, float *microglia, float *tCytotoxic, float *antibody, float *conventionalDC, float  *activatedDC, float time);
#endif