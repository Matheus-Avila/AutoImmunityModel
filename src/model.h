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

    float alphaTHelper;
    float alphaTCytotoxic;
    float alphaB;
    float alphaP;
    
    float bTHelper;
    float bTCytotoxic;
    float bRho;
    float bRhoB;
    
    float rhoTHelper;
    float rhoTCytotoxic;
    float rhoB;
    float rhoP;
    float bRhoP;
    float rhoAntibody;
    float estableTHelper;
    float estableB;
    float estableP;
    float estableTCytotoxic;
    int V_LV;
    int V_BV;
    int V_LN;
}structParameters;

typedef struct structModel
{
    float ht;
    float hx;

    int tFinal;
    int xFinal;

    int numFiguras;

    // int *timeVec;
    // int *spaceVec;

    int timeLen;
    int spaceLen;

    float **thetaBV;
    float **thetaPV;

    float activatedDCTissueVessels;
    float antibodyTissueVessels;
    float tCytotoxicTissueVessels;
    
    float ***microglia;
    float ***oligodendrocyte;
    float ***conventionalDc;
    float ***activatedDc;
    float ***antibody;
    float ***tCytotoxic;

    float dendriticLymphNode[2];
    float tHelperLymphNode[2];
    float tCytotoxicLymphNode[2];
    float bCellLymphNode[2];
    float plasmaCellLymphNode[2];
    float antibodyLymphNode[2];

    structParameters parametersModel;

}structModel;

int VerifyCFL(structParameters parametersModel);

float AdvectionTerm(float populationPoint, float avgValue);

float UpDownWind(float frontIPoint, float ijPoint, float avgValue);

float CalculateChemottaxis(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual,\
 float valor_medio, float gradiente_odc_i, float gradiente_odc_j, float hx);

float CalculateDiffusion(float ponto_posterior_j, float ponto_anterior_j, float ponto_posterior_i, float ponto_anterior_i, float ponto_atual, float hx);

float fFunc(float valuePopulation, float avgPopulation);

float* EquationsLymphNode(structModel model, float* populationLN, int stepPos);

void SolverLymphNode(structModel *model, int stepPos);

structModel ModelInitialize(structParameters params, int dt, int dx, int tFinal, int xFinal, int numFiguras);

void DefineBVPV(structModel *model);

void InitialConditionTissueMicroglia(structModel* model);

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN);

void RunModel(structModel *model);