typedef struct structParameters
{
    float mu_m;
    float r_m;

    float chi;
    float d_mic;
    float d_dc;
    float d_da;
    float d_t_cit;
    float d_anti;
    
    float lamb_f_m;
    float b_d;
    float r_dc;
    float r_t;
    float microgliaClearance;
    float mu_dc;

    float gamma_D;
    float gamma_F;
    float gamma_T;

    float avgT;
    float avgDc;
    float avgMic;
    float avgOdc;
    
    float c_dc;
    float c_da;
    float c_dl;

    float alpha_T_h;
    float alpha_T_c;
    float alpha_B;
    float alpha_P;
    
    float b_T;
    float b_Tc;
    float b_rho;
    float b_rho_b;
    
    float rho_T;
    float rho_Tc;
    float rho_B;
    float rho_P;
    float b_rho_p;
    float rho_F;
    float estable_T_h;
    float estable_B;
    float estable_P;
    float estable_T_c;
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

float* EquationsLymphNode(float* populationLN, float step);

void SolverLymphNode(structModel *model, float step);

structModel ModelInitialize(structParameters params, int dt, int dx, int tFinal, int xFinal, int numFiguras);

void DefineBVPV(structModel *model);

void InitialConditionTissueMicroglia(structModel* model);

void InitialConditionLymphNode(structModel* model, float dendriticLN, float thelperLN, float tcytotoxicLN, float bcellLN, float plasmacellLN, float antibodyLN);

void RunModel(structModel *model);