#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "model.h"

#define DT 0.0002
#define DX 0.5
#define TFINAL 28
#define XFINAL 20
#define NUMFIGURAS 28

structParameters ParametersInitialize(){
    structParameters params;
    params.aDcDiffusion = 60*24/pow(2.5,2)*6.6*pow(10,-5);
    params.antibodyDiffusion = 60*24/pow(2.5,2)*6.6*pow(10,-5);
    params.cDcDiffusion = 60*24/pow(2.5,2)*6.6*pow(10,-5);
    params.aDcDiffusion = 60*24/pow(2.5,2)*6.6*pow(10,-5);
    params.tCytoDiffusion = 60*24/pow(2.5,2)*6.6*pow(10,-5);
    params.chi = 0.298*60*2;
    
    params.muCDc = 60*24*3*pow(10,-4);
    params.muMic = 60*24*3*pow(10,-6);
    params.rM = 60*24*3.96*pow(10,-6);
    params.rT = 0.1;
    params.lambAntMic = 5.702*pow(10,-3);
    params.bD = 0.001;
    
    params.gammaD = 0.01;
    params.gammaAntibody = 0.3;
    params.gammaT = 0.2;

    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;

    params.cMic = 0.1;
    params.cCDc = 0.1;
    params.cADc = 0.1;
    params.cDl = 0.1;
    params.alphaTHelper = 0.1;
    params.alphaTCytotoxic = 0.1;
    params.alphaB = 0.1;
    params.alphaP = 0.1;
    params.bTHelper = 0.17;
    params.bTCytotoxic = 0.017;
    params.bRho = 0.6;
    params.bRhoB = 3.02;
    params.rhoTHelper = 2;
    params.rhoTCytotoxic = 2;
    params.rhoB = 11;
    params.rhoP = 3;
    params.bRhoP = 1.02;
    params.rhoAntibody = pow(5.1,-2);
    params.estableTHelper = 84;
    params.estableTCytotoxic = 40;
    params.estableB = 25;
    params.estableP = 2.5;
    params.V_LN = 160;

    return params;
}


int main(int argc, char* argv){
    structParameters parameters = ParametersInitialize();

    structModel model = ModelInitialize(parameters, DT, DX, TFINAL, XFINAL, NUMFIGURAS);
    RunModel(&model);

    return 0;
}