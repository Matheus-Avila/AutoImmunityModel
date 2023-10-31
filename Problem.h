#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include "model.h"
#include <fstream>
#include <string>
#include <float.h>
#include "math.h"

structParameters ParametersInitialize(){
    structParameters params;
    params.micDiffusion = 0.015206;
    params.antibodyDiffusion = 0.15206;
    params.cDcDiffusion = 0.015206;
    params.aDcDiffusion = 0.015206;
    params.tCytoDiffusion = 0.015206;
    params.chi = 0.03;
    
    params.muCDc = 60*24*3*pow(10,-5);
    params.muMic = 60*24*3*pow(10,-6);
    params.rM = 60*24*6*pow(10,-7);
    params.rT = 0.001;
    params.lambAntMic = 5.702*pow(10,-3);
    params.bD = 0.001;
    
    params.gammaD = 0.1;
    params.gammaAntibody = 0.3;
    params.gammaT = 0.1;

    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;

    params.cMic = 0.1;
    params.cCDc = 1;
    params.cADc = 1;
    params.cDl = 0.1;
    params.cF = 0.1;
    params.alphaTHelper = 0.1;
    params.alphaTCytotoxic = 0.00995612;
    params.alphaB = 0.1;
    params.alphaP = 1;
    params.bTHelper = 0.17;
    params.bTCytotoxic = 0.00197612;
    params.bRho = 0.6;
    params.bRhoB = 3.02;
    params.bRhoP = 1.02;
    params.rhoTHelper = 2;
    params.rhoTCytotoxic = 2;
    params.rhoB = 11;
    params.rhoP = 3;
    params.rhoAntibody = 5.1*pow(10,-2);
    params.stableTHelper = 70;
    params.stableTCytotoxic = 28.4;
    params.stableB = 25;
    params.stableP = 2.5;
    params.V_LN = 40;
    params.V_BV = 0;
    params.V_PV = 0;
    params.epslon_x = 0.457328;

    return params;
}

float sumVector(float vector[], int sizeX, int sizeY) {
    float sum = 0.0;
    
  
    for (int i = 0; i < sizeX; i++) {
        for(int j = 0; j < sizeY; j++) {
            sum += vector[i * sizeX + j];
        }
    }
  
    return sum;
}

int appendOnFile(const char* filename, double value) {
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir o arquivo." << std::endl;
        return 1;
    }

    file << value << std::endl;

    file.close();
    return 0;
}

float sumOneDimensionVector(float vector[], int si) {
    float sum = 0.0;
    for(int i = 0; i < si; i++) {
      sum += vector[i];
    }
    return sum;
}

using namespace pagmo;

class MSProblemTCytoParams : public pagmo::problem {
    public:
        MSProblemTCytoParams() : pagmo::problem() {}
        std::string get_name() const {
            return "Multiple Sclerosis Problem with Tcyto alpha and beta";
        }

        vector_double::size_type get_nx() const {
            return 0;
        }

        vector_double::size_type get_nec() const {
            return 0;
        }

        vector_double::size_type get_nic() const {
            return 0;
        }

        vector_double::size_type get_ncx() const {
            return 2;
        }

        std::pair<vector_double, vector_double> get_bounds() const
        {
            // 
            return {{0.0001, 0.0001, 0.1, 0.05}, {0.5, 0.5, 0.99, 0.8}};
        }

        vector_double fitness(const vector_double& variables) const {
            double alpha = variables[0];
            double beta = variables[1];
            double epslon = variables[2];
            double gammaT = variables[3];
            float ht = 0.0002, hx = 0.5;
            int numFigs = 28, numPointsLN = 1000, time = 56, space = 20, numStepsLN = 100, saveFigs = 0;
            structParameters parameters = ParametersInitialize();
            parameters.alphaTHelper = alpha;
            parameters.bTHelper = beta;
            parameters.epslon_x = epslon;
            parameters.gammaT = gammaT;
            structModel model = ModelInitialize(parameters, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);
            int size = 3;
            int save_times[3] = {70000, 140000 - 1, 280000 - 1};
            float* result_ = RunModel(&model, save_times, size);
            if(isnan(result_[0])) result_[0] = DBL_MAX;
            if(isnan(result_[1])) result_[1] = DBL_MAX;
            if(isnan(result_[2])) result_[2] = DBL_MAX;
            //std::cout << std::endl << "Result 1: " << result_[0] << " Result 2: " << result_[1] << " Result 3: " << result_[2] << std::endl;
            double result = (double)(pow(1.883 - result_[0], 2) + pow(3.35 - result_[1], 2) + pow(1.44 - result_[2], 2));
            result = sqrt(result);
            std::cout << "Gamma T " << gammaT << " Epslon: " << epslon << " Alpha: " << alpha << " Beta: " << beta << " Result: " << result << std::endl;
            return {result};
        }   

            
};