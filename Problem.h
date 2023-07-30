#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include "model.h"
#include <fstream>
#include <string>
#include <float.h>

structParameters ParametersInitialize(){
    structParameters params;
    params.micDiffusion = 9.6*24*6.6*pow(10,-5);
    params.antibodyDiffusion = 9.6*24*6.6*pow(10,-4);
    params.cDcDiffusion = 9.6*24*6.6*pow(10,-5);
    params.aDcDiffusion = 9.6*24*6.6*pow(10,-5);
    params.tCytoDiffusion = 9.6*24*6.6*pow(10,-5);
    params.chi = 0.298*60*2;
    params.epslon_x = 0.212607;
    params.muCDc = 60*24*3*pow(10,-4);
    params.muMic = 60*24*3*pow(10,-6);
    params.rM = 60*24*3.96*pow(10,-6);
    params.rT = 0.1;
    params.lambAntMic = 5.702*pow(10,-3);
    params.bD = 0.001;
    
    params.gammaD = 0.01;
    params.gammaAntibody = 0.3;
    params.gammaT = 2;

    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;

    params.cMic = 0.1;
    params.cCDc = 0.1;
    params.cADc = 0.1;
    params.cDl = 0.1;
    params.cF = 0.1;
    params.alphaTHelper = 0.1;
    params.alphaTCytotoxic = 0.00274764;
    params.alphaB = 0.1;
    params.alphaP = 1;
    params.bTHelper = 0.17;
    params.bTCytotoxic = 0.179413;
    params.bRho = 0.6;
    params.bRhoB = 3.02;
    params.bRhoP = 1.02;
    params.rhoTHelper = 2;
    params.rhoTCytotoxic = 2;
    params.rhoB = 11;
    params.rhoP = 3;
    params.rhoAntibody = 5.1*pow(10,-2);
    params.estableTHelper = 84;
    params.estableTCytotoxic = 40;
    params.estableB = 25;
    params.estableP = 2.5;
    params.V_LN = 160;
    params.V_BV = 0;
    params.V_PV = 0;

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
            return {{0.001, 0.001}, {0.5, 0.5}};
        }

        vector_double fitness(const vector_double& variables) const {
            double alpha = variables[0];
            double beta = variables[1];
            float ht = 0.0002, hx = 0.5;
            int numFigs = 28, numPointsLN = 1000, time = 28, space = 20, numStepsLN = 100, saveFigs = 0;
            structParameters parameters = ParametersInitialize();
            parameters.alphaTHelper = alpha;
            parameters.bTHelper = beta;
            structModel model = ModelInitialize(parameters, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);
            int size = 2;
            int save_times[2] = {70000, 140000 - 1};
            float* result_ = RunModel(&model, save_times, size);
            if(isnan(result_[0])) result_[0] = DBL_MAX;
            if(isnan(result_[1])) result_[1] = DBL_MAX;
            double result = (double)(fabs(6.81 - result_[0]) + fabs(9.93 - result_[1]));
            std::cout << alpha << ',' << beta << ',' << result_[0] << ',' << result_[1] << std::endl;
            //appendOnFile("iterations_tcyto.txt", result);
            return {result};
        }   

            
};