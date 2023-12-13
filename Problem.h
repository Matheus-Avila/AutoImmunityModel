#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include "model.h"
#include <fstream>
#include <string>
#include <float.h>
#include <vector>
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
    params.alphaTCytotoxic = 0.1;
    params.alphaB = 0.1;
    params.alphaP = 1;
    params.bTHelper = 0.17;
    params.bTCytotoxic = 0.001;
    params.bRho = 0.6;
    params.bRhoB = 3.02;
    params.bRhoP = 1.02;
    params.rhoTHelper = 2;
    params.rhoTCytotoxic = 2;
    params.rhoB = 11;
    params.rhoP = 3;
    params.rhoAntibody = 5.1*pow(10,-2);
    params.stableTHelper = 58.9;
    params.stableTCytotoxic = 28.4;
    params.stableB = 25;
    params.stableP = 2.5;
    params.V_LN = 40;
    params.V_BV = 0;
    params.V_PV = 0;

    params.epslon_x = 0.1;

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
            return {{0.0001, 0.01, 0.005}, {1.0, 0.99, 1.0}};
        }

        std::vector<double> fitness(const vector_double& variables) const {
            double alpha = variables[0];
            double epslon = variables[1];
            double gammaT = variables[2];
            float ht = 0.0002, hx = 0.5;
            int numFigs = 56, numPointsLN = 1000, time = 56, space = 20, numStepsLN = 100, saveFigs = 0;
            structParameters parameters = ParametersInitialize();
            parameters.alphaTCytotoxic = alpha;
            parameters.epslon_x = epslon;
            parameters.gammaT = gammaT;
            structModel model = ModelInitialize(parameters, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);
            int size = 3;
            int save_times[3] = {70000, 140000, 280000};
            
            std::cout << "Gamma T " << gammaT << " Epslon: " << epslon << " Alpha: " << alpha << std::endl;
            float points[size] = {6.7, 9.7, 9.5};
            float error = RunModel(&model, save_times, size, points);
            vector_double _error = (vector_double) error;
            std::vector<double> v;
            v.resize(1);
            std::cout << "Error: " << error << std::endl;
            v[0] = (double)error;
            return v;
        }   

            
};

