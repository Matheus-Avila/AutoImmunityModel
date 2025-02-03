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
    params.gammaT = 0.9138;//0.1;

    params.avgT = 37;
    params.avgDc = 33;
    params.avgMic = 350;
    params.avgOdc = 400;

    params.ct = 0.0678;//0.0852465;//0.078673;//0.01;
    params.cMic = 0.1;
    params.cCDc = 1;
    params.cADc = 1;
    params.cDl = 0.1;
    params.cF = 0.1;
    params.alphaTHelper = 0.1;
    params.alphaTCytotoxic = 0.6836;//0.1;
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
    params.stableTHelper = 0.44;//58.9; //ponto inicial
    params.stableTCytotoxic =  17.5586;//28.4;
    params.stableB = 11.43;//25;
    params.stableP = 2.5;
    params.V_LN = 40;
    params.V_BV = 0;
    params.V_PV = 0;
    
    params.epslon_x =  0.8836; //0.99

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
        MSProblemTCytoParams() : pagmo::problem() {} //A classe MSProblemTCytoParams herda de pagmo::problem, 
                                                    //que é a classe base para definir problemas personalizados em PaGMO.
        std::string get_name() const { //retorna o nome do problema
            return "Multiple Sclerosis Problem with Tcyto alpha and beta";
        }

        vector_double::size_type get_nx() const { //Retorna o número de variáveis de decisão do problema. 
                                                //Aqui, retorna 0, o que indica que esse método pode estar incompleto ou que a função não precisa dessa informação.
            return 0;
        }

        vector_double::size_type get_nec() const {
            return 0;
        }

        vector_double::size_type get_nic() const {
            return 0;
        }
        //Retornam o número de restrições de igualdade e desigualdade, 
        //respectivamente. Ambos retornam 0, indicando que não há restrições.

        vector_double::size_type get_ncx() const { //Retorna o número de variáveis de controle, que é 2 neste caso. 
                                                //Isso indica que há duas variáveis que o algoritmo de otimização deve ajustar.
            return 2;
        }

        //Define os limites das variáveis de decisão. Aqui, os limites são {0.0001, 0.55, 0.005, 0.001} como limites inferiores 
        //e {1.0, 0.99, 1.0, 0.1} como limites superiores. 
        //Esses limites são aplicáveis às variáveis epslon, alpha, gammaT e ct, respectivamente.
        std::pair<vector_double, vector_double> get_bounds() const
        {
            return {{0.0001, 0.55, 0.005,0.001}, {1.0, 0.99, 1.0, 0.1}};//0.001
            //return {{0.1}, {27.6}};
        }

        std::vector<double> fitness(const vector_double& variables) const {
            //std::vector<double> epslon ;
            double epslon = variables[0];
            double alpha = variables[1];
            double gammaT = variables[2];
            double ct = variables[3];
            //double tcyto = variables[0];
            
            float ht = 0.0002, hx = 0.5;
            int numFigs = 7, numPointsLN = 1000, time = 360, space = 20, numStepsLN = 100, saveFigs = 1 ;
            structParameters parameters = ParametersInitialize();
            // parameters.alphaTCytotoxic = alpha;
            // parameters.epslon_x = epslon;
            // parameters.gammaT = gammaT;
            // parameters.ct = ct;
            //parameters.stableTCytotoxic = tcyto;
            structModel model = ModelInitialize(parameters, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);
            // int size = 5;
            // int save_times[size] = {0,14,30,60,90}; //14 e 28

            int targetSize = 7;
            int targetDays[targetSize] = {0,14,30,60,90,180,360};
            
            //std::cout << " Epslon: " << parameters.eps_new << std::endl;
            std::cout << std::endl;
            std::cout << "Gamma T: " << parameters.gammaT << " Epslon: " << parameters.epslon_x << " Alpha: " << parameters.alphaTCytotoxic << " Ct: " << parameters.ct << std::endl;
            //std::cout << "stableTCytotoxic: " << parameters.stableTCytotoxic << std::endl;
            //float points[targetSize] = {27.6, 5.4, 8.6, 7.8, 8.2, 5.1}; 
            //float points[targetSize] = {22.81, 0.68, 3.83, 3.11, 1.34, 2.18, 1.44}; //paciente 7
            //float points[targetSize] = {11.89, 0.66, 4.27, 0.87, 1.02, 0.89, 1.04}; //paciente 2
            //float points[targetSize] = {43.19, 1.59, 2.01, 1.78, 5.19, 3.94, 3.27}; //paciente 16
            //float points[targetSize] = {35.07, 6.79, 24.89, 22.19, 15.29, 8.86, 10.55}; //paciente 5
            //float points[targetSize] = {19.76, 2.28, 2.26, 1.32, 2.24, 1.32, 1.11}; //paciente 20
            //float points[targetSize] = {20.59, 14.32, 8.88, 22.82, 15.65, 10.80, 8.51}; //paciente 21
            //float points[targetSize] = {40.10, 24.50, 19.19, 19.72, 0.17, 6.16, 3.66}; //paciente 23
            //float points[targetSize] = {24.38, 2.06, 4.23, 1.28, 3.47, 2.96, 2.65}; //paciente 18
            //float points[targetSize] = {11.57, 1.48, 2.01, 3.98, 3.15, 7.41, 2.56}; //paciente 14
            //float points[targetSize] = {43.70, 4.43, 4.33, 3.40, 1.68, 1.61, 1.44}; //paciente 13
            float points[targetSize] = {21.39, 7.17, 13.45, 10.42, 5.03, 5.74, 5.59}; //paciente 11

            float error = RunModel(&model, points, targetDays, targetSize);
            //std::cout << "Error: (voltou runmodel) " << error << std::endl;
            //vector_double _error = (vector_double) error;
            std::vector<double> v;
            v.resize(1);
            std::cout << "Error: " << error << std::endl;
            v[0] = (double)error;
            //exit(0);
            return v;
        }   

            
};



// structParameters ParametersInitialize() {
//     structParameters params;
//     params.micDiffusion = 0.015206;
//     params.antibodyDiffusion = 0.15206;
//     params.cDcDiffusion = 0.015206;
//     params.aDcDiffusion = 0.015206;
//     params.tCytoDiffusion = 0.015206;
//     params.chi = 0.03;

//     params.muCDc = 60 * 24 * 3 * std::pow(10, -5);
//     params.muMic = 60 * 24 * 3 * std::pow(10, -6);
//     params.rM = 60 * 24 * 6 * std::pow(10, -7);
//     params.rT = 0.001;
//     params.lambAntMic = 5.702 * std::pow(10, -3);
//     params.bD = 0.001;

//     params.gammaD = 0.1;
//     params.gammaAntibody = 0.3;
//     params.gammaT = 0.8430;//0.947922; //0.9589417
    
//     params.avgT = 37;
//     params.avgDc = 33;
//     params.avgMic = 350;
//     params.avgOdc = 400;

//     params.ct = 0.1;//0.0678;//0.128673;//0.0999403;//0.078673; //0.0786403
//     params.cMic = 0.1;
//     params.cCDc = 1;
//     params.cADc = 1;
//     params.cDl = 0.1;
//     params.cF = 0.1;
//     params.alphaTHelper = 0.1;
//     params.alphaTCytotoxic =  0.6062;//0.565958; //0.6168163
//     params.alphaB = 0.1;
//     params.alphaP = 1;
//     params.bTHelper = 0.17;
//     params.bTCytotoxic = 0.001; 
//     params.bRho = 0.6;
//     params.bRhoB = 3.02;
//     params.bRhoP = 1.02;
//     params.rhoTHelper = 2;
//     params.rhoTCytotoxic = 2;
//     params.rhoB = 11;
//     params.rhoP = 3;
//     params.rhoAntibody = 5.1 * std::pow(10, -2);
//     params.stableTHelper = 0.44;
//     params.stableTCytotoxic = 27.6;
//     params.stableB = 11.43;
//     params.stableP = 2.5;
//     params.V_LN = 40;
//     params.V_BV = 0;
//     params.V_PV = 0;
    
//     params.epslon_x = 0.99;//0.8127;//0.9354767;//0.55;//0.9354767; //0.971843;//0.55;

//     return params;
// }

// float sumVector(const std::vector<float>& vector, int sizeX, int sizeY) {
//     float sum = 0.0;
//     for (int i = 0; i < sizeX; ++i) {
//         for (int j = 0; j < sizeY; ++j) {
//             sum += vector[i * sizeX + j];
//         }
//     }
//     return sum;
// }

// int appendOnFile(const std::string& filename, double value) {
//     std::ofstream file(filename, std::ios::app);
//     if (!file.is_open()) {
//         std::cerr << "Erro ao abrir o arquivo." << std::endl;
//         return 1;
//     }

//     file << value << std::endl;
//     return 0;
// }

// float sumOneDimensionVector(const std::vector<float>& vector) {
//     float sum = 0.0;
//     for (float v : vector) {
//         sum += v;
//     }
//     return sum;
// }

// using namespace pagmo;

// class MSProblemTCytoParams : public pagmo::problem {
// public:
//     MSProblemTCytoParams() : pagmo::problem() {}

//     std::string get_name() const {
//         return "Multiple Sclerosis Problem with Tcyto alpha and beta";
//     }

//     vector_double::size_type get_nx() const  {
//         return 0;
//     }

//     vector_double::size_type get_nec() const  {
//         return 0;
//     }

//     vector_double::size_type get_nic() const  {
//         return 0;
//     }

//     vector_double::size_type get_ncx() const  {
//         return 2;
//     }

//     std::pair<vector_double, vector_double> get_bounds() const  {
//         return {{0.0001, 0.55, 0.005, 0.001}, {1.0, 0.99, 1.0, 0.1}};
//     }

//     std::vector<double> fitness(const vector_double& variables) const  {
//         double epslon = variables[0];
//         double alpha = variables[1];
//         double gammaT = variables[2];
//         double ct = variables[3];
        
//         float ht = 0.0002, hx = 0.5;
//         int numFigs = 7, numPointsLN = 1000, time = 360, space = 20, numStepsLN = 100, saveFigs = 0;
//         structParameters parameters = ParametersInitialize();
//         parameters.alphaTCytotoxic = alpha;
//         parameters.epslon_x = epslon;
//         parameters.gammaT = gammaT;
//         parameters.ct = ct;

//         structModel model = ModelInitialize(parameters, ht, hx, time, space, numFigs, numPointsLN, numStepsLN, saveFigs);

//         std::vector<int> targetDays = {0, 14, 30, 60, 90, 180, 360};
//         //std::vector<float> points = {27.6, 5.4, 8.6, 7.8, 8.2, 5.1, 4.7};
//         //std::vector<float> points = {40.1, 24.5, 19.2, 19.7, 0.2, 6.2, 3.7}; //paciente 23
//         std::vector<float> points = {24.38, 2.06, 4.23, 1.28, 3.47, 2.96, 2.65}; //paciente 18
//         std::cout << "\nGamma T " << parameters.gammaT << " Epslon: " << parameters.epslon_x << " Alpha: " << parameters.alphaTCytotoxic << " Ct: " << parameters.ct << std::endl;
//         float error = RunModel(&model, points.data(), targetDays.data(), targetDays.size());

//         std::vector<double> v = {static_cast<double>(error)};
//         v.resize(1);
//         std::cout << "Error: " << error << std::endl;
//         v[0] = (double)error;
//         return v;
//     }
// };

