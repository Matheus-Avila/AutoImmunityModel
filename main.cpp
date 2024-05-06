#include "model.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Problem.h"
#include <chrono>
#include <iostream>
#include <pagmo/algorithm.hpp>

#include <pagmo/algorithms/de1220.hpp>

#include <pagmo/archipelago.hpp>


void WriteTime(float ExecTime){
    FILE *fileTime;
    fileTime = fopen("./ExecsTimes.txt", "a");
    if(fileTime != NULL){
    fprintf(fileTime, "%f\n", ExecTime);
    fclose(fileTime);
    }else{
        printf("Error execution time file\n");
        exit(0);
    }
}

void clearPhgTxt(){
    system("find ./result/ -name '*.png' -type f -delete");
    //system("find ./result/ -name '*.txt' -type f -delete");
    system("mkdir result");
    system("mkdir result/matrix");
    system("mkdir result/odc");
    system("mkdir result/mic");
    system("mkdir result/tke");
    system("mkdir result/ant");
    system("mkdir result/da");
    system("mkdir result/dc");
}

// int main(){
//     clearPhgTxt();
//     printf("Finishing computing\n");
//     return 0;
// }

int main(){
    
    clearPhgTxt();

    std::vector<double> epslon_values = {0, 0.55, 0.99};
    
    for (double epslon : epslon_values) {
        MSProblemTCytoParams pr; // modelo
        problem prob{pr}; // configuracoes da evolucao

        std::cout << prob << std::endl;

        //esse algoritmo eh da pagmo DE1220 = pDE
        // Self-adaptive Differential Evolution
        algorithm algo{de1220(50)}; // numero de geracoes? 10
        // eh para imprimir na tela???
        algo.set_verbosity(1);
        // cria arquipelago ... explicar ... 
        archipelago archi{8u, algo, prob, 10u};
        /* Evolves the population for a maximum number of generations, 
        until one of tolerances set on the population flatness 
        (x_tol, f_tol) are met.*/ 
        archi.evolve(20); // population 5

        // espera por que?
        archi.wait_check();

        // Print to screen the fitness of the
        // best solution in the new population.
        // std::cout << "Fitness of the best solution: "
        //           << new_pop.champion_f()[0] << '\n';


        //comentar a funcao! 
        for (const auto &isl : archi) {
            std::cout << isl.get_population().champion_f()[0] << '\n';
        }
    }
    return 0;
}