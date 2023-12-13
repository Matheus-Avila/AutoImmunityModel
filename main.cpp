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
    system("find ./result/ -name '*.txt' -type f -delete");
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
    MSProblemTCytoParams pr;
    problem prob{pr};

    std::cout << prob << std::endl;

    algorithm algo{de1220(10)};
    algo.set_verbosity(1);
    archipelago archi{8u, algo, prob, 10u};
    archi.evolve(5);


    archi.wait_check();

    for (const auto &isl : archi) {

        std::cout << isl.get_population().champion_f()[0] << '\n';

    }
    
    return 0;
}