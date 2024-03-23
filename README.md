# AutoImmunityModel

C version of [NeuroImmunoModel](https://github.com/quintelabm/NeuroImmunoModel) - Mathematical Model of the Activation of the Immune Response in the Central Nervous System.

## Publication

Graphical abstract of the paper that describes the model:

![Modelo Matematico](https://ars.els-cdn.com/content/image/1-s2.0-S0377042723001073-ga1_lrg.jpg)



DE PAULA, M. A. M.; QUINTELA, B. de M.; LOBOSCO, M. On the use of a coupled
mathematical model for understanding the dynamics of multiple sclerosis. Journal of
Computational and Applied Mathematics, v. 428, p. 115163, 2023. ISSN 0377-0427
DOI: [https://doi.org/10.1016/j.cam.2023.115163](https://doi.org/10.1016/j.cam.2023.115163)

## Branch main

Stable serial version of the code. 

Compile: 

ˋˋˋ
gcc -O3 -o main main.c model.c
ˋˋˋ

## omp

Stable parallel version of the code in OMP. 

Compile: 
ˋˋˋ
gcc -O3 -fopenmp -o main mainOMP.c modelOMP.c -lm
ˋˋˋ

Execute: 

ˋˋˋ
./main [ThreadsNum]
ˋˋˋ

## mpi

Stable parallel version of the code in MPI. 

Compile with mpicc.[mpich or openmpi] modelMPI.c mainMPI.c -O3 -o main 

Execute: 
ˋˋˋ
mpiexec.[mpich or openmpi] main -n [ProcNum] ./main
ˋˋˋ

## cuda

Stable parallel version of the code in CUDA. CUDA version GPU Solving Tissue and CPU solving Lymph node

Compile: 
ˋˋˋ
nvcc -O3 -use_fast_math -o main main.c model.c
ˋˋˋ

## cudaHtDinamico

CUDA version GPU using dynamic time step.

Compile: 
ˋˋˋ
nvcc -O3 -use_fast_math -o main main.c model.c
ˋˋˋ

## Stream

CUDA version GPU Solving Tissue and CPU solving Lymph node. Assync copies.

Compile: 

ˋˋˋ
nvcc -O3 -use_fast_math -o main main.c model.c
ˋˋˋ

## blockGroup

CUDA version solving Tissue and Lymph Node in GPU using cooperativeGroup for reduction in tissue.

Compile: 

ˋˋˋ
nvcc -O3 -use_fast_math -o main main.c model.c
ˋˋˋ

## lymphNodeGPU

CUDA version solving Tissue and Lymph Node in GPU. CPU only for reduction in tissue.

Compile: 

ˋˋˋ
nvcc -O3 -use_fast_math -o main main.c model.c
ˋˋˋ

# Authors:

* Matheus Avila Moreira de Paula 
* Barbara de M Quintela
* Marcelo Lobosco

