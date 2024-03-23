# AutoImmunityModel

C version of [NeuroImmunoModel](https://github.com/quintelabm/NeuroImmunoModel) - Mathematical Model of the Activation of the Immune Response in the Central Nervous System.

![Modelo Matematico](https://ars.els-cdn.com/content/image/1-s2.0-S0377042723001073-ga1_lrg.jpg)

## Branch main

Stable serial version of the code. 
Compile: gcc -O3 -o main main.c model.c

## omp

Stable parallel version of the code in OMP. 
Compile: gcc -O3 -fopenmp -o main mainOMP.c modelOMP.c -lm
Execute: ./main [ThreadsNum]

## mpi

Stable parallel version of the code in MPI. Compile with mpicc.[mpich or openmpi] modelMPI.c mainMPI.c -O3 -o main 
Execute: mpiexec.[mpich or openmpi] main -n [ProcNum] ./main

## cuda

Stable parallel version of the code in CUDA. CUDA version GPU Solving Tissue and CPU solving Lymph node
Compile: nvcc -O3 -use_fast_math -o main main.c model.c

## cudaHtDinamico

CUDA version GPU using dynamic time step.
Compile: nvcc -O3 -use_fast_math -o main main.c model.c

## Stream

CUDA version GPU Solving Tissue and CPU solving Lymph node. Assync copies.
Compile: nvcc -O3 -use_fast_math -o main main.c model.c

## blockGroup

CUDA version solving Tissue and Lymph Node in GPU using cooperativeGroup for reduction in tissue.
Compile: nvcc -O3 -use_fast_math -o main main.c model.c

## lymphNodeGPU

CUDA version solving Tissue and Lymph Node in GPU. CPU only for reduction in tissue.
Compile: nvcc -O3 -use_fast_math -o main main.c model.c

# Authors:

* Matheus Avila Moreira de Paula 
* Barbara de M Quintela
* Marcelo Lobosco

