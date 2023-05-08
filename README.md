# AutoImmunityModel

C version of NeuroImmunoModel

# Branch main

Stable serial version of the code. 
Compile: gcc -O3 -o main main.c model.c

# omp

Stable parallel version of the code in OMP. 
Compile: gcc -O3 -fopenmp -o main main.c model.c
Execute: ./main [ThreadsNum]

# mpi

Stable parallel version of the code in MPI. Compile with mpicc 
Execute: mpiexec main -n [ProcNum]

# cuda

Stable parallel version of the code in CUDA. CUDA version GPU Solving Tissue and CPU solving Lymph node
Compile: nvcc -O3 -o main main.c model.c

# Stream

CUDA version GPU Solving Tissue and CPU solving Lymph node. Assync copies.

# blockGroup

CUDA version solving Tissue and Lymph Node in GPU using cooperativeGroup for reduction in tissue.

# lymphNodeGPU

CUDA version solving Tissue and Lymph Node in GPU. CPU only for reduction in tissue.
