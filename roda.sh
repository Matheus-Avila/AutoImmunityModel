#!/bin/bash
for counter in $(seq 1 12)
do
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_ALL 22 100 100
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_CPU 22 100 100
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_GPU 22 100 100
done
