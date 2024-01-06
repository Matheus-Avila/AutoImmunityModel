#!/bin/bash
for counter in $(seq 1 12)
do
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_ALL 11 10 2000
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_CPU 11 10 2000
  mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_GPU 11 10 2000
done