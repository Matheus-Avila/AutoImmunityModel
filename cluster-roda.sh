#!/bin/bash

for counter in $(seq 1 8)
do
#  echo "estatico all"
#  /usr/bin/time -a -o timeALL-ESTATICO.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_ALL 20 20 600 | tee -a timeALL-ESTATICO.txt
  echo "estatico CPU"
  /usr/bin/time -a -o timeCPU-ESTATICO.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_CPU 20 20 600 | tee -a timeCPU-ESTATICO.txt
#  echo "estatico GPU"
#  /usr/bin/time -a -o timeGPU-ESTATICO.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_GPU 20 20 600 | tee -a timeGPU-ESTATICO.txt
  
  echo "dinamico all"
  /usr/bin/time -a -o timeALL-DINAMICO.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL 20 20 600 | tee -a timeALL-DINAMICO.txt
  echo "dinamico GPU"
  /usr/bin/time -a -o timeGPU-DINAMICO.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU 20 20 600 | tee -a timeGPU-DINAMICO.txt

  echo "dinamico all"
  /usr/bin/time -a -o timeALL-DINAMICO3.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL 20 20 300 | tee -a timeALL-DINAMICO3.txt
  echo "dinamico GPU"
  /usr/bin/time -a -o timeGPU-DINAMICO3.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU 20 20 300 | tee -a timeGPU-DINAMICO3.txt

done
