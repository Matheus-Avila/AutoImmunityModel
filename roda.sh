!/bin/bash

for counter in $(seq 1 10)
do
  /usr/bin/time -a -o timeESTALL-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_ALL 20 20 600 | tee -a timeESTALL-EST.txt
  /usr/bin/time -a -o timeESTCPU-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_CPU 20 20 600 | tee -a timeESTCPU-EST.txt
  /usr/bin/time -a -o timeESTGPU-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_GPU 20 20 600 | tee -a timeESTGPU-EST.txt

  /usr/bin/time -a -o timeUnBALL-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_ALL 20 20 600 | tee -a timeUnBALL-UnB.txt
  /usr/bin/time -a -o timeUnBCPU-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_CPU 20 20 600 | tee -a timeUnBCPU-UnB.txt
  /usr/bin/time -a -o timeUnBGPU-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_GPU 20 20 600 | tee -a timeUnBGPU-UnB.txt
  
  /usr/bin/time -a -o timeDIMALL-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-85 20 20 600 | tee -a timeDIMALL-85.txt
  /usr/bin/time -a -o timeDIMCPU-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-85 20 20 600 | tee -a timeDIMCPU-85.txt
  /usr/bin/time -a -o timeDIMGPU-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-85 20 20 600 | tee -a timeDIMGPU-85.txt 
  
  /usr/bin/time -a -o timeDIMALL-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-17 20 20 600 | tee -a timeDIMALL-17.txt
  /usr/bin/time -a -o timeDIMCPU-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-17 20 20 600 | tee -a timeDIMCPU-17.txt
  /usr/bin/time -a -o timeDIMGPU-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-17 20 20 600 | tee -a timeDIMGPU-17.txt
  
  /usr/bin/time -a -o timeDIMALL-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-6 20 20 600 | tee -a timeDIMALL-6.txt
  /usr/bin/time -a -o timeDIMCPU-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-6 20 20 600 | tee -a timeDIMCPU-6.txt
  /usr/bin/time -a -o timeDIMGPU-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-6 20 20 600 | tee -a timeDIMGPU-6.txt

done


for counter in $(seq 1 10)
do
  /usr/bin/time -a -o timeESTALL-300-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_ALL 20 20 300 | tee -a timeESTALL-300-EST.txt
  /usr/bin/time -a -o timeESTCPU-300-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_CPU 20 20 300 | tee -a timeESTCPU-300-EST.txt
  /usr/bin/time -a -o timeESTGPU-300-EST.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Estatico_GPU 20 20 300 | tee -a timeESTGPU-300-EST.txt

  /usr/bin/time -a -o timeUnBALL-300-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_ALL 20 20 300 | tee -a timeUnBALL-300-UnB.txt
  /usr/bin/time -a -o timeUnBCPU-300-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_CPU 20 20 300 | tee -a timeUnBCPU-300-UnB.txt
  /usr/bin/time -a -o timeUnBGPU-300-UnB.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_UnB_GPU 20 20 300 | tee -a timeUnBGPU-300-UnB.txt
  
  /usr/bin/time -a -o timeDIMALL-300-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-85 20 20 300 | tee -a timeDIMALL-300-85.txt
  /usr/bin/time -a -o timeDIMCPU-300-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-85 20 20 300 | tee -a timeDIMCPU-300-85.txt
  /usr/bin/time -a -o timeDIMGPU-300-85.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-85 20 20 300 | tee -a timeDIMGPU-300-85.txt 
  
  /usr/bin/time -a -o timeDIMALL-300-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-17 20 20 300 | tee -a timeDIMALL-300-17.txt
  /usr/bin/time -a -o timeDIMCPU-300-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-17 20 20 300 | tee -a timeDIMCPU-300-17.txt
  /usr/bin/time -a -o timeDIMGPU-300-17.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-17 20 20 300 | tee -a timeDIMGPU-300-17.txt
  
  /usr/bin/time -a -o timeDIMALL-300-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_ALL-6 20 20 300 | tee -a timeDIMALL-300-6.txt
  /usr/bin/time -a -o timeDIMCPU-300-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_CPU-6 20 20 300 | tee -a timeDIMCPU-300-6.txt
  /usr/bin/time -a -o timeDIMGPU-300-6.txt mpiexec -machinefile ./maquinas -n 1 ./HIS_Dinamico_GPU-6 20 20 300 | tee -a timeDIMGPU-300-6.txt

done