export LD_LIBRARY_PATH=/share/apps/amd-opencl/AMDAPPSDK-3.0/lib/x86_64/sdk:$LD_LIBRARY_PATH
/opt/mpich2/gnu/bin/mpic++ -L/share/apps/amd-opencl/AMDAPPSDK-3.0/lib/x86_64/sdk -D MAX_NUMBER_OF_DEVICES_PER_PLATFORM=3 -D CPU_DEVICES  -O3 -o HIS_Estatico_CPU Balanceador.cpp OpenCLWrapper.cpp main.cpp -lOpenCL -lrt
/opt/mpich2/gnu/bin/mpic++ -L/share/apps/amd-opencl/AMDAPPSDK-3.0/lib/x86_64/sdk -D MAX_NUMBER_OF_DEVICES_PER_PLATFORM=3 -D GPU_DEVICES  -O3 -o HIS_Estatico_GPU Balanceador.cpp OpenCLWrapper.cpp main.cpp -lOpenCL -lrt
/opt/mpich2/gnu/bin/mpic++ -L/share/apps/amd-opencl/AMDAPPSDK-3.0/lib/x86_64/sdk -D MAX_NUMBER_OF_DEVICES_PER_PLATFORM=3 -D ALL_DEVICES  -O3 -o HIS_Estatico_ALL Balanceador.cpp OpenCLWrapper.cpp main.cpp -lOpenCL -lrt
