#ifndef OPENCLWRAPPER_H
#define OPENCLWRAPPER_H

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#endif

int InitParallelProcessor();	//Returns number of devices.
void FinishParallelProcessor();

int CreateKernel(int devicePosition, const char *source, const char *kernelName);				//Return ID of the kernel.
void SetKernelAttribute(int devicePosition, int kernelID, int attribute, int memoryObjectID);
bool RemoveKernel(int devicePosition, int kernelID);

int CreateMemoryObject(int devicePosition, int size, cl_mem_flags memoryType, void *hostMemory);		//Return ID of the cl memory object.
int WriteToMemoryObject(int devicePosition, int memoryObjectID, const char *data, int offset, int size);	//Return position of the event generated.
int ReadFromMemoryObject(int devicePosition, int memoryObjectID, char *data, int offset, int size);		//Return position of the event generated.
bool RemoveMemoryObject(int devicePosition, int memoryObjectID);

int RunKernel(int devicePosition, int kernelID, int parallelDataOffset, int parallelData, int workGroupSize);	//Return position of the event generated.
void SynchronizeCommandQueue(int devicePosition);

void SynchronizeEvent(int devicePosition, int eventPosition);
long int GetEventTaskOverheadTicks(int devicePosition, int eventPosition);
long int GetEventTaskTicks(int devicePosition, int eventPosition);

cl_device_type GetDeviceType(int devicePosition);
size_t GetDeviceMaxWorkItemsPerWorkGroup(int devicePosition);
cl_uint GetDeviceComputeUnits(int devicePosition);

bool isDeviceCPU(int devicePosition);

#endif
