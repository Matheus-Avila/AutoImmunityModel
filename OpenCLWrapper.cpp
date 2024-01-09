#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "OpenCLWrapper.h"

#define MAX_SOURCE_BUFFER_LENGTH	1000000
#define MAX_NUMBER_OF_PLATFORMS		10
#define MAX_NUMBER_OF_DEVICES		10

#define MAX_MEMORY_OBJECTS	100
#define MAX_KERNELS		100
#define MAX_EVENTS		1000

cl_platform_id platformIDs[MAX_NUMBER_OF_PLATFORMS];

struct Device
{
	cl_device_id deviceID;
	cl_device_type deviceType;
	cl_context context;
	cl_command_queue kernelCommandQueue;
	cl_command_queue dataCommandQueue;
	cl_program program;

	cl_mem *memoryObjects;
	cl_kernel *kernels;

	//The position of the ID vector points to the position of the correspondent OpenCL ID position.
	int *memoryObjectID;
	int *kernelID;

	//Events don't need ID, they are piled up on top of each other and discarded when synchronization occurs.
	cl_event *events;

	size_t deviceMaxWorkItemsPerWorkGroup;	//Maximum number of work items (run on each processing element) per work group (run on each compute unit). (Compute units have many processing elements)
	cl_uint deviceComputeUnits;		//Compute units in the device.

	int numberOfMemoryObjects;
	int numberOfKernels;
	int numberOfEvents;
};

cl_uint numberOfDevices;
Device *devices;

int automaticNumber = 0;

int Maximum(int a, int b)
{
	return (a > b) ? a : b;
}

int GetMemoryObjectPosition(int devicePosition, int memoryObjectID)
{
	for(int count = 0; count < devices[devicePosition].numberOfMemoryObjects; count++)
	{
		if(devices[devicePosition].memoryObjectID[count] == memoryObjectID)
		{
			return count;
		}
	}
	return -1;
}

int GetKernelPosition(int devicePosition, int kernelID)
{
	for(int count = 0; count < devices[devicePosition].numberOfKernels; count++)
	{
		if(devices[devicePosition].kernelID[count] == kernelID)
		{
			return count;
		}
	}
	return -1;
}

int InitParallelProcessor()
{
	cl_int state;

	//Get platforms.
	cl_uint numberOfPlatforms = 0;
	state = clGetPlatformIDs(MAX_NUMBER_OF_PLATFORMS, platformIDs, &numberOfPlatforms);
	if(state != CL_SUCCESS)
	{
		printf("OpenCL Error: Platform couldn't be found.\n");
	}
	printf("%i platform(s) found.\n", numberOfPlatforms);
	
	cl_device_id deviceList[MAX_NUMBER_OF_DEVICES];
	devices = new Device[MAX_NUMBER_OF_DEVICES];
	numberOfDevices = 0;
	for(int platform = 0; platform < numberOfPlatforms; platform++)
	{
		//Get devices.
		cl_uint numberOfDevicesOfPlatform;
		#if defined(ALL_DEVICES)
		state = clGetDeviceIDs(platformIDs[platform], CL_DEVICE_TYPE_CPU, MAX_NUMBER_OF_DEVICES_PER_PLATFORM, deviceList+numberOfDevices, &numberOfDevicesOfPlatform);
		if(state != CL_SUCCESS)
			state = clGetDeviceIDs(platformIDs[platform], CL_DEVICE_TYPE_GPU, MAX_NUMBER_OF_DEVICES_PER_PLATFORM, deviceList+numberOfDevices, &numberOfDevicesOfPlatform);
		#elif defined(GPU_DEVICES)
		state = clGetDeviceIDs(platformIDs[platform], CL_DEVICE_TYPE_GPU, MAX_NUMBER_OF_DEVICES_PER_PLATFORM, deviceList+numberOfDevices, &numberOfDevicesOfPlatform);
		#elif defined(CPU_DEVICES)
		state = clGetDeviceIDs(platformIDs[platform], CL_DEVICE_TYPE_CPU, MAX_NUMBER_OF_DEVICES_PER_PLATFORM, deviceList+numberOfDevices, &numberOfDevicesOfPlatform);
		#endif
		if(state != CL_SUCCESS)
		{
			printf("OpenCL Error: Devices couldn't be resolved.\n");
		}
		else
		{
			if(numberOfDevicesOfPlatform > MAX_NUMBER_OF_DEVICES_PER_PLATFORM)
			{
				numberOfDevicesOfPlatform = MAX_NUMBER_OF_DEVICES_PER_PLATFORM;
			}
			printf("%i device(s) found on platform %i.\n", numberOfDevicesOfPlatform, platform);
		
		//Set devices.
		for(int count = numberOfDevices; count < numberOfDevices+numberOfDevicesOfPlatform; count++)
		{
			//Get ID.
			devices[count].deviceID = deviceList[count];

			//Get type.
			clGetDeviceInfo(devices[count].deviceID, CL_DEVICE_TYPE, sizeof(cl_device_type), &devices[count].deviceType, NULL);
			clGetDeviceInfo(devices[count].deviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_device_type), &devices[count].deviceMaxWorkItemsPerWorkGroup, NULL);
			clGetDeviceInfo(devices[count].deviceID, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_device_type), &devices[count].deviceComputeUnits, NULL);

			if(devices[count].deviceType == CL_DEVICE_TYPE_GPU)
			{
				printf("Device (%i) type: GPU\n", count);
			}
			else if(devices[count].deviceType == CL_DEVICE_TYPE_CPU)
			{
				printf("Device (%i) type: CPU\n", count);
			}
			printf("CL_DEVICE_MAX_COMPUTE_UNITS %d\n", devices[count].deviceComputeUnits);

			
			//Create context.
			devices[count].context = clCreateContext(NULL, 1, &devices[count].deviceID, NULL, NULL, &state);
			if(state != CL_SUCCESS)
			{
				printf("OpenCL Error: Context couldn't be created.\n");
			}

			//Create command queue.
			devices[count].kernelCommandQueue = clCreateCommandQueue(devices[count].context, devices[count].deviceID, CL_QUEUE_PROFILING_ENABLE, &state);
			if(state != CL_SUCCESS)
			{
				printf("OpenCL Error: Kernel message queue couldn't be created.\n");
			}

			//Create command queue.
			devices[count].dataCommandQueue = clCreateCommandQueue(devices[count].context, devices[count].deviceID, CL_QUEUE_PROFILING_ENABLE, &state);
			if(state != CL_SUCCESS)
			{
				printf("OpenCL Error: Data message queue couldn't be created.\n");
			}

			//Initialize memory objects, kernel and events.
			devices[count].numberOfMemoryObjects = 0;
			devices[count].numberOfKernels = 0;
			devices[count].numberOfEvents = 0;

			devices[count].memoryObjects = new cl_mem[MAX_MEMORY_OBJECTS];
			devices[count].kernels = new cl_kernel[MAX_KERNELS];
			memset(devices[count].memoryObjects, 0, sizeof(cl_mem)*MAX_MEMORY_OBJECTS);
			memset(devices[count].kernels, 0, sizeof(cl_kernel)*MAX_KERNELS);

			devices[count].memoryObjectID = new int[MAX_MEMORY_OBJECTS];
			devices[count].kernelID = new int[MAX_KERNELS];
			memset(devices[count].memoryObjectID, 0, sizeof(int)*MAX_MEMORY_OBJECTS);
			memset(devices[count].kernelID, 0, sizeof(int)*MAX_KERNELS);

			devices[count].events = new cl_event[MAX_EVENTS];
			memset(devices[count].events, 0, sizeof(cl_kernel)*MAX_EVENTS);

			devices[count].program = 0;
		}
		numberOfDevices += numberOfDevicesOfPlatform;
	}
	}	
	return numberOfDevices;
}

void FinishParallelProcessor()
{
	for(int count = 0; count < numberOfDevices; count++)
	{
		clFlush(devices[count].kernelCommandQueue);
		clFinish(devices[count].kernelCommandQueue);

		clFlush(devices[count].dataCommandQueue);
		clFinish(devices[count].dataCommandQueue);

		for(int count2 = 0; count2 < devices[count].numberOfKernels; count2++)
		{
			clReleaseKernel(devices[count].kernels[count2]);
		}
		clReleaseProgram(devices[count].program);
		for(int count2 = 0; count2 < devices[count].numberOfMemoryObjects; count2++)
		{
			clReleaseMemObject(devices[count].memoryObjects[count2]);
		}
		clReleaseCommandQueue(devices[count].kernelCommandQueue);
		clReleaseCommandQueue(devices[count].dataCommandQueue);

		delete [] devices[count].memoryObjects;
		delete [] devices[count].kernels;
		devices[count].memoryObjects = NULL;
		devices[count].kernels = NULL;

		delete [] devices[count].memoryObjectID;
		delete [] devices[count].kernelID;
		devices[count].memoryObjectID = NULL;
		devices[count].kernelID = NULL;

		delete [] devices[count].events;
		devices[count].events = NULL;

		devices[count].numberOfMemoryObjects = 0;
		devices[count].numberOfKernels = 0;
		devices[count].numberOfEvents = 0;
	}
	delete [] devices;
	devices = NULL;
}

int CreateKernel(int devicePosition, const char *source, const char *kernelName)
{
	if(devices[devicePosition].program != 0)
	{
		clReleaseProgram(devices[devicePosition].program);
	}
	devices[devicePosition].program = 0;

	cl_int state;

	//Read kernel file.
	FILE *fileHandle;
	char *sourceBuffer = (char *)malloc(sizeof(char)*MAX_SOURCE_BUFFER_LENGTH);
	if((fileHandle = fopen(source, "r")) == NULL)
	{
		printf("Error reading %s\n!", source);
		return -1;
	}
	size_t sourceBufferLength = fread(sourceBuffer, 1, sizeof(char)*MAX_SOURCE_BUFFER_LENGTH, fileHandle);

	//Create program.
	devices[devicePosition].program = clCreateProgramWithSource(devices[devicePosition].context, 1, (const char **)&sourceBuffer, (const size_t *)&sourceBufferLength, &state);

	//Close kernel file.
	fclose(fileHandle);
	fileHandle = NULL;
	free(sourceBuffer);
	sourceBuffer = NULL;

	//Program created?
	if(state != CL_SUCCESS)
	{
		printf("Error creating program!\n");
		return -1;
	}
	cl_build_status status;
        size_t size;
        char* log;
	//Compile program.
	state = clBuildProgram(devices[devicePosition].program, 1, &devices[devicePosition].deviceID, NULL, NULL, NULL);
	clGetProgramBuildInfo(devices[devicePosition].program, devices[devicePosition].deviceID, CL_PROGRAM_BUILD_LOG, 0, NULL, &size );
	if(status != CL_BUILD_SUCCESS || state != CL_SUCCESS)
	{       
	        clGetProgramBuildInfo(devices[devicePosition].program, devices[devicePosition].deviceID, CL_PROGRAM_BUILD_LOG, 0, NULL, &size);
	        log = (char*)malloc(size);
	        clGetProgramBuildInfo(devices[devicePosition].program, devices[devicePosition].deviceID, CL_PROGRAM_BUILD_LOG, size, log, NULL);
	        
		printf("Error compiling program!\n%s\n", log);
		return -1;
	}

	//Create kernel.
	devices[devicePosition].kernels[devices[devicePosition].numberOfKernels] = clCreateKernel(devices[devicePosition].program, kernelName, &state);
	if(state != CL_SUCCESS)
	{
		printf("Error creating kernel! \n");
		return -1;
	}
	devices[devicePosition].kernelID[devices[devicePosition].numberOfKernels] = automaticNumber;
	devices[devicePosition].numberOfKernels += 1;
	automaticNumber += 1;
	return automaticNumber-1;
}

void SetKernelAttribute(int devicePosition, int kernelID, int attribute, int memoryObjectID)
{
	int kernelPosition = GetKernelPosition(devicePosition, kernelID);
	int memoryObjectPosition = GetMemoryObjectPosition(devicePosition, memoryObjectID);
	if(kernelPosition != -1 && memoryObjectPosition != -1)
	{
		cl_int state = clSetKernelArg(devices[devicePosition].kernels[kernelPosition], attribute, sizeof(cl_mem), (void *)&devices[devicePosition].memoryObjects[memoryObjectPosition]);
		if(state != CL_SUCCESS)
		{
			printf("Error setting kernel argument!\n");
		}
	}
	else
	{
		printf("Error setting kernel argument: Either kernel ID=%i or MemOBJ=%i don't exist!\n", kernelID, memoryObjectID);
	}
}

bool RemoveKernel(int devicePosition, int kernelID)
{
	int kernelPosition = GetKernelPosition(devicePosition, kernelID);
	if(kernelPosition != -1)
	{
		clReleaseKernel(devices[devicePosition].kernels[kernelPosition]);
		memcpy(devices[devicePosition].kernels+kernelPosition, devices[devicePosition].kernels+kernelPosition+1, sizeof(cl_kernel)*(devices[devicePosition].numberOfKernels-1));
		memcpy(devices[devicePosition].kernelID+kernelPosition, devices[devicePosition].kernelID+kernelPosition+1, sizeof(cl_kernel)*(devices[devicePosition].numberOfKernels-1));
		devices[devicePosition].numberOfKernels -= 1;
		return true;
	}
	return false;
}

int CreateMemoryObject(int devicePosition, int size, cl_mem_flags memoryType, void *hostMemory)
{
	cl_int state;
	if(devices[devicePosition].numberOfMemoryObjects < MAX_MEMORY_OBJECTS)
	{
		devices[devicePosition].memoryObjects[devices[devicePosition].numberOfMemoryObjects] = clCreateBuffer(devices[devicePosition].context, memoryType, size, hostMemory, &state);
		if(state != CL_SUCCESS)
		{
			printf("Error creating memory object!\n");
			return -1;
		}
		else
		{
			devices[devicePosition].memoryObjectID[devices[devicePosition].numberOfMemoryObjects] = automaticNumber;
			devices[devicePosition].numberOfMemoryObjects += 1;
		}
		automaticNumber += 1;
		return automaticNumber-1;
	}
	printf("Error creating memory object, limit exceeded!");
	return -1;
}

int WriteToMemoryObject(int devicePosition, int memoryObjectID, const char *data, int offset, int size)
{
	cl_int state;
	int memoryObjectPosition = GetMemoryObjectPosition(devicePosition, memoryObjectID);
	if(memoryObjectPosition != -1 && devices[devicePosition].numberOfEvents < MAX_EVENTS)
	{
		state = clEnqueueWriteBuffer(devices[devicePosition].dataCommandQueue, devices[devicePosition].memoryObjects[memoryObjectPosition], CL_FALSE, offset, size, data, 0, NULL, &devices[devicePosition].events[devices[devicePosition].numberOfEvents]);
		if(state != CL_SUCCESS)
		{
			printf("Error writing to memory object %i.\n", state);
			exit(1);
		}
		else
		{
			clFlush(devices[devicePosition].dataCommandQueue);
			devices[devicePosition].numberOfEvents += 1;
			return devices[devicePosition].numberOfEvents-1;
		}
	}


	printf("Error! Couldn't find memory object position %i or number of events %i exceeded limit.\n", memoryObjectPosition, devices[devicePosition].numberOfEvents);
	return -1;
}

int ReadFromMemoryObject(int devicePosition, int memoryObjectID, char *data, int offset, int size)
{
	cl_int state;
	int memoryObjectPosition = GetMemoryObjectPosition(devicePosition, memoryObjectID);
	if(memoryObjectPosition != -1 && devices[devicePosition].numberOfEvents < MAX_EVENTS)
	{
		state = clEnqueueReadBuffer(devices[devicePosition].dataCommandQueue, devices[devicePosition].memoryObjects[memoryObjectPosition], CL_FALSE, offset, size, data, 0, NULL, &devices[devicePosition].events[devices[devicePosition].numberOfEvents]);
		if(state != CL_SUCCESS)
		{
			printf("Error reading from memory object %i.\n", state);
			return -1;
		}
		else
		{
			clFlush(devices[devicePosition].dataCommandQueue);
			devices[devicePosition].numberOfEvents += 1;
			return devices[devicePosition].numberOfEvents-1;
		}
	}

	printf("Error! Couldn't find memory object position %i or number of events %i exceeded limit.\n", memoryObjectPosition, devices[devicePosition].numberOfEvents);
	return -1;	
}

bool RemoveMemoryObject(int devicePosition, int memoryObjectID)
{
	int memoryObjectPosition = GetMemoryObjectPosition(devicePosition, memoryObjectID);
	if(memoryObjectPosition != -1)
	{
		clReleaseMemObject(devices[devicePosition].memoryObjects[memoryObjectPosition]);
		memcpy(devices[devicePosition].memoryObjects+memoryObjectPosition, devices[devicePosition].memoryObjects+memoryObjectPosition+1, sizeof(cl_mem)*(devices[devicePosition].numberOfMemoryObjects-1));
		memcpy(devices[devicePosition].memoryObjectID+memoryObjectPosition, devices[devicePosition].memoryObjectID+memoryObjectPosition+1, sizeof(cl_mem)*(devices[devicePosition].numberOfMemoryObjects-1));
		devices[devicePosition].numberOfMemoryObjects -= 1;
		return true;
	}
	return false;
}

int RunKernel(int devicePosition, int kernelID, int parallelDataOffset, int parallelData, int workGroupSize)
{
	//----------------------------------------------------------------------------
	//Obs: OpenGL data used in this kernel must be synchronized before proceeding.
	//----------------------------------------------------------------------------

	int kernelPosition = GetKernelPosition(devicePosition, kernelID);
	if(kernelPosition != -1 && devices[devicePosition].numberOfEvents < MAX_EVENTS)
	{
		//Make sure parallelData is a power of 2.
		size_t globalItemsOffset = Maximum(parallelDataOffset, 0);
		size_t globalItems = parallelData;
		size_t mask = 0;
		globalItems = Maximum(workGroupSize, parallelData + workGroupSize - (parallelData%workGroupSize));
		
		cl_int state;
		size_t localItems = workGroupSize;

		state = clEnqueueNDRangeKernel(devices[devicePosition].kernelCommandQueue, devices[devicePosition].kernels[kernelPosition], 1, &globalItemsOffset, &globalItems, &localItems, 0, NULL, &devices[devicePosition].events[devices[devicePosition].numberOfEvents]);
		if(state != CL_SUCCESS)
		{
			printf("Error queueing task! %i\n", state);
			return -1;
		}
		else
		{
			clFlush(devices[devicePosition].kernelCommandQueue);

			devices[devicePosition].numberOfEvents += 1;
			return devices[devicePosition].numberOfEvents-1;
		}
	}

	printf("Error! Couldn't find kernel position %i or number of events %i exceeded limit.\n", kernelPosition, devices[devicePosition].numberOfEvents);
	return -1;
}

void SynchronizeCommandQueue(int devicePosition)
{
	clFinish(devices[devicePosition].kernelCommandQueue);
	clFinish(devices[devicePosition].dataCommandQueue);
	devices[devicePosition].numberOfEvents = 0;
}

void SynchronizeEvent(int devicePosition, int eventPosition)
{
	clWaitForEvents(1, &devices[devicePosition].events[eventPosition]);
}

long int GetEventTaskOverheadTicks(int devicePosition, int eventPosition)
{
	long int ticksStart;
	long int ticksEnd;

	clGetEventProfilingInfo(devices[devicePosition].events[eventPosition], CL_PROFILING_COMMAND_QUEUED, sizeof(long int), &ticksStart, NULL);
	clGetEventProfilingInfo(devices[devicePosition].events[eventPosition], CL_PROFILING_COMMAND_START, sizeof(long int), &ticksEnd, NULL);
	return (ticksEnd - ticksStart);
}

long int GetEventTaskTicks(int devicePosition, int eventPosition)
{
	long int ticksStart;
	long int ticksEnd;

	clGetEventProfilingInfo(devices[devicePosition].events[eventPosition], CL_PROFILING_COMMAND_START, sizeof(long int), &ticksStart, NULL);
	clGetEventProfilingInfo(devices[devicePosition].events[eventPosition], CL_PROFILING_COMMAND_END, sizeof(long int), &ticksEnd, NULL);
	return (ticksEnd - ticksStart);
}

cl_device_type GetDeviceType(int devicePosition)
{
	return devices[devicePosition].deviceType;
}

size_t GetDeviceMaxWorkItemsPerWorkGroup(int devicePosition)
{
	return devices[devicePosition].deviceMaxWorkItemsPerWorkGroup;
}

cl_uint GetDeviceComputeUnits(int devicePosition)
{
	return devices[devicePosition].deviceComputeUnits;
}

bool isDeviceCPU(int devicePosition)
{
	return devices[devicePosition].deviceType == CL_DEVICE_TYPE_CPU ? true : false;
}

