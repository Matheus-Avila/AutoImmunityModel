#include <stdio.h>
#include <mpi.h>

#include <sys/time.h>
#include <stdio.h>
#include <string.h>

#include "Balanceador.h"
#include "OpenCLWrapper.h"

//********************************************************************************************************************************
//Esta versao nao tem sobreposicao de comunicacao com computacao, provavelmente existe um problema no driver do OpenCL, verificar.
//********************************************************************************************************************************

// #define PRINT
#define CPU_WORK_GROUP_SIZE	6
#define GPU_WORK_GROUP_SIZE	48
#define SIMULACOES		1/0.0002
#define INTERVALO_BALANCEAMENTO	1000
#define BALANCEAMENTO_THRESHOLD	0.000025
#define PRECISAO_BALANCEAMENTO	10

//Tipos de celulas.
#define MICROGLIA		0
#define OLIGODENDROCYTES		1
#define CONVENTIONAL_DC		2
#define ACTIVATED_DC		3
#define T_CYTOTOXIC		4
#define ANTIBODY		5
#define MALHA_TOTAL_CELULAS	6

//Informacoes de acesso à estrutura "parametrosMalha".
#define OFFSET_COMPUTACAO               0
#define LENGTH_COMPUTACAO               1
#define COMPRIMENTO_GLOBAL_X            2
#define COMPRIMENTO_GLOBAL_Y            3
#define COMPRIMENTO_GLOBAL_Z            4
#define MALHA_DIMENSAO_POSICAO_Z        5
#define MALHA_DIMENSAO_POSICAO_Y        6
#define MALHA_DIMENSAO_POSICAO_X        7
#define MALHA_DIMENSAO_CELULAS          8
#define NUMERO_PARAMETROS_MALHA         9

//Habilitar.
#define HABILITAR_ESTATICO
#define HABILITAR_BALANCEAMENTO
#define HABILITAR_BENCHMARK

void InicializarParametrosMalhaHIS(int **parametrosMalha, unsigned int offsetComputacao, unsigned int lengthComputacao, unsigned int xMalhaLength, unsigned int yMalhaLength, unsigned int zMalhaLength)
{
	*parametrosMalha = new int[NUMERO_PARAMETROS_MALHA];

	(*parametrosMalha)[OFFSET_COMPUTACAO] = offsetComputacao;
	(*parametrosMalha)[LENGTH_COMPUTACAO] = lengthComputacao;
	(*parametrosMalha)[COMPRIMENTO_GLOBAL_X] = xMalhaLength;
	(*parametrosMalha)[COMPRIMENTO_GLOBAL_Y] = yMalhaLength;
	(*parametrosMalha)[COMPRIMENTO_GLOBAL_Z] = zMalhaLength;
	(*parametrosMalha)[MALHA_DIMENSAO_POSICAO_Z] = yMalhaLength*xMalhaLength*MALHA_TOTAL_CELULAS;
	(*parametrosMalha)[MALHA_DIMENSAO_POSICAO_Y] = xMalhaLength*MALHA_TOTAL_CELULAS;
	(*parametrosMalha)[MALHA_DIMENSAO_POSICAO_X] = MALHA_TOTAL_CELULAS;
	(*parametrosMalha)[MALHA_DIMENSAO_CELULAS] = 1;
}

void InicializarPontosHIS(float **malha, const int *parametrosMalha)
{
	*malha = new float[parametrosMalha[COMPRIMENTO_GLOBAL_X]*parametrosMalha[COMPRIMENTO_GLOBAL_Y]*parametrosMalha[COMPRIMENTO_GLOBAL_Z]*MALHA_TOTAL_CELULAS*sizeof(float)]; // esse sizeof(float) está errado
	for(unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++)
	{	
		for(unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++)
		{
			for(unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++)
			{
				if(x >= (0.40f*parametrosMalha[COMPRIMENTO_GLOBAL_X] && x <= (0.60f*parametrosMalha[COMPRIMENTO_GLOBAL_X])) && y >= (0.40f*parametrosMalha[COMPRIMENTO_GLOBAL_Y] && y <= (0.60f*parametrosMalha[COMPRIMENTO_GLOBAL_Y])))
				{
					(*malha)[(MICROGLIA * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 100.0f;
				}
				else
				{
					(*malha)[(MICROGLIA * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
				}

				(*malha)[(OLIGODENDROCYTES * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
				(*malha)[(CONVENTIONAL_DC * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
				(*malha)[(ACTIVATED_DC * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
				(*malha)[(T_CYTOTOXIC * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
				(*malha)[(ANTIBODY * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])] = 0.0f;
			}
		}
	}
}

void LerPontosHIS(const float *malha, const int *parametrosMalha)
{
	for(unsigned int x = 0; x < parametrosMalha[COMPRIMENTO_GLOBAL_X]; x++)
	{
		for(unsigned int y = 0; y < parametrosMalha[COMPRIMENTO_GLOBAL_Y]; y++)
		{
			for(unsigned int z = 0; z < parametrosMalha[COMPRIMENTO_GLOBAL_Z]; z++)
			{
				if((OLIGODENDROCYTES * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X]) >= parametrosMalha[OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS && (OLIGODENDROCYTES * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X]) < (parametrosMalha[OFFSET_COMPUTACAO]+parametrosMalha[LENGTH_COMPUTACAO])*MALHA_TOTAL_CELULAS)
				{
					printf("%f ", malha[(OLIGODENDROCYTES * parametrosMalha[MALHA_DIMENSAO_CELULAS]) + (z * parametrosMalha[MALHA_DIMENSAO_POSICAO_Z]) + (y * parametrosMalha[MALHA_DIMENSAO_POSICAO_Y]) + (x *parametrosMalha[MALHA_DIMENSAO_POSICAO_X])]);
				}
				else
				{
					printf("%f ", 0.0f);
				}
			}
			printf("\n");
		}
	}
}

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("Por favor, digite as dimensões da malha.\n");
		return 0;
	}

	int xMalhaLength = atoi(argv[1]);
	int yMalhaLength = atoi(argv[2]);
	int zMalhaLength = atoi(argv[3]);

	//***************
	//*Inicializacao.
	//***************

	#ifdef HABILITAR_BENCHMARK
	double tempoInicio, tempoFim;
	double tempoComputacaoInterna = 0.0;
	double tempoTrocaBorda = 0.0;
	double tempoComputacaoBorda = 0.0;
	double tempoBalanceamento = 0.0;
	#endif

	int world_size;
	int world_rank;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	MPI_Request sendRequest, receiveRequest;

	int dispositivos = InitParallelProcessor();

	int dispositivosLocal[world_size];
	int dispositivosWorld[world_size];
	memset(dispositivosLocal, 0, sizeof(int)*world_size);
	dispositivosLocal[world_rank] = dispositivos;
	MPI_Reduce(dispositivosLocal, dispositivosWorld, world_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(dispositivosWorld, world_size, MPI_INT, 0, MPI_COMM_WORLD);

	int todosDispositivos = 0;
	int meusDispositivosOffset;
	int meusDispositivosLength;
	for(int count = 0; count < world_size; count++)
	{
		if(count == world_rank)
		{
			meusDispositivosOffset = todosDispositivos;
			meusDispositivosLength = dispositivosWorld[count];
		}
		todosDispositivos += dispositivosWorld[count];
	}

	int *parametrosMalha[todosDispositivos];
	float *malhaSwapBuffer[2];
	long int tempos[todosDispositivos];
	float cargasNovas[todosDispositivos];
	float cargasAntigas[todosDispositivos];

	int parametrosMalhaDispositivo[todosDispositivos];
	int malhaSwapBufferDispositivo[todosDispositivos][2];
	int kernelDispositivo[todosDispositivos];
	int dataEventoDispositivo[todosDispositivos];
	int kernelEventoDispositivo[todosDispositivos];

	int offsetComputacao = 0;
	int lengthComputacao = (xMalhaLength*yMalhaLength*zMalhaLength)/todosDispositivos;
	for(int count = 0; count < todosDispositivos; count++)
	{
		if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
		{	
			InicializarParametrosMalhaHIS(&parametrosMalha[count], offsetComputacao, (count+1 == todosDispositivos) ? (xMalhaLength*yMalhaLength*zMalhaLength)-offsetComputacao : lengthComputacao, xMalhaLength, yMalhaLength, zMalhaLength);
			
			printf("offset %d legth %d final %d \n", parametrosMalha[count][OFFSET_COMPUTACAO], parametrosMalha[count][LENGTH_COMPUTACAO], parametrosMalha[count][OFFSET_COMPUTACAO] + parametrosMalha[count][LENGTH_COMPUTACAO]-1);

			if(count == meusDispositivosOffset)
			{
				InicializarPontosHIS(&malhaSwapBuffer[0], parametrosMalha[count]);
				InicializarPontosHIS(&malhaSwapBuffer[1], parametrosMalha[count]);
			}

			parametrosMalhaDispositivo[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(int)*NUMERO_PARAMETROS_MALHA, CL_MEM_READ_ONLY, NULL);
			malhaSwapBufferDispositivo[count][0] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS), CL_MEM_READ_WRITE, NULL);
			malhaSwapBufferDispositivo[count][1] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS), CL_MEM_READ_WRITE, NULL);
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosMalhaDispositivo[count], (char *)parametrosMalha[count], 0, sizeof(int)*NUMERO_PARAMETROS_MALHA);
			WriteToMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)malhaSwapBuffer[0], 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS));
			WriteToMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][1], (char *)malhaSwapBuffer[1], 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS));
			SynchronizeCommandQueue(count-meusDispositivosOffset);

			kernelDispositivo[count] = CreateKernel(count-meusDispositivosOffset, "kernels.cl", "ProcessarPontos");
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][0]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][1]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 2, parametrosMalhaDispositivo[count]);
		}
		offsetComputacao += lengthComputacao;
	}

	for(int count = 0; count < todosDispositivos; count++)
	{
		cargasNovas[count] = ((float)(count+1))*(1.0f/((float)todosDispositivos));
		cargasAntigas[count] = cargasNovas[count];
		tempos[count] = 1;
	}

	//*******
	//Tempos.
	//*******

	struct timeval timeStart, timeEnd;
	gettimeofday(&timeStart, NULL);

	//**********
	//Simulaço.
	//**********

	bool balanceamento = false;
	#ifdef HABILITAR_BALANCEAMENTO
	balanceamento = true;
	#endif

	for(int simulacao = 0; simulacao < SIMULACOES; simulacao++)
	{
		//Balanceamento de carga.
		if(balanceamento && ((simulacao == 0) || (simulacao == 1) || (simulacao % INTERVALO_BALANCEAMENTO == 0)))
		{
			
			#ifdef HABILITAR_ESTATICO
			if(simulacao > 1)
			{
				balanceamento = false;
			}
			#endif
			
			#ifdef HABILITAR_BENCHMARK
			tempoInicio = MPI_Wtime();
			#endif

			//Precisao do balanceamento.
			memset(tempos, 0, sizeof(long int)*todosDispositivos);
			for(int precisao = 0; precisao < PRECISAO_BALANCEAMENTO; precisao++)
			{
				//Computação.
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{
						if((simulacao%2)==0)
						{
							SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][0]);
							SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][1]);
						}
						else
						{
							SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][1]);
							SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][0]);
						}
	
						kernelEventoDispositivo[count] = RunKernel(count-meusDispositivosOffset, kernelDispositivo[count], parametrosMalha[count][OFFSET_COMPUTACAO], parametrosMalha[count][LENGTH_COMPUTACAO], isDeviceCPU(count-meusDispositivosOffset) ? CPU_WORK_GROUP_SIZE :  GPU_WORK_GROUP_SIZE);
					}
				}

				//Tempos.
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{
						SynchronizeCommandQueue(count-meusDispositivosOffset);
						tempos[count] += GetEventTaskTicks(count-meusDispositivosOffset, kernelEventoDispositivo[count]);
					}
				}
			}

			//Reduzir tempos.
			long int temposRoot[todosDispositivos];
			MPI_Reduce(tempos, temposRoot, todosDispositivos, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(temposRoot, todosDispositivos, MPI_LONG, 0, MPI_COMM_WORLD);
			memcpy(tempos, temposRoot, sizeof(long int)*todosDispositivos);
			ComputarCargas(tempos, cargasAntigas, cargasNovas, todosDispositivos);

			//Computar novas cargas.
			if(ComputarNorma(cargasAntigas, cargasNovas, todosDispositivos) > BALANCEAMENTO_THRESHOLD)
			{
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{
						// int overlapNovoOffset = ((int)(((count == 0) ? 0.0f : cargasNovas[count-1])*((float)(xMalhaLength*yMalhaLength*zMalhaLength))));
						int aux = 0;
						if(count == todosDispositivos-1)
							aux = (xMalhaLength*yMalhaLength*zMalhaLength) - (parametrosMalha[count-1][OFFSET_COMPUTACAO] + parametrosMalha[count-1][LENGTH_COMPUTACAO]);
						int overlapNovoOffset = (count == 0) ? 0 : parametrosMalha[count-1][OFFSET_COMPUTACAO] + parametrosMalha[count-1][LENGTH_COMPUTACAO];
						int overlapNovoLength = ((int)(((count == 0) ? cargasNovas[count]-0.0f : cargasNovas[count]-cargasNovas[count-1])*((float)(xMalhaLength*yMalhaLength*zMalhaLength))));
						if(count == todosDispositivos -1 ){
							overlapNovoLength = aux;
						}
						for(int count2 = 0; count2 < todosDispositivos; count2++)
						{
							if(count > count2)
							{
								//Atender requisicoes de outros processos.
								if(RecuperarPosicaoHistograma(dispositivosWorld, world_size, count) != RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2))
								{
									int overlap[2];
									int alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2);
									float *malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
									int malhaDevice = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
									MPI_Recv(overlap, 2, MPI_INT, alvo, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								
									//Podem ocorrer requisicoes vazias.
									if(overlap[1] > 0)
									{
										ReadFromMemoryObject(count-meusDispositivosOffset, malhaDevice, (char *)(malha+(overlap[0]*MALHA_TOTAL_CELULAS)), overlap[0]*MALHA_TOTAL_CELULAS*sizeof(float), overlap[1]*MALHA_TOTAL_CELULAS*sizeof(float));
										SynchronizeCommandQueue(count-meusDispositivosOffset);
										MPI_Send(malha+(overlap[0]*MALHA_TOTAL_CELULAS), overlap[1]*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD);
									}
								}
							}
							else if(count < count2)
							{
								//Fazer requisicoes a outros processos.
								int overlapAntigoOffset = ((int)(((count2 == 0) ? 0 : cargasAntigas[count2-1])*(xMalhaLength*yMalhaLength*zMalhaLength)));
								int overlapAntigoLength = ((int)(((count2 == 0) ? cargasAntigas[count2]-0.0f : cargasAntigas[count2]-cargasAntigas[count2-1])*(xMalhaLength*yMalhaLength*zMalhaLength)));

								int intersecaoOffset;
								int intersecaoLength;

								if(	((overlapAntigoOffset <= overlapNovoOffset-(xMalhaLength*yMalhaLength)) && ComputarIntersecao(overlapAntigoOffset, overlapAntigoLength, overlapNovoOffset-(xMalhaLength*yMalhaLength), overlapNovoLength+(xMalhaLength*yMalhaLength), &intersecaoOffset, &intersecaoLength)) ||
									((overlapAntigoOffset > overlapNovoOffset-(xMalhaLength*yMalhaLength)) && ComputarIntersecao(overlapNovoOffset-(xMalhaLength*yMalhaLength), overlapNovoLength+(xMalhaLength*yMalhaLength), overlapAntigoOffset, overlapAntigoLength, &intersecaoOffset, &intersecaoLength)))
								{
									if(count2 >= meusDispositivosOffset && count2 < meusDispositivosOffset+meusDispositivosLength)
									{
										float *malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];

										int malhaDevice[2] = {	((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1],
													((simulacao%2)==0) ? malhaSwapBufferDispositivo[count2][0] : malhaSwapBufferDispositivo[count2][1]};

										ReadFromMemoryObject(count2-meusDispositivosOffset, malhaDevice[1], (char *)(malha+(intersecaoOffset*MALHA_TOTAL_CELULAS)), intersecaoOffset*MALHA_TOTAL_CELULAS*sizeof(float), intersecaoLength*MALHA_TOTAL_CELULAS*sizeof(float));
										SynchronizeCommandQueue(count2-meusDispositivosOffset);
										
										WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(intersecaoOffset*MALHA_TOTAL_CELULAS)), intersecaoOffset*MALHA_TOTAL_CELULAS*sizeof(float), intersecaoLength*MALHA_TOTAL_CELULAS*sizeof(float));
										SynchronizeCommandQueue(count-meusDispositivosOffset);
									}
									else
									{
										//Fazer uma requisicao.
										if(RecuperarPosicaoHistograma(dispositivosWorld, world_size, count) != RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2))
										{
											int overlap[2] = {intersecaoOffset, intersecaoLength};
											int alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2);
											float *malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
											int malhaDevice = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
											MPI_Send(overlap, 2, MPI_INT, alvo, 0, MPI_COMM_WORLD);
											MPI_Recv(malha+(overlap[0]*MALHA_TOTAL_CELULAS), overlap[1]*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
											
											WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice, (char *)(malha+(overlap[0]*MALHA_TOTAL_CELULAS)), overlap[0]*MALHA_TOTAL_CELULAS*sizeof(float), overlap[1]*MALHA_TOTAL_CELULAS*sizeof(float));
											SynchronizeCommandQueue(count-meusDispositivosOffset);
										}
									}
								}
								else
								{
									//Fazer uma requisicao vazia.
									if(RecuperarPosicaoHistograma(dispositivosWorld, world_size, count) != RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2))
									{
										int overlap[2] = {0, 0};
										int alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count2);
										float *malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
										MPI_Send(overlap, 2, MPI_INT, alvo, 0, MPI_COMM_WORLD);
									}
								}
							}
						}
						parametrosMalha[count][OFFSET_COMPUTACAO] = overlapNovoOffset;
						parametrosMalha[count][LENGTH_COMPUTACAO] = overlapNovoLength;

						WriteToMemoryObject(count-meusDispositivosOffset, parametrosMalhaDispositivo[count], (char *)parametrosMalha[count], 0, sizeof(int)*NUMERO_PARAMETROS_MALHA);
						SynchronizeCommandQueue(count-meusDispositivosOffset);
					}
				}
				memcpy(cargasAntigas, cargasNovas, sizeof(float)*todosDispositivos);
			}

			#ifdef HABILITAR_BENCHMARK
			MPI_Barrier(MPI_COMM_WORLD);
			tempoFim = MPI_Wtime();
			tempoBalanceamento += tempoFim-tempoInicio;
			#endif
		}
		else
		{
			#ifdef HABILITAR_BENCHMARK
			tempoInicio = MPI_Wtime();
			#endif

			//Computação interna.
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					RunKernel(count-meusDispositivosOffset, kernelDispositivo[count], parametrosMalha[count][OFFSET_COMPUTACAO]+(xMalhaLength*yMalhaLength), parametrosMalha[count][LENGTH_COMPUTACAO]-(xMalhaLength*yMalhaLength), isDeviceCPU(count-meusDispositivosOffset) ? CPU_WORK_GROUP_SIZE :  GPU_WORK_GROUP_SIZE);
				}
			}

			//Sincronizacao da computação interna.
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					SynchronizeCommandQueue(count-meusDispositivosOffset);
				}
			}

			#ifdef HABILITAR_BENCHMARK
			MPI_Barrier(MPI_COMM_WORLD);
			tempoFim = MPI_Wtime();
			tempoComputacaoInterna += tempoFim-tempoInicio;
			tempoInicio = MPI_Wtime();
			#endif

			//Transferencia de bordas, feita em quatro passos.
			for(int passo = 0; passo < 4; passo++)
			{
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{
						int tamanhoBorda = xMalhaLength*yMalhaLength;
						float *malha;
						int malhaDevice[2];
						int borda[2];
						int alvo;
					
						//Entre processos diferentes, no quarto passo.
						if(passo == 3)
						{
							if(count == meusDispositivosOffset && count > 0)
							{
								malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
								malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
								borda[0] = parametrosMalha[count][OFFSET_COMPUTACAO]-(tamanhoBorda);
								borda[0] = (borda[0] < 0) ? 0 : borda[0];
								borda[1] = parametrosMalha[count][OFFSET_COMPUTACAO];
								alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count-1);

								if(alvo%2 == 0)
								{
									MPI_Irecv(malha+(borda[0]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &receiveRequest);

									dataEventoDispositivo[count] = ReadFromMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);
									MPI_Isend(malha+(borda[1]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &sendRequest);
									MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
									MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
									
									WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);

								}
							}
							if(count == meusDispositivosOffset+meusDispositivosLength-1 && count < todosDispositivos-1)
							{
								malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
								malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
								borda[0] = (parametrosMalha[count][OFFSET_COMPUTACAO]+parametrosMalha[count][LENGTH_COMPUTACAO])-(tamanhoBorda);
								borda[0] = (borda[0] < 0) ? 0 : borda[0];
								borda[1] = (parametrosMalha[count][OFFSET_COMPUTACAO]+parametrosMalha[count][LENGTH_COMPUTACAO]);
								alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count+1);

								if(alvo%2 == 1)
								{
									dataEventoDispositivo[count] = ReadFromMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);
									MPI_Isend(malha+(borda[1]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &sendRequest);
									MPI_Irecv(malha+(borda[0]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &receiveRequest);
									MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
									MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
									WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);

								}
							}
						}

						//Entre processos diferentes, no terceiro passo.
						if(passo == 2)
						{
							if(count == meusDispositivosOffset && count > 0)
							{
								malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
								malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
								borda[0] = parametrosMalha[count][OFFSET_COMPUTACAO]-(tamanhoBorda);
								borda[0] = (borda[0] < 0) ? 0 : borda[0];
								borda[1] = parametrosMalha[count][OFFSET_COMPUTACAO];
								alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count-1);

								if(alvo%2 == 1)
								{
									MPI_Irecv(malha+(borda[0]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &receiveRequest);

									dataEventoDispositivo[count] = ReadFromMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);
									MPI_Isend(malha+(borda[1]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &sendRequest);
									MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
									MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);

									WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);

								}
							}
							if(count == meusDispositivosOffset+meusDispositivosLength-1 && count < todosDispositivos-1)
							{
								malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
								malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count][0] : malhaSwapBufferDispositivo[count][1];
								borda[0] = (parametrosMalha[count][OFFSET_COMPUTACAO]+parametrosMalha[count][LENGTH_COMPUTACAO])-(tamanhoBorda);
								borda[0] = (borda[0] < 0) ? 0 : borda[0];
								borda[1] = (parametrosMalha[count][OFFSET_COMPUTACAO]+parametrosMalha[count][LENGTH_COMPUTACAO]);
								alvo = RecuperarPosicaoHistograma(dispositivosWorld, world_size, count+1);

								if(alvo%2 == 0)
								{
									dataEventoDispositivo[count] = ReadFromMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);
									MPI_Isend(malha+(borda[1]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &sendRequest);
									MPI_Irecv(malha+(borda[0]*MALHA_TOTAL_CELULAS), tamanhoBorda*MALHA_TOTAL_CELULAS, MPI_FLOAT, alvo, 0, MPI_COMM_WORLD, &receiveRequest);
									MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
									MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);

									WriteToMemoryObject(count-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
									SynchronizeCommandQueue(count-meusDispositivosOffset);

								}
							}
						}

						//No mesmo processo, no primeiro passo.
						if(passo == 0 && count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength-1)
						{
							malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
							malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count+0][0] : malhaSwapBufferDispositivo[count+0][1];
							malhaDevice[1] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count+1][0] : malhaSwapBufferDispositivo[count+1][1];
							borda[0] = parametrosMalha[count+1][OFFSET_COMPUTACAO]-(tamanhoBorda);
							borda[0] = (borda[0] < 0) ? 0 : borda[0];
							borda[1] = parametrosMalha[count+1][OFFSET_COMPUTACAO];

							dataEventoDispositivo[count+0] = ReadFromMemoryObject(count+0-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
							SynchronizeCommandQueue(count+0-meusDispositivosOffset);

							dataEventoDispositivo[count+1] = ReadFromMemoryObject(count+1-meusDispositivosOffset, malhaDevice[1], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
							SynchronizeCommandQueue(count+1-meusDispositivosOffset);

						}

						//No mesmo processo, no segundo passo.
						if(passo == 1 && count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength-1)
						{
							malha = ((simulacao%2)==0) ? malhaSwapBuffer[0] : malhaSwapBuffer[1];
							malhaDevice[0] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count+0][0] : malhaSwapBufferDispositivo[count+0][1];
							malhaDevice[1] = ((simulacao%2)==0) ? malhaSwapBufferDispositivo[count+1][0] : malhaSwapBufferDispositivo[count+1][1];
							borda[0] = parametrosMalha[count+1][OFFSET_COMPUTACAO]-(tamanhoBorda);
							borda[0] = (borda[0] < 0) ? 0 : borda[0];
							borda[1] = parametrosMalha[count+1][OFFSET_COMPUTACAO];

							WriteToMemoryObject(count+0-meusDispositivosOffset, malhaDevice[0], (char *)(malha+(borda[1]*MALHA_TOTAL_CELULAS)), borda[1]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
							SynchronizeCommandQueue(count+0-meusDispositivosOffset);

							WriteToMemoryObject(count+1-meusDispositivosOffset, malhaDevice[1], (char *)(malha+(borda[0]*MALHA_TOTAL_CELULAS)), borda[0]*MALHA_TOTAL_CELULAS*sizeof(float), tamanhoBorda*MALHA_TOTAL_CELULAS*sizeof(float));
							SynchronizeCommandQueue(count+1-meusDispositivosOffset);
						}
					}
				}
			}

			//Sincronizacao da comunicacao.
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					SynchronizeCommandQueue(count-meusDispositivosOffset);
				}
			}

			#ifdef HABILITAR_BENCHMARK
			MPI_Barrier(MPI_COMM_WORLD);
			tempoFim = MPI_Wtime();
			tempoTrocaBorda += tempoFim-tempoInicio;
			tempoInicio = MPI_Wtime();
			#endif

			//Computação das bordas.
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					if((simulacao%2)==0)
					{
						SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][0]);
						SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][1]);
					}
					else
					{
						SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][1]);
						SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][0]);
					}

					RunKernel(count-meusDispositivosOffset, kernelDispositivo[count], parametrosMalha[count][OFFSET_COMPUTACAO], xMalhaLength*yMalhaLength, isDeviceCPU(count-meusDispositivosOffset) ? CPU_WORK_GROUP_SIZE :  GPU_WORK_GROUP_SIZE);
					RunKernel(count-meusDispositivosOffset, kernelDispositivo[count], parametrosMalha[count][OFFSET_COMPUTACAO]+parametrosMalha[count][LENGTH_COMPUTACAO]-(xMalhaLength*yMalhaLength), xMalhaLength*yMalhaLength, isDeviceCPU(count-meusDispositivosOffset) ? CPU_WORK_GROUP_SIZE :  GPU_WORK_GROUP_SIZE);
				}
			}

			//Sincronizacao da computação das borda
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					SynchronizeCommandQueue(count-meusDispositivosOffset);
				}
			}

			#ifdef HABILITAR_BENCHMARK
			MPI_Barrier(MPI_COMM_WORLD);
			tempoFim = MPI_Wtime();
			tempoComputacaoBorda += tempoFim-tempoInicio;
			#endif
		}
	}

	//*******
	//Tempos.
	//*******

	MPI_Barrier(MPI_COMM_WORLD);
	int contadorDePontos = 0;

	/******
	******** Saving figures 
	*******/
	// Provavelmente vai funcionar so com uma thread do MPI. nao sei como criar uma regiao critica e so um processo escrver no arquivo por vez
	FILE *file;
	char filename[30];
	sprintf(filename, "results.vtk");
	file = fopen(filename, "w");
	fclose(file);
	for(int count2 = 0; count2 < world_size; count2++)
	{	
		if(count2 == world_rank)
		{
			printf("Malha do processo %i\n", world_rank);
			for(int count = 0; count < todosDispositivos; count++)
			{	
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{	
					printf("Malha do dispositivo %i\n", count);
					ReadFromMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)(malhaSwapBuffer[0]+(parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS)), parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float), parametrosMalha[count][LENGTH_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float));
					SynchronizeCommandQueue(count-meusDispositivosOffset);
					sprintf(filename, "results.vtk");
					file = fopen(filename, "a");
					if(count2 == 0 && count == meusDispositivosOffset){
						fprintf(file, "# vtk DataFile Version 3.0\n");
						fprintf(file, "results.vtk\n");
						fprintf(file, "ASCII\n");
						fprintf(file, "DATASET RECTILINEAR_GRID\n");
						fprintf(file, "DIMENSIONS %d %d %d\n", xMalhaLength, yMalhaLength, zMalhaLength);

						fprintf(file, "X_COORDINATES %d double\n", xMalhaLength);
						for(int i=0; i<xMalhaLength; i++)
						{
							fprintf(file, "%f ", i*(0.1/xMalhaLength));
						}
						fprintf(file, "\n");

						fprintf(file, "Y_COORDINATES %d double\n", yMalhaLength);
						for(int j=0; j<yMalhaLength; j++)
						{
							fprintf(file, "%f ", j*(0.1/yMalhaLength));
						}
						fprintf(file, "\n");

						fprintf(file, "Z_COORDINATES %d double\n", zMalhaLength);
						for(int m=0; m<zMalhaLength; m++)
						{
							fprintf(file, "%f ", m*(0.1/zMalhaLength));
						}
						fprintf(file, "\n");

						fprintf(file, "POINT_DATA %d \n", xMalhaLength*yMalhaLength*zMalhaLength);
						fprintf(file, "FIELD FieldData 1 \n");
						fprintf(file, "Temperatura 1 %d double \n", xMalhaLength * yMalhaLength * zMalhaLength);
					}
					const float *malha = malhaSwapBuffer[0];
					const int *param = parametrosMalha[count];
					printf("count %d: inicio %d -- length %d final: %d \n", count, param[OFFSET_COMPUTACAO], param[LENGTH_COMPUTACAO], param[OFFSET_COMPUTACAO] + param[LENGTH_COMPUTACAO]-1);
					for(unsigned int x = 0; x <  param[COMPRIMENTO_GLOBAL_X]; x++)
					{
						for(unsigned int y = 0; y <  param[COMPRIMENTO_GLOBAL_Y]; y++)
						{
							for(unsigned int z = 0; z <  param[COMPRIMENTO_GLOBAL_Z]; z++)
							{
								if((OLIGODENDROCYTES * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) >= param[OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS && (OLIGODENDROCYTES * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) <  ((param[OFFSET_COMPUTACAO]+param[LENGTH_COMPUTACAO])*MALHA_TOTAL_CELULAS) )
								{	contadorDePontos++;
									fprintf(file, "%f \n", malha[(OLIGODENDROCYTES * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X])]);
								}
							}
						}
					}
					fclose(file);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	printf("contadorDePontos %d\n", contadorDePontos);
	// sprintf(filename, "results.vtk");
	// file = fopen(filename, "a");
	// fprintf(file, "\n");
	// fclose(file);

	if(world_rank == 0)
	{
		gettimeofday(&timeEnd, NULL);
		printf("Overall ticks (1tick->1ms): %lu\n", (timeEnd.tv_sec - timeStart.tv_sec)*1000000 + (timeEnd.tv_usec - timeStart.tv_usec));

		#ifdef HABILITAR_BENCHMARK
		printf("Internal computation (s): %f\n", tempoComputacaoInterna);
		printf("Border swap time (s): %f\n", tempoTrocaBorda);
		printf("Border computation time (s): %f\n", tempoComputacaoBorda);
		printf("Balancing time (s): %f\n", tempoBalanceamento);
		#endif

		for(int count = 0; count < todosDispositivos; count++)
		{
			printf("Tempo dispositivo (1tick->1nanosegundo) %i: %li\n", count, tempos[count]);
		}
	}

	//************
	//Finalização.
	//************

	#ifdef PRINT
	for(int count2 = 0; count2 < world_size; count2++)
	{
		if(count2 == world_rank)
		{
			printf("Malha do processo %i\n", world_rank);
			for(int count = 0; count < todosDispositivos; count++)
			{
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{
					printf("Malha do dispositivo %i\n", count);
					ReadFromMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)(malhaSwapBuffer[0]+(parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS)), parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float), parametrosMalha[count][LENGTH_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float));
					SynchronizeCommandQueue(count-meusDispositivosOffset);
					LerPontosHIS(malhaSwapBuffer[0], parametrosMalha[count]);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	#endif
	
	for(int count = 0; count < todosDispositivos; count++)
	{
		if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
		{
			delete [] parametrosMalha[count];
			parametrosMalha[count] = NULL;
		}
	}
	delete [] malhaSwapBuffer[0];
	delete [] malhaSwapBuffer[1];
	malhaSwapBuffer[0] = NULL;
	malhaSwapBuffer[1] = NULL;

	FinishParallelProcessor();
	MPI_Finalize();
	return 0;
}

