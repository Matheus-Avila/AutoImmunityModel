#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>

#include "Balanceador.h"
#include "OpenCLWrapper.h"

//********************************************************************************************************************************
//Esta versao nao tem sobreposicao de comunicacao com computacao, provavelmente existe um problema no driver do OpenCL, verificar.
//********************************************************************************************************************************

# define NUM_POINTSLN 1000
// #define PRINT
#define CPU_WORK_GROUP_SIZE	6
#define GPU_WORK_GROUP_SIZE	48
#define SIMULACOES		30/0.0002
#define INTERVALO_BALANCEAMENTO	1000
#define BALANCEAMENTO_THRESHOLD	0.000125 //000125: quarto de face; 00025: metade de face; 0005: face inteira 
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
// #define HABILITAR_BENCHMARK

void WritePopulationLymphNode(float *population, char* fileName){
	FILE *file;
	file = fopen(fileName, "w");
	if(file != NULL){
		for(int i=0;i<NUM_POINTSLN;i++){
			fprintf(file, "%f\n", population[i]);
		}
		fclose(file);
	}else{
		printf("Error matrix file\n");
		exit(0);
	}
}

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
				// if((x >= (0.3f*parametrosMalha[COMPRIMENTO_GLOBAL_X])) && (x <= (0.5f*parametrosMalha[COMPRIMENTO_GLOBAL_X])))
				if( pow((x - 0.5f*parametrosMalha[COMPRIMENTO_GLOBAL_X]),2) + pow((y - 0.5f*parametrosMalha[COMPRIMENTO_GLOBAL_Y]),2) < 5)
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

// Provavelmente vai funcionar so com uma thread do MPI. nao sei como criar uma regiao critica e so um processo escrver no arquivo por vez
//ADAPTAR COM MPI SE RODAR EM 2 MAQUINAS -- ATENCAO
void SaveFigure(int malhaSwapBufferDispositivo[][2], float **malhaSwapBuffer, int **parametrosMalha, int xMalhaLength, int yMalhaLength, int zMalhaLength, int meusDispositivosOffset, int meusDispositivosLength, float time, int world_size, int world_rank, int todosDispositivos)
{
	FILE *file;
	char filename[30];
	char charTime[4];
	sprintf(filename, "result/results");
	snprintf(charTime, sizeof(charTime), "%f", time);
	strcat(filename, charTime);
	strcat(filename, ".vtk");
	file = fopen(filename, "w");
	fclose(file);
	for(int count2 = 0; count2 < world_size; count2++)
	{	
		if(count2 == world_rank)
		{
			// printf("Malha do processo %i\n", world_rank);
			for(int count = 0; count < todosDispositivos; count++)
			{	
				if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
				{	
					// printf("Malha do dispositivo %i\n", count);
					ReadFromMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)(malhaSwapBuffer[0]+(parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS)), parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float), parametrosMalha[count][LENGTH_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float));
					SynchronizeCommandQueue(count-meusDispositivosOffset);
					file = fopen(filename, "a");
					if(count2 == 0 && count == meusDispositivosOffset){
						fprintf(file, "# vtk DataFile Version 2.0\nReally cool data\n");
						fprintf(file, "ASCII\n");
						fprintf(file, "DATASET STRUCTURED_GRID\n");
						fprintf(file, "DIMENSIONS %d %d %d\n", xMalhaLength, yMalhaLength, zMalhaLength);
						fprintf(file, "POINTS %d float\n", xMalhaLength*yMalhaLength*zMalhaLength);
						for(unsigned int x = 0; x <  parametrosMalha[count][COMPRIMENTO_GLOBAL_X]; x++)
						{
							for(unsigned int y = 0; y <  parametrosMalha[count][COMPRIMENTO_GLOBAL_Y]; y++)
							{
								for(unsigned int z = 0; z <  parametrosMalha[count][COMPRIMENTO_GLOBAL_Z]; z++)
								{
									fprintf(file, "%d %d %d\n", x, y, z);
								}
							}
						}
						fprintf(file, "POINT_DATA %d \n", xMalhaLength*yMalhaLength*zMalhaLength);
						fprintf(file, "SCALARS volume_scalars float 1\n");
						fprintf(file, "LOOKUP_TABLE default\n");
					}
					const float *malha = malhaSwapBuffer[0];
					const int *param = parametrosMalha[count];
					// printf("count %d: inicio %d -- length %d final: %d \n", count, param[OFFSET_COMPUTACAO], param[LENGTH_COMPUTACAO], param[OFFSET_COMPUTACAO] + param[LENGTH_COMPUTACAO]-1);
					
					for(unsigned int z = 0; z < param[COMPRIMENTO_GLOBAL_Z]; z++)
					{
						for(unsigned int y = 0; y < param[COMPRIMENTO_GLOBAL_Y]; y++)
						{	
							for(unsigned int x = 0; x < param[COMPRIMENTO_GLOBAL_X]; x++)
							{
								if((OLIGODENDROCYTES * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) >= param[OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS && (OLIGODENDROCYTES * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) <  ((param[OFFSET_COMPUTACAO]+param[LENGTH_COMPUTACAO])*MALHA_TOTAL_CELULAS) )
								{	
									fprintf(file, "%f ", malha[(T_CYTOTOXIC * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X])]);
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
	double tempoInicio, tempoFim;
	#ifdef HABILITAR_BENCHMARK
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

	float cDl = 0.1;
	float cF = 0.1;
	float alphaTHelper = 0.1;
	float alphaTCytotoxic = 0.1;
	float alphaB = 0.1;
	float alphaP = 1;
	float bTHelper = 0.17;
	float bTCytotoxic = 0.001;
	float bRho = 0.6;
	float bRhoB = 3.02;
	float bRhoP = 1.02;
	float rhoTHelper = 2;
	float rhoTCytotoxic = 2;
	float rhoB = 11;
	float rhoP = 3;
	float rhoAntibody = 0.051;
	float stableTHelper = 70;
	float stableTCytotoxic = 40;
	float stableB = 25;
	float stableP = 2.5;
	float deltaT = 0.0002;
	float gammaD = 0.1;
	float gammaAntibody = 0.3;
	float gammaT = 0.1;
	float V_LN = 160;
	int V_BV = 0;
	int V_PV = 0;

	int tamanhoLN[world_size];
	int vetBV[xMalhaLength*yMalhaLength*zMalhaLength];
	int vetPV[xMalhaLength*yMalhaLength*zMalhaLength];

	//ADAPTAR COM MPI SE RODAR EM 2 MAQUINAS -- ATENCAO
	//inicializa BV PV
	int randomVal;
    for(int k = 1; k < xMalhaLength*yMalhaLength*zMalhaLength - 1; k++){
        randomVal = rand() % 100;
        if(randomVal <10){
            V_BV++;
            V_PV++;
            vetBV[k] = 1;
			vetPV[k+1] = 1;
        }else{
			vetBV[k] = 0;
			vetPV[k+1] = 0;
		}
    }

    V_BV = V_BV * 0.5 * 0.5 * 0.5;
    V_PV = V_PV * 0.5 * 0.5 * 0.5;

	float integralTecido[3] = {0, 0, 0};
	float tecidoIntegraisPontoAPonto[xMalhaLength*yMalhaLength*zMalhaLength];
	float tecidoInteiro[xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS];

	float populacoesLinfonodo[2][6];
	populacoesLinfonodo[0][0] = populacoesLinfonodo[1][0] = 0;
	populacoesLinfonodo[0][1] = populacoesLinfonodo[1][1] = stableTHelper;
	populacoesLinfonodo[0][2] = populacoesLinfonodo[1][2] = stableTCytotoxic;
	populacoesLinfonodo[0][3] = populacoesLinfonodo[1][3] = stableB;
	populacoesLinfonodo[0][4] = populacoesLinfonodo[1][4] = 0;
	populacoesLinfonodo[0][5] = populacoesLinfonodo[1][5] = 0;

	int *parametrosMalha[todosDispositivos];
	float *malhaSwapBuffer[2];
	long int tempos[todosDispositivos];
	float cargasNovas[todosDispositivos];
	float cargasAntigas[todosDispositivos];

	int parametrosMalhaDispositivo[todosDispositivos];
	int malhaSwapBufferDispositivo[todosDispositivos][2];

	int parametrosPopulacoesLinfonodo[todosDispositivos];
	int parametrosDendriticaIntegralTecidoPontoAPonto[todosDispositivos];
	int parametrosCytotoxicaIntegralTecidoPontoAPonto[todosDispositivos];
	int parametrosAnticorpoIntegralTecidoPontoAPonto[todosDispositivos];
	int parametrosPosicaoBVTecido[todosDispositivos];
	int parametrosPosicaoPVTecido[todosDispositivos];

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
			parametrosPosicaoBVTecido[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(int)*(xMalhaLength*yMalhaLength*zMalhaLength), CL_MEM_READ_WRITE, NULL);
			parametrosPosicaoPVTecido[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(int)*(xMalhaLength*yMalhaLength*zMalhaLength), CL_MEM_READ_WRITE, NULL);
			parametrosPopulacoesLinfonodo[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(6), CL_MEM_READ_WRITE, NULL);
			parametrosDendriticaIntegralTecidoPontoAPonto[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength), CL_MEM_READ_WRITE, NULL);
			parametrosCytotoxicaIntegralTecidoPontoAPonto[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength), CL_MEM_READ_WRITE, NULL);
			parametrosAnticorpoIntegralTecidoPontoAPonto[count] = CreateMemoryObject(count-meusDispositivosOffset, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength), CL_MEM_READ_WRITE, NULL);


			WriteToMemoryObject(count-meusDispositivosOffset, parametrosMalhaDispositivo[count], (char *)parametrosMalha[count], 0, sizeof(int)*NUMERO_PARAMETROS_MALHA);
			WriteToMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)malhaSwapBuffer[0], 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS));
			WriteToMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][1], (char *)malhaSwapBuffer[1], 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength*MALHA_TOTAL_CELULAS));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosPosicaoBVTecido[count], (char *)vetBV, 0, sizeof(int)*(xMalhaLength*yMalhaLength*zMalhaLength));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosPosicaoPVTecido[count], (char *)vetPV, 0, sizeof(int)*(xMalhaLength*yMalhaLength*zMalhaLength));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosPopulacoesLinfonodo[count], (char *)populacoesLinfonodo[0], 0, sizeof(float)*(6));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosDendriticaIntegralTecidoPontoAPonto[count], (char *)tecidoIntegraisPontoAPonto, 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosCytotoxicaIntegralTecidoPontoAPonto[count], (char *)tecidoIntegraisPontoAPonto, 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength));
			WriteToMemoryObject(count-meusDispositivosOffset, parametrosAnticorpoIntegralTecidoPontoAPonto[count], (char *)tecidoIntegraisPontoAPonto, 0, sizeof(float)*(xMalhaLength*yMalhaLength*zMalhaLength));

			SynchronizeCommandQueue(count-meusDispositivosOffset);

			kernelDispositivo[count] = CreateKernel(count-meusDispositivosOffset, "kernels.cl", "ProcessarPontos");
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 0, malhaSwapBufferDispositivo[count][0]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 1, malhaSwapBufferDispositivo[count][1]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 2, parametrosMalhaDispositivo[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 3, parametrosPosicaoBVTecido[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 4, parametrosPosicaoPVTecido[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 5, parametrosPopulacoesLinfonodo[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 6, parametrosDendriticaIntegralTecidoPontoAPonto[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 7, parametrosCytotoxicaIntegralTecidoPontoAPonto[count]);
			SetKernelAttribute(count-meusDispositivosOffset, kernelDispositivo[count], 8, parametrosAnticorpoIntegralTecidoPontoAPonto[count]);

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

	float dcLN         = 0;
    float tCytoLN      = stableTHelper;
    float tHelperLN    = stableTCytotoxic;
    float bCellLN      = stableB;
    float plasmaCellLN = 0;
    float antibodyLN   = 0;

	float activatedDcMigration = 0;
	float activatedDcClearance = 0;
	float tCytoActivation = 0;
	float tCytoHomeostasis = 0;
	float tCytoMigration = 0;
	float tHelperActivation = 0;
	float tHelperHomeostasis = 0;
	float tHelperDispendure = 0;
	float bCellActivation = 0;
	float bcellHomeostasis = 0;
	float plasmaActivation = 0;
	float plasmaHomeostasis = 0;
	float antibodyProduction = 0;
	float antibodyDecayment = 0;
	float antibodyMigration = 0;
	float dendriticLymphNodeSavedPoints[NUM_POINTSLN];
	float tCytotoxicLymphNodeSavedPoints[NUM_POINTSLN];
	float tHelperLymphNodeSavedPoints[NUM_POINTSLN];
	float bCellLymphNodeSavedPoints[NUM_POINTSLN];
	float plasmaCellLymphNodeSavedPoints[NUM_POINTSLN];
	float antibodyLymphNodeSavedPoints[NUM_POINTSLN];

	float resultadoParcialEdo[6];
	
	int intervalPoints = (int)(SIMULACOES/NUM_POINTSLN);
	int stepKPlus, stepKMinus;

	tempoInicio = MPI_Wtime();

	for(int simulacao = 0; simulacao < SIMULACOES; simulacao++)
	{	
		/*********************SOLUCAO DO LINFONODO ***********************/
		stepKPlus = simulacao%2;
		stepKMinus = !(stepKPlus && 1);

		dcLN = populacoesLinfonodo[stepKMinus][0];		
		tCytoLN = populacoesLinfonodo[stepKMinus][1];
		tHelperLN = populacoesLinfonodo[stepKMinus][2];
		bCellLN = populacoesLinfonodo[stepKMinus][3];
		plasmaCellLN = populacoesLinfonodo[stepKMinus][4];
		antibodyLN = populacoesLinfonodo[stepKMinus][5];

		//Dendritic cell
		activatedDcMigration = gammaD * (integralTecido[0] - dcLN) * (float)(V_PV/V_LN);
		activatedDcClearance = cDl * dcLN;
		resultadoParcialEdo[0] = activatedDcMigration - activatedDcClearance;

		//T Cytotoxic
		tCytoActivation = bTCytotoxic * (rhoTCytotoxic*tCytoLN*dcLN - tCytoLN*dcLN);
		tCytoHomeostasis = alphaTCytotoxic * (stableTCytotoxic - tCytoLN);
		tCytoMigration = gammaT * (tCytoLN - integralTecido[1]) * (float)(V_BV/V_LN);
		resultadoParcialEdo[1] = tCytoActivation + tCytoHomeostasis - tCytoMigration;

		//T Helper
		tHelperActivation = bTHelper * (rhoTHelper * tHelperLN * dcLN - tHelperLN * dcLN);
		tHelperHomeostasis = alphaTHelper * (stableTHelper - tHelperLN);
		tHelperDispendure = bRho * dcLN * tHelperLN * bCellLN;
		resultadoParcialEdo[2] = tHelperActivation + tHelperHomeostasis - tHelperDispendure;

		//B Cell
		bCellActivation = bRhoB * (rhoB * tHelperLN * dcLN - tHelperLN * dcLN * bCellLN);
		bcellHomeostasis = alphaB * (stableB - bCellLN);
		resultadoParcialEdo[3] = bcellHomeostasis + bCellActivation;

		//Plasma Cells
		plasmaActivation = bRhoP * (rhoP * tHelperLN * dcLN * bCellLN);
		plasmaHomeostasis = alphaP * (stableP - plasmaCellLN);
		resultadoParcialEdo[4] = plasmaHomeostasis + plasmaActivation;

		//Antibody
		antibodyProduction = rhoAntibody * plasmaCellLN;
		antibodyDecayment = cF * antibodyLN;
		antibodyMigration = gammaAntibody * (antibodyLN - integralTecido[2]) * (float)(V_BV/V_LN);
		resultadoParcialEdo[5] = antibodyProduction - antibodyMigration - antibodyDecayment;

		populacoesLinfonodo[stepKPlus][0] = populacoesLinfonodo[stepKMinus][0] + deltaT * resultadoParcialEdo[0];
		populacoesLinfonodo[stepKPlus][1] = populacoesLinfonodo[stepKMinus][1] + deltaT * resultadoParcialEdo[1];
		populacoesLinfonodo[stepKPlus][2] = populacoesLinfonodo[stepKMinus][2] + deltaT * resultadoParcialEdo[2];
		populacoesLinfonodo[stepKPlus][3] = populacoesLinfonodo[stepKMinus][3] + deltaT * resultadoParcialEdo[3];
		populacoesLinfonodo[stepKPlus][4] = populacoesLinfonodo[stepKMinus][4] + deltaT * resultadoParcialEdo[4];
		populacoesLinfonodo[stepKPlus][5] = populacoesLinfonodo[stepKMinus][5] + deltaT * resultadoParcialEdo[5];

		if(simulacao%intervalPoints){
			int posSave = simulacao/intervalPoints;
			dendriticLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][0];
			tCytotoxicLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][1];
			tHelperLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][2];
			bCellLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][3];
			plasmaCellLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][4];
			antibodyLymphNodeSavedPoints[posSave] = populacoesLinfonodo[stepKPlus][5];
		}

		// /*ENVIAR POPULACAO LN PARA KERNELS*/
		for(int count2 = 0; count2 < world_size; count2++)
		{
			if(count2 == world_rank)
			{
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{	
						WriteToMemoryObject(count-meusDispositivosOffset, parametrosPopulacoesLinfonodo[count], (char *)populacoesLinfonodo[stepKPlus], 0, sizeof(float)*(6));
						
						SynchronizeCommandQueue(count-meusDispositivosOffset);
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		#ifdef SAVEFIGURES
		if(simulacao == 0){
			SaveFigure(malhaSwapBufferDispositivo, malhaSwapBuffer, parametrosMalha, xMalhaLength, yMalhaLength, zMalhaLength, meusDispositivosOffset, meusDispositivosLength, simulacao*0.0002, world_size, world_rank, todosDispositivos);
		}
		#endif
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

		integralTecido[0] = 0.0;
		integralTecido[1] = 0.0;
		integralTecido[2] = 0.0;

		for(int count2 = 0; count2 < world_size; count2++)
		{
			if(count2 == world_rank)
			{
				for(int count = 0; count < todosDispositivos; count++)
				{
					if(count >= meusDispositivosOffset && count < meusDispositivosOffset+meusDispositivosLength)
					{
						ReadFromMemoryObject(count-meusDispositivosOffset, malhaSwapBufferDispositivo[count][0], (char *)(tecidoInteiro +(parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS)), parametrosMalha[count][OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float), parametrosMalha[count][LENGTH_COMPUTACAO]*MALHA_TOTAL_CELULAS*sizeof(float));
						SynchronizeCommandQueue(count-meusDispositivosOffset);
						const float *malha = tecidoInteiro;
						const int *param = parametrosMalha[count];
						// printf("count %d: inicio %d -- length %d final: %d \n", count, param[OFFSET_COMPUTACAO], param[LENGTH_COMPUTACAO], param[OFFSET_COMPUTACAO] + param[LENGTH_COMPUTACAO]-1);
						
						for(unsigned int z = 0; z < param[COMPRIMENTO_GLOBAL_Z]; z++)
						{
							for(unsigned int y = 0; y < param[COMPRIMENTO_GLOBAL_Y]; y++)
							{	
								for(unsigned int x = 0; x < param[COMPRIMENTO_GLOBAL_X]; x++)
								{
									if((MICROGLIA * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) >= param[OFFSET_COMPUTACAO]*MALHA_TOTAL_CELULAS && (MICROGLIA * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X]) <  ((param[OFFSET_COMPUTACAO]+param[LENGTH_COMPUTACAO])*MALHA_TOTAL_CELULAS) )
									{	
										integralTecido[0] += malha[(MICROGLIA * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X])] * vetPV[(z * param[MALHA_DIMENSAO_POSICAO_Z] / 6) + (y * param[MALHA_DIMENSAO_POSICAO_Y] / 6) + (x *param[MALHA_DIMENSAO_POSICAO_X] / 6)];
										integralTecido[1] += malha[(MICROGLIA * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X])] * vetBV[(z * param[MALHA_DIMENSAO_POSICAO_Z] / 6) + (y * param[MALHA_DIMENSAO_POSICAO_Y] / 6) + (x *param[MALHA_DIMENSAO_POSICAO_X] / 6)];
										integralTecido[2] += malha[(MICROGLIA * param[MALHA_DIMENSAO_CELULAS]) + (z * param[MALHA_DIMENSAO_POSICAO_Z]) + (y * param[MALHA_DIMENSAO_POSICAO_Y]) + (x *param[MALHA_DIMENSAO_POSICAO_X])] * vetBV[(z * param[MALHA_DIMENSAO_POSICAO_Z] / 6) + (y * param[MALHA_DIMENSAO_POSICAO_Y] / 6) + (x *param[MALHA_DIMENSAO_POSICAO_X] / 6)];
									}
								}
							}
						}
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Percorrer tecidoInteiro e colocar o resultado nas variáveis de integral do tecido
		#ifdef SAVEFIGURES
		if(simulacao != 0 && simulacao%((int)(1/0.0002)) == 0){
			printf("simulacao %d ** total %d\n",simulacao, SIMULACOES);
			SaveFigure(malhaSwapBufferDispositivo, malhaSwapBuffer, parametrosMalha, xMalhaLength, yMalhaLength, zMalhaLength, meusDispositivosOffset, meusDispositivosLength, simulacao*0.0002, world_size, world_rank, todosDispositivos);
		}
		#endif
	}

	//*******
	//Tempos.
	//*******
	/******
	******** Saving figures 
	*******/
	#ifdef SAVEFIGURES
	SaveFigure(malhaSwapBufferDispositivo, malhaSwapBuffer, parametrosMalha, xMalhaLength, yMalhaLength, zMalhaLength, meusDispositivosOffset, meusDispositivosLength, SIMULACOES*0.0002, world_size, world_rank, todosDispositivos);
	#endif
	MPI_Barrier(MPI_COMM_WORLD);
	tempoFim = MPI_Wtime();
	
	#ifdef GPU_DEVICES
	#ifdef HABILITAR_BALANCEAMENTO
	#ifdef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/GPU_DEVICES-HABILITAR_ESTATICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#ifndef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/GPU_DEVICES-DINAMICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#endif
	#endif

	#ifdef ALL_DEVICES
	#ifdef HABILITAR_BALANCEAMENTO
	#ifdef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/ALL_DEVICES-HABILITAR_ESTATICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#ifndef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/ALL_DEVICES-DINAMICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#endif
	#endif

	#ifdef CPU_DEVICES
	#ifdef HABILITAR_BALANCEAMENTO
	#ifdef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/CPU_DEVICES-HABILITAR_ESTATICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#ifndef HABILITAR_ESTATICO
	char filename[200];
	char charThreshold[10];
	sprintf(filename, "tempo/CPU_DEVICES-DINAMICO-BALANCEAMENTO_THRESHOLD:");
	snprintf(charThreshold, sizeof(charThreshold), "%f", BALANCEAMENTO_THRESHOLD);
	strcat(filename, charThreshold);
	strcat(filename, ".txt");
	FILE *file;
	file = fopen(filename, "a");
	fprintf(file, "%lf\n",tempoFim-tempoInicio);
	fclose(file);
	#endif
	#endif
	#endif



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
	#ifdef SAVEFIGURES
	WritePopulationLymphNode(dendriticLymphNodeSavedPoints, "./result/dendritic.txt");
    WritePopulationLymphNode(tCytotoxicLymphNodeSavedPoints, "./result/tHelper.txt");
    WritePopulationLymphNode(tHelperLymphNodeSavedPoints, "./result/tCyto.txt");
    WritePopulationLymphNode(bCellLymphNodeSavedPoints, "./result/bCell.txt");
    WritePopulationLymphNode(plasmaCellLymphNodeSavedPoints, "./result/plasmaCell.txt");
    WritePopulationLymphNode(antibodyLymphNodeSavedPoints, "./result/antibody.txt");
	char buffer[10];
    char command[40] = {};
    strcat(command, "python3 plotLymphNode.py ");
    snprintf(buffer, sizeof(buffer), "%d", SIMULACOES*deltaT);
    strcat(command, buffer);
    strcat(command, " ");
    snprintf(buffer, sizeof(buffer), "%f", deltaT);
    strcat(command, buffer);
    // system(command);
    #endif
	
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

