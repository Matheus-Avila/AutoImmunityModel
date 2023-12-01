//*****************
//Orientacao Malha.
//*****************
//
//	 (Z)
//     	  |  
//	  | 
//	  |____(X)
//       /
//      /
//     /
//    (Y)
//
//*****************

//Tipos de celulas.
#define MALHA_TOTAL_CELULAS	6 // quantidade de células na malha
#define MICROGLIA		0 // célula 1
#define OLIGODENDROCYTES		1 // célula 2
#define CONVENTIONAL_DC		2
#define ACTIVATED_DC		3
#define T_CYTOTOXIC		4
#define ANTIBODY		5

//Estrutura de celulas.
#define CELULAS_SIZEOF			48 // qt de células * sizeof(double)
#define CELULAS_SINGLE_CELL_SIZEOF	8	//Pode ser usado junto com as definicoes de tipos de celulas para substituir os valores abaixo.
#define MICROGLIA_OFFSET		0
#define OLIGODENDROCYTES_OFFSET		8
#define CONVENTIONAL_DC_OFFSET		16
#define ACTIVATED_DC_OFFSET		24
#define T_CYTOTOXIC_OFFSET		32
#define ANTIBODY_OFFSET		40
//Nao vou usar os 2 ultimos!!
#define CELULAS_G_OFFSET		48
#define CELULAS_CA_OFFSET		56
//copiando as vizinhanças (OR-> original, XP -> X Proximo, XM-> X anterior...)
#define CELULAS_POSICAO_OR_OFFSET	0
#define CELULAS_POSICAO_XP_OFFSET	1
#define CELULAS_POSICAO_XM_OFFSET	2
#define CELULAS_POSICAO_YP_OFFSET	3
#define CELULAS_POSICAO_YM_OFFSET	4
#define CELULAS_POSICAO_ZP_OFFSET	5
#define CELULAS_POSICAO_ZM_OFFSET	6
#define CELULAS_NOVO_VALOR_OFFSET	7

//Informacoes de acesso à estrutura "parametrosMalhaGPU".
#define OFFSET_COMPUTACAO		0
#define LENGTH_COMPUTACAO		1
#define MALHA_DIMENSAO_X		2
#define MALHA_DIMENSAO_Y		3
#define MALHA_DIMENSAO_Z		4
#define MALHA_DIMENSAO_POSICAO_Z	5
#define MALHA_DIMENSAO_POSICAO_Y	6
#define MALHA_DIMENSAO_POSICAO_X	7
#define MALHA_DIMENSAO_CELULAS		8
#define NUMERO_PARAMETROS_MALHA		9

#define days 1
__constant float micDiffusion = 0.015206;
__constant float antibodyDiffusion = 0.15206;
__constant float cDcDiffusion = 0.015206;
__constant float aDcDiffusion = 0.015206;
__constant float tCytoDiffusion = 0.015206;
__constant float chi = 0.03;

__constant float muCDc = 0.0432;
__constant float muMic = 0.00432;
__constant float rM = 0.000432;
__constant float rT = 0.001;
__constant float lambAntMic = 0.005702;
__constant float bD = 0.001;

__constant float gammaD = 0.1;
__constant float gammaAntibody = 0.3;
__constant float gammaT = 0.1;

__constant float avgT = 37;
__constant float avgDc = 33;
__constant float avgMic = 350;
__constant float avgOdc = 400;

__constant float cMic = 0.1;
__constant float cCDc = 1;
__constant float cADc = 1;
__constant float cDl = 0.1;
__constant float cF = 0.1;
__constant float alphaTHelper = 0.1;
__constant float alphaTCytotoxic = 0.1;
__constant float alphaB = 0.1;
__constant float alphaP = 1;
__constant float bTHelper = 0.17;
__constant float bTCytotoxic = 0.001;
__constant float bRho = 0.6;
__constant float bRhoB = 3.02;
__constant float bRhoP = 1.02;
__constant float rhoTHelper = 2;
__constant float rhoTCytotoxic = 2;
__constant float rhoB = 11;
__constant float rhoP = 3;
__constant float rhoAntibody = 0.051;
__constant float stableTHelper = 70;
__constant float stableTCytotoxic = 40;
__constant float stableB = 25;
__constant float stableP = 2.5;
__constant float V_LN = 40;
__constant float V_BV = 0;
__constant float V_PV = 0;
__constant float deltaX = 0.1;
__constant float deltaY = 0.1;
__constant float deltaZ = 0.1;

float fFunc(float valuePopulation, float avgPopulation){
    return valuePopulation*valuePopulation/(float)(valuePopulation + avgPopulation);
}

float Laplaciano(int celulaOffset, float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
	return ((xPosicaoGlobal > 0 && xPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET])/(deltaX*deltaX) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)-xPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaX*deltaX)) + (xPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_X]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaX*deltaX)))) + ((yPosicaoGlobal > 0 && yPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET])/(deltaY*deltaY) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)-yPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaY*deltaY)) + (yPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_Y]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaY*deltaY)))) + ((zPosicaoGlobal > 0 && zPosicaoGlobal < (parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) ? (celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] -2 * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] + celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET])/(deltaZ*deltaZ) : ((((parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)-zPosicaoGlobal )/(parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaZ*deltaZ)) + (zPosicaoGlobal /(parametrosMalhaGPU[MALHA_DIMENSAO_Z]-1)) * ((celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET] - celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET])/(deltaZ*deltaZ))));
}

float Quimiotaxia(int celulaOffset, float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
	return (((((xPosicaoGlobal > 0) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XM_OFFSET]) > 0) ? (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_XM_OFFSET] / deltaX) : (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaX)) : 0) + ((xPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_X] - 1) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaX) : ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_XP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_XP_OFFSET] / deltaX)) : 0))/deltaX) + ((((yPosicaoGlobal > 0) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YM_OFFSET]) > 0) ? (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_YM_OFFSET] / deltaY) : (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaY)) : 0) + ((yPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_Y] - 1) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaY) : ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_YP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_YP_OFFSET] / deltaY)) : 0))/deltaY) + ((((zPosicaoGlobal > 0) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) > 0) ? (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_ZM_OFFSET] / deltaZ) : (-(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZM_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaZ)) : 0) + ((zPosicaoGlobal < parametrosMalhaGPU[MALHA_DIMENSAO_Z] - 1) ? (((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) > 0) ? ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_OR_OFFSET] / deltaZ) : ((celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_ZP_OFFSET] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]) * celulas[celulaOffset + CELULAS_POSICAO_ZP_OFFSET] / deltaZ)) : 0))/deltaZ));
}

void CalcularPontos(float *celulas, int xPosicaoGlobal, int yPosicaoGlobal, int zPosicaoGlobal, __constant int *parametrosMalhaGPU)
{
	//**************************************************Novo************************

	float microgliaDiffusion = micDiffusion * Laplaciano(MICROGLIA_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float conventionalDcDiffusion = cDcDiffusion * Laplaciano(CONVENTIONAL_DC_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float tCytotoxicDiffusion = tCytoDiffusion * Laplaciano(T_CYTOTOXIC_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float activatedDCDiffusion = aDcDiffusion * Laplaciano(ACTIVATED_DC_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float antibodyDiffusion = antibodyDiffusion * Laplaciano(ANTIBODY_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	
	float microgliaChemotaxis = chi * Quimiotaxia(MICROGLIA_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float conventionalDcChemotaxis = chi * Quimiotaxia(CONVENTIONAL_DC_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	float tCytotoxicChemotaxis = chi * Quimiotaxia(T_CYTOTOXIC_OFFSET, celulas, xPosicaoGlobal, yPosicaoGlobal, zPosicaoGlobal, parametrosMalhaGPU);
	//Microglia update

	float microgliaReaction = muMic*celulas[MICROGLIA_OFFSET + CELULAS_POSICAO_OR_OFFSET]*(avgMic - celulas[MICROGLIA_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	float microgliaClearance = cMic*celulas[MICROGLIA_OFFSET + CELULAS_POSICAO_OR_OFFSET];

	celulas[MICROGLIA_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = microgliaDiffusion - microgliaChemotaxis + microgliaReaction - microgliaClearance;

	//Conventional DC update
	float conventionalDcReaction = muCDc*celulas[OLIGODENDROCYTES_OFFSET + CELULAS_POSICAO_OR_OFFSET]*(avgDc - celulas[CONVENTIONAL_DC_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	float conventionalDcActivation = bD*celulas[CONVENTIONAL_DC_OFFSET + CELULAS_POSICAO_OR_OFFSET]*celulas[OLIGODENDROCYTES_OFFSET + CELULAS_POSICAO_OR_OFFSET];
	float conventionalDcClearance = cCDc*celulas[CONVENTIONAL_DC_OFFSET + CELULAS_POSICAO_OR_OFFSET];

	celulas[CONVENTIONAL_DC_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = conventionalDcDiffusion - conventionalDcChemotaxis - conventionalDcClearance + conventionalDcReaction - conventionalDcActivation;

	//Activated DC update
	float activatedDcClearance = cADc*celulas[ACTIVATED_DC_OFFSET + CELULAS_POSICAO_OR_OFFSET];
	float activatedDcMigration = 0;//model->thetaPV[kPos]*gammaD*(model->dendriticLymphNode[stepKPlus] - celulas[ACTIVATED_DC_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	
	celulas[ACTIVATED_DC_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = activatedDCDiffusion + conventionalDcActivation + activatedDcMigration - activatedDcClearance;

	//CD8 T update
	float tCytotoxicMigration = 0;//model->thetaBV[kPos]*gammaT*(model->tCytotoxicLymphNode[stepKPlus] - celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	
	celulas[T_CYTOTOXIC_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = tCytotoxicDiffusion - tCytotoxicChemotaxis + tCytotoxicMigration;
	
	//Antibody update
	float odcAntibodyMicrogliaFagocitosis = lambAntMic*celulas[ANTIBODY_OFFSET + CELULAS_POSICAO_OR_OFFSET]*(avgOdc - celulas[OLIGODENDROCYTES_OFFSET + CELULAS_POSICAO_OR_OFFSET])*fFunc(celulas[MICROGLIA_OFFSET + CELULAS_POSICAO_OR_OFFSET], avgMic);
	float antibodyMigration = 0;//model->thetaBV[kPos]*gammaAntibody*(model->antibodyLymphNode[stepKPlus] - celulas[ANTIBODY_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	
	celulas[ANTIBODY_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = antibodyDiffusion + antibodyMigration - odcAntibodyMicrogliaFagocitosis;
	
	//Oligodendrocytes update
	float odcMicrogliaFagocitosis = rM*fFunc(celulas[MICROGLIA_OFFSET + CELULAS_POSICAO_OR_OFFSET], avgMic)*(avgOdc - celulas[OLIGODENDROCYTES_OFFSET + CELULAS_POSICAO_OR_OFFSET]);
	float odcTCytotoxicApoptosis = rT*fFunc(celulas[T_CYTOTOXIC_OFFSET + CELULAS_POSICAO_OR_OFFSET], avgT)*(avgOdc - celulas[OLIGODENDROCYTES_OFFSET + CELULAS_POSICAO_OR_OFFSET]);

	celulas[OLIGODENDROCYTES_OFFSET + CELULAS_NOVO_VALOR_OFFSET] = odcAntibodyMicrogliaFagocitosis + odcMicrogliaFagocitosis + odcTCytotoxicApoptosis;
	
}

__kernel void ProcessarPontos(__global float *malhaPrincipalAtual, __global float *malhaPrincipalAnterior, __constant int *parametrosMalhaGPU)
{
	int globalThreadID = get_global_id(0);

	float celulas[CELULAS_SIZEOF];

	//Descobrir posicao 3D local na malha.
	int posicaoGlobalZ = (globalThreadID/(parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X]));
	int posicaoGlobalY = (globalThreadID%(parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X])) / parametrosMalhaGPU[MALHA_DIMENSAO_X];
	int posicaoGlobalX = (globalThreadID%(parametrosMalhaGPU[MALHA_DIMENSAO_Y] * parametrosMalhaGPU[MALHA_DIMENSAO_X])) % parametrosMalhaGPU[MALHA_DIMENSAO_X];

	if(posicaoGlobalZ >= parametrosMalhaGPU[MALHA_DIMENSAO_Z])
	{
		return;
	}

	//**************************************
	//Preencher celulas para calcular EDO's.
	//**************************************

	int malhaIndex = ((posicaoGlobalZ) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z]) + ((posicaoGlobalY) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y]) + ((posicaoGlobalX) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X]);

	//Loop por todas as celulas.
	for(int count = 0; count < MALHA_TOTAL_CELULAS; count++)
	{
		//Origem.
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_OR_OFFSET] = malhaPrincipalAnterior[malhaIndex + ((MICROGLIA  + count) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS])];

		//Vizinhança.
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_ZP_OFFSET] = ((posicaoGlobalZ + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_Z]))	? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((+1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z])] : 0.0f;
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_ZM_OFFSET] = ((posicaoGlobalZ - 1 >= 0))						? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((-1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Z])] : 0.0f;
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_YP_OFFSET] = ((posicaoGlobalY + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_Y]))	? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((+1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y])] : 0.0f;
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_YM_OFFSET] = ((posicaoGlobalY - 1 >= 0))						? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((-1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_Y])] : 0.0f;
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_XP_OFFSET] = ((posicaoGlobalX + 1 < parametrosMalhaGPU[MALHA_DIMENSAO_X]))	? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((+1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X])] : 0.0f;
		celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_POSICAO_XM_OFFSET] = ((posicaoGlobalX - 1 >= 0))						? malhaPrincipalAnterior[malhaIndex + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS]) + ((-1) * parametrosMalhaGPU[MALHA_DIMENSAO_POSICAO_X])] : 0.0f;
	}

	//**************************************
	//Atualizar malha com pontos calculados.
	//**************************************
	
	CalcularPontos(celulas, posicaoGlobalX, posicaoGlobalY, posicaoGlobalZ, parametrosMalhaGPU);

	//Loop por todas as celulas.
	for(int count = 0; count < MALHA_TOTAL_CELULAS; count++)
	{
		malhaPrincipalAtual[malhaIndex  + (((MICROGLIA  + count)) * parametrosMalhaGPU[MALHA_DIMENSAO_CELULAS])] = celulas[((MICROGLIA  + count) * CELULAS_SINGLE_CELL_SIZEOF) + CELULAS_NOVO_VALOR_OFFSET];
	}
}

