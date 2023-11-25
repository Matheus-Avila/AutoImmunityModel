#include "Balanceador.h"

#include <stdio.h>
#include <math.h>

void ComputarCargas(const long int *tempos, const float *cargasAntigas, float *cargas, int participantes)
{
	if(participantes == 1)
        {
                cargas[0] = 1.0f;
                return;
        }

        float cargaTotal = 0.0f;
        for(int count = 0; count < participantes; count++)
        {
                cargaTotal += ((count == 0) ? (cargasAntigas[count]-0.0f) : (cargasAntigas[count]-cargasAntigas[count-1])) * ((count == 0) ? 1.0f : ((float)tempos[0])/((float)tempos[count]));
        }

        for(int count = 0; count < participantes; count++)
        {
                float cargaNova = (((count == 0) ? (cargasAntigas[count]-0.0f) : (cargasAntigas[count]-cargasAntigas[count-1])) * ((count == 0) ? 1.0f : ((float)tempos[0])/((float)tempos[count]))) / cargaTotal;
                cargas[count] = ((count == 0) ? cargaNova : cargas[count-1]+cargaNova);
        }
}

//Offset1 deve ser sempre menor ou igual que offset2.
bool ComputarIntersecao(int offset1, int length1, int offset2, int length2, int *intersecaoOffset, int *intersecaoLength)
{
	if(offset1+length1 <= offset2)
	{
		return false;
	}

	if(offset1+length1 > offset2+length2)
	{
		*intersecaoOffset = offset2;
		*intersecaoLength = length2;
	}
	else
	{
		*intersecaoOffset = offset2;
		*intersecaoLength = (offset1+length1)-offset2;
	}
	return true;
}

int RecuperarPosicaoHistograma(int *histograma, int tamanho, int indice)
{
	int offset = 0;
	for(int count = 0; count < tamanho; count++)
	{
		if(indice >= offset && indice < offset+histograma[count])
		{
			return count;
		}
		offset += histograma[count];
	}
	return -1;
}

float ComputarDesvioPadraoPercentual(const long int *tempos, int participantes)
{
	double media = 0.0;
	for(int count = 0; count < participantes; count++)
	{
		media += (double)tempos[count];
	}
	media /= (double)participantes;
	
	double variancia = 0.0;
	for(int count = 0; count < participantes; count++)
	{
		variancia += ((double)tempos[count]-media)*((double)tempos[count]-media);
	}
	variancia /= (double)participantes;
	return sqrt(variancia)/media;
}

float ComputarNorma(const float *cargasAntigas, const float *cargasNovas, int participantes)
{
	float retorno = 0.0;
	for(int count = 0; count < participantes; count++)
	{
		retorno += (cargasAntigas[count]-cargasNovas[count])*(cargasAntigas[count]-cargasNovas[count]);
	}
	return sqrt(retorno);
}

