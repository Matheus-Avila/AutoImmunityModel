#ifndef BALANCEADOR_H
#define BALANCEADOR_H

void ComputarCargas(const long int *tempos, const float *cargasAntigas, float *cargas, int participantes);
bool ComputarIntersecao(int offset1, int length1, int offset2, int length2, int *intersecaoOffset, int *intersecaoLength);
int RecuperarPosicaoHistograma(int *histograma, int tamanho, int indice);
float ComputarDesvioPadraoPercentual(const long int *tempos, int participantes);
float ComputarNorma(const float *cargasAntigas, const float *cargaNovas, int participantes);

#endif
