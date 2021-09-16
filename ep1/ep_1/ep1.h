#ifndef _EP1_H
#define _EP1_H

struct complexo
{
	double parteReal;
	double parteImaginaria;
};
typedef struct complexo Complexo;


/* PARTE 1 */
double modulo (double a);
double pontoFixo(double x_0, int nMax, int i);

/* PARTE 2 */

int procuraResposta(Complexo *respostas, Complexo res, int numResp);

/* Funções de operadores de complexos */
double moduloComplexo(Complexo complexo);
Complexo soma (Complexo num, double valor);
Complexo adiconaParteImaginaria(Complexo num, double valor);
Complexo multiplica(Complexo complexo, double constante);
Complexo somaComplexos (Complexo complexo_1, Complexo complexo_2);
Complexo subtraiComplexos (Complexo complexo_1, Complexo complexo_2);
Complexo multiplicaComplexos (Complexo complexo_1, Complexo complexo_2);
Complexo divideComplexos (Complexo complexo_1, Complexo complexo_2);
Complexo powComplexo(Complexo complexo, int n);


/* Funções obrigatórias */
Complexo newton (Complexo x_0, int nMax);
Complexo evalf(Complexo x);
Complexo evalDf(Complexo x);
void newton_basins (double l, double u, int p);

#endif