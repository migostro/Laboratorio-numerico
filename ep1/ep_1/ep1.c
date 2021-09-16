#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <unistd.h>

#include "ep1.h"

#define epsilon 0.0000000000000001
#define linha 2
#define coluna 2
#define pixels 2000

// func = 0 (x⁴ - 1)
// func = 1 (x³ - 1)
// func = 2 (x⁶ - 1)
// func = 3 (x⁵ - 1)
// func = 4 (x² + 1)
// func = 5 (x⁷ - 1)
// func = 6 (x⁸ - 1)
#define func 6

int main()
{
	/*********************************************************** PARTE 1 **************************************************************************/

	double x_0[3], aux;

	// pontos que foram escolhidos para cada função
	// foram escolhida por serem razoavelmente proximos de cada raiz
	x_0[0] = 0;
	x_0[1] = 1;
	x_0[2] = 2.5;	

	printf("PARTE 1\n");

	for (int i = 0; i < 3; ++i)
	{
		aux = pontoFixo(x_0[i], 100, i);
		printf("Raíz %d: %lf \n", i+1, aux);
	}

	/* PARTE 2 */

	printf("PARTE 2\n");



	newton_basins(linha, coluna, pixels);

	printf("acabou o programa\n");


	return 0;
}

/* PARTE 1 */

double modulo (double a)
{
	if (a < 0)
	{
		return -a;
	}
	return a;
}

double pontoFixo(double x_0, int nMax, int i)
{
	int cont = 0;
	double x_anterior = 1000000, x_atual = x_0;

	while (modulo(x_atual - x_anterior) > epsilon && cont < nMax)
	{
		// g(x) = -((e^x)/2)^(1/2)
		if (i == 0)
		{
			x_anterior = x_atual;

			x_atual =  -sqrt(exp(x_atual)/2);
		}
		// g(x) = ((e^x)/2)^(1/2)
		else if (i == 1)
		{
			x_anterior = x_atual;

			x_atual =  sqrt(exp(x_atual)/2);
		}
		// g(x) = ln(2) + 2ln(x)
		else if (i == 2)
		{
			x_anterior = x_atual;

			x_atual = log(2) + 2*log(x_atual);
		}
		
		cont++;
	}

	return x_atual;
	
}

/*********************************************************** PARTE 2 **************************************************************************/
Complexo newton(Complexo x_0, int nMax)
{
	int cont = 0;
	Complexo x_atual, dx;

	x_atual = x_0;

	dx = divideComplexos(evalf(x_atual), evalDf(x_atual));
	//
	//
	while ( moduloComplexo(dx) > epsilon && moduloComplexo(evalf(x_atual)) > epsilon && cont < nMax)
	{
		//x_anterior = x_atual;
		// aplica o metodo de newton
		// x_k+1 = x_k - f(x_k)/f'(x_k)
		x_atual = subtraiComplexos(x_atual, dx);

		dx = divideComplexos(evalf(x_atual), evalDf(x_atual));

		cont++;
	}

	//esse é o possivel root
	return x_atual;
}

Complexo evalf(Complexo x)
{
	// x⁴ - 1
	if (func == 0) return (soma(powComplexo(x, 4), -1));
	// x³ - 1
	else if (func == 1) return (soma(powComplexo(x, 3), -1));
	// x⁶ - 1
	else if (func == 2) return (soma(powComplexo(x,6), -1));
	// x⁵ - 1
	else if (func == 3)	return (soma(powComplexo(x,5), -1));
	// x² + 1
	else if (func == 4)	return (soma(powComplexo(x,2), 1));
	// x⁷ - 1
	else if (func == 5)	return (soma(powComplexo(x,7), -1));
	// x⁸ - 1
	else if (func == 6)	return (soma(powComplexo(x,8), -1));
}

Complexo evalDf(Complexo x)
{
	// x⁴ - 1
	if (func == 0) return (multiplica(powComplexo(x, 3), 4));
	// x³ - 1
	else if (func == 1) return (multiplica(powComplexo(x, 2), 3));
	// x⁶ - 1
	else if (func == 2) return (multiplica(powComplexo(x, 5), 6));
	// x⁵ - 1
	else if (func == 3)	return (multiplica(powComplexo(x, 4), 5));
	// x² + 1
	else if (func == 4)	return (multiplica(powComplexo(x, 1), 2));
	// x⁷ - 1
	else if (func == 5)	return (multiplica(powComplexo(x, 6), 7));
	// x⁸ - 1
	else if (func == 6)	return (multiplica(powComplexo(x, 7), 8));
}

void newton_basins(double l, double u, int p)
{
	int achou, **matriz, numResp = 0, maxRaizes = 2000;
	Complexo x_0, res, *raizes, novo;
	double eps_x, eps_y;
	FILE *pont_arq;

	matriz = malloc(p*sizeof(int*));
	raizes = malloc(maxRaizes*sizeof(Complexo));
	
	for (int i = 0; i < p; ++i)
	{
		matriz[i] = malloc(p*sizeof(int));
	}

	if(pont_arq == NULL)
	{
		printf("Erro na abertura do arquivo!");
	 	return 0;
	}

	eps_x = (2*u)/p;
	eps_y = (2*l)/p;

	x_0.parteImaginaria = l;
	for (int i = 0; i < p; i++)
	{

		x_0.parteReal = u;
		for (int j = 0; j < p; j++)
		{
			// o segundo argumento é o maximo de vezes que o metodo pode rodar
			res = newton(x_0, 60);

			// verifica se não divergiu
			// x⁵ - 1 ... 3
			// x⁶ - 1 ... 5.5
			// para a penultima função 5.565
			// 6 .. 6.75
			if (moduloComplexo(evalf(res)) < epsilon*6.75)
			{

				achou = 0;
				// verifica se essa raiz já foi encontrada
				// Se não foi encontrada ela é adicionada e seta a posição do vetor que ela pertence na matriz
				// Se ela já existe guarda na posição da matriz a posição relativa da raiz
				for (int k = 0; k < numResp && achou == 0; ++k)
				{
					// verifica uma tolerância nos resultados
					if (comparaComplexo(raizes[k], res))
					{
						// adiciona a "cor" da raiz encontrada
						matriz[i][j] = k;
						achou = 1;
					}
				}

				// adiciona nova raíz caso não exista ainda
				if (!achou)
				{
					// adiciona a nova raiz na coleção de raizes descobertas
					raizes[numResp] = res;

					// adiciona a cor da nova raíz encontrada
					matriz[i][j] = numResp;
					
					numResp++;


				}
			}
			// caso divergiu
			else
			{
				matriz[i][j] = -1;
			}

			x_0.parteReal -= eps_x;
		}
		x_0.parteImaginaria -= eps_y;
	}

	//abrindo o arquivo com tipo de abertura w
	pont_arq = fopen("matriz_simulacao_6.txt", "w");
	// guarda a matriz no txt para poder colocar no grafico
	for (int i = 0; i < p; ++i)// linha
	{
		for (int j = 0; j < p; ++j)// coluna
		{
			fprintf(pont_arq, "%d ", matriz[j][i]);
		}
		fprintf(pont_arq, "\n");
	}
	close(pont_arq);

	printf("fim\n");

	for (int i = 0; i < p; ++i)
	{
		free(matriz[i]);
	}
	free(matriz);

	free(raizes);
}

int procuraResposta(Complexo *respostas, Complexo res, int numResp)
{
	for (int i = 0; i < numResp; ++i)
	{
		if (comparaComplexo(respostas[i], res))
		{
			return i;
		}
	}
	return -1;
}

int comparaComplexo(Complexo complexo_1, Complexo complexo_2)
{
	double a, b, c, d;

	a = complexo_1.parteReal;
	b = complexo_1.parteImaginaria;
	c = complexo_2.parteReal;
	d = complexo_2.parteImaginaria;

	if (modulo(a - c) < 0.00001 && modulo(b - d) < 0.00001)
	{
		return 1;
	}

	return 0;
}

/* FUNÇÕES ULTILIZADAS PARA SIMULAR O PLANO IMAGINARIO */

double moduloComplexo(Complexo complexo)
{
	double a, b;

	a = complexo.parteReal;
	b = complexo.parteImaginaria;

	return sqrt(pow(a, 2) + pow(b, 2));
}

// Soma a parte real a uma constante
Complexo soma(Complexo num, double valor)
{
	Complexo novo;

	novo.parteReal = num.parteReal + valor;
	novo.parteImaginaria = num.parteImaginaria;

	return novo;
}

// Soma a parte imaginaria a uma constante
Complexo adiconaParteImaginaria(Complexo num, double valor)
{
	Complexo novo;

	novo = num;

	novo.parteReal = num.parteReal;
	novo.parteImaginaria += valor;

	return novo;
}

// Multiplica uma constante real aos complexos
// devolve a soma
Complexo multiplica(Complexo complexo, double constante)
{
	Complexo novo;

	novo.parteReal = complexo.parteReal * constante;
	novo.parteImaginaria = complexo.parteImaginaria * constante;

	return novo;
}


// soma dois complexos e devolve sua soma
// soma na forma (a + c) + (b + d)i
Complexo somaComplexos (Complexo complexo_1, Complexo complexo_2)
{
	Complexo novo;
	double a, b, c, d;

	a = complexo_1.parteReal;
	b = complexo_1.parteImaginaria;

	c = complexo_2.parteReal;
	d = complexo_2.parteImaginaria;

	

	novo.parteReal = a + c;//soma(complexo_1, complexo_2.parteReal)->parteReal;
	novo.parteImaginaria = b + d;//adiconaParteImaginaria(complexo_1, complexo_2.parteImaginaria)->parteImaginaria;

	return novo;
}

// subtrai dois complexos e retorna o resultado
// faz complexo_1 - complexo_2
Complexo subtraiComplexos (Complexo complexo_1, Complexo complexo_2)
{
	Complexo novo;
	double a, b, c, d;


	a = complexo_1.parteReal;
	b = complexo_1.parteImaginaria;

	c = complexo_2.parteReal;
	d = complexo_2.parteImaginaria;

	novo.parteReal = a - c;//soma(complexo_1, -complexo_2.parteReal)->parteReal;
	novo.parteImaginaria = b - d;//adiconaParteImaginaria(complexo_1, -complexo_2.parteImaginaria)->parteImaginaria;

	return novo;
}

// Multiplica os dois complexos respectivamente na forma a + bi e c + di
// a multiplicação deixa na forma (ac - bd) + (ad + bc)i
// onde i é o numero imaginario
Complexo multiplicaComplexos (Complexo complexo_1, Complexo complexo_2)
{
	Complexo novo;

	double a, b, c, d;

	a = complexo_1.parteReal;
	b = complexo_1.parteImaginaria;

	c = complexo_2.parteReal;
	d = complexo_2.parteImaginaria;

	novo.parteReal = a*c - b*d;
	novo.parteImaginaria = a*d + b*c;

	return novo;
}

// divide os complexos e devolve o resultado
// divide na forma ( (ac + bd) + (bc - ad)i )/ (c² + d²)
// onde i é a parte imaginaria
// faz complexo_1/complexo_2
Complexo divideComplexos (Complexo complexo_1, Complexo complexo_2)
{
	Complexo novo;

	double a, b, c, d;

	a = complexo_1.parteReal;
	b = complexo_1.parteImaginaria;
	
	c = complexo_2.parteReal;
	d = complexo_2.parteImaginaria;

	novo.parteReal = (a*c + b*d) / (c*c + d*d);
	novo.parteImaginaria = (b*c - a*d) / (c*c + d*d);

	return novo;
}

Complexo powComplexo(Complexo complexo, int n)
{
	Complexo novo, aux;

	if (n == 0)
	{
		novo.parteReal = 1;
		novo.parteImaginaria = 0;

		return novo;
	}

	novo = complexo;

	for (int i = 1; i < n; ++i)
	{
		novo = multiplicaComplexos(novo, complexo);
	}

	return novo;

}

