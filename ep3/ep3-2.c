#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FUNC 0 // para os valores 0 calcula sin(x), 1 calcula x³, 2 calcula exp(-x) e para 3 calcula o valor de pi

double f(double x);
double aleatorio(double a,double b);
double MonteCarlo(double a,double b,int n,int dim);

int main()
{
    // numero de pontos que sera calculado
    int n = 1000000;

    // numero de dimensões da integral
    int dim = 2;
    
    // intervalo que será calculado a integral [a,b]
    double a = 0;
    double b = 1;

    srand((unsigned)time(NULL));
    
    printf("volume calculado: %lf\n", MonteCarlo(a,b,n,dim));
}

double f(double x)
{
    if (FUNC == 0)
        return  sin(x);
    else if (FUNC == 1)
        return pow(4*x+3,3)*4;
    else if (FUNC == 2)
        return exp(1-(1/x))/pow(x,2);
    else if (FUNC == 3)
        return pow(x,2);
}

double aleatorio(double a, double b)
{
    return (double)rand()/(double)(RAND_MAX)*(b-a) + a;
}

double MonteCarlo(double a,double b,int n,int dim)
{
    // valor da integral
    double I = 0;
    if (FUNC < 3)
    {
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < n; j++)
            {
                I += f(aleatorio(a,b));
            }
            I /= n;
        }
    }
    else // estima um quarto da area do circulo
    {
        for (int i = 0; i < n; i++)
        {
            double soma = 0;
            for (int j = 0; j < dim; j++)
            {
                soma += f(aleatorio(a,b));
            }

            if (soma <= 1)
            {
                I++;
            }
        }
        I = 4*I/n;
    }

    return I;
}