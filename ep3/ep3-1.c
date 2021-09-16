#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Pn (double x, double xi[7] ,double a[7]);
double trapezio (double x[7],double a[7],int m);
double Simpson (double x[7],double a[7],int m);
double Romberg (double x[7],double a[7],int m);
double * Newton(double x[7],double y[7]);

#define MAX 6
#define MIN 0

int main()
{
    // numero de iterações que será rodado na potencia de 2
    // roda 2^m
    int m = 25;

    double areaTrap, areaSimp, areaRomb;
    double * a ; // guarda as constantes do polinomiode newton
    double x[7] = {0, 5,      10,     15,     20,     25,     30}; // intervalos dados de x
    double y[7] = {0, 1.5297, 9.5120, 8.7025, 2.8087, 1.0881, 0.3537}; //valores das forças da função dada na tabela

    // guarda as constantes do polinomio de newton
    a = Newton(x,y);

    areaTrap = trapezio(x,a,m);
    printf("TRAPEZIO %.20lf\n",areaTrap);
    areaSimp = Simpson(x,a,m);
    printf("SIMPSON  %.20lf\n",areaSimp);
    areaRomb = Romberg(x,a,m);
    printf("ROMBERG  %.20lf\n",areaRomb);
    
    free(a);
}

// calcula o valor no ponto dado com o polinomi de newton
double Pn (double x, double xi[7] ,double * a)
{
    double res = a[MAX];
    
    int i = MAX-1;
    for (int i = MAX-1; i >= 0;i--)
    {
        // por manipulação algebrica deixa esse jeito mais eficiente por
        // não fazer a mesma conta varias vezes
        res = res*(x-xi[i]) + a[i];
    }
    return res;
}

// calcula a area do trapezio no intervalo [x1,x2] sendo x1 < x2
double trapezio (double x[7],double * a,int m)
{
    double h = x[MAX] - x[MIN];
    double area, n;

    // T(f;hm)
    double T = h*(Pn(x[MIN],x,a) + Pn(x[MAX],x,a))/2;
    for (int i = 1; i <= m; i++)
    {
        // calcula o novo h
        h = h/2;

        area = 0;
        n = pow(2,i-1)-1;
        for (int j = 0; j <= n; j++)
        {
            area = area + Pn(x[MIN]+(2*j+1)*h, x, a);
        }
        T = T/2 + h*area;
    }
    return T;
}

double Simpson (double x[7],double * a,int m)
{
    double h = (x[MAX] - x[MIN])/pow(2,m+1);
    double T = 0;

    double somaImpar;
    double somaPar;

    somaImpar = 0;

    int n = pow(2,m);

    // vai de x1 até xn
    // calcula a soma parcial dos pontos impares
    for (int i = 1; i <= n; i++)
    {
        somaImpar = somaImpar + Pn((2*i+1)*h, x, a);
    }

    somaPar = 0;
    // calcula a soma parcial dos pontos pares
    for (int i = 1; i <= n-1; i++)
    {
        somaPar = somaPar + Pn(2*i*h, x, a);
    }

    // LEI DE SIMPOSON COMPOSTA
    T = (h/3)*(Pn(x[MIN], x, a) + Pn(x[MAX], x, a) + 4*somaImpar + 2*somaPar);
}

double Romberg (double x[7],double *a,int m)
{
    double h = x[MAX]-x[MIN];
    double **R;
    double aux, area;
    int n;

    R = malloc((m + 1) * sizeof(double));
    for (int i = 0; i < (m + 1); i++)
        R[i] = malloc((m + 1)* sizeof(double));;
    
    R[0][0] = h*(Pn(x[MIN], x, a) + Pn(x[MAX], x, a))/2;
    // constroi a coluna de R para rodar o método
    for (int i = 1; i <= m; i++)
    {
        // calcula o novo h
        h = h/2;

        area = 0;
        n = pow(2,i-1)-1;
        for (int j = 0; j <= n; j++)
        {
            area = area + Pn(x[MIN]+(2*j+1)*h, x, a);
        }
        R[i][0] = R[i-1][0]/2 + h*area;
    }

    // constrói a matriz para achar os ai's
    // começa da segunda coluna (i) já que a primeira já esta feita
    // as linhas (j) começa em m e vai até i sengdo i o numero da coluna
    for (int i = 1; i <= m+1 ;i++)
    {
        for (int j = m; j >= i; j--)
        {
            R[j][i] = R[j][i-1] + (R[j][i-1] - R[j-1][i-1])/(pow(2,2*i)-1);
        }
    }

    aux = R[m][m];

    for (int i = 0; i < m+1; i++)
    {
        free(R[i]);
    }
    free(R);
    
    return aux;
}

// encontra as constrantes ai do polinomio de newton pela interpolação de newton
double * Newton(double x[7],double y[7])
{
    int m = MAX;
    double A[7][7];
    double * a = malloc(7*sizeof(double));
    double h = x[1]-x[0];

    for (int i = 0; i < 7; i++)
    {
        A[i][0] = y[i];
    }

    // constrói a matriz para achar os ai's
    for (int i = 0; i < MAX; i++)
    {
        for (int j = 0; j < m-i; j++)
        {
            A[j][i+1] = (A[j+1][i] - A[j][i])/(x[j+i+1]-x[j]);
        }
    }

    for (int i = 0; i < MAX+1; i++)
    {
        a[i] = A[0][i];
    }

    return a;
}