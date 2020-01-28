#define _XOPEN_SOURCE



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define eps 1E-6
#define alpha 1E-4
#define beta 1E-3
#define toler 1E-8

/* Structure definition*/
typedef struct{
	double modularity;
	int *community;
}comstruct;

double smprng(void);
comstruct uni(double **a, int n);
double bisection(double **B, int *g, int *s, int n);
double eigen(double **B, double *v, int n);
double abmax(double *v,int n);
double KLtuning(double **B, int *s, int n);
double Finaltuning(double **B, int **S, int n, int *ng);
double Agglomeration(double **B, int n, int *ng, int *gs, int **G);

void trans(int **S, int ng, int *gs, int **G);
void invertrans(int **S, int n, int *ng, int *gs, int **G);
double Modularity(int n,int *s,double **b);
double Addmatrix(double **b,int n);
double dQA(int n1, int n2, int *g1, int **G, double **b);
void con_to_G(int *con, int n, int *ng, int *gs, int **G);

static double m;
static double Q;
static int **G;		  

			  
