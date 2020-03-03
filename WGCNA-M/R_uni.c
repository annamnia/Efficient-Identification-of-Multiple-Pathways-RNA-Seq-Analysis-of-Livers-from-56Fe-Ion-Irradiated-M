#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "R_uni.h"

//this function takes data passed from R and applies the Modularity algorithm
void uni_ComDet(int *size, double *net, int *iter, double *MODULARITY, int *PARTITION)
{

	int i,j, n = *size, run = *iter;
	double **A;
	comstruct result;

	A = malloc(n * sizeof(double *));
	A[0] = malloc(n * n * sizeof(double));

	for (i = 1; i < n; i++) A[i] = A[i-1] + n;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = net[i * n + j];
			//if (A[i][j] > 1) A[i][j] = 1;
		}
		A[i][i] = 0;
	}

	double Qmax = 0;
	for (i = 0; i < run; i++) {
		result = uni(A, n);
		if (result.modularity > Qmax) {
			Qmax = result.modularity;
			for (j = 0; j < n; j++)
				PARTITION[j] = result.community[j];
		}
		//free(result.community);
		//printf("current Q is %g\n", result.modularity);
	}
	MODULARITY[0] = Qmax;

	free(A[0]);
	free(A);
	return;
}
