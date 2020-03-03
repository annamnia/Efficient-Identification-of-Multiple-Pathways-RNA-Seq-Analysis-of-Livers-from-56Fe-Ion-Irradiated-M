#define _XOPEN_SOURCE

/* Structure definition*/
typedef struct{
	double modularity;
	int *community;
}comstruct;

/* Prototypes */
comstruct uni(double **a, int n);
