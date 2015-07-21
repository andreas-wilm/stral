#include <stdio.h>

void ludcmp(float **a, int n, int *indx, float *d);
void lubksb(float **a, int n, int *indx, float b[]);
float **matinverse(float **a,int n);
float **matmult(float **a,float **b,float **c,int i, int k, int j);
