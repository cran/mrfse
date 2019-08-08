#ifndef UTIL_H
#define UTIL_H

int* malloc_int(int n);
int **matrixINT(int n, int m);
void matrixINTcpy(int **A, int **B, int n, int m);
void free_matrixINT(int **A, int n);

extern const double INF;

#endif
