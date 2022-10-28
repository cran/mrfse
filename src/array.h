#ifndef ARRAY_H
#define ARRAY_H

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

typedef struct {
    int* array2;
    int size;
}array2;

array2* array2_arange(int n);
array2* array2_zeros(int n);
array2* array2_copy(array2* a);
int array2_equals(array2* a, array2* b);
void array2_destroy(array2* a);
array2** array2_matrix(int nrow);
void array2_reverse(array2* a);
void array2_matrix_destroy(array2** a, int nrow);
array2* array2_erase(array2* a, int e);
void array2_print(array2* a);
int array2_contains(array2* a, int e);
array2* array2_remove(array2* a, array2* b);
array2* array2_sub(array2* a, int i);
array2* array2_random_x(int x, int n);
vector<int> array_to_vec(array2 *a);

#endif
