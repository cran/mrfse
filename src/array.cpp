#include <Rcpp.h>
#include "array.h"
#include "util.h"

using namespace Rcpp;
using namespace std;

array2* malloc_array2(int n);
array2** malloc_matrix(int n);

array2* array2_arange(int n) {
    array2* a = malloc_array2(1);
    a->array2 = malloc_int(n);
    a->size = n;
    for (int i = 0; i < n; i++) a->array2[i] = i;
    return a;    
}

array2* array2_zeros(int n) {
    array2* a = malloc_array2(1);
    a->array2 = malloc_int(n);
    a->size = n;
    for (int i = 0; i < n; i++) a->array2[i] = 0;
    return a;
}

array2* array2_copy(array2* a) {
    array2* b = malloc_array2(1);
    b->array2 = malloc_int(a->size);
    b->size = a->size;
    for (int i = 0; i < b->size; i++) b->array2[i] = a->array2[i];
    return b;
}

int array2_equals(array2 *a, array2 *b) {
    if (a->size != b->size) return 0;
    for (int i = 0; i < a->size; i++)
        if (a->array2[i] != b->array2[i]) return 0;
    return 1;    
}

array2** array2_matrix(int nrow) {
    array2** a = malloc_matrix(nrow);
    return a;
}

array2* array2_erase(array2* a, int e) {
    int index = -1;
    for (int i = 0; i < a->size; i++) {
        if (a->array2[i] == e) {
            index = i;
            break;
        }
    }
    if (index != -1) {
        array2 *b = array2_zeros(a->size -1);
        for (int i = 0, j = 0; i < a->size; i++) {
            if (i != index) b->array2[j++] = a->array2[i];                
        }
        return b;
    } else
        return a;
}

void array2_reverse(array2* a) {
    int mid = a->size / 2;
    for (int i = 0; i < mid; i++) {
        int tmp = a->array2[i];
        a->array2[i] = a->array2[a->size - i - 1];
        a->array2[a->size - i - 1] = tmp;
    }
}

void array2_destroy(array2* a) {
    free(a->array2);
    a->array2 = NULL;
    free(a);
    a = NULL;
}

void array2_matrix_destroy(array2** a, int nrow) {
    for (int i = 0; i < nrow; i++) {
        array2_destroy(a[i]);            
    }
    free(a);
    a = NULL;
}

int array2_contains(array2 *a, int e) {
    for (int i = 0; i < a->size; i++)
        if (a->array2[i] == e) return 1;

    return 0;
}

array2* array2_remove(array2* a, array2* b) {
    int n = a->size;
    for (int i = 0; i < b->size; i++) {
        if (array2_contains(a, b->array2[i])) n--;
    }
    array2 *r = array2_zeros(n);
    for (int i = 0, j = 0; i < a->size; i++)
        if (!array2_contains(b, a->array2[i])) r->array2[j++] = a->array2[i];

    return r;
}

array2* array2_sub(array2* a, int i) {
    array2* r = array2_zeros(i + 1);
    for (int j = 0; j <= i; j++)
        r->array2[j] = a->array2[j];

    return r;    
}

void array2_print(array2* a) {
    /* printf("size = %d\n", a->size); */
    for (int i = 0; i < a->size; i++) {
	// printf("%d ",  a->array2[i]);
    }
}

 array2* malloc_array2(int n) {
    array2* a = (array2*) malloc(n * sizeof(array2));
    if (a == NULL) {
        // error("malloc returned NULL!\n");
    }
    return a;
}

array2** malloc_matrix(int n) {
    array2** v = (array2**) malloc(n * sizeof(array2*));
    if (v == NULL) {
        // error("malloc returned NULL!\n");
    }
    return v;
}

array2* array2_random_x(int x, int n) {
    array2* a = array2_zeros(n);
    for (int i = 0; i < n; i++) {
	a->array2[i] = int_unif(x) % x;
    }
    return a;
}

vector<int> array_to_vec(array2* a) {
    int n = a->size;
    vector<int> result(n);
    for (int i = 0; i < n; i++) {
	result[i] = a->array2[i];
    }
    array2_destroy(a);
    return result;
}
