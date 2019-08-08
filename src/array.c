#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include "array.h"
#include "util.h"

static array* malloc_array(int n);
static array** malloc_matrix(int n);

array* array_arange(int n) {
    array* a = malloc_array(1);
    a->array = malloc_int(n);
    a->size = n;
    for (int i = 0; i < n; i++) a->array[i] = i;
    return a;    
}

array* array_zeros(int n) {
    array* a = malloc_array(1);
    a->array = malloc_int(n);
    a->size = n;
    for (int i = 0; i < n; i++) a->array[i] = 0;
    return a;
}

array* array_copy(array* a) {
    array* b = malloc_array(1);
    b->array = malloc_int(a->size);
    b->size = a->size;
    for (int i = 0; i < b->size; i++) b->array[i] = a->array[i];
    return b;
}

int array_equals(array *a, array *b) {
    if (a->size != b->size) return 0;
    for (int i = 0; i < a->size; i++)
        if (a->array[i] != b->array[i]) return 0;
    return 1;    
}

array** array_matrix(int nrow) {
    array** a = malloc_matrix(nrow);
    return a;
}

array* array_erase(array* a, int e) {
    int index = -1;
    for (int i = 0; i < a->size; i++) {
        if (a->array[i] == e) {
            index = i;
            break;
        }
    }
    if (index != -1) {
        array *b = array_zeros(a->size -1);
        for (int i = 0, j = 0; i < a->size; i++) {
            if (i != index) b->array[j++] = a->array[i];                
        }
        return b;
    } else
        return a;
}

void array_reverse(array* a) {
    int mid = a->size / 2;
    for (int i = 0; i < mid; i++) {
        int tmp = a->array[i];
        a->array[i] = a->array[a->size - i - 1];
        a->array[a->size - i - 1] = tmp;
    }
}

void array_destroy(array* a) {
    free(a->array);
    a->array = NULL;
    free(a);
    a = NULL;
}

void array_matrix_destroy(array** a, int nrow) {
    for (int i = 0; i < nrow; i++) {
        array_destroy(a[i]);            
    }
    free(a);
    a = NULL;
}

int array_contains(array *a, int e) {
    for (int i = 0; i < a->size; i++)
        if (a->array[i] == e) return 1;

    return 0;
}

array* array_remove(array* a, array* b) {
    int n = a->size;
    for (int i = 0; i < b->size; i++) {
        if (array_contains(a, b->array[i])) n--;
    }
    array *r = array_zeros(n);
    for (int i = 0, j = 0; i < a->size; i++)
        if (!array_contains(b, a->array[i])) r->array[j++] = a->array[i];

    return r;
}

array* array_sub(array* a, int i) {
    array* r = array_zeros(i + 1);
    for (int j = 0; j <= i; j++)
        r->array[j] = a->array[j];

    return r;    
}

static array* malloc_array(int n) {
    array* a = (array*) malloc(n * sizeof(array));
    if (a == NULL) {
        error("malloc returned NULL!\n");
    }
    return a;
}

static array** malloc_matrix(int n) {
    array** v = (array**) malloc(n * sizeof(array*));
    if (v == NULL) {
        error("malloc returned NULL!\n");
    }
    return v;
}
