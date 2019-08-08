#ifndef ARRAY_H
#define ARRAY_H

typedef struct {
    int* array;
    int size;
}array;

array* array_arange(int n);
array* array_zeros(int n);
array* array_copy(array* a);
int array_equals(array* a, array* b);
void array_destroy(array* a);
array** array_matrix(int nrow);
void array_reverse(array* a);
void array_matrix_destroy(array** a, int nrow);
array* array_erase(array* a, int e);
void array_print(array* a);
int array_contains(array* a, int e);
array* array_remove(array* a, array* b);
array* array_sub(array* a, int i);

#endif
