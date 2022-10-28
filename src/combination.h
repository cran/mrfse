#ifndef COMBINATION_H
#define COMBINATION_H

#include "array.h"

typedef struct{
    int n;
    int k;
    int END;
    int *comb;
    array2* data;
}combination;

combination* combination_init(array2* a, int size);
array2* combination_next(combination* c);
int combination_has_next(combination* c);
void combination_finish(combination* c);
#endif










