#ifndef COMBINATION_H
#define COMBINATION_H

#include "array.h"

typedef struct{
    int n;
    int k;
    int END;
    int *comb;
    array* data;
}combination;

combination* combination_init(array* a, int size);
array* combination_next(combination* c);
int combination_has_next(combination* c);
void combination_finish(combination* c);
#endif










