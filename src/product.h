#ifndef PRODUCT_H
#define PRODUCT_H

#include "array.h"

typedef struct{
    int counter;
    array2* pointer;
    array2* A;
    int END;
}product;

product* product_init(int a, int repeat);
array2* product_next(product* p);
int product_has_next(product* p);
void product_finish(product* p);

#endif
