#ifndef PRODUCT_H
#define PRODUCT_H

#include "array.h"

typedef struct{
    int counter;
    array* pointer;
    array* A;
    int END;
}product;

product* product_init(array* a, int repeat);
array* product_next(product* p);
int product_has_next(product* p);
void product_finish(product* p);

#endif
