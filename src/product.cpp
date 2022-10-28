#include <stdio.h>
#include <R.h>
#include <stdlib.h>
#include "array.h"
#include "util.h"
#include "product.h"


static void product_step_pointer(product* p);
static product* malloc_product();

static int pow2(int a, int b) {
    int result = 1;
    for (int i = 0; i < b; i++) result *= a;
    return result;
}

product* product_init(int a, int repeat) {
    product* p = malloc_product();
    p->pointer = array2_zeros(repeat);
    // printf("AAA = %d\n", a);
    p->A = array2_arange(a);
    p->counter = 0;
    p->END = pow2(p->A->size, p->pointer->size);
    return p;

}

int product_has_next(product* p) {
    if (p->counter == p->END) return 0;
    return 1;
}

static void product_step_pointer(product* p) {
    if (p->pointer->size > 0) {
        p->pointer->array2[0]++;
        for (int i = 0; i < p->pointer->size-1 ; i++)
            if (p->pointer->array2[i] % p->A->size == 0) {
                p->pointer->array2[i] = 0;
                p->pointer->array2[i+1]++;
            } else
                break;
    }
}

void product_finish(product* p) {
    array2_destroy(p->A);
    array2_destroy(p->pointer);
    free(p);
    p = NULL;
}

array2* product_next(product* p) {
    array2* result = array2_zeros(p->pointer->size);
    for (int i = 0; i < p->pointer->size; i++)
        result->array2[i] = p->A->array2[p->pointer->array2[i]];
    product_step_pointer(p);
    p->counter++;
    array2_reverse(result);
    return result;
}

static product* malloc_product() {
    product* ptr = (product*) malloc(sizeof(product));
    if (ptr == NULL) {
        error("malloc returned NULL!\n");
    }
    return ptr;
}

