#include <stdio.h>
#include <R.h>
#include <stdlib.h>
#include "array.h"
#include "combination.h"
#include "util.h"

static combination* malloc_combination();

combination* combination_init(array* a, int size) {
    combination* c = malloc_combination();
    c->n = a->size;
    c->k = size;
    c->comb = (int*) malloc_int(c->k);
    c->data = array_copy(a);
    if (size == 0)
	c->END = 1;
    else	    
	c->END = 0;
    for (int i = 0; i < size; i++)
	c->comb[i] = i;
    
    return c;
}

static int combination_step(combination* c) {
    int n = c->n;
    int k = c->k;
    int i = k - 1;
    while (i > 0 && c->comb[i] == n - k + i)
	i--;
    if (i == 0 && c->comb[i] == n - k)
	return 0;
    c->comb[i]++;
    for (; i < k - 1; i++)
	c->comb[i + 1] = c->comb[i] + 1;
    return 1;    
}



array* combination_next(combination* c) {
    array* result = array_zeros(c->k);
    for (int i = 0; i < c->k; i++) {
	result->array[i] = c->data->array[c->comb[i]];
    }
    if(!combination_step(c))
	c->END = 1;	
    return result;    
}

int combination_has_next(combination* c) {
    return !c->END;
}

void combination_finish(combination* c) {
    array_destroy(c->data);
    free(c->comb);
    free(c);
    c = NULL;
}

static combination* malloc_combination() {
    combination* ptr = (combination*) malloc(sizeof(combination));
    if (ptr == NULL) {
        error("malloc returned NULL!\n");
    }
    return ptr;
}


