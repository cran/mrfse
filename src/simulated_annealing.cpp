#include "simulated_annealing.h"
#include "array.h"
#include "util.h"
#include "mrfse.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <R.h>

static array2* add_neigh(int v, array2* W, int n) {
    int m = W->size + 1;
    int e;
    array2* new_neigh = array2_zeros(m);
    for (int i = 0; i < W->size; i++)
	new_neigh->array2[i] = W->array2[i];
    new_neigh->array2[m-1] = -1;
    e = int_unif(n) % n;
    while (array2_contains(new_neigh, e) || e == v)
	e = int_unif(n) % n;
    
    new_neigh->array2[m-1] = e;
    return new_neigh;
}

static array2* remove_neigh(int v, array2* W, int n) {
    int m = W->size - 1;
    array2* new_neigh = array2_zeros(m);
    int e = int_unif(n) % n;
    while (!array2_contains(W, e) || e == v)
	e = int_unif(n) % n;

    for (int i = 0, j = 0; i < m+1; i++) {
	if (W->array2[i] != e)
	    new_neigh->array2[j++] = W->array2[i];
    }
    return new_neigh;    
}

array2* sa_next_neigh(int v, array2* W, int n, int max_degree) {
    if (W->size == 0) return add_neigh(v, W, n);
    if (W->size == n-1 || W->size == max_degree) return remove_neigh(v, W, n);

    float p = unif_rand();
    if (p < 0.5) return add_neigh(v, W, n);
    return remove_neigh(v, W, n);    
}

static long double sa_temp(int v, vector<int>& sa_t) {
    long double alpha = 0.99;
    long double result = sa_t[v];
    sa_t[v] = sa_t[v] * alpha;
    return result;
}

long double sa_boltzmann_factor(long double actual, long double neigh, int v,
			   vector<int>& sa_t) {
    return expl((neigh - actual) / sa_temp(v, sa_t));    
}

array2* sa_initial_neigh(int v, int n)  {
    /* int m = (n - 1) / 2; */
    int m = 0;
    array2* initial_neigh = array2_zeros(m);

    for (int i = 0; i < m; i++)
    	initial_neigh->array2[i] = v;

    for (int i = 0; i < m; i++) {
    	int e = int_unif(n) % n;
    	while (array2_contains(initial_neigh, e) || e == v) {
    	    e = int_unif(n) % n;
    	}
    	initial_neigh->array2[i] = e;
    }
    return initial_neigh;
}
