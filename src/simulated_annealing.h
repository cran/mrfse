#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H
#include "array.h"
/* #include "mrfse.h" */
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

array2* sa_next_neigh(int v, array2* W, int n, int max_degree);
array2* sa_initial_neigh(int v, int n);
long double sa_boltzmann_factor(long double actual, long double neigh, int k, vector<int>& sa_t);

#endif
