#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <Rcpp.h>
#include <omp.h>
#include "mrfse.h"
#include "product.h"
#include "combination.h"
#include "array.h"
#include "util.h"
#include "simulated_annealing.h"

using namespace Rcpp;
using namespace std;


static vector<vector<int>> mysample;
static vector<int> sa_t;
static int A;
static int max_degree;
static int iterations;
static double t0;
static double c;
static int p;
static int n;

void init_data(int a, IntegerMatrix sample, double cons, double temp, int iter, int md) {

    A          = a;
    p          = sample.ncol();
    n          = sample.nrow();
    t0         = temp;
    iterations = iter;
    max_degree = md;
    c          = cons;
    
    sa_t.resize(p);
    mysample.resize(n);
    for (int i = 0; i < n; i++) {
	mysample[i].resize(p);
	for (int j = 0; j < p; j++) mysample[i][j] = sample(i, j);
    }
}

void count_in_sample(int v, array2 *W, array2* a,
			    array2* aW, int *N_W, int *N_v_W) {
    *N_W = 0, *N_v_W = 0;
    int m = W->size;
    array2* x_W = array2_zeros(m);    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
	    x_W->array2[j] = mysample[i][W->array2[j]];
        if (array2_equals(x_W, aW)) {
            if (mysample[i][v] == a->array2[0])
		*N_v_W = *N_v_W + 1;
            *N_W = *N_W + 1;
        }
    }
    array2_destroy(x_W);        
}

long double penalized_factor(int W) {
    return (c * pow(A, W) * (log(n) / log(A)));
}

long double likelihood(int v, array2* W, array2* a, array2* aW) {

    int N_W, N_v_W;
    long double p_hat = 0.0;
    count_in_sample(v, W, a, aW, &N_W, &N_v_W);
    if (N_W == 0)
	p_hat = 1.0 / A;
    else 
	p_hat = (long double) N_v_W / N_W;
    if (N_v_W == 0)
	return 0.0;    
    return (long double) N_v_W * log(p_hat);
}

long double L_vertex(int v, array2* W) {
    int m = W->size;
    long double L_value = 0.0;
    product* p = product_init(A, 1);
    while (product_has_next(p)) {
        array2* a = product_next(p);
        product* pW = product_init(A, m);
        while (product_has_next(pW)) {
            array2* aW = product_next(pW);
	    // array2_print(a);
	    // printf("\n");
	    // array2_print(aW);
	    // printf("\n");
            L_value += likelihood(v, W, a, aW);
            array2_destroy(aW);            
        }
        array2_destroy(a);
        product_finish(pW);
    }
    product_finish(p);
    return L_value - penalized_factor(W->size);
    /* return -(log(fabsl(L_value - penalized_factor(W->size, data)))); */
}


static vector<int> estimate_neighborhood(int v) {
    double best_value = -1 * INF;
    array2* best_neighborhood = array2_zeros(0);
    array2* V = array2_arange(p);
    V = array2_erase(V, v);    
    for (int i = 0; i <= max_degree; i++) {
        combination* c = combination_init(V, i);
        while (combination_has_next(c)) {

            array2* W = combination_next(c);
            double PL_value = L_vertex(v, W);
            if (PL_value > best_value) {
                best_value = PL_value;
                array2_destroy(best_neighborhood);
                best_neighborhood = W;
            } else 
                array2_destroy(W);
        }

        combination_finish(c);
    }

    array2_destroy(V);
    // array2_print(best_neighborhood);
    // printf("best = %lf\n", best_value);
    // Printf("\n");
    return array_to_vec(best_neighborhood);
}

vector<int> estimate_neighborhood_sa(int v) {

    array2* current_neighborhood = sa_initial_neigh(v, p);
    array2* best_neighborhood = array2_copy(current_neighborhood);
    long double current_value = L_vertex(v, current_neighborhood);
    long double best_value = L_vertex(v, best_neighborhood);
    int accepted = 0;

    for (int i = 0; i < iterations; i++) {
	array2* W = sa_next_neigh(v, current_neighborhood, p, max_degree);
    	long double PL_value = L_vertex(v, W);
	long double boltz_factor = sa_boltzmann_factor(current_value, PL_value, v, sa_t);
	if (PL_value > best_value) {
	    array2_destroy(best_neighborhood);
	    best_neighborhood = array2_copy(W);
	    best_value = PL_value;
	}	
	if (PL_value > current_value || unif_rand() < boltz_factor) {
    	    array2_destroy(current_neighborhood);
    	    current_neighborhood = W;
	    current_value = PL_value;
	    accepted++;
    	}
	else
	    array2_destroy(W);
    }
    array2_destroy(current_neighborhood);
    return array_to_vec(best_neighborhood);
}

// static void estimate_graph() {
//     #pragma omp parallel for
//     for (int v = 0; v < p; v++)
// 	data->adj[v] = estimate_neighborhood(v, data);
// }

// static void estimate_graph_sa(struct mrfse_data* data) {
//     data->sa_t = (long double*) malloc(data->V_size * sizeof(long double));
//     #pragma omp parallel for
//     for (int v = 0; v < data->V_size; v++) {
// 	data->sa_t[v] = data->sa_t0;
// 	data->adj[v] = estimate_neighborhood_sa(v, data);
//     }    
// }

// static void cv_blocs(struct mrfse_data *data) {
//     int q = data->sample_size / data->k;
//     int r = data->sample_size % data->k;
//     data->fold_bloc->array[0] = -1;
//     for (int i = 1; i <= r; i++) {
//         data->fold_bloc->array[i] = data->fold_bloc->array[i-1] + q + 1;
//     }
//     for (int i = r+1; i <= data->k; i++)
//         data->fold_bloc->array[i] = data->fold_bloc->array[i-1] + q;
// }

// static void get_fold(int k, struct mrfse_data *data) {
//     int a = data->fold_bloc->array[k];
//     int b = data->fold_bloc->array[k-1];
//     data->fold_size = a - b;
//     for (int i = b+1, j = 0; i <= a; i++) {
// 	for (int k = 0; k < data->V_size; k++)
// 	    data->fold[j][k] = data->sample[i][k];
// 	j++;
//     }
// }

// static void get_out_fold(int k, struct mrfse_data *data) {
//     int a = data->fold_bloc->array[k];
//     int b = data->fold_bloc->array[k-1];
//     data->out_fold_size = data->sample_size - (a - b);
//     for (int i = 0, j = 0; i < data->sample_size; i++) {
// 	if (!(i > b && i <= a)) {
// 	    for (k = 0; k < data->V_size; k++)
// 		data->out_fold[j][k] = data->sample[i][k];
// 	    j++;
// 	}
//     }
// }

// static void sample_cv(struct mrfse_data *data) {
//     matrixINTcpy(data->sample, data->out_fold,
// 		 data->out_fold_size, data->V_size);
//     data->sample_size = data->out_fold_size;    
// }

// static void un_sample_cv(int **tmp, struct mrfse_data *data) {
//     data->sample_size = data->out_fold_size + data->fold_size;
//     matrixINTcpy(data->sample, tmp, data->sample_size, data->V_size);
// }

// static double L_vertex_cv(int v, struct mrfse_data *data) {
//     double value = 0.0;
//     product *p = product_init(data->A, 1);
//     array *W = data->adj[v];
//     while (product_has_next(p)) {
// 	array *a = product_next(p);
// 	product *pW = product_init(data->A, W->size);
// 	while(product_has_next(pW)) {
// 	    array *aW = product_next(pW);
// 	    value += likelihood_cv(v, W, a, aW, data);
// 	    array_destroy(aW);
// 	}
// 	product_finish(pW);
// 	array_destroy(a);
//     }
//     array_destroy(W);
//     product_finish(p);
//     return value;    
// }

// static double cv_value(struct mrfse_data *data) {
//     double value = 0.0;    
//     int **tmp = matrixINT(data->sample_size, data->V_size);
//     matrixINTcpy(tmp, data->sample, data->sample_size, data->V_size);
//     for (int i = 1; i <= data->k; i++) {
// 	get_fold(i, data);
// 	get_out_fold(i, data);
// 	sample_cv(data);
// 	estimate_graph(data); /* Estimate graph with K - 1 folds */
//     	for (int v = 0; v < data->V_size; v++) {
// 	    value += L_vertex_cv(v, data);
//     	}
//     	un_sample_cv(tmp, data);
//     }    
//     free_matrixINT(tmp, data->sample_size);
//     return value / data->k;
// }

// void mrfse_cv(struct mrfse_data *data) {
//     double best_value = -INF, best_c = 0.0;
//     cv_blocs(data);    
//     for (int i = 0 ; i < data->c_values_size; i++) {
//     	data->c = data->c_values[i];
//     	double value = cv_value(data);
//     	if (value > best_value) {
//     	    best_value = value;
//     	    best_c = data->c_values[i];
//     	}
//     }
//     data->c = best_c;
// }


static vector<vector<int>> permutations(int A, int r) {
    int N = (int) pow((double) (A + 1), r);
    vector<vector<int>> result;
    result.resize(N);
    for (int i = 0; i < result.size(); i++)
	result[i].resize(r, 0);

    int i = 0, j = 0;
    while (i < N) {
	int j = i, k = 0;
	while(j > 0) {
	    int a = j % (A + 1);
	    j /= (A + 1);
	    result[i][k++] = a;
	}
	i++;
    }
    return result;
}

// [[Rcpp::export]]
List mrfse(int A, IntegerMatrix sample, double c, int max_degree) {

    // A--;
    init_data(A, sample, c, -1, -1, max_degree);
    List result(p);
    vector<vector<int>> aa = permutations(A, 1);
    int n = aa.size();
    vector<vector<int>> result_vect(p);
    // #pragma omp parallel for shared(mysample, A, c, p, n, max_degree,	\
    				    result)
    for (int v = 0; v < p; v++) {
	result_vect[v] = estimate_neighborhood(v);
    }

    for (int v = 0; v < p; v++) {
    	IntegerVector S = wrap(result_vect[v]);
    	result[v] = S+1;
    }
    return result;
}

// [[Rcpp::export]]
List mrfse_sa(int A, IntegerMatrix sample, double c, double t0, int iterations, int max_degree) {
    // srand(12345);
    // estimate_graph_sa(data);
    init_data(A, sample, c, t0, iterations, max_degree);
    // A--;
    List result(p);

    vector<vector<int>> aa = permutations(A, 1);
    int n = aa.size();

    vector<vector<int>> result_vect(p);
    #pragma omp parallel for shared(mysample, A, c, t0, iterations, \
    p, n, max_degree, result)
    for (int v = 0; v < p; v++) {
	sa_t[v] = t0;
	result_vect[v] = estimate_neighborhood_sa(v);
    }
    for (int v = 0; v < p; v++) {
    	IntegerVector S = wrap(result_vect[v]);
    	result[v] = S+1;
    }
    return result;
}
