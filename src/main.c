#include <R.h>
#include <Rinternals.h>
#include "mrfse.h"
#include "array.h"
#include "util.h"

void flatten_to_matrix(int **M, int n, int m, int *V) {
    for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
	    M[i][j] = V[n*j + i];
}

void matrix_to_flatten(int *V, int **M, int n, int m) {
    for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
	    V[n*j + i] = M[i][j];
}

void setUp(struct mrfse_data *data, SEXP a_size,
		 SEXP sample, SEXP c, SEXP max_neigh, SEXP k) {
    data->V_size = ncols(sample);
    data->A_size = asInteger(a_size);
    data->sample_size = nrows(sample);
    data->adj = array_matrix(data->V_size);
    data->sample = matrixINT(data->sample_size, data->V_size);
    flatten_to_matrix(data->sample, data->sample_size,
    		      data->V_size, INTEGER(sample));
    data->V = array_arange(data->V_size);
    data->A = array_arange(data->A_size);
    if (data->cv_enable == 1) {
	data->k = asInteger(k);
	data->fold_bloc = array_zeros(data->k + 1);
	data->fold = matrixINT(data->sample_size, data->V_size);
	data->out_fold = matrixINT(data->sample_size, data->V_size);
	data->c_values = REAL(c);
	data->c_values_size = length(c);
    } else 
	data->c = asReal(c);
    
    if (isNull(max_neigh))
	data->max_neigh = data->V_size - 1;
    else
	data->max_neigh = asInteger(max_neigh);
}

void setDown(struct mrfse_data *data) {
    if (data->cv_enable == 1) {
	free_matrixINT(data->fold, data->sample_size);
	free_matrixINT(data->out_fold, data->sample_size);
	array_destroy(data->fold_bloc);
    } else {
	for (int v = 0; v < data->V_size; v++)
	    array_destroy(data->adj[v]);
    }
    free(data->adj);
    array_destroy(data->V);
    array_destroy(data->A);
    free_matrixINT(data->sample, data->sample_size);
}

static SEXP array_to_vector(int i, struct mrfse_data *data) {
    SEXP v = PROTECT(allocVector(INTSXP, data->adj[i]->size));
    for (int j = 0; j < data->adj[i]->size; j++)
	INTEGER(v)[j] = data->adj[i]->array[j]+1;
    return v;
}

static void check_A(SEXP a_size) {    
    if (!isNumeric(a_size) || asInteger(a_size) <= 0)
	error("a_size argument must be a scalar positive integer");    
}

static void check_sample(SEXP sample, int a_size) {
    
    if (!isNumeric(sample) || !isMatrix(sample))
	error("sample argument must be a integer-entry matrix");

    PROTECT(sample = coerceVector(sample, INTSXP));
    for (int i = 0; i < length(sample); i++)
    	if (INTEGER(sample)[i] < 0 || INTEGER(sample)[i] >= a_size)
    	    error("sample elements must be belong range 0 and a_size - 1");
    UNPROTECT(1);
}

static void check_c(SEXP c) {
    if (!isNumeric(c))
	error("c argument must be double");
}

static void check_max_neigh(SEXP max_neigh, int V_size) {
    if (!isNull(max_neigh) && (!isNumeric(max_neigh) ||		       
			       V_size <= asInteger(max_neigh)))
	error("max_neigh, if used, must be a non-negative scalar "
	      "integer and be less than ncols(sample)");
}

static void check_k(SEXP k) {
    if (k != NULL && length(k) != 1 && !isNumeric(k))
	error("k argument must be a scalar integer");    
}

void input_checking(SEXP a_size, SEXP sample, SEXP c, SEXP max_neigh,
		    SEXP k) {
    check_A(a_size);
    check_sample(sample, asInteger(a_size));
    check_c(c);
    check_max_neigh(max_neigh, ncols(sample));
    check_k(k);
}

SEXP Rmrfse(SEXP a_size, SEXP sample, SEXP c, SEXP max_neigh) {
    input_checking(a_size, sample, c, max_neigh, NULL);
    SEXP ans;
    ans = PROTECT(allocVector(VECSXP, ncols(sample)));
    PROTECT(sample = coerceVector(sample, INTSXP));
    PROTECT(c = coerceVector(c, REALSXP));
    struct mrfse_data *data = (struct mrfse_data *)
	R_alloc(1, sizeof(struct mrfse_data));
    data->cv_enable = 0;
    setUp(data, a_size, sample, c, max_neigh, NULL);
    mrfse(data);
    for (int i = 0; i < ncols(sample); i++)
	SET_VECTOR_ELT(ans, i, array_to_vector(i, data));
    setDown(data);
    UNPROTECT(ncols(sample) + 3);
    return ans;
}

SEXP Rmrfse_cv(SEXP a_size, SEXP sample, SEXP c, SEXP k,
	      SEXP max_neigh) {
    input_checking(a_size, sample, c, max_neigh, k);
    SEXP ans = PROTECT(allocVector(REALSXP, 1));
    PROTECT(sample = coerceVector(sample, INTSXP));
    PROTECT(c = coerceVector(c, REALSXP));
    struct mrfse_data *data = (struct mrfse_data *)
	R_alloc(1, sizeof(struct mrfse_data));
    data->cv_enable = 1;
    setUp(data, a_size, sample, c, max_neigh, k);
    mrfse_cv(data);
    REAL(ans)[0] = data->c;
    setDown(data);
    UNPROTECT(3);
    return ans;
}

