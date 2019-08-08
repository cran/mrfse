#include "array.h"

#ifndef MRFSE_H
#define MRFSE_H

struct mrfse_data {
    int V_size;
    array *V;
    int A_size;
    array *A;
    double c;
    int max_neigh;
    int max_edges;
    array **adj;
    int **sample;
    int sample_size;
    int num_graphs;
    /* For Cross Validation use */
    int k;
    array *fold_bloc;
    int **fold;
    int fold_size;
    int **out_fold;
    int out_fold_size;
    int cv_enable;
    double *c_values;
    int c_values_size;
};

void mrfse(struct mrfse_data *data);
void mrfse_cv(struct mrfse_data *data);

#endif
