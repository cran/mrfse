#ifndef MRFSE_H
#define MRFSE_H

#include <Rcpp.h>

using namespace Rcpp;

struct mrfse_data{
    double tau;
    IntegerMatrix sample;
    int p;
    int n;
    int A;    
};

#endif
