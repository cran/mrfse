#include <Rcpp.h>
#include <stdlib.h>
#include <omp.h>
#include "mrfse.h"
#include "util.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

// Global variables

static vector<vector<int>> mysample;
static double tau;
static int p;
static int n;
static int A;
static int max_degree;
static bool use_kl;

// Enviroment
Environment gtools = Environment::namespace_env("gtools");


// R native functions

Function asMatrix("as.matrix");
Function asVector("as.vector");
// Function sda("sd");
Function expand_grid("expand.grid");
// Function fpermutations = gtools["permutations"];
Function unlist("unlist");

template<class C, class T>
auto contains(const C& v, const T& x)
-> decltype(end(v), true)
{
    return end(v) != find(begin(v), end(v), x);
}


void init_data(int a, IntegerMatrix sample, double t, int md) {

    tau        = t;
    p          = sample.ncol();
    n          = sample.nrow();
    A          = a;
    max_degree = md;

    mysample.resize(n);
    for (int i = 0; i < n; i++) {
	mysample[i].resize(p);
	for (int j = 0; j < p; j++) mysample[i][j] = sample(i, j);
    }
}

// vector<int> subvector(int i, vector<int> &S) {
//     vector<int> a;
//     a.resize(S.size());
//     for (int j = 0; i < a.size(); i++) a[i] = mysample[S[i]];
//     return a;
    
// }
// void printVector(IntegerVector a) {
//     for (int i = 0; i < a.size(); i++)
// 	printf("%d ", a[i]);
//     printf("\n");
// }

bool equals(IntegerVector a, IntegerVector b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++)
	if (a[i] != b[i]) return false;
    return true;
}

double pviS(int v, int i, vector<int> &S, int xv, int xi,
	    vector<int> &xS) {
    int N_iS = 0, N_v_iS = 0;
    vector<int> sam_xS(S.size());
    for (int j = 0; j < n; j++) {
	for (int k = 0; k < S.size(); k++) sam_xS[k] = mysample[j][S[k]];
	if (mysample[j][i] == xi && sam_xS == xS ) {
	    N_iS++;
	    if (mysample[j][v] == xv) N_v_iS++;
	}	
    }
    if (N_iS == 0) return (double) 1 / (A + 1);
    return (double) N_v_iS / N_iS;
}

double pvS(int v, int xv, vector<int> &S, vector<int> &xS) {
    int NvS = 0, NS = 0;
    vector<int> sam_xS(S.size());
    for (int i = 0; i < n; i++) {
	// IntegerVector mysamplei = mysample(i, _);
	for (int j = 0; j < S.size(); j++) sam_xS[j] = mysample[i][S[j]];
	if (sam_xS == xS) {
	    if (mysample[i][v] == xv) NvS++;
	    NS++;
	}
    }
    if (NS == 0) return (double) 1 / (A + 1);
    return (double) NvS / NS;
}

double pS(vector<int> &S, vector<int> &xS) {
    int N = 0;
    vector<int> sam_S(S.size());
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < S.size(); j++) sam_S[j] = mysample[i][S[j]];
	// IntegerVector mysamplei = mysample(i, _);
	if (sam_S == xS) N++;
    }
    return (double) N / n;    
}

double weigth(int i, vector<int> &S, vector<int> &xS) {
    // return 1;
    double result = 1;
    for (int a = 0; a < A; a++) {
	double p = pvS(i, a, S, xS);
	result *= p * (1 - p);
    }
    return (A + 1)  * result;
}

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

double average_eta_2(int v, int i, vector<int> &S) {
    vector<vector<int>> xS, xv, xi;
    double average = 0.0;
    bool empty_flag = true;	
    xv = permutations(A, 1), xi = permutations(A, 2);
    xS = permutations(A, S.size());
    for (int j = 0; j < xS.size() || empty_flag; j++) {
	empty_flag = false;
    	double eta = 0.0;
    	for (int k = 0; k < xv.size(); k++) {
	    for (int l = 0; l < xi.size(); l++) {
		if(xi[l][0] != xi[l][1]) {
		    double p1 = pviS(v, i, S, xv[k][0], xi[l][0], xS[j]);
		    double p2 = pviS(v, i, S, xv[k][0], xi[l][1], xS[j]);
		    eta += abs(p1 - p2);
		}
	    }
	}
    	average += eta * weigth(i, S, xS[j]) * pS(S, xS[j]);
    }
    return average;
}

double kullback(int v, int i, vector<int> &S) {
    vector<vector<int>> xS, xv, xi, Si;
    double average = 0.0;
    bool empty_flag = true;	
    xv = permutations(A, 1), xi = permutations(A, 1);
    xS = permutations(A, S.size());
    for (int j = 0; j < xS.size() || empty_flag; j++) {
	empty_flag = false;
	for (int l = 0; l < xi.size(); l++) {
	    bool zero_flag = false;
	    double eta = 0.0;
	    for (int k = 0; k < xv.size(); k++) {
		double pvS_var = pvS(v, xv[k][0], S, xS[j]);
		if (pvS_var == 0.0) zero_flag = true;
	    }

	    if (zero_flag) continue;		

	    for (int k = 0; k < xv.size(); k++) {
		double pvS_var = pvS(v, xv[k][0], S, xS[j]);
		double pviS_var = pviS(v, i, S, xv[k][0], xi[l][0], xS[j]);
		if (pviS_var > 0.0) {
		    eta += pviS_var * log(pviS_var / pvS_var);
		}
	    }
	    S.push_back(i);
	    xS[j].push_back(xi[l][0]);
	    average += eta * pS(S, xS[j]);
	    xS[j].pop_back();
	    S.pop_back();
	}
    }
    return average;
    
}

double average_eta(int v, int i, vector<int> &S) {
    vector<vector<int>> xS, xv, xi;
    double average = 0.0;
    bool empty_flag = true;	
    xv = permutations(A, 1), xi = permutations(A, 1);
    xS = permutations(A, S.size());
    for (int j = 0; j < xS.size() || empty_flag; j++) {
	empty_flag = false;
    	double eta = 0.0;
    	for (int k = 0; k < xv.size(); k++) {
	    double pvS_var = pvS(v, xv[k][0], S, xS[j]);
	    for (int l = 0; l < xi.size(); l++) {
		double p = pviS(v, i, S, xv[k][0], xi[l][0], xS[j]);
		eta += abs(p - pvS_var);
	    }
	}
    	average += eta * weigth(i, S, xS[j]) * pS(S, xS[j]);
    }
    return average;
}

void best_average_eta(int* i, double* eta, int v, vector<int> &S) {
    *i   = -1;
    *eta = 0;       
    for (int j = 0; j < p; j++) {
    	if (j != v && !contains(S, j)) {
	    double e;
	    // if (use_kl) e = kullback(v, j, S);
	    // else        e = average_eta_2(v, j, S);
	    e = kullback(v, j, S);
	    if (e > *eta) {
		*eta = e;
		*i   = j;
	    }	    
    	}
    }
}

vector<int> pseudo_neighborhood(int v) {
    vector<int> S;
    int i;
    double eta;

    for(;;) {
    	best_average_eta(&i, &eta, v, S);
    	if(eta < tau || S.size() >= max_degree || i == -1) break;
    	S.push_back(i);
    }
    return S;
}

vector<int> pruning(int v, vector<int> &S) {
    list<int> able;
    vector<int> vectAble;
    for(int j = 0; j < S.size(); j++) {
	int i = S[j];
	S.erase(S.begin() + j);
	double eta;
	// if (use_kl) eta = kullback(v, i, S);
	// else        eta = average_eta_2(v, i, S);
	eta = kullback(v, i, S);
	S.insert(S.begin() + j, i);
	if (eta > tau) able.push_back(i);	
    }    
    vectAble = vector<int>(begin(able), end(able));    
    return vectAble;
}

static vector<int> estimate_neighborhood(int v) {
    vector<int> S = pseudo_neighborhood(v);
    S = pruning(v, S);
    return S;
}



vector<int> bla2(int v) {
    vector<int> a = {v};
    return a;
}

void aaa(vector<vector<int>> &b) {
    for (int i = 0; i < b.size(); i++)
	for (int j = 0; j < b[i].size(); j++)
	    int e = 3 + b[i][j];
}

IntegerMatrix cc;

// void bbb(IntegerMatrix bb) {
//     for (int i = 0; i < bb.nrow(); i++)
// 	for (int j = 0; j < bb.ncol(); j++)
// 	    int e = bb(i, j) + 34;
// }

void bbb(vector<vector<int>> &c) {
    for (int i = 0; i < c.size(); i++)
	for (int j = 0; j < c[i].size(); j++)
	    int e = c[i][j] + 34;
}


// [[Rcpp::export]]
List mrfse_ci(int A, IntegerMatrix sample, double tau, int max_degree) {

    A--;
    init_data(A, sample, tau, max_degree);
    List result(p);
    vector<vector<int>> aa = permutations(A, 1);
    int n = aa.size();
    vector<vector<int>> result_vect(p);
    #pragma omp parallel for shared(mysample, A, tau, p, n, max_degree,	\
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


