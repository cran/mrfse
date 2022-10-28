#source("RcppExports.R")

"mrfse.exact" <- function(a_size, sample, c, max_neigh=ncol(sample)-1) {
    ## if (max_neigh == NULL)
    ##     max_neigh = ncol(sample)
    return (mrfse(a_size, sample, c, max_neigh))
}

## "mrfse.exact.cv" <- function(a_size, sample, can, k=10, max_neigh=NULL) {
##     return (.Call('Rmrfse_cv', a_size, sample, can, k, max_neigh))
## }

"mrfse.exact.con" <- function(a_size, sample, c, max_neigh=ncol(sample)-1) {
    ## if (max_neigh == NULL)
        ## max_neigh = ncol(sample)
    list.adj <- mrfse(a_size, sample, c, max_neigh)
    n <- length(list.adj)
    matrix.adj <- matrix(rep(0, n**2), ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (is.element(j, list.adj[[i]]) &&
                is.element(i, list.adj[[j]])) {
                matrix.adj[i, j] <- 1
                matrix.adj[j, i] <- 1
            }
        }
    }
    return (matrix.adj)
}

"mrfse.exact.ncon" <- function(a_size, sample, c, max_neigh=ncol(sample)-1) {
    ## if (max_neigh == NULL)
    ##     max_neigh = ncol(sample)
    list.adj <- mrfse(a_size, sample, c, max_neigh)
    n <- length(list.adj)
    matrix.adj <- matrix(rep(0, n**2), ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (is.element(j, list.adj[[i]])) {
                matrix.adj[i, j] <- 1
                matrix.adj[j, i] <- 1
            }
        }
    }
    return (matrix.adj)
}

"mrfse.sa" <- function(a_size, sample, c, t0, iterations=1000,
                       max_neigh=ncol(sample)-1) {
    ## if (max_neigh == NULL)
    ##     max_neigh = ncol(sample)
    return (mrfse_sa(a_size, sample, c, t0, iterations, max_neigh))
}

"mrfse.sa.ncon" <- function(a_size, sample, c, t0, iterations=1000,
                            max_neigh=ncol(sample)-1) {
    ## if (max_neigh == NULL)
    ##     max_neigh = ncol(sample)
    list.adj <- mrfse_sa (a_size, sample, c, t0, iterations, max_neigh)
    n <- length(list.adj)
    matrix.adj <- matrix(rep(0, n**2), ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (is.element(j, list.adj[[i]])) {
                matrix.adj[i, j] <- 1
                matrix.adj[j, i] <- 1
            }
        }
    }
    return (matrix.adj)
}


"mrfse.sa.con" <- function(a_size, sample, c, t0, iterations=1000,
                           max_neigh=ncol(sample)-1) {
    list.adj <- mrfse_sa(a_size, sample, c, t0, iterations, max_neigh)
    n <- length(list.adj)
    matrix.adj <- matrix(rep(0, n**2), ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (is.element(j, list.adj[[i]]) &&
                is.element(i, list.adj[[j]])) {
                matrix.adj[i, j] <- 1
                matrix.adj[j, i] <- 1
            }
        }
    }
    return (matrix.adj)
}

"mrfse.ci" <- function(a_size, sample, tau, max_degree=ncol(sample)-1) {
    return (mrfse_ci(a_size, sample, tau, max_degree))
}

"mrfse.ci.ncon" <- function(a_size, sample, tau, max_degree=ncol(sample)-1) {
    list.adj <- mrfse_ci(a_size, sample, tau, max_degree)
    matrix.adj <- list_to_ncon(list.adj)
    return (matrix.adj)
}

"mrfse.ci.con" <- function(a_size, sample, tau, max_degree=ncol(sample)-1) {
    list.adj <- mrfse_ci(a_size, sample, tau, max_degree)
    matrix.adj <- list_to_con(list.adj)
    return (matrix.adj)
}

list_to_con <- function(adj) {
    n <- length(adj)
    mat <- matrix(0, ncol=n, nrow=n)
    for (i in 1:n)
        for (j in 1:n)
            if (is.element(i, adj[[j]]) && is.element(j, adj[[i]])) {
                mat[i, j] <- 1
                mat[j, i] <- 1
            }
    mat
}

list_to_ncon <- function(adj) {
    n <- length(adj)
    mat <- matrix(0, ncol=n, nrow=n)
    for (i in 1:n)
        for (j in 1:n)
            if (is.element(i, adj[[j]]) ) {
                mat[i, j] <- 1
                mat[j, i] <- 1
            }
    mat
}
