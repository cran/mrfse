"mrfse" <- function(a_size, sample, c, max_neigh=NULL) {
     return (.Call('Rmrfse', a_size, sample, c, max_neigh))
}

"cv.mrfse" <- function(a_size, sample, can, k=10, max_neigh=NULL) {
    return (.Call('Rmrfse_cv', a_size, sample, can, k, max_neigh))
}

"con.mrfse" <- function(a_size, sample, c, max_neigh=NULL) {
    list.adj <- .Call('Rmrfse', a_size, sample, c, max_neigh)
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

"ncon.mrfse" <- function(a_size, sample, c, max_neigh=NULL) {
    list.adj <- .Call('Rmrfse', a_size, sample, c, max_neigh)
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
