library(Rfast)
library(gtools)

moral_graph <- function(dag_adj) {
    moral_adj <- dag_adj
    for (i in 1:ncol(moral_adj)) {
        same_child <- which(dag_adj[, i] == 1)
        if (length(same_child) >= 2) {
            all_comb <- combn(same_child, 2)
            for (j in 1:ncol(all_comb)) {
                x <- all_comb[1, j]
                y <- all_comb[2, j]
                moral_adj[x, y] <- 1
                moral_adj[y, x] <- 1
            }
        }
    }
    for (i in 1:nrow(moral_adj)) {
        for (j in 1:ncol(moral_adj)) {
            if (moral_adj[i, j] == 1)
                moral_adj[j, i] = 1
        }
    }
    moral_adj
}

mrfse.create.sampler <- function(dag_adj, A) {
    sampler <- NULL

    for (i in 1:ncol(dag_adj)) {
        neigh <- which(dag_adj[,i] > 0)
        sampler$neigh[[i]] <- neigh
        sampler$probs[[i]] <- matrix(ncol=A, nrow=A^length(neigh))
        ## sampler$probs[[i]] <- rdirichlet(A^length(neigh), rep(1, A))
        for (j in 1:A^length(neigh)) {
            sampler$probs[[i]][j,] <- rdirichlet(1, rep(1, A))
            while(min(sampler$probs[[i]][j,]) < 1/(3 * A))
                sampler$probs[[i]][j,] <- rdirichlet(1, rep(1, A))
        }
    }

    sampler$moral_adj <- moral_graph(dag_adj)

    sampler$topol_sort <- topological_sort(dag_adj)

    sampler$num_nodes <- ncol(dag_adj)
    sampler$A <- A
    sampler$max_degree <- max_degree_2(sampler$moral_adj)
    sampler
}

idx_to_num <- function(v, A) {
    n = length(v)
    result = 1
    for (i in 1:n) {
        result = result + ifelse(v[i] > 0, (v[i] * A^(i-1)), 0)
    }
    if (n == 0)
        result = 1
    result
}

mrfse.sample <- function(sampler, n) {
    sample <- matrix(ncol=sampler$num_nodes, nrow=n)
    A <- sampler$A
    for (i in 1:n) {
        for (v in sampler$topol_sort) {
            neigh <- sampler$neigh[[v]]
            config <- sample[i, neigh]
            prob <- sampler$probs[[v]][idx_to_num(config, A),]
            sample[i, v] <- sample(0:(A - 1), 1, prob=prob)
        }
    }
    sample
}

random_DAG <- function(n, prob=0.2) {
    adj <- matrix(0, ncol = n, nrow = n)
    for (i in 1:n) {
        j <- 1
        while (j < i) {
            adj[i, j] <- ifelse(prob > runif(1), 1, 0)
            j <- j + 1
        }

    }
    adj
}

max_degree_2 <- function(adj) {
    result = 0
    me <- 0
    for (i in 1:ncol(adj)) {
        result = max(result, sum(adj[i, ]))
        me <- me + sum(adj[i, ])
    }
    me <- me / ncol(adj)
    return (result)
}
