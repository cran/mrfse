\name{mrfse.create.sampler}
\encoding{utf-8}
\alias{mrfse.create.sampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a sampler for Markov random field.
}
\description{
  Create a sampler for Markov random field from a DAG
}
\usage{
mrfse.create.sampler(dag_adj, A)
}

\arguments{
  \item{dag_adj}{
	An direct acyclic graph adjacency matrix
}
  \item{A}{
	Size of alphabet

}
}
\value{
  A list filled with the following components:

neigh: A list of neighborhood. For each i, neigh[[i]] is a markovian neighborhood of vertex i

probs: A list of probabilities. For each i, probs[[i]] is matrix of probabilities of vertex i given your markovian neighborhood. Those probabilites will be used to generate a sample.

moral_adj: moral graph of adj_dag

topol_sort: topological sort of adj_dag

num_nodes: number of nodes de adj_dag

A: alphabet size
}

\author{
  Rodrigo Carvalho
}

\examples{
library(mrfse)
adj = matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0), byrow=TRUE, ncol=3)
mrfse.create.sampler(adj, 3)
}


