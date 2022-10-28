\name{mrfse.ci}
\encoding{utf-8}
\alias{mrfse.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bresler's non-binary Markov random field structure estimator
}
\description{
  A greedy algorithm to estimate Markovian neighborhoods.
}
\usage{
mrfse.ci(a_size, sample, tau, max_degree=ncol(sample)-1)
}

\arguments{
  \item{a_size}{
    Size of the alphabet.
}
  \item{sample}{
    A integer-valued matrix. Each value must belong range \code{0}
    and \code{a_size - 1}. Matrix has dimension \code{n} x \code{V}, where
    \code{n} is number of samples and \code{V} is number of nodes.
}
  \item{tau}{
    A hyperparameter. See references.
  }
  \item{max_degree}{
    The maximum length of a candidate Markovian neighborhood. Must be
  non-negative and less than \code{ncol(sample)}.
}

}

\value{
  A list filled with estimated Markov neighborhood for each graph vertex
}
\references{
  Guy Bresler. 2015. Efficiently Learning Ising Models on Arbitrary Graphs. In Proceedings of the forty-seventh annual ACM symposium on Theory of Computing (STOC '15). Association for Computing Machinery, New York, NY, USA, 771–782. DOI:https://doi.org/10.1145/2746539.2746631
}

\author{
  Rodrigo Carvalho
}

\examples{
library(mrfse)
a_size = c(0, 1)
s = matrix(sample(a_size, size=1000, replace=TRUE), ncol=5)
mrfse.ci(length(a_size), s, 0.2)
}