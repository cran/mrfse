\name{mrfse.sa.ncon}
\encoding{utf-8}
\alias{mrfse.sa.ncon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Non-conservative approach for Frondana's mrfse using simulated annealing
}
\description{
  A penalized likelihood BIC-based to estimate Markovian neighborhoods.
}
\usage{
mrfse.sa.ncon(a_size, sample, c, t0, iterations=1000, max_neigh=ncol(sample)-1)
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
  \item{c}{
    The penalization constant. Must be positive.
  }
  \item{t0}{
    Inital temperature
  }

  \item{iterations}{
    Number of simulated annealing iterations
}
  \item{max_neigh}{
    The maximum length of a candidate Markovian neighborhood. Must be
  non-negative and less than \code{ncol(sample)}. 
}

}

\value{
  A adjacency matrix of the estimated Markov random field graph.
}
\references{
  FRONDANA, Iara Moreira. \emph{Model selection for discrete Markov random
  fields on graphs}. São Paulo : Instituto de Matemática e Estatística,
  University of São Paulo, 2016. Doctoral Thesis in Estatística.
  <doi:10.11606/T.45.2018.tde-02022018-151123>
  \url{http://www.teses.usp.br/teses/disponiveis/45/45133/tde-02022018-151123/publico/tese_Iara_Frondana.pdf}
}
\author{
  Rodrigo Carvalho
}

\examples{
library(mrfse)
a_size = c(0, 1)
s = matrix(sample(a_size, size=1000, replace=TRUE), ncol=5)
mrfse.sa.ncon(length(a_size), s, 1.0, 500, 1000)
}