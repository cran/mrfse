\name{mrfse.exact.con}
\encoding{utf-8}
\alias{mrfse.exact.con}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Conservative approach for Frondana's mrfse
}
\description{
  Conservative construction of the estimated Markov random
  field graph.
}
\usage{
mrfse.exact.con(a_size, sample, c, max_neigh = ncol(sample) - 1)

}
%- maybe also 'usage' for other objects documented here.
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
a = c(0, 1)
s = matrix(sample(a, size=1000, replace=TRUE), ncol=5)
mrfse.exact.con(length(a), s, 1.0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
