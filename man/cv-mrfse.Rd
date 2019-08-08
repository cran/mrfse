\name{cv.mrfse}
\encoding{UTF-8}
\alias{cv.mrfse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cross-validation for mrfse
}
\description{
  K-fold Cross-validation to estimate a penalization constant.
}
\usage{
cv.mrfse(a_size, sample, can, k = 10, max_neigh = NULL)
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
  \item{can}{
    Candidates values for penalization constant.
}
  \item{k}{
    Size of fold.
}
  \item{max_neigh}{
    The maximum length of a candidate Markovian neighborhood. Must be
  non-negative and less than \code{ncol(sample)}. If
  \code{NULL} then \code{max_neigh} will be set as \code{ncol(sample) - 1}.
}
}
\value{
  Return the best fitted for k-fold-cv in \code{can}.
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
can = seq(0.1, 1.0, 0.1)
s = matrix(sample(a_size, size=1000, replace=TRUE), ncol=5)
cv.mrfse(length(a_size), s, can, 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
