\name{mrfse.sample}
\encoding{utf-8}
\alias{mrfse.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a independent sample of a Markov random field
}
\description{
Generate a independent sample of a Markov random field according to the probilities of the sampler.

}
\usage{
mrfse.sample(sampler, n)
}

\arguments{
  \item{sampler}{
	A sampler created by \code{mrfse.create.sampler} function
}
  \item{n}{
	Size of sample
}
}

\value{
  A matrix whose number of columns is the number of nodes. Each line is a single independent sample of Markov random field given by the probabilites of sampler.
}
\author{
  Rodrigo Carvalho
}

\examples{
library(mrfse)
adj = matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0), byrow=TRUE, ncol=3)
sampler = mrfse.create.sampler(adj, 3)
mrfse.sample(sampler, 3000)
}


