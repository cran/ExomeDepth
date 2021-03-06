% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{viterbi.hmm}
\alias{viterbi.hmm}
\title{Computes the Viterbi path for a hidden markov model}
\usage{
viterbi.hmm(transitions, loglikelihood, positions, expected.CNV.length)
}
\arguments{
\item{transitions}{Transition matrix}

\item{loglikelihood}{numeric matrix containing the loglikelihood of the data
under the possible states}

\item{positions}{Positions of the exons}

\item{expected.CNV.length}{Expected length of CNV calls, which has an impact
on the transition matrix between CNV states.}
}
\value{
\item{comp1 }{Description of 'comp1'} \item{comp2 }{Description of
'comp2'}
}
\description{
Estimates the most likely path for a hidden Markov Chain using the maximum
likelihood Viterbi algorithm.
}
\details{
Standard forward-backward Viterbi algorithm using a precomputed matrix of
likelihoods.
}
