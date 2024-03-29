% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_definition.R
\name{initialize,ExomeDepth-method}
\alias{initialize,ExomeDepth-method}
\title{ExomeDepth initialization tool}
\usage{
\S4method{initialize}{ExomeDepth}(
  .Object,
  data = NULL,
  test,
  reference,
  formula = "cbind(test, reference) ~ 1",
  phi.bins = 1,
  prop.tumor = 1,
  subset.for.speed = NULL,
  positions = GenomicRanges::GRanges(),
  verbose = TRUE
)
}
\arguments{
\item{.Object}{ExomeDepth object}

\item{data}{Data frame containing potential covariates.}

\item{test}{Numeric, vector of counts for the test sample.}

\item{reference}{Numeric, vector of counts for the reference sample.}

\item{formula}{Linear model to be used when fitting the data.}

\item{phi.bins}{Numeric, defaults to 1. Number of discrete bins for the over-dispersion parameter phi, depending on read depth.
Do not modify this parameter for the standard use of ExomeDepth.}

\item{prop.tumor}{Numeric, defaults to 1. For the somatic variant calling, this assesses the proportion of the test sample data originating from the tumour.
Do not modify this parameter for the standard use of ExomeDepth.}

\item{subset.for.speed}{Numeric, defaults to NULL. If non-null, this sets the number of data points to be used for an accelerated fit of the data.}

\item{positions}{Optional GRanges argument specifying the positions of the exons (or DNA regions) where the reads were counted for test and reference.}

\item{verbose}{Logical, controls the output level.}
}
\value{
An ExomeDepth object, which contains the CNV calls after running a Viterbi algorithm.
}
\description{
Builds an exomeDepth object from test and reference vectors
}
\examples{

data(ExomeCount)  #pick an example count file
small_count <- ExomeCount[1:100, ]  #reduce the size for speedy computations

## remove exons without data below
small_count <- small_count[ small_count$Exome2 + small_count$Exome3 > 0, ]

example_object <- new('ExomeDepth', test = small_count$Exome2,
                                    reference = small_count$Exome3,
                                    formula = 'cbind(test, reference) ~ 1')
print(example_object)
print( mean(example_object@expected)) ## proportion of reads expected to match the test set
}
