\name{countBamInGRanges.exomeDepth}
\alias{countBamInGRanges.exomeDepth}
\title{
  Compute read count data from BAM files.
}
\description{
  Parses a BAM file and count reads that are located within a target
  region defined by a GenomicRanges object.
}
\usage{
countBamInGRanges.exomeDepth(bam.file, granges, min.mapq = 1, read.width = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam.file}{BAM file to be parsed}
  \item{granges}{Genomic ranges object defining the bins}
  \item{min.mapq}{Minimum mapping quality}
  \item{read.width}{Window around the target region to count a read
    within the bin}
}
\details{
  Largely derived from its equivalent function in the exomeCopy package.
}
\value{
  A GRanges object with count data.
}
