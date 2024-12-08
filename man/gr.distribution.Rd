% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{gr.distribution}
\alias{gr.distribution}
\title{Visualizing the distribution of CMs across samples}
\usage{
gr.distribution(nmf_res, ...)
}
\arguments{
\item{nmf_res}{The result of function nmf which includes the coefficient matrix.}

\item{meta}{The groupings of each sample(not required). If provided, this file should contain at least two columns:"sampleID" and the groupings category.}

\item{group}{The column name for groupings in meta.}
}
\value{
The distribution of CMs across samples.
}
\description{
Visualizing the distribution of CMs across samples
}
\examples{
gr.distribution(nmf_res)
gr.distribution(nmf_res,meta=meta,group="tissue")
}
