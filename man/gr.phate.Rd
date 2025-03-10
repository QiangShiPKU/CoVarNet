% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{gr.phate}
\alias{gr.phate}
\title{Scatter plots in PHATE basis}
\usage{
gr.phate(data, feature, ...)
}
\arguments{
\item{data}{An AnnData object after PHATE reduction.}

\item{feature}{Input vector of features.}

\item{color}{Colors to use for continous variables or categorical groups.}
}
\value{
Scatter plots in PHATE basis
}
\description{
Scatter plots in PHATE basis
}
\examples{
gr.phate(data, feature)
}
