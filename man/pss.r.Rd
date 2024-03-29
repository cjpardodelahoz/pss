% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pss_functions.R
\name{pss.r}
\alias{pss.r}
\title{PSS for rows of a matrix}
\usage{
pss.r(samp, phyl, q, runs = 999)
}
\arguments{
\item{samp}{Required. A matrix or dataframe with focal taxa (interaction networks) or sites/communities (CDM) as rows, and partner taxa as columns. Elements of the matrix are interaction frequencies (for interaction matrices), or site-specific species abundances (CDM)}

\item{phyl}{Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the columns of samp. This tree may have taxa not included in samp.}

\item{q}{Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the columns of samp. Default is to estimate these values from interaction frequencies or relative species abundances summed across the matrix.}

\item{runs}{Optional. The number of randomizations to estimate the null distribution of klMPD values used to calculate the SES. Default is 999 runs.}
}
\value{
A dataframe with results for each species or site/community in the rows of samp.
}
\description{
PSS index for the taxa or sites/communities in the rows of an interaction matrix or community data matrix (CDM). This is a standardized effect size of klMPD.
}
\author{
Carlos J. Pardo De la Hoz
}
