% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pss_functions.R
\name{kl.mpd}
\alias{kl.mpd}
\title{KL-Weighted MPD (klMPD)}
\usage{
kl.mpd(samp, dis, q)
}
\arguments{
\item{samp}{Required. A matrix or dataframe with focal taxa (interaction networks) or sites/communities (CDM) as rows, and partner taxa as columns. Elements of the matrix are interaction frequencies (for interaction matrices), or site-specific species abundances (CDM)}

\item{dis}{Required. A symmetric matrix with the phylogentic distances among the taxa in the columns of samp. If you have a phylogenetic tree, you can generate this with cophenetic(your.tree). All taxa in the columns should be part of the matrix, but the matrix may have taxa not included in samp.}

\item{q}{Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the columns of samp. Default is to estimate these values from interaction frequencies or relative species abundances summed across the matrix.}
}
\value{
A named numerical vector with the klMPD estimated for each taxon or community/site in the rows of samp.
}
\description{
KL-weighted version of the mean pairwise phylogenetic distance (klMPD). This is a metric of phylogenetic diversity that is corrected for partner availability in the case of interaction networks, and species regional availability in the case of community data matrices (CDM).
}
\author{
Carlos J. Pardo De la Hoz
}
