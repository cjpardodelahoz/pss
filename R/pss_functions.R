###### KL-Weighted MPD (klMPD) for rows of interaction matrix ######
#' KL-Weighted MPD (klMPD)
#'
#' KL-weighted version of the mean pairwise phylogenetic distance (klMPD). This is a metric of phylogenetic diversity that is corrected for partner availability in the case of interaction networks, and species regional availability in the case of community data matrices (CDM).
#'
#' @param samp Required. A matrix or dataframe with focal taxa (interaction networks) or sites/communities (CDM) as rows, and partner taxa as columns. Elements of the matrix are interaction frequencies (for interaction matrices), or site-specific species abundances (CDM)
#' @param dis Required. A symmetric matrix with the phylogentic distances among the taxa in the columns of samp. If you have a phylogenetic tree, you can generate this with cophenetic(your.tree). All taxa in the columns should be part of the matrix, but the matrix may have taxa not included in samp.
#' @param q Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the columns of samp. Default is to estimate these values from interaction frequencies or relative species abundances summed across the matrix.
#' @return A named numerical vector with the klMPD estimated for each taxon or community/site in the rows of samp.
#' @author Carlos J. Pardo De la Hoz
#' @export
kl.mpd <- function (samp, dis, q)
{
  if (missing(q)) {
    q <- colSums(samp)/sum(samp)
  }
  dist.mat <- as.matrix(dis)
  temp.function <- function(samp)
  {
    ##Get names of partners of a single species i
    com.names <- names(samp[samp > 0])
    ##get interaction frequences of pecies i with all its partners
    com <- t(as.matrix(samp[samp > 0]))
    ##Calculate P'ik by dividing the com vector by the com rowSums
    Pmn <- as.vector(com/rowSums(com))
    ##Calculate qk for the species that associate with m
    com.q <- q[com.names]
    ##Calculate KL factos for the species that associate with i
    KL <- as.matrix(Pmn*log(Pmn/com.q))
    ##Subet the KL matrix to include only partners with KL >0
    KL.included <- as.matrix(KL[KL[,1] > 0, 1])
    ##Get the names of the partners with KL > 0
    KL.included.names <- names(KL[KL[,1] > 0, 1])
    ##Make phylogenetic distance matrix for partners of focal species
    com.dist.mat <- dist.mat[KL.included.names, KL.included.names]
    ##Calculate the product of KL factors of all partners of the focal species
    KL.products <-  KL.included%*%t(KL.included)
    ##calculate mean of the community phylogenetic distance matrix weighted by the KL products
    weighted.mean(com.dist.mat, KL.products)
  }
  ##apply to the entire community dataset
  apply(samp, MARGIN = 1, temp.function)
}

###### SES KL-weighted MPD function (PSS) for rows of interaction matrix ######
#' PSS for rows of a matrix
#'
#' PSS index for the taxa or sites/communities in the rows of an interaction matrix or community data matrix (CDM). This is a standardized effect size of klMPD.
#'
#' @param samp Required. A matrix or dataframe with focal taxa (interaction networks) or sites/communities (CDM) as rows, and partner taxa as columns. Elements of the matrix are interaction frequencies (for interaction matrices), or site-specific species abundances (CDM)
#' @param phyl Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the columns of samp. This tree may have taxa not included in samp.
#' @param q Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the columns of samp. Default is to estimate these values from interaction frequencies or relative species abundances summed across the matrix.
#' @param runs Optional. The number of randomizations to estimate the null distribution of klMPD values used to calculate the SES. Default is 999 runs.
#' @return A dataframe with results for each species or site/community in the rows of samp.
#' @author Carlos J. Pardo De la Hoz
#' @export
pss.r <- function (samp, phyl, q, runs = 999)
{
  ##if q is not supplied by the user, calculate it from interaction freqs
  if (missing(q)) {
    q <- colSums(samp)/sum(samp)
  }
  ##make sure the interaction matrix is of class matrix
  samp <- as.matrix(samp)
  ##get phylogenetic distance matrix
  dis <- as.matrix(cophenetic(phyl))
  ##get the observed klMPD values
  kl.mpd.obs <- kl.mpd(samp, dis, q)
  ##get the null distribution of klMPD by suffling taxa
  kl.mpd.rand <- t(replicate(runs, kl.mpd(samp, picante::taxaShuffle(dis), q)))
  ## get the mean and sd of the null distribution of klMPD values
  kl.mpd.rand.mean <- apply(X = kl.mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
  kl.mpd.rand.sd <- apply(X = kl.mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
  ##get the names of the species for which the null mena is 0.
  spp.one.partner <- setdiff(names(kl.mpd.rand.mean), names(kl.mpd.rand.mean[kl.mpd.rand.mean > 0]))
  ##if there are species with only one partner, then implement the following coercion
  if (length(spp.one.partner) >= 1) {
    for (sp in spp.one.partner) {
      ##get the name of the sole partner of the species with one partner
      row <- samp[sp,]
      partner <- names(row[row == max(row)])
      partner <- partner[1]
      ##get the partner column divided by 2
      new.col <- samp[,partner]/2
      ##add the new column to the new matrix
      temp.samp <- cbind(samp, new.col)
      ##replace the original partner column for the one divided by 2
      temp.samp[,partner] <- samp[,partner]/2
      ##duplicate q vector
      q.partner <- q
      ##divide the q of the partner by 2
      q.partner[partner] <- q.partner[partner]/2
      ##add the new column to the q vector as the q of the original partner divided by 2
      new.col <- q.partner[partner]
      names(new.col) <- c("new.col")
      q.partner <- c(q.partner, new.col)
      ##get the partner column from the distance matrix
      new.col <- dis[partner,]
      ##add new.col to the distance matrix
      temp.dis <- cbind(dis, new.col)
      ##replace distance by minimum divided by 2
      temp.dis[partner, "new.col"] <- min(dis[dis > 0])/2
      ##add new colum as a row with 0 as distance to itself
      new.col <- c(temp.dis[,"new.col"], temp.dis[partner,partner])
      temp.dis <- rbind(temp.dis, new.col)
      ##get new obs values with temp data
      temp.obs <- kl.mpd(temp.samp, temp.dis, q.partner)
      kl.mpd.obs[sp] <- temp.obs[sp]
      ##get null distribution witht the temp data
      temp.kl.mpd.rand <- t(replicate(999, kl.mpd(temp.samp, picante::taxaShuffle(temp.dis), q.partner)))
      temp.kl.mpd.rand.mean <- apply(X = temp.kl.mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
      temp.kl.mpd.rand.sd <- apply(X = temp.kl.mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
      ##replace the mean and sd values on the vectors that will go into the z-score
      kl.mpd.rand.mean[sp] <- temp.kl.mpd.rand.mean[sp]
      kl.mpd.rand.sd[sp] <- temp.kl.mpd.rand.sd[sp]
    }
  }
  pss <- (kl.mpd.obs - kl.mpd.rand.mean)/kl.mpd.rand.sd
  kl.mpd.obs.rank <- apply(X = rbind(kl.mpd.obs, kl.mpd.rand), MARGIN = 2, FUN = rank)[1, ]
  kl.mpd.obs.rank <- ifelse(is.na(kl.mpd.rand.mean), NA, kl.mpd.obs.rank)
  data.frame(ntaxa = specnumber(samp), kl.mpd.obs, kl.mpd.rand.mean, kl.mpd.rand.sd, kl.mpd.obs.rank, pss, pss.p = kl.mpd.obs.rank/(runs +
                                                                                                                                 1), runs = runs, row.names = row.names(samp))
}

###### SES KL-weighted MPD (PSS) for rows and columns of interaction matrix ######
#' PSS for rows and columns
#'
#' PSS for the rows and columns of an interaction matrix. PSS is a standardized effect size (SES) of klMPD.
#'
#' @param samp Required. A matrix. Elements of the matrix are interaction frequencies.
#' @param phyl.rows Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the rows of samp. This tree may have taxa not included in samp.
#' @param phyl.cols Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the columns of samp. This tree may have taxa not included in samp.
#' @param q.rows Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the rows of samp. Default is to estimate these values from interaction frequencies.
#' @param q.cols Optional. A named numerical vector with the relative availabilities (i.e., sum(q) = 1) of the taxa in the columns of samp. Default is to estimate these values from interaction frequencies.
#' @return A list with two dataframes (rows and columns) with results for each species in samp.
#' @author Carlos J. Pardo De la Hoz
#' @export
pss.rc <- function(samp, phyl.rows, phyl.cols, q.rows, q.cols) {
  ##if q.cols is not supplied by the user, calculate it from interaction freqs
  if (missing(q.cols)) {
    q.cols <- colSums(samp)/sum(samp)
  }
  rows <- pss.r(samp, phyl.cols, q.cols)
  samp.cols <- as.data.frame(t(samp))
  ##if q.rows is not supplied by the user, calculate it from interaction freqs
  if (missing(q.rows)) {
    q.rows <- colSums(samp.cols)/sum(samp.cols)
  }
  cols <- pss.r(samp.cols, phyl.rows, q.rows)
  out <- list(rows, cols)
  out
}

###### SES abundance weighted MPD (wMPD) for rows and columns of interaction matrix ######
#' Standardized Effect Size of wMPD for rows and columns
#'
#' Standardized Effect Size of the interaction-frequency-weighted mean pairwise phylogenetic distance (wMPD) for the taxa in the rows and columns of an interaction matrix. This is simply an extension of picante's function ses.mpd
#'
#' @param samp Required. A matrix. Elements of the matrix are interaction frequencies.
#' @param phyl.rows Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the rows of samp. This tree may have taxa not included in samp.
#' @param phyl.cols Required. An object of class phylo. This is a phylogenetic tree that includes all the taxa in the columns of samp. This tree may have taxa not included in samp.
#' @return A list with two dataframes (rows and columns) with results for each taxon in samp. See ?ses.mpd in Picante for more details.
#' @author Carlos J. Pardo De la Hoz
#' @export
ses.mpd.rc <- function(samp, phyl.rows, phyl.cols) {
  dis.rows <- cophenetic(phyl.rows)
  dis.cols <- cophenetic(phyl.cols)
  rows <- picante::ses.mpd(samp, dis.cols, abundance.weighted = TRUE, null.model = "taxa.labels")
  samp.cols <- as.data.frame(t(samp))
  cols <- picante::ses.mpd(samp.cols, dis.rows, abundance.weighted = TRUE, null.model = "taxa.labels")
  out <- list(rows, cols)
  out
}
