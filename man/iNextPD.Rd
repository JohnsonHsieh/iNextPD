\name{iNextPD}
\alias{iNextPD}
\title{iNterpolation and EXTrapolation of Fiath's PD}
\usage{
  iNextPD(abun, lables, phy,
    size = NULL, endpoint = NULL, knots = 40, se = TRUE,
    nboot = 50)
}
\arguments{
  \item{abun}{a vector/data.frame/list of species abundance.}

  \item{size}{an integer vector of
  rarefaction/extrapolation sample size (number of
  individuals or sampling units), default is NULL. If
  \code{size} is not be specified, it would compute
  rarefaction/extrapolation by endpoint and knots.}

  \item{endpoint}{a integer of sample size that is the
  endpoint for rarefaction/extrapolation. Default is double
  reference sample size.}

  \item{knots}{a number of knots (say K, default is 40)
  specifying separate sample size between 1 and
  \code{endpoint}. If \code{endpoint} is smaller than
  reference sample size, then \code{iNEXT()} compute
  rarefaction part only and divided by approximately equall
  spaced K. If \code{endpoint} is larger than reference
  sample size, then \code{iNEXT()} will compute
  approximately K/2 equally spaced for rarefaction part and
  divided extrapolation part as approximately K/2 equally
  spaced between reference sample size and
  \code{endpoint}.}

  \item{se}{calculate bootstrap standard error and show 95%
  confidence interval; default is TRUE.}

  \item{nboot}{the number of bootstrap resampling times,
  default is 50.}
}
\value{
  a list of interpolation and extrapolation Faith's PD
  with specific sample size (m) and sample coverage (SC)
}
\description{
  \code{iNextPD} Estimation of interpolation and
  extrapolation of Faith's PD
}
\examples{

}

