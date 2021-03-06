% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R
\name{ggiNEXT}
\alias{ggiNEXT}
\alias{ggiNEXT.default}
\alias{ggiNEXT.iNextPD}
\title{ggplot2 extension for an iNextPD object}
\usage{
ggiNEXT(x, type = 1, se = TRUE, facet.var = "none", color.var = "site",
  grey = FALSE)

\method{ggiNEXT}{iNextPD}(x, type = 1, se = TRUE, facet.var = "none",
  color.var = "site", grey = FALSE)

\method{ggiNEXT}{default}(x, ...)
}
\arguments{
\item{x}{an \code{iNextPD} object computed by \code{\link{iNextPD}}.}

\item{type}{three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).}

\item{se}{a logical variable to display confidence interval around the estimated sampling curve.}

\item{facet.var}{create a separate plot for each value of a specified variable: 
no separation \cr (\code{facet.var="none"}); 
a separate plot for each diversity order (\code{facet.var="order"}); 
a separate plot for each site (\code{facet.var="site"}); 
a separate plot for each combination of order x site (\code{facet.var="both"}).}

\item{color.var}{create curves in different colors for values of a specified variable:
all curves are in the same color (\code{color.var="none"}); 
use different colors for diversity orders (\code{color.var="order"}); 
use different colors for sites (\code{color.var="site"}); 
use different colors for combinations of order x site (\code{color.var="both"}).}

\item{grey}{a logical variable to display grey and white ggplot2 theme.}

\item{...}{other arguments passed on to methods. Not currently used.}
}
\value{
a ggplot2 object
}
\description{
\code{ggiNEXT}: the \code{ggplot} extension for \code{\link{iNextPD}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
}
\examples{
# single-assemblage abundance data
data(bird)
bird.phy <- ade4::newick2phylog(bird$tre)
bird.lab <- rownames(bird$abun)
out1 <- iNextPD(bird$abun$North.site, bird.lab, bird.phy, 
        q=1, datatype="abundance", endpoint=400, se=TRUE)
ggiNEXT(x=out1, type=1)
ggiNEXT(x=out1, type=2)
ggiNEXT(x=out1, type=3)

\dontrun{
# single-assemblage incidence data with three orders q
out2 <- iNextPD(bird$inci$North.site, bird.lab, bird.phy, 
        q=c(0,1,2), datatype="incidence_raw", endpoint=25)
ggiNEXT(out2, se=FALSE, color.var="order")

# multiple-assemblage abundance data with three orders q
out3 <-  iNextPD(bird$abun, bird.lab, bird.phy, 
        q=c(0,1,2), datatype="abundance", endpoint=400)
ggiNEXT(out3, facet.var="site", color.var="order")
ggiNEXT(out3, facet.var="both", color.var="both")
}
}

