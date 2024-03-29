% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spBwplot.R
\name{spBwplot}
\alias{spBwplot}
\alias{panelSpBwplot}
\alias{getBwplotStats}
\alias{prepBwplotStats.data.frame}
\alias{prepBwplotStats.default}
\alias{prepBwplotStats}
\title{Weighted box plots}
\usage{
spBwplot(
  inp,
  x,
  cond = NULL,
  horizontal = TRUE,
  coef = 1.5,
  zeros = TRUE,
  minRatio = NULL,
  do.out = FALSE,
  ...
)

panelSpBwplot(x, y, coef = 1.5, zeros = TRUE, ratio, outliers, subscripts, ...)

getBwplotStats(x, weights = NULL, cond = NULL, data, ..., name = "")

prepBwplotStats(x, w, ..., name = "")

\method{prepBwplotStats}{data.frame}(x, w, ..., name = "")

\method{prepBwplotStats}{default}(x, w, ..., name = "")
}
\arguments{
\item{inp}{an object of class \code{\linkS4class{simPopObj}} containing
survey sample and synthetic population data.}

\item{x}{a character vector specifying the columns of data available in the
sample and the population (specified in input object 'inp') to be plotted.}

\item{cond}{an optional character vector (of length 1, if used) specifying
the conditioning variable.}

\item{horizontal}{a logical indicating whether the boxes should be
horizontal or vertical.}

\item{coef}{a numeric value that determines the extension of the whiskers.}

\item{zeros}{a logical indicating whether the variables specified by
\code{x} are semi-continuous, i.e., contain a considerable amount of zeros.
If \code{TRUE}, the box widths correspond to the proportion of non-zero data
points and the (weighted) box plot statistics are computed for these
non-zero data points only.}

\item{minRatio}{a numeric value in \eqn{(0,1]}; if \code{zeros} is
\code{TRUE}, the boxes may become unreadable for a large proportion of
zeros. In such a case, this can be used to specify a minimum ratio for the
box widths. Variable box widths for semi-continuous variables can be
suppressed by setting this value to 1.}

\item{do.out}{a logical indicating whether data points that lie beyond the
extremes of the whiskers should be plotted. Note that this is \code{FALSE}
by default.}

\item{\dots}{further arguments to be passed to
\code{\link[lattice:xyplot]{bwplot}}.}
}
\value{
An object of class \code{"trellis"}, as returned by
\code{\link[lattice:xyplot]{bwplot}}.
}
\description{
Produce box-and-whisker plots of continuous or semi-continuous variables,
possibly broken down according to conditioning variables and taking into
account sample weights.
}
\details{
Missing values are ignored for producing box plots and weights are directly
extracted from the input object \code{inp}.
}
\examples{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
\donttest{
## approx. 20 seconds computation time
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
  strata="db040", weight="db090")
simPop <- simStructure(data=inp, method="direct",
  basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))

# multinomial model with random draws
eusilcM <- simContinuous(simPop, additional="netIncome",
  regModel  = ~rb090+hsize+pl030+pb220a+hsize,
  upper=200000, equidist=FALSE, nr_cpus=1)

# plot results
spBwplot(eusilcM, x="netIncome", cond=NULL)
spBwplot(eusilcM, x="netIncome", cond="rb090", layout=c(1,2))
}

}
\references{
A. Alfons, M. Templ (2011) Simulation of close-to-reality population data for household surveys with application to EU-SILC. 
\emph{Statistical Methods & Applications}, \strong{20} (3), 383--407. \doi{10.1007/s10260-011-0163-2}
}
\seealso{
\code{\link{spBwplotStats}}, \code{\link[lattice:xyplot]{bwplot}}
}
\author{
Andreas Alfons and Bernhard Meindl
}
\keyword{hplot}
\keyword{internal}
