% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0classes.R
\docType{class}
\name{simPopObj-class}
\alias{simPopObj-class}
\alias{show,simPopObj-method}
\title{Class \code{"simPopObj"}}
\description{
An object that is used throughout the package containing information on the
sample (in slot \code{sample}), the population (slot \code{pop}) and
optionally some margins in form of a table (slot \code{table}).
}
\section{Objects from the Class}{
 Objects are automatically created in
function \code{\link{simStructure}}.
}

\examples{

showClass("simPopObj")

## show method: generate an object of class simPop first
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
class(eusilcP)
## shows some basic information:
eusilcP

}
\seealso{
\code{\linkS4class{dataObj}}
}
\author{
Bernhard Meindl and Matthias Templ
}
\keyword{classes}
