#' Construct a matrix of binary variables for calibration
#'
#' Construct a matrix of binary variables for calibration of sample weights
#' according to known marginal population totals. The following methods
#' are implemented:
#' \itemize{
#' \item \code{calibVars.default(x)}
#' \item \code{calibVars.matrix(x)}
#' \item \code{calibVars.matrix(x)}
#' \item \code{calibVars.data.frame(x)}
#' }
#'
#' @name calibVars
#' @param x a vector that can be interpreted as factor, or a matrix or
#' \code{data.frame} consisting of such variables.
#' @return A matrix of binary variables that indicate membership to the
#' corresponding factor levels.
#' @author Bernhard Meindl and Andreas Alfons
#' @references
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. doi: 10.18637/jss.v079.i10
#' @export calibVars
#' @seealso \code{\link{calibSample}}
#' @keywords survey
#' @examples
#' data(eusilcS)
#' # default method
#' \dontrun{
#' aux <- calibVars(eusilcS$rb090)
#' head(aux)
#' # data.frame method
#' aux <- calibVars(eusilcS[, c("db040", "rb090")])
#' head(aux)
#' }
calibVars <- function(x) UseMethod("calibVars")

#' @export
calibVars.default <- function(x) {
  if(length(x) == 0) matrix(integer(), 0, 0)
  x <- as.factor(x)
  levs <- levels(x)
  res <- sapply(levels(x), function(l) as.integer(x == l))
  rownames(res) <- names(x)
  res
}

#' @export
calibVars.matrix <- function(x) calibVars(as.data.frame(x))

#' @export
calibVars.data.frame <- function(x) {
  res <- lapply(x, calibVars)  # list of matrices for each variable
  res <- mapply(function(x, nam) {
    colnames(x) <- paste(nam, colnames(x), sep=".")
    x
  }, res, names(x), SIMPLIFY=FALSE)
  res <- do.call("cbind", res)
  rownames(res) <- row.names(x)
  res
}
