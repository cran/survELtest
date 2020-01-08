#' Simulated survival data with crossing hazard functions from the piecewise exponential distribution
#'
#' @description The data frame \code{hazardcross} is obtained as follows. The survival time is generated 
#' from the piecewise exponential model displayed in the left column of Figure 1 in Chang and McKeague (2016). 
#' The censoring distributions (the same in each arm) is specified to
#' be uniform with administrative censoring at \eqn{t = 10}, and a censoring rate of \eqn{25\%} in the first group.
#' @format The \code{hazardcross} is a data frame with 100 observations of 3 variables,
#' and has the following columns:
#' \itemize{
#'    \item \code{time} the observed times to first remission and censoring times
#'    \item \code{censor} the censoring indicator
#'    \item \code{group} the grouping variable
#' }
#' @references H. Chang, I.W. McKeague, "Empirical likelihood based tests for stochastic ordering under right censorship,"
#'  \emph{Electronic Journal of Statistics}, Vol. 10, No. 2, pp. 2511-2536, (2016).

#' @seealso \code{\link{nocrossings}}, \code{\link{ptwiseELtest}}, \code{\link{supELtest}}

"hazardcross"