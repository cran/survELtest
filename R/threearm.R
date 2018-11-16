#' Time to first remission data
#'
#' @description The data frame \code{threearm} is obtained by perturbing the time-to-remission data 
#' from patients in a three-arm randomized clinical trial for the treatment of major depression. 
#' See \code{\link{nocrossings}}, \code{\link{ptwiseELtest}} and \code{\link{supELtest}} for the application.
#' @format The \code{threearm} is a data frame with 664 observations of 3 variables,
#' and has the following columns:
#' \itemize{
#'    \item \code{time} the observed times to first remission and censoring times
#'    \item \code{censor} the censoring indicator
#'    \item \code{group} the grouping variable
#' }
#' @seealso \code{\link{nocrossings}}, \code{\link{ptwiseELtest}}, \code{\link{supELtest}}

"threearm"
