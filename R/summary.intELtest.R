#' Summary function for intELtest object
#' 
#' @description Returns a list containing the integrated EL statistics, the critical value based on bootstrap, 
#' and the p-value of the test.
#' @name summary.intELtest
#' @param object the result of a call to the \code{intELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @return \code{summary.intELtest} returns a list with following components:
#' \itemize{
#'    \item \code{call} the statement used to create the \code{intELtest} object
#'    \item \code{teststat} the resulting integrated EL statistics
#'    \item \code{critval} the critical value based on bootstrap
#'    \item \code{pvalue} the p-value of the test
#'    \item \code{sided} the value of the input argument of intELtest
#'    \item \code{alpha} the value of the input argument of intELtest
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{print.intELtest}}
#' @examples
#' library(survELtest)
#' result = intELtest(survival::Surv(hepatitis$time, hepatitis$censor) ~ hepatitis$group)
#' summary(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## intELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group)
#' ## 
#' ## Two-sided integrated EL test statistic = 1.42, p = 0.007,
#' ## critical value based on bootstrap = 0.875 at a significance level of 0.05
#' @export
## generic function to summary an S3 object of class "intELtest"
summary.intELtest = function(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(object$sided == 1){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nOne-sided integrated EL test statistic = ", format(object$teststat, digits = digits), ", p = ", format(object$pvalue, digits = digits), ",", sep = "")
      cat("\ncritical value based on bootstrap = ", format(object$critval, digits = digits), " at a significance level of ", format(object$alpha, digits = digits), "\n", sep = "")
    }
    result = list(call = object$call, teststat = object$teststat, critval = object$critval, pvalue = object$pvalue, sided = object$sided, alpha = object$alpha)
  }else if(object$sided == 2){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nTwo-sided integrated EL test statistic = ", format(object$teststat, digits = digits), ", p = ", format(object$pvalue, digits = digits), ",", sep = "")
      cat("\ncritical value based on bootstrap = ", format(object$critval, digits = digits),  " at a significance level of ", format(object$alpha, digits = digits), "\n", sep = "")
    }
    result = list(call = object$call, teststat = object$teststat, critval = object$critval, pvalue = object$pvalue, sided = object$sided, alpha = object$alpha)
  }
  invisible(result)
}

