#' Summary function for supELtest object
#' 
#' @description Returns a list containing the maximally selected EL statistics, the critical value based on bootstrap, 
#' and the p-value of the test.
#' @name summary.supELtest
#' @param object the result of a call to the \code{supELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @return \code{summary.supELtest} returns a list with following components:
#' \itemize{
#'    \item \code{call} the statement used to create the \code{supELtest} object
#'    \item \code{teststat} the resulting integrated EL statistics
#'    \item \code{critval} the critical value based on bootstrap
#'    \item \code{pvalue} the p-value of the test
#'    \item \code{sided} the value of the input argument of supELtest
#'    \item \code{alpha} the value of the input argument of supELtest
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{print.supELtest}}
#' @examples
#' library(survELtest)
#' nocrossings(survival::Surv(hepatitis$time, hepatitis$censor)~
#'     hepatitis$group, sided = 1)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## nocrossings(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## Decision = 1
#' 
#' ## A decision value of 1 means the case of crossings or alternative orderings among the 
#' ## survival functions is excluded. Thus, we can proceed to the one-sided test.
#' 
#' result = supELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
#'              hepatitis$group, sided = 1)
#' summary(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## supELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## One-sided maximally selected EL test statistic = 10.36, p = 0.006,
#' ## critical value based on bootstrap = 6.289 at a significance level of 0.05
#' @export
## generic function to summary an S3 object of class "supELtest"
summary.supELtest = function(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(object$sided == 1){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nOne-sided maximally selected EL test statistic = ", format(object$teststat, digits = digits), ", p = ", format(object$pvalue, digits = digits), ",", sep = "")
      cat("\ncritical value based on bootstrap = ", format(object$critval, digits = digits), " at a significance level of ", format(object$alpha, digits = digits), "\n", sep = "")
    }
    result = list(call = object$call, teststat = object$teststat, critval = object$critval, pvalue = object$pvalue, sided = object$sided, alpha = object$alpha)
  }else if(object$sided == 2){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nTwo-sided maximally selected EL test statistic = ", format(object$teststat, digits = digits), ", p = ", format(object$pvalue, digits = digits), ",", sep = "")
      cat("\ncritical value based on bootstrap = ", format(object$critval, digits = digits),  " at a significance level of ", format(object$alpha, digits = digits), "\n", sep = "")
    }
    result = list(call = object$call, teststat = object$teststat, critval = object$critval, pvalue = object$pvalue, sided = object$sided, alpha = object$alpha)
  }
  invisible(result)
}
