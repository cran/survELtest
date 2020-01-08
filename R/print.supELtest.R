#' Print a supELtest object
#' 
#' @description Print the maximally selected EL statistics and the p-value of the test.
#' @name print.supELtest
#' @param x the result of a call to the \code{supELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @seealso \code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{summary.supELtest}}
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
#' print(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## supELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~  
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## One-sided maximally selected EL test statistic = 10.36, p = 0.006
#' @export
## generic function to print an S3 object of class "supELtest"
print.supELtest = function(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(x$sided == 1){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nOne-sided maximally selected EL test statistic = ", format(x$teststat, digits = digits), ", p = ", format(x$pvalue, digits = digits), "\n", sep = "")
    }
    result = list(call = x$call, teststat = x$teststat, pvalue = x$pvalue, sided = x$sided)
  }else if(x$sided == 2){
    if(quiet == FALSE){
      cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\nTwo-sided maximally selected EL test statistic = ", format(x$teststat, digits = digits), ", p = ", format(x$pvalue, digits = digits), "\n", sep = "")
    }
    result = list(call = x$call, teststat = x$teststat, pvalue = x$pvalue, sided = x$sided)
  }
  invisible(result)
}