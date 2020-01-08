#' Print an intELtest object
#' 
#' @description Print the integrated EL statistics and the p-value of the test.
#' @name print.intELtest
#' @param x the result of a call to the \code{intELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{summary.intELtest}}
#' @examples
#' library(survELtest)
#' result = intELtest(survival::Surv(hepatitis$time, hepatitis$censor) ~ hepatitis$group)
#' print(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## intELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group)
#' ## 
#' ## Two-sided integrated EL test statistic = 1.42, p = 0.007
#' @export
## generic function to print an S3 object of class "intELtest"
print.intELtest = function(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
    if(x$sided == 1){
        if(quiet == FALSE){
            cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
            cat("\nOne-sided integrated EL test statistic = ", format(x$teststat, digits = digits), ", p = ", format(x$pvalue, digits = digits), "\n", sep = "")
        }
        result = list(call = x$call, teststat = x$teststat, pvalue = x$pvalue, sided = x$sided)
    }else if(x$sided == 2){
        if(quiet == FALSE){
            cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
            cat("\nTwo-sided integrated EL test statistic = ", format(x$teststat, digits = digits), ", p = ", format(x$pvalue, digits = digits), "\n", sep = "")
        }
        result = list(call = x$call, teststat = x$teststat, pvalue = x$pvalue, sided = x$sided)
    }
    invisible(result)
}

