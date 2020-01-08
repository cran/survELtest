#' Print a nocrossings object
#' 
#' @description Returns the decision for rejection of the null hypothesis that there are crossings or alternative orderings among the survival functions.
#' @name print.nocrossings
#' @param x the result of a call to the \code{nocrossings} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @seealso \code{\link{hepatitis}}, \code{\link{nocrossings}}, \code{\link{summary.nocrossings}}
#' @examples
#' library(survELtest)
#' result = nocrossings(survival::Surv(hepatitis$time, hepatitis$censor)~
#'              hepatitis$group, sided = 1)
#' print(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## nocrossings(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## Decision = 1
#' @export
## generic function to print an S3 object of class "nocrossings"
print.nocrossings = function(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(quiet == FALSE){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nDecision = ", x$decision, "\n", sep = "")
  }
  result = list(call = x$call, decision = x$decision)
  invisible(result)
}

