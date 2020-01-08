#' Summary function for nocrossings object
#' 
#' @description Returns the decision for rejection of the null hypothesis that there are crossings or alternative orderings among the survival functions.
#' @name summary.nocrossings
#' @param object the result of a call to the \code{nocrossings} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @return \code{summary.nocrossings} returns a list with following components:
#' \itemize{
#'    \item \code{call} the statement used to create the \code{nocrossings} object
#'    \item \code{decision} \eqn{1} for rejection of the null hypothesis that there are crossings or alternative orderings among the survival functions, and \eqn{0} otherwise
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{nocrossings}}, \code{\link{print.nocrossings}}
#' @examples
#' library(survELtest)
#' result = nocrossings(survival::Surv(hepatitis$time, hepatitis$censor)~
#'           hepatitis$group, sided = 1)
#' summary(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## nocrossings(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## Decision = 1
#' @export
## generic function to summary an S3 object of class "nocrossings"
summary.nocrossings = function(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
    if(quiet == FALSE){
        cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
        cat("\nDecision = ", object$decision, "\n", sep = "")
    }
    result = list(call = object$call, decision = object$decision)
    invisible(result)
}

