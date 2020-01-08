#' Print a ptwiseELtest object
#' 
#' @description Print some summary statistics for the observed uncensored time points, and the decisions, 
#' statistics, and critical values of the pointwise EL tests at those time points.
#' @name print.ptwiseELtest
#' @param x the result of a call to the \code{ptwiseELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @seealso \code{\link{hepatitis}}, \code{\link{ptwiseELtest}}, \code{\link{summary.ptwiseELtest}}
#' @examples
#' library(survELtest)
#' result = ptwiseELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
#'              hepatitis$group, sided = 1)
#' print(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## ptwiseELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## Range of time_pts is from 5.2 to 153.1
#' ## 30 out of 45 decisions are 1, the other 15 decisions are 0
#' ## -----
#' ## Summary of stat_ptwise:
#' ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' ##   0.000   2.293   3.694   4.263   6.288  10.360 
#' ## -----
#' ## Summary of critval_ptwise:
#' ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' ##   2.117   2.346   2.483   2.509   2.669   2.951
#' @export
## generic function to print an S3 object of class "ptwiseELtest"
print.ptwiseELtest = function(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(quiet == FALSE){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nRange of time_pts is from ", format(min(x$result_dataframe$time_pts), digits = digits)," to ", format(max(x$result_dataframe$time_pts), digits = digits), sep = "")
    cat("\n", sum(x$result_dataframe$decision == 1), " out of ", nrow(x$result_dataframe), " decisions are 1, the other ", sum(x$result_dataframe$decision == 0), " decisions are 0", sep = "")
    cat("\n-----\nSummary of stat_ptwise:\n")
    print(summary(x$result_dataframe$stat_ptwise, digits = digits), digits = digits)
    cat("-----\nSummary of critval_ptwise:\n")
    print(summary(x$result_dataframe$critval_ptwise, digits = digits), digits = digits)
  }
  result = list(call = x$call, time_pts = x$result_dataframe$time_pts, decision = x$result_dataframe$decision, stat_ptwise = x$result_dataframe$stat_ptwise, critval = x$result_dataframe$critval)
  invisible(result)
}

