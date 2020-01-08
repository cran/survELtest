#' The pointwise EL testing
#' 
#' @description \code{ptwiseELtest} gives pointwise EL testing to compare the survival curves at
#' each time point.
#' @name ptwiseELtest
#' @param formula a formula object with a \code{Surv} object as the response on the left of the \code{~} operator and 
#' the grouping variable as the term on the right. The \code{Surv} object involves two variables: the observed survival 
#' and censoring times, and the censoring indicator, which takes a value of \eqn{1} if the observed time is uncensored 
#' and \eqn{0} otherwise. The grouping variable takes different values for different groups. 
#' @param data an optional data frame containing the variables in the \code{formula}: the observed survival and censoring times,
#'  the censoring indicator, and the grouping variable. If not found in \code{data}, the variables in the \code{formula} should 
#'  be already defined by the user or in attached \code{R} objects. The default is the data frame with three columns of 
#'  variables taken from the formula: column 1 contains the observed survival and censoring times, column 2 the censoring 
#'  indicator, and column 3 the grouping variable.
#' @param group_order a \eqn{k}-vector containing the values of the grouping variable, with the \eqn{j}-th element being the group
#' hypothesized to have the \eqn{j}-th highest survival rates, \eqn{j=1,\ldots,k}. The default is the vector of sorted grouping variables.
#' @param t1  the first endpoint of a prespecified time interval, if any, to which the comparison of the survival functions is restricted. 
#' The default value is \eqn{0}.
#' @param t2 the second endpoint of a prespecified time interval, if any, to which the comparison of the survival 
#' functions is restricted. The default value is \eqn{\infty}.
#' @param sided \eqn{2} if two-sided test, and \eqn{1} if one-sided test. The default value is \eqn{2}.
#' @param nboot the number of bootstrap replications in calculating critical values for the tests. 
#' The default value is \eqn{1000}.
#' @param alpha the pre-specified significance level of the tests. The default value is \eqn{0.05}.
#' @param seed  the seed for the random number generator in \code{R}, for generating bootstrap samples needed 
#' to calculate the critical values for the tests. The default value is \eqn{1011}.
#' @param nlimit a number used to calculate \code{nsplit=} \eqn{m}\code{/nlimit}, the number of 
#' parts into which the calculation of the \code{nboot} bootstrap replications is split. 
#' The use of this variable can make computation faster when the number of time points \eqn{m} is large. 
#' The default value for \code{nlimit} is 200.
#' @return \code{ptwiseELtest} returns a \code{ptwiseELtest} object, a list with 12 elements:
#' \itemize{
#'    \item \code{call} the function call
#'    \item \code{result_dataframe} a dataframe with \code{time_pts} in the first column, \code{decision} 
#'    in the second column, \code{stat_ptwise} in the third column and \code{critval_ptwise} in the fourth column.
#'    \itemize{
#'          \item \code{time_pts} a vector containing the observed uncensored time points at which the
#'          Kaplan—Meier estimate is positive and less than \eqn{1} for each sample.
#'          \item \code{decision} a vector containing the decisions of the pointwise EL tests at \code{time_pts}.
#'          The decision at each of \code{time_pts} is \eqn{1} for rejection of the null hypothesis that the survival
#'          functions are the same at the specific time point, and \eqn{0} otherwise.
#'          \item \code{stat_ptwise} a vector containing the pointwise EL statistics at \code{time_pts}.
#'          \item \code{critval_ptwise} a vector containing the critical values for pointwise EL testing
#'          at \code{time_pts}.
#'    }
#'    \item \code{formula} the value of the input argument of ptwiseELtest
#'    \item \code{data} the value of the input argument of ptwiseELtest
#'    \item \code{group_order} the value of the input argument of ptwiseELtest
#'    \item \code{t1} the value of the input argument of ptwiseELtest
#'    \item \code{t2} the value of the input argument of ptwiseELtest
#'    \item \code{sided} the value of the input argument of ptwiseELtest
#'    \item \code{nboot} the value of the input argument of ptwiseELtest
#'    \item \code{alpha} the value of the input argument of ptwiseELtest
#'    \item \code{seed} the value of the input argument of ptwiseELtest
#'    \item \code{nlimit} the value of the input argument of ptwiseELtest
#' }
#' Methods defined for \code{ptwiseELtest} objects are provided for \code{print} and \code{summary}.
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{supELtest}}, \code{\link{nocrossings}}, \code{\link{print.ptwiseELtest}}, \code{\link{summary.ptwiseELtest}}
#' @examples
#' library(survELtest)
#' ptwiseELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
#'     hepatitis$group, sided = 1)
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
ptwiseELtest = function(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200) {
    #check and process inputs
    call = match.call()
    options(warn = -1)
    checkInput_formula_data("ptwiseELtest", formula, data)
    group_order = processInput_group_order(formula, group_order)
    formula = processInput_formula(formula, data)
    data = processInput_data(formula, data)
    checkInputs(data, group_order, t1, t2, sided, nboot, wt = NULL, alpha, seed, nlimit)
    options(warn = 0)
    
    k = length(unique(data$group))
    if (k == 2){
        at_ts = neg2ELratio(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
        boot_ptw = apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 2, quantile, 1 - alpha)
    }else{
        at_ts = teststat(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
        boot_ptw = apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts), 2, quantile, 1 - alpha)
    }
    
    result_dataframe = data.frame(time_pts       = at_ts$ts,
                                  decision       = as.numeric(at_ts$neg2ELratio_at_ts > boot_ptw),
                                  stat_ptwise    = at_ts$neg2ELratio_at_ts,
                                  critval_ptwise = boot_ptw)

    result = list(call = call, result_dataframe = result_dataframe, formula = formula, data = data, group_order = group_order, t1 = t1, t2 = t2, sided = sided, nboot = nboot, alpha = alpha, seed = seed, nlimit = nlimit)
    attr(result, "class") = "ptwiseELtest"
    return(result)
}
