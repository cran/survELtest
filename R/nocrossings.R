#' The test that excludes the possibility of crossings or alternative orderings among the survival functions
#' 
#' @description The test \code{nocrossings} should be used before one-sided testing via \code{\link{intELtest}}
#' or \code{\link{supELtest}} to exclude the possibility of crossings or alternative orderings among the survival functions.
#' @name nocrossings
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
#' @return \code{nocrossings} returns a \code{nocrossings} object, a list with 12 elements:
#' \itemize{
#'    \item \code{call} the function call
#'    \item \code{decision} \eqn{1} for rejection of the null hypothesis that there are crossings or alternative orderings among the survival functions, and \eqn{0} otherwise
#'    \item \code{formula} the value of the input argument of nocrossings
#'    \item \code{data} the value of the input argument of nocrossings
#'    \item \code{group_order} the value of the input argument of nocrossings
#'    \item \code{t1} the value of the input argument of nocrossings
#'    \item \code{t2} the value of the input argument of nocrossings
#'    \item \code{sided} the value of the input argument of nocrossings
#'    \item \code{nboot} the value of the input argument of nocrossings
#'    \item \code{alpha} the value of the input argument of nocrossings
#'    \item \code{seed} the value of the input argument of nocrossings
#'    \item \code{nlimit} the value of the input argument of nocrossings
#' }
#' Methods defined for \code{nocrossings} objects are provided for \code{print} and \code{summary}.
#' @references
#' \itemize{
#'    \item H. Chang, I.W. McKeague, "Empirical likelihood based tests
#'      for stochastic ordering under right censorship," \emph{Electronic Journal of Statistics},
#'      Vol. 10, No. 2, pp. 2511-2536 (2016).
#'    \item H. Chang, I.W. McKeague, "Nonparametric testing for multiple survival functions with 
#'    non-inferiority margins," \emph{Annals of Statistics}, Vol. 47, No. 1, pp. 205-232, (2019).
#' }
#' @seealso  \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{supELtest}}, \code{\link{ptwiseELtest}}, \code{\link{print.nocrossings}}, \code{\link{summary.nocrossings}}
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
#' @export
nocrossings = function(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200) {
    #check and process inputs
    call = match.call()
    options(warn = -1)
    checkInput_formula_data("nocrossings", formula, data)
    group_order = processInput_group_order(formula, group_order)
    formula = processInput_formula(formula, data)
    data = processInput_data(formula, data)
    checkInputs(data, group_order, t1, t2, sided, nboot, wt = NULL, alpha, seed, nlimit)
    options(warn = 0)
    
    k = length(unique(data$group))
    if (k == 2){
        at_ts = neg2ELratio(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }else{
        at_ts = teststat(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }
    
    if (is.null(at_ts)){
      decision = NULL
    }else{
      decision = at_ts$test_nocross
    }
    
    result = list(call = call, decision = decision, formula = formula, data = data, group_order = group_order, t1 = t1, t2 = t2, sided = sided, nboot = nboot, alpha = alpha, seed = seed, nlimit = nlimit)
    attr(result, "class") = "nocrossings"
    return(result)
}
