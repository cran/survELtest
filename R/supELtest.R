#' The maximally selected EL test
#' 
#' @description \code{supELtest} provides the maximally selected EL statistics that
#' is better adapted at detecting local differences:
#' \deqn{\sup_{i=1,\ldots,m}\{-2\log R(t_i)\},}
#' where \eqn{R(t)} is the EL ratio that compares the survival functions at each given time \eqn{t},
#' and \eqn{ 0<t_1<\ldots<t_m<\infty} are the (ordered) observed uncensored times at which the
#' Kaplan--Meier estimate is positive and less than 1 for each sample.
#' @name supELtest
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
#' @return \code{supELtest} returns a \code{supELtest} object, a list with 14 elements:
#' \itemize{
#'    \item \code{call} the function call
#'    \item \code{teststat} the resulting integrated EL statistics
#'    \item \code{critval} the critical value based on bootstrap
#'    \item \code{pvalue} the p-value of the test
#'    \item \code{formula} the value of the input argument of supELtest
#'    \item \code{data} the value of the input argument of supELtest
#'    \item \code{group_order} the value of the input argument of supELtest
#'    \item \code{t1} the value of the input argument of supELtest
#'    \item \code{t2} the value of the input argument of supELtest
#'    \item \code{sided} the value of the input argument of supELtest
#'    \item \code{nboot} the value of the input argument of supELtest
#'    \item \code{alpha} the value of the input argument of supELtest
#'    \item \code{seed} the value of the input argument of supELtest
#'    \item \code{nlimit} the value of the input argument of supELtest
#' }
#' Methods defined for \code{supELtest} objects are provided for \code{print} and \code{summary}.
#' @references
#' \itemize{
#'    \item H. Chang, I.W. McKeague, "Empirical likelihood based tests
#'      for stochastic ordering under right censorship," \emph{Electronic Journal of Statistics},
#'      Vol. 10, No. 2, pp. 2511-2536 (2016).
#'    \item H. Chang, I.W. McKeague, "Nonparametric testing for multiple survival functions 
#'    with non-inferiority margins," \emph{Annals of Statistics}, Vol. 47, No. 1, pp. 205-232, (2019).
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{ptwiseELtest}}, \code{\link{nocrossings}}, \code{\link{print.supELtest}}, \code{\link{summary.supELtest}}
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
#' supELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
#'     hepatitis$group, sided = 1)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## supELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~   
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ## One-sided maximally selected EL test statistic = 10.36, p = 0.006
#' @export
supELtest = function(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200) {
  #check and process inputs
  call = match.call()
  options(warn = -1)
  checkInput_formula_data("supELtest", formula, data)
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
  
  if (is.null(at_ts)) return (NULL)
  critval  = at_ts$EL_SOcrit
  teststat = at_ts$suptest
  pvalue   = at_ts$p_value_suptest
  
  result = list(call = call, teststat = teststat, critval = critval, pvalue = pvalue, formula = formula, data = data, group_order = group_order, t1 = t1, t2 = t2, sided = sided, nboot = nboot, alpha = alpha, seed = seed, nlimit = nlimit)
  attr(result, "class") = "supELtest"
  return(result)
}
