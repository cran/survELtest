#' The test that excludes the possibility of crossings or alternative orderings
#' 
#' @description The test \code{nocrossings} should be used before one-sided testing via \code{\link{intELtest}}
#' or \code{\link{supELtest}} to exclude the possibility of crossings or alternative orderings.
#' @name nocrossings
#' @param data a data frame/matrix with 3 columns: column 1 contains the observed survival
#' and censoring times, column 2 the censoring indicator, and column 3 the grouping variable.
#' This is a compulsory input.
#' @param group_order a \eqn{k}-vector containing the values of the grouping variable
#' (in column 3 of the data frame/matrix), with the \eqn{j}-th element being the group
#' hypothesized to have the \eqn{j}-th highest survival rates, \eqn{j=1,\ldots,k}.
#' The default is the vector of sorted grouping variables.
#' @param t1 the first endpoint of a prespecified time interval, if any,
#' during which the comparison of the survival functions is restricted to.
#' The default value is \eqn{0}.
#' @param t2 the second endpoint of a prespecified time interval, if any,
#' during which the comparison of the survival functions is restricted to.
#' The default value is \eqn{\infty}.
#' @param sided \eqn{2} if two-sided test, and \eqn{1} if one-sided test. The default value is \eqn{2}.
#' @param nboot the number of bootstrap replications in calculating critical values for the tests.
#' The defualt value is \eqn{1000}.
#' @param alpha the pre-specified significance level of the tests. The default value is \eqn{0.05}.
#' @param seed the seed of random number generation in \code{R} for generating bootstrap samples
#' needed to calculate critical values for the tests. The default value is \eqn{1011}.
#' @param nlimit a number used to calculate \code{nsplit=} \eqn{m}\code{/nlimit}, the number of parts
#' we split the calculation of the \code{nboot} bootstrap replications  into.
#' This can make computation faster when the number of time points \eqn{m} is too large.
#' The default value for \code{nlimit} is \eqn{200}.
#' @return
#' \itemize{
#'    \item \code{decision} \eqn{1} for rejection of the null hypothesis that there are crossings or alternative orderings, and \eqn{0} otherwise
#' }
#' @references
#' \itemize{
#'    \item H. Chang, I.W. McKeague, "Empirical likelihood based tests
#'      for stochastic ordering under right censorship," \emph{Electronic Journal of Statistics},
#'      Vol. 10, No. 2, pp. 2511-2536 (2016).
#'    \item H. Chang, I.W. McKeague, "Nonparametric testing for multiple survival functions with
#'      non-inferiority margins," \emph{Annals of Statistics}, accepted (2018).
#' }
#' @seealso \code{\link{intELtest}}, \code{\link{supELtest}}, \code{\link{threearm}}
#' @examples
#' library(survELtest)
#' nocrossings(threearm[1:30,],group_order=c(3,2,1),sided=1)
#' ## OUTPUT:
#' ## $decision
#' ## [1] 1
#' @export
#' @importFrom stats quantile

nocrossings = function(data, group_order = sort(unique(data[,3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200) {
    k = length(unique(data[,3]))
    
    if (k != length(unique(group_order))){
        stop("Parameter \"group_order\" doesn't match the actual number of groups in your input data.")
    }
    if (k == 2){
        at_ts = neg2ELratio(data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }else{
        at_ts = teststat(data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }
    
    if (is.null(at_ts)) return (NULL)
    return(list(decision = at_ts$test_nocross))
}

