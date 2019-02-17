#' The pointwise EL testing
#' 
#' @description \code{ptwiseELtest} gives pointwise EL testing to compare the survival curves at
#' each time point.
#' @name ptwiseELtest
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
#' @return \code{ptwiseELtest} returns a list with four elements:
#' \itemize{
#'    \item \code{time_pts} a vector containing the observed uncensored time points at which the
#'    Kaplan—Meier estimate is positive and less than \eqn{1} for each sample.
#'    \item \code{decision} a vector containing the decisions of the pointwise EL tests at \code{time_pts}.
#'    The decision at each of \code{time_pts} is \eqn{1} for rejection of the null hypothesis that the survival
#'    functions are the same at the specific time point, and \eqn{0} otherwise.
#'    \item \code{stat_ptwise} a vector containing the pointwise EL statistics at \code{time_pts}.
#'    \item \code{critval_ptwise} a vector containing the critical values for pointwise EL testing
#'    at \code{time_pts}.
#' }
#' @seealso \code{\link{intELtest}}, \code{\link{supELtest}}, \code{\link{threearm}}
#' @examples
#' \dontrun{
#' library(survELtest)
#' ptwiseELtest(threearm[1:30,],group_order=c(3,2,1),sided=1)
#' ## It produces the estimates at 2 observed uncensored time points.
#' }
#' @export
#' @importFrom stats quantile

ptwiseELtest = function(data, group_order = sort(unique(data[,3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200) {
      k = length(unique(data[,3]))
    
    if (k != length(unique(group_order))){
        stop("Parameter \"group_order\" doesn't match the actual number of groups in your input data.")
    }
    if (k == 2){
        at_ts = neg2ELratio(data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
        boot_ptw = apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 2, quantile, 1 - alpha)
    }else{
        at_ts = teststat(data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
        boot_ptw = apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts), 2, quantile, 1 - alpha)
    }
    
    return(list(time_pts       = at_ts$ts,
    decision       = as.numeric(at_ts$neg2ELratio_at_ts > boot_ptw),
    stat_ptwise    = at_ts$neg2ELratio_at_ts,
    critval_ptwise = boot_ptw))
}
