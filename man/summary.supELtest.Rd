% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.supELtest.R
\name{summary.supELtest}
\alias{summary.supELtest}
\title{Summary function for supELtest object}
\usage{
\method{summary}{supELtest}(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)
}
\arguments{
\item{object}{the result of a call to the \code{supELtest} function}

\item{digits}{significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}}

\item{quiet}{a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}}

\item{...}{for future method}
}
\value{
\code{summary.supELtest} returns a list with following components:
\itemize{
   \item \code{call} the statement used to create the \code{supELtest} object
   \item \code{teststat} the resulting integrated EL statistics
   \item \code{critval} the critical value based on bootstrap
   \item \code{pvalue} the p-value of the test
   \item \code{sided} the value of the input argument of supELtest
   \item \code{alpha} the value of the input argument of supELtest
}
}
\description{
Returns a list containing the maximally selected EL statistics, the critical value based on bootstrap, 
and the p-value of the test.
}
\examples{
library(survELtest)
nocrossings(survival::Surv(hepatitis$time, hepatitis$censor)~
    hepatitis$group, sided = 1)

## OUTPUT:
## Call:
## nocrossings(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
##     hepatitis$group, sided = 1)
## 
## Decision = 1

## A decision value of 1 means the case of crossings or alternative orderings among the 
## survival functions is excluded. Thus, we can proceed to the one-sided test.

result = supELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
             hepatitis$group, sided = 1)
summary(result)

## OUTPUT:
## Call:
## supELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
##     hepatitis$group, sided = 1)
## 
## One-sided maximally selected EL test statistic = 10.36, p = 0.006,
## critical value based on bootstrap = 6.289 at a significance level of 0.05
}
\seealso{
\code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{print.supELtest}}
}
