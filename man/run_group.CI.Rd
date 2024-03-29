% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities_statistics.R
\name{run_group.CI}
\alias{run_group.CI}
\title{Group Confidence Interval}
\usage{
run_group.CI(x, data, ci = 0.95)
}
\arguments{
\item{x}{(Required). an \code{aggregate} compatible formula.}

\item{data}{(Required). a data frame (or list) from which the variables in formula should be taken.}

\item{ci}{(optional). numeric. confidence interval to be calculated.}
}
\value{
A data frame consisting of one column for each grouping factor plus
three columns for the upper bound, mean and lower bound of the confidence interval
for each level of the grouping factor
}
\description{
Calculates the confidence interval of grouped data
}
\examples{
run_group.CI(Sepal.Length~Species, iris, 0.95)

}
\references{
https://github.com/cran/Rmisc/
}
\author{
Created by Hua Zou (5/19/2022 Shenzhen China)
}
