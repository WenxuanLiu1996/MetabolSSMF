% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{bootstrap}
\alias{bootstrap}
\title{bootstrap}
\usage{
bootstrap(data, k, btimes, H)
}
\arguments{
\item{data}{Numeric matrix of data.}

\item{k}{The number of prototypes/clusters.}

\item{btimes}{Integer, number of bootstrap samples.}

\item{H}{Matrix, input \eqn{H} matrix to start the algorithm.}
}
\value{
\code{W.est} The \eqn{W} matrix estimated by bootstrap.

\code{lower} Lower bound of confidence interval.

\code{upper} Upper bound of confidence interval.
}
\description{
Bootstrap resampling approach to estimated the confidence interval.
}
\details{
Bootstrap resampling approach is designed to estimate the variance, confidence intervals and other properties of statistics.
The basic idea is given by the \eqn{\frac{\alpha}{2}} and \eqn{1-\frac{\alpha}{2}} percentiles of the \eqn{B} bootstrap estimated values, where \eqn{B} is the times of bootstrapping.

Bootstrap samples of size \eqn{n} is created by applying a random number generator and substituting observations. The \eqn{b^{th}} bootstrap sample is denoted as
\deqn{X^{*(b)}=(x_1^{*b}, x_2^{*b},...,x_n^{*b}),}

where each \eqn{x_i^{*b}} is a random selection from the original observation.
Then, apply the SSMF algorithm to each bootstrap sample and calculate the \eqn{b^{th}} bootstrap replicate of the prototypes matrix,
which is denoted as \eqn{W^{*(b)}}

The estimate variance of \eqn{B} bootstrap replicates can be calculated by

\deqn{sd(W^*) = \sqrt {\frac{1}{B-1} \sum_{b=1}^{B} [W^{*(b)}-\bar{W}^{*}]^2}}

where \eqn{\bar{W}^{*}=\frac{1}{B} \sum_{b=1}^{B} W^{*(b)}}. Therefore, the 95\% CIs for the prototypes can be calculated by

\deqn{(W^*-t_{(0.025, B-1)} \times sd, W^*+t_{(0.975, B-1)} \times sd),}
where \eqn{t_{(0.975, n-1)}} is the quantile of student \eqn{t} distribution with 95\% significance and \eqn{(B-1)} degree freedom.
}
\references{
STINE, R. (1989). An Introduction to Bootstrap Methods: Examples and Ideas. Sociological Methods & Research, 18(2-3), 243-291. \url{https://doi.org/10.1177/0049124189018002003}
}
\author{
Wenxuan Liu
}
