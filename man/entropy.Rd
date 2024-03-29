% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{entropy}
\alias{entropy}
\title{entropy}
\usage{
entropy(x)
}
\arguments{
\item{x}{Vector of memberships}
}
\value{
A numeric value of entropy.
}
\description{
Calculate the entropy of the memberships of each observation. The base of the logarithm is 2.
}
\details{
Given the membership vector of the \eqn{i} observation \eqn{h_i}, the entropy is defined as
\deqn{Entropy(h_i) = -\sum_{r=1}^k h_{ir} \cdot log_b h_{ir},}
where b is the base of the logarithm used.
}
\examples{
membership <- c(0.1, 0.2, 0.3, 0.4)
entropy(membership)
}
\author{
Wenxuan Liu
}
