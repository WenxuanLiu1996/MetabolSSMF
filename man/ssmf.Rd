% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{ssmf}
\alias{ssmf}
\title{ssmf(Simplex-structure matrix factorisation)}
\usage{
ssmf(
  data,
  k,
  H = NULL,
  W = NULL,
  meth = c("dirichlet", "uniform", "kmeans", "nmf"),
  lr = 0.01,
  nruns = 50
)
}
\arguments{
\item{data}{Numeric matrix of data.}

\item{k}{The number of prototypes/clusters.}

\item{H}{Matrix, user input \eqn{H} matrix to start the algorithm. If input is empty, the function will initialise \eqn{H} matrix automatically.}

\item{W}{Matrix, user input \eqn{W} matrix to start the algorithm. If input is empty, the function will initialise \eqn{W} matrix automatically.}

\item{meth}{Specification of method to initialise the \eqn{W} and \eqn{H} matrix, see 'method' in \code{init()}.}

\item{lr}{Optimisation learning rate.}

\item{nruns}{The maximum times of running the algorithm.}
}
\value{
\code{W} The optimised \eqn{W} matrix, containing the values of prototypes.

\code{H} The optimised \eqn{H} matrix, containing the values of soft clustering memberships.

\code{SSE} The sum of square errors/residuals sum of square.
}
\description{
This function is running ssmf model.
}
\details{
Let \eqn{X \in R^{n \times p}} be the data set consisting of \eqn{p} variables and \eqn{n} samples.
Given an integer \eqn{k \ll \text{min}(n,p)}, the rank of factorisation,
SSMF finds a matrix \eqn{W_{k \times p}} containing the basis vectors in each row and a matrix \eqn{H_{n \times k}},
whose rows are the mixing weights of the observations, such that \deqn{X \approx HW,}
where \eqn{H \geq 0} and each row of \eqn{H} belongs to the unit simplex, that is \deqn{\sum_r^k h_{ir}=1, \text{for}\ i=1,2,...,n.}

The optimisation is described as gradient descent that alternates between determining the optimal \eqn{W} for a given \eqn{H} and determining the optimal
\eqn{H} for a given \eqn{W}. The loss function is defined as residual sum of squares (RSS), that is \deqn{\mathrm{RSS} = \| X-HW \|^2 = \sum_{i,j} (X_{ij}-(HW)_{ij})^2,}
\eqn{\text{for}\ i=1,2,...,n; j=1,2,...,p}. The algorithm will stop until RSS reduction is sufficiently small or the number of maximum iterations is reached.
}
\examples{
library(MetabolSSMF)
fit1 <- ssmf(data, k = 4, lr = 0.001, H = H)
fit2 <- ssmf(data, k = 4, H = H)
fit3 <- ssmf(data, k = 4, W = W)
fit4 <- ssmf(data, k = 4, meth = 'dirichlet')
fit5 <- ssmf(data, k = 4)
}
\references{
Abdolali, Maryam & Gillis, Nicolas. (2020). Simplex-Structured Matrix Factorization: Sparsity-based Identifiability and Provably Correct Algorithms. \url{https://doi.org/10.1137/20M1354982}
}
\author{
Wenxuan Liu
}
