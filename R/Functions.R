# library(tidyr)
# library(lsei)
# library(LaplacesDemon)
# library(iterators)
# library(foreach)
# library(doParallel)


#' @title init
#' @description This function initialise the \eqn{H_{n \times k} = (h_1, \dots, h_i, \dots, h_n)} matrix
#' and \eqn{W_{k \times p} = (w_1, \dots, w_j, \dots, w_p)^T} matrix to start the SSMF model.
#' This function is often used in conjunction with the function ssmf().
#' @usage init(k, n, data = NULL, method = c('dirichlet', 'uniform', 'kmeans', 'nmf'))
#' @param k The number of prototypes/clusters.
#' @param n The number of observations.
#' @param data Numeric matrix of data.
#' @param method character: 'dirichlet', 'uniform', 'kmeans' or 'nmf'.
#'
#' @details
#' 'dirichlet': sampling \eqn{H} matrix from dirichlet distribution, grouping the observations with their maximum memberships in \eqn{h_i} and calculate the means in each group as prototype vectors of \eqn{W} matrix;
#'
#' 'uniform': sampling \eqn{H} matrix from uniform distribution, \eqn{W} matrix is created in the same way of 'dirichlet';
#'
#' 'kmeans': applying the centres of the kmeans clustering to create \eqn{W} matrix and converts the classification into a \eqn{H} matrix;
#'
#' 'nmf': using the matrix of basis components and the coefficient matrix of an NMF model to create \eqn{W} matrix and \eqn{H} matrix.
#'
#' @import mclust
#' @import NMF
#' @import tidyr
#' @import LaplacesDemon
#'
#' @return Initialised \eqn{H}, \eqn{W} matrix
#' @author Wenxuan Liu
#' @export
init <- function(k, n, data = NULL, method = c('dirichlet', 'uniform', 'kmeans', 'nmf')) {
  res <- list()

  if (length(method) > 1) {
    m <- method[1]            #default method is uniform
  }
  else {
    m <- method
  }

  if (m == 'kmeans') {
    km <- kmeans(data, centers = k)
    h <- mclust::unmap(km$cluster)
    w <- km$centers
  }

  else if (m == 'uniform') {
    h <- matrix(1, nrow = n, ncol = k)
    w <- t(apply(data, 2, mean))
    if (k > 1) {
      h <- matrix(runif(n * k, 0, 1), nrow = n, ncol = k)
      h <- t(apply(h, 1, function(x) x / sum(x)))
      w <- aggregate(. ~ group, mean, data = as.data.frame(cbind(data, 'group' = apply(h, 1, which.max))))[,-1]
    }
  }
  else if (m == 'dirichlet') {
    h <- rdirichlet(n, sample(seq(0.1, 2, 0.1), k))
    w <- aggregate( . ~ group, mean, data = as.data.frame(cbind(data, 'group' = apply(h, 1, which.max))))[,-1]
  }
  else if (m == 'nmf') {
    if (is.null(data) == F) {
      h <- matrix(1, nrow = n, ncol = k)
      w <- t(apply(data, 2, mean))
      if (k > 1) {
        nmf <- NMF::nmf(t(data), rank = k)
        w <- t(NMF::basis(nmf))
        h <- NMF::coef(nmf)
        h <- t(apply(h, 2, function(x) x / sum(x)))
      }
    }
    else {return('Data should be included while using NMF')}
  }
  else {return('No definded method')}

  res$H <- h
  res$W <- w

  return(res)
}


#' @title ssmf(Simplex-structure matrix factorisation)
#' @description This function is running ssmf model.
#' @param data Numeric matrix of data.
#' @param k The number of prototypes/clusters.
#' @param H Matrix, user input \eqn{H} matrix to start the algorithm. If input is empty, the function will initialise \eqn{H} matrix automatically.
#' @param W Matrix, user input \eqn{W} matrix to start the algorithm. If input is empty, the function will initialise \eqn{W} matrix automatically.
#' @param meth Specification of method to initialise the \eqn{W} and \eqn{H} matrix, see 'method' in \code{init()}.
#' @param lr Optimisation learning rate.
#' @param nruns The maximum times of running the algorithm.

#' @details
#' Let \eqn{X \in R^{n \times p}} be the data set consisting of \eqn{p} variables and \eqn{n} samples.
#' Given an integer \eqn{k \ll \text{min}(n,p)}, the rank of factorisation,
#' SSMF finds a matrix \eqn{W_{k \times p}} containing the basis vectors in each row and a matrix \eqn{H_{n \times k}},
#' whose rows are the mixing weights of the observations, such that \deqn{X \approx HW,}
#' where \eqn{H \geq 0} and each row of \eqn{H} belongs to the unit simplex, that is \deqn{\sum_r^k h_{ir}=1, \text{for}\ i=1,2,...,n.}
#'
#' The optimisation is described as gradient descent that alternates between determining the optimal \eqn{W} for a given \eqn{H} and determining the optimal
#' \eqn{H} for a given \eqn{W}. The loss function is defined as residual sum of squares (RSS), that is \deqn{\mathrm{RSS} = \| X-HW \|^2 = \sum_{i,j} (X_{ij}-(HW)_{ij})^2,}
#' \eqn{\text{for}\ i=1,2,...,n; j=1,2,...,p}. The algorithm will stop until RSS reduction is sufficiently small or the number of maximum iterations is reached.
#'
#' @import lsei
#'
#' @return \code{W} The optimised \eqn{W} matrix, containing the values of prototypes.
#' @return \code{H} The optimised \eqn{H} matrix, containing the values of soft clustering memberships.
#' @return \code{SSE} The sum of square errors/residuals sum of square.
#' @export
#' @references Abdolali, Maryam & Gillis, Nicolas. (2020). Simplex-Structured Matrix Factorization: Sparsity-based Identifiability and Provably Correct Algorithms. \url{https://doi.org/10.1137/20M1354982}
#' @author Wenxuan Liu
#' @examples
#' library(MetabolSSMF)
#' fit1 <- ssmf(data, k = 4, lr = 0.001, H = H)
#' fit2 <- ssmf(data, k = 4, H = H)
#' fit3 <- ssmf(data, k = 4, W = W)
#' fit4 <- ssmf(data, k = 4, meth = 'dirichlet')
#' fit5 <- ssmf(data, k = 4)

ssmf <- function(data, k, H = NULL, W = NULL, meth = c('dirichlet', 'uniform', 'kmeans', 'nmf'), lr = 0.01, nruns = 50) {
  res <- list()
  itr <- 1
  rssflow <- c()

  critold <- Inf
  crit <- 0

  if (is.null(H) == T & is.null(W) == T) {
    initialise <- init( k = k, n = nrow(data), method = meth, data = data)
    H <- initialise$H

    while ((critold - crit) > lr & itr <= nruns) {
#     print(itr)

      W <- coef(lm(data ~ H - 1))

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      for (i in 1:nrow(data)) {
        h <- pnnls(a = t(W), b = t(data[i,]), sum = 1)$x
        H[i,] <- h
      }

      critold <- crit

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      rssflow <- c(rssflow, critold, crit)

      itr <- itr + 1
    }
  }



  else if (is.null(H) == T & is.null(W) == F) {
    H <- matrix(NA, nrow = nrow(data), ncol = k)
    while ((critold - crit) > lr & itr <= nruns) {
#      print(itr)

      for (i in 1:nrow(data)) {
        h <- pnnls(a = t(W), b = t(data[i,]), sum = 1)$x
        H[i,] <- h
      }

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      W <- coef(lm(data ~ H - 1))

      critold <- crit

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      rssflow <- c(rssflow, critold, crit)

      itr <- itr + 1
    }

  }
  else if (is.null(H) == F & is.null(W) == T) {
    while ((critold - crit) > lr & itr <= nruns) {
#      print(itr)

      W <- coef(lm(data ~ H - 1))

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      for (i in 1:nrow((data))) {
        h <- pnnls(a = t(W), b = t(data[i,]), sum = 1)$x
        H[i,] <- h
      }

      critold <- crit

      crit <- sum((data - H %*% W)^2)

#      print(crit)

      rssflow <- c(rssflow, critold, crit)

      itr <- itr + 1
    }
  }

  colnames(W) <- colnames(data)
  rownames(H) <- rownames(data)
  rownames(W) <- 1:k
  colnames(H) <- 1:k

  res$H <- H
  res$W <- W
  res$SSE <- crit
#  res$rssflow <- rssflow
#  res$det <- det

  return(res)

}

#Entropy
#' @title entropy
#' @description Calculate the entropy of the memberships of each observation. The base of the logarithm is 2.
#' @param x Vector of memberships
#'
#' @details Given the membership vector of the \eqn{i} observation \eqn{h_i}, the entropy is defined as
#' \deqn{Entropy(h_i) = -\sum_{r=1}^k h_{ir} \cdot log_b h_{ir},}
#' where b is the base of the logarithm used.
#'
#' @return A numeric value of entropy.
#' @export
#' @author Wenxuan Liu
#' @examples
#' membership <- c(0.1, 0.2, 0.3, 0.4)
#' entropy(membership)
entropy <- function(x) {
  e <- rep(NA, length(x))
  for (i in 1:length(x)) {
    if (x[i] == 0) {
      e[i] <- 0
    }
    else {
      e[i] <- x[i] * log2(x[i])
    }
  }
  return(-sum(e))
}

#Gap
#' @title gap
#' @description Estimating the number of prototypes/clusters in a data set via the gap statistics.
#' @param data Numeric matrix of data.
#' @param k.vector A vector of the number of prototypes/clusters.
#' @param rss Numeric vector, sum of square errors/residual sum of square of of ssmf model
#' @param meth Character, specification of method to initialise the \eqn{W} and \eqn{H} matrix, see 'method' in init().
#' @param itr Integer, number of Monte Carlo samples.
#'
#' @details
#' This gap statistics technique is selecting the biggest difference between the original RSS and the RSS under an appropriate null reference distribution of the data, which is defined to be
#' \deqn{\mathrm{Gap}(k) = \frac{1}{B} \sum_{b} \log(\mathrm{RSS}^*_{kb}) - \log(\mathrm{RSS}_{k})}
#'
#' where \eqn{B, (b=1, 2, ..., B)} is the number of sampled null reference distributions of the data;
#' \eqn{RSS^*_{kb}} is the residual sum of squares of the \eqn{b^\mathrm{th}} null reference distributions of the data that fitted in the SSMF model;
#' \eqn{RSS_{k}} is the residual sum of squares of the original data \eqn{X}. The estimated gap suggests the number of prototypes via
#'
#' \deqn{\hat{k} = \mathrm{smallest} \ k \  \mathrm{such \ that} \ \mathrm{Gap}(k) \geq \mathrm{Gap}(k+1) - s_{k+1}}
#'
#' where \eqn{s_{k+1}} is standard error that is defined as
#'
#' \deqn{s_{k+1}=sd_k \sqrt{1+\frac{1}{B}}}
#'
#' and \eqn{sd_k} is the standard deviation:
#'
#' \deqn{sd_k=\{ \frac{1}{B} \sum_{b} [\log(RSS^*_{kb})-\frac{1}{B} \sum_{b} \log(RSS^8_{kb})]^2 \}^\frac{1}{2}}.
#'
#' @import iterators
#' @import foreach
#' @import doParallel
#'
#' @return \code{gap} Gap values.
#' @return \code{optimal.k} The optimal number of prototypes/clusters that gap statistics suggested.
#' @return \code{standard.error} Standard error
#' @references Tibshirani, R., Walther, G., & Hastie, T. (2001). Estimating the Number of Clusters in a Data Set via the Gap Statistic. Journal of the Royal Statistical Society. Series B (Statistical Methodology), 63(2), 411–423. http://www.jstor.org/stable/2680607
#' @author Wenxuan Liu
#' @export
#'

gap <- function(data, k.vector, rss, meth = 'kmeans', itr = 50) {
  result <- list()
  registerDoParallel()
  rss_gap <- foreach(icount(itr), .combine = rbind) %dopar% {
    rss_b <- rep(NA, length(k.vector))
    data_b <- apply(data, 2, function(x) runif(length(x)))
    for (k in k.vector) {
      fit_gap <- ssmf(as.matrix(data_b), k = k, meth = meth)
      rss_b[k] <- fit_gap$SSE
    }
    rss_b
  }
  gap <- apply(log(rss_gap), 2, mean) - log(rss)

  l <- apply(log(rss_gap), 2, mean)
  sdk <- sqrt(apply((log(rss_gap) - apply(log(rss_gap), 2, mean))^2, 2, mean))
  sk <- sqrt(1 + (1 / itr)) * sdk

  opt.k <- which( gap[-length(k.vector)] - (gap - sk)[-1] >= 0)[1]

  result$gap <- gap
  result$optimal.k <- opt.k
  result$standard.error <- sk
  return(result)
}


#Bootstrap
#' @title bootstrap
#' @description Bootstrap resampling approach to estimated the confidence interval.
#' @param data Numeric matrix of data.
#' @param k The number of prototypes/clusters.
#' @param btimes Integer, number of bootstrap samples.
#' @param H Matrix, input \eqn{H} matrix to start the algorithm.
#'
#' @details
#' Bootstrap resampling approach is designed to estimate the variance, confidence intervals and other properties of statistics.
#' The basic idea is given by the \eqn{\frac{\alpha}{2}} and \eqn{1-\frac{\alpha}{2}} percentiles of the \eqn{B} bootstrap estimated values, where \eqn{B} is the times of bootstrapping.
#'
#' Bootstrap samples of size \eqn{n} is created by applying a random number generator and substituting observations. The \eqn{b^{th}} bootstrap sample is denoted as
#' \deqn{X^{*(b)}=(x_1^{*b}, x_2^{*b},...,x_n^{*b}),}
#'
#' where each \eqn{x_i^{*b}} is a random selection from the original observation.
#' Then, apply the SSMF algorithm to each bootstrap sample and calculate the \eqn{b^{th}} bootstrap replicate of the prototypes matrix,
#' which is denoted as \eqn{W^{*(b)}}
#'
#' The estimate variance of \eqn{B} bootstrap replicates can be calculated by
#'
#' \deqn{sd(W^*) = \sqrt {\frac{1}{B-1} \sum_{b=1}^{B} [W^{*(b)}-\bar{W}^{*}]^2}}
#'
#' where \eqn{\bar{W}^{*}=\frac{1}{B} \sum_{b=1}^{B} W^{*(b)}}. Therefore, the 95\% CIs for the prototypes can be calculated by
#'
#' \deqn{(W^*-t_{(0.025, B-1)} \times sd, W^*+t_{(0.975, B-1)} \times sd),}
#
#' where \eqn{t_{(0.975, n-1)}} is the quantile of student \eqn{t} distribution with 95\% significance and \eqn{(B-1)} degree freedom.
#'
#' @import iterators
#' @import foreach
#' @import doParallel
#'
#' @return \code{W.est} The \eqn{W} matrix estimated by bootstrap.
#' @return \code{lower} Lower bound of confidence interval.
#' @return \code{upper} Upper bound of confidence interval.
#' @references STINE, R. (1989). An Introduction to Bootstrap Methods: Examples and Ideas. Sociological Methods & Research, 18(2-3), 243-291. \url{https://doi.org/10.1177/0049124189018002003}
#' @author Wenxuan Liu
#' @export
#'

bootstrap <- function(data, k, btimes, H) {
  result <- list()
  n <- dim(data)[1]
  registerDoParallel()
  out <- foreach(icount(btimes)) %dopar% {
    ind <- sample(1:n, replace = TRUE, size = n)
    data_boot <- data[ind, ]
    H_boot <- H[ind,]
    res_boot <- ssmf(as.matrix(data_boot), k = k, H = H_boot, lr = 0.001)
    res_boot
  }
  W_list <- lapply(out, function(x) x$W)

  W_b <- matrix(0, k, dim(data)[2])
  for (i in 1:btimes) {
    W_b <- W_b + W_list[[i]]
  }
  W_bmean <- W_b / btimes


  W_bvar <- matrix(0, k, dim(data)[2])
  for (b in 1:btimes) {
    W_bvar <- W_bvar + (W_list[[b]] - W_bmean)^2
  }
  W_bse <- sqrt((W_bvar) / (btimes-1))

  upper <- W_bmean + 1.96 * W_bse
  lower <- W_bmean - 1.96 * W_bse

  result$W.est <- W_bmean
  result$lower <- lower
  result$upper <- upper
#  result$se <- W_bse

  return(result)
}


#Soft adjust rand index
#' @title sARI
#' @description Soft adjust rand index, a soft agreement measure for class partitions incorporating assignment probabilities
#' @param partition1 Numeric matrix of the (posterior) probabilities of assignment of observations in partition 1, memberships matrix.
#' @param partition2 Numeric matrix of the (posterior) probabilities of assignment of observations in partition 2, memberships matrix.
#'
#' @details
#' A soft adjusted Rand index (sARI) is proposed to compare the two soft partitions (e.g., the true and estimated fuzzy memberships).
#' Let \eqn{p_{rci} = h_{ri} \cdot h_{ci}} be the posterior probabilities of assignment of observation \eqn{i} to \eqn{r^{th}}
#' cluster in the true partition and to \eqn{c^{th}} cluster in the estimated partition, that is \eqn{h_{ri}} and \eqn{h_{ci}} are the given true and estimated fuzzy membership of sample \eqn{i}.
#' we have the sum of the posterior probabilities of assignment to \eqn{r^{th}} cluster in the true partition and to \eqn{c^{th}} cluster in the estimated partition over all observations.
#'
#' \deqn{p_{rc.} = \sum_{i=1}^n p_{rci};}
#'
#' the sum of the posterior probabilities for belonging to assignment to \eqn{r^{th}} cluster in the true partition over all observations.
#' \deqn{p_{r..} = \sum_{c=1}^C \sum_{i=1}^n p_{rci};}
#'
#' the sum of the posterior probabilities for belonging to assignment to \eqn{c^{th}} cluster in the estimated partition over all observations.
#' \deqn{p_{.c.} = \sum_{r=1}^R \sum_{i=1}^n p_{rci}.}
#'
#' The sARI for true and estimated soft partitions is defined as
#' \deqn{\mathrm{sARI}(H_{true}, H_{est}) = \frac{\sum_{r,c} \frac{\Gamma(p_{rc.}+1)}{\Gamma(p_{rc.}-1)} - \frac{1}{n(n-1)} \Lambda_{rc}}{\frac{1}{2} \Lambda_{rc} - \frac{1}{n(n-1)} \Lambda_{rc}}}
#'
#' where \eqn{\Lambda_{rc} = \sum_{r} \frac{\Gamma(p_{r..}+1)}{\Gamma(p_{r..}-1)} + \sum_{c} \frac{\Gamma(p_{.c.}+1)}{\Gamma(p_{.c.}-1)}}. A larger sARI illustrates the more similar between these two partitions.
#'
#'
#'
#' @return Soft adjust rand index.
#' @references Flynt, A., Dean, N. & Nugent, R. sARI: a soft agreement measure for class partitions incorporating assignment probabilities. Adv Data Anal Classif 13, 303–323 (2019). \url{https://doi.org/10.1007/s11634-018-0346-x}
#' @author Wenxuan Liu
#' @export
#'
sARI <- function(partition1, partition2) {
  N <- nrow(partition1)
  cftab <- t(as.matrix(partition1)) %*% as.matrix(partition2)
  L <- sum(gamma(colSums(cftab) + 1) / gamma(colSums(cftab) - 1)) + sum(gamma(rowSums(cftab) + 1) / gamma(rowSums(cftab) - 1))
  result <- (sum(gamma(cftab + 1) / gamma(cftab - 1)) - (1 / (N * (N - 1))) * L) / (0.5 * L - (1 / (N * (N - 1))) * L)
  return(result)
}
