# library(tidyr)
# library(lsei)
# library(LaplacesDemon)
# library(iterators)
# library(foreach)
# library(doParallel)


#' @title Initialise the membership matrix \eqn{H} or prototype matrix \eqn{W}.
#' @description This function initialises the \eqn{H_{n \times k}} matrix
#' or the \eqn{W_{k \times p}} matrix to start the SSMF model.
#' This function is often used in conjunction with the function ssmf( ). Also, the code can be run separately from the function
#' ssmf( ). This function returns to simplex-structured soft membership matrix \eqn{H} and prototype matrix \eqn{W}.
#'
#' @param data Data matrix or data frame.
#' @param k The number of prototypes/clusters.
#' @param method Character: 'kmeans', 'uniform', 'dirichlet' or 'nmf'. If there are more than one method,
#' the default is selecting the first method in the vector.
#'
#' @details
#' 'kmeans': create the \eqn{W} matrix using the centres of the kmeans output; create the \eqn{H} matrix by converting the classification into a binary matrix.
#'
#' 'uniform': create the \eqn{H} matrix by sampling the values from uniform distribution and making the rows of the matrix lie in the unit simplex; group the observations with their maximum memberships
#' and create the \eqn{W} matrix by combining the mean vector in each group.
#'
#' 'dirichlet': create the \eqn{H} matrix by sampling the values from Dirichlet distribution; group the observations with their maximum memberships
#' and create the \eqn{W} matrix by combining the mean vector in each group.
#'
#' 'nmf': create the \eqn{W} matrix using the matrix of basic components from NMF model; the coefficient matrix is acquired from NMF model,
#' then the \eqn{H} is created by making the rows of the coefficient matrix lie in the unit simplex.
#'
#' @import mclust
#' @import NMF
#' @import tidyr
#' @import LaplacesDemon
#'
#' @return Initialised \eqn{H}, \eqn{W} matrix.
#' @author Wenxuan Liu
#' @export
#' @examples
#' # example code
#' \donttest{
#' init(data = SimulatedDataset, k = 4, method = 'kmeans')
#' }

init <- function(data = NULL, k = NULL, method = c('kmeans', 'uniform', 'dirichlet', 'nmf')) {
  if(is.null(k) | is.null(data)){
    stop('Please don\'t leave k or data empty')
  }

  n <- nrow(data)

  res <- list()

  if (length(method) > 1) {
    m <- method[1]            #default method is kmeans
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


#' @title Simplex-structured matrix factorisation algorithm (SSMF).
#' @description This function implements on SSMF on a data matrix or data frame.
#' @param data Data matrix or data frame.
#' @param k The number of prototypes/clusters.
#' @param H Matrix, user input \eqn{H} matrix to start the algorithm. If input is empty, the function will initialise \eqn{H} matrix automatically.
#' @param W Matrix, user input \eqn{W} matrix to start the algorithm. If input is empty, the function will initialise \eqn{W} matrix automatically.
#' @param meth Specification of method to initialise the \eqn{W} and \eqn{H} matrix, see 'method' in \code{init()}.
#' @param lr Optimisation learning rate.
#' @param nruns The maximum times of running the algorithm.

#' @details
#' Let \eqn{X \in R^{n \times p}} be the data set  with \eqn{n} observations and \eqn{p} variables.
#' Given an integer \eqn{k \ll \text{min}(n,p)},
#' the data set is clustered by simplex-structured matrix factorisation (SSMF), which aims to process soft clustering
#' and partition the observations into \eqn{k} fuzzy clusters such that the sum of squares from observations to the
#' assigned cluster prototypes is minimised.

#' SSMF finds \eqn{H_{n \times k}} and \eqn{W_{k \times p}},
#' such that \deqn{X \approx HW,}

#' A cluster prototype refers to a vector that represent the characteristics of a particular cluster,
#' denoted by \eqn{w_r \in \mathbb{R}^{p}} , where \eqn{r} is the \eqn{r^{th}} cluster.
#' A cluster membership vector \eqn{h_i \in \mathbb{R}^{k}} describes the proportion of the cluster prototypes
#' of the \eqn{i^{th}} observation. \eqn{W} is the prototype matrix where each row is the cluster prototype and
#' \eqn{H} is the soft membership matrix where each row gives the soft cluster membership of each observation.

#' The problem of finding the approximate matrix factorisation is solved by minising residual sum of squares (RSS), that is
#' \deqn{\mathrm{RSS} = \| X-HW \|^2 = \sum_{i=1}^{n}\sum_{j=1}^{p} \left\{ X_{ij}-(HW)_{ij}\right\}^2,}
#' such that \eqn{\sum_{r=1}^k h_{ir}=1}  and  \eqn{h_{ir}\geq 0}.
#'
#' @import lsei
#'
#' @return \code{W} The optimised \eqn{W} matrix, containing the values of prototypes.
#' @return \code{H} The optimised \eqn{H} matrix, containing the values of soft memberships.
#' @return \code{SSE} The residuals sum of square.
#' @export
#' @references Abdolali, Maryam & Gillis, Nicolas. (2020). Simplex-Structured Matrix Factorization: Sparsity-based Identifiability and Provably Correct Algorithms. <doi:10.1137/20M1354982>
#' @author Wenxuan Liu
#' @examples
#'
#' \donttest{
#' library(MetabolSSMF)
#'
#' # Initialisation by user
#' data <- SimulatedDataset
#' k <- 4
#'
#' ## Initialised by kmeans
#' fit.km <- kmeans(data, centers = k)
#'
#' H <- mclust::unmap(fit.km$cluster)
#' W <- fit.km$centers
#'
#' fit1 <- ssmf(data, k = k, H = H) #start the algorithm from H
#' fit2 <- ssmf(data, k = k, W = W) #start the algorithm from W
#'
#' # Initialisation inside the function
#' fit3 <- ssmf(data, k = 4, meth = 'dirichlet')
#' fit4 <- ssmf(data, k = 4)
#' }

ssmf <- function(data, k, H = NULL, W = NULL, meth = c('kmeans', 'uniform', 'dirichlet', 'nmf'), lr = 0.01, nruns = 50) {

  if(is.data.frame(data)){
    data <- as.matrix(data)
  }


  res <- list()
  itr <- 1
  rssflow <- c()

  critold <- Inf
  crit <- 0

  if (is.null(H) == T & is.null(W) == T) {
    initialise <- init( k = k, method = meth, data = data)
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

#Shannon diversity index
#' @title Shannon diversity index
#' @description Calculate the Shannon diversity index of the memberships of an observation. The base of the logarithm is 2.
#' @param x A membership vector.
#' @param two.power Logical, whether return to the value of \eqn{2^{\mathrm{E}(h_{i})}}.
#'
#' @details Given a membership vector of the \eqn{i^{th}} observation \eqn{h_i}, the Shannon diversity index is defined as
#' \deqn{\mathrm{E}(h_{i}) = -\sum_{r=1}^k h_{ir} \mathrm{log}_2 (h_{ir}).}
#' Specifically, in the case of \eqn{h_{ir}=0}, the value of \eqn{h_{ir} \mathrm{log}_2 (h_{ir})} is taken to be 0.
#'
#'
#' @return A numeric value of Shannon diversity index \eqn{\mathrm{E}(h_{i})} or \eqn{2^{\mathrm{E}(h_{i})}}.
#' @export
#' @author Wenxuan Liu
#' @examples
#' # Memberships vector
#' membership1 <- c(0.1, 0.2, 0.3, 0.4)
#' diversity(membership1)
#' diversity(membership1, two.power = TRUE)
#'
#' # Memberships matrix
#' membership2 <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.4, 0.1, 0.2, 0.3, 0.1, 0.4),
#'                       nrow=3, ncol=4, byrow=TRUE)
#'
#' E <- rep(NA, nrow(membership2))
#' for(i in 1:nrow(membership2)){
#'   E[i] <- diversity(membership2[i,])
#' }
#' E

diversity <- function(x, two.power=FALSE) {
  e <- rep(NA, length(x))
  for (i in 1:length(x)) {
    if (x[i] == 0) {
      e[i] <- 0
    }
    else {
      e[i] <- x[i] * log2(x[i])
    }
  }
  if(two.power==TRUE){
    return(2^(-sum(e)))
  }
  else{
    return(-sum(e))
  }
}

#Gap
#' @title Gap statistic algorithm.
#' @description Estimating the number of prototypes/clusters in a data set using the gap statistic.
#' @param data Data matrix or data frame.
#' @param rss Numeric vector, residual sum of squares from ssmf model using the number of clusters \eqn{1,2, \ldots, k}.
#' @param meth Character, specification of method to initialise the \eqn{W} and \eqn{H} matrix, see 'method' in init( ).
#' @param itr Integer, number of Monte Carlo samples.
#' @param ncore The number of cores to use for parallel execution.
#'
#' @details
#' This gap statistic selects the biggest difference between the original residual sum of squares (RSS) and the RSS under an appropriate null reference distribution of the data, which is defined to be
#' \deqn{\mathrm{Gap}(k) = \frac{1}{B} \sum_{b=1}^{B} \log(\mathrm{RSS}^*_{kb}) - \log(\mathrm{RSS}_{k}),}
#'
#' where \eqn{B} is the number of samples from the reference distribution;
#' \eqn{\mathrm{RSS}^*_{kb}} is the residual sum of squares of the \eqn{b^th} sample from the reference distribution fitted in the SSMF model model using \eqn{k} clusters;
#' \eqn{RSS_{k}} is the residual sum of squares for the original data \eqn{X} fitted the model using the same \eqn{k}.
#' The estimated gap suggests the number of prototypes/clusters (\eqn{\hat{k}}) using
#'
#' \deqn{\hat{k} = \mathrm{smallest} \ k \  \mathrm{such \ that} \ \mathrm{Gap}(k) \geq \mathrm{Gap}(k+1) - s_{k+1},}
#'
#' where \eqn{s_{k+1}} is standard error that is defined as
#'
#' \deqn{s_{k+1}=sd_k \sqrt{1+\frac{1}{B}},}
#'
#' and \eqn{sd_k} is the standard deviation:
#'
#' \deqn{sd_k=\sqrt{ \frac{1}{B} \sum_{b} [\log(\mathrm{RSS}^*_{kb})-\frac{1}{B} \sum_{b} \log(\mathrm{RSS}^*_{kb})]^2}.}
#'
#' @import iterators
#' @import foreach
#' @import doParallel
#'
#' @return \code{gap} Gap value vector.
#' @return \code{optimal.k} The optimal number of prototypes/clusters.
#' @return \code{standard.error} Standard error vector.
#' @references Tibshirani, R., Walther, G., & Hastie, T. (2001). Estimating the Number of Clusters in a Data Set via the Gap Statistic. Journal of the Royal Statistical Society. Series B (Statistical Methodology), 63(2), 411–423. <doi:10.1111/1467-9868.00293>
#' @author Wenxuan Liu
#' @export
#' @examples
#' # example code
#' \donttest{
#' data <- SimulatedDataset
#'
#' k <- 6
#'
#' rss <- rep(NA, k)
#' for(i in 1:k){
#'   rss[i] <- ssmf(data = data, k = i)$SSE
#' }
#'
#' gap(data = data, rss = rss)
#' }

gap <- function(data, rss, meth = c('kmeans', 'uniform', 'dirichlet', 'nmf'), itr = 50, ncore = 2) {

  k.vector <- 1:length(rss)

  result <- list()
  registerDoParallel(cores = ncore)
  rss_gap <- foreach(icount(itr), .combine = rbind) %dopar% {
    rss_b <- rep(NA, length(k.vector))
    data_b <- apply(data, 2, function(x) runif(length(x)))
    for (k in k.vector) {
      fit_gap <- ssmf(data_b, k = k, meth = meth)
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
#' @title Bootstrap algorithm function.
#' @description Bootstrap resampling approach to estimate the confidence intervals for the cluster prototypes.
#' @param data Data matrix or data frame.
#' @param k The number of prototypes/clusters.
#' @param H Matrix, input \eqn{H} matrix to start the algorithm. Usually the \eqn{H} matrix is the output of the function ssmf( ).
#' If \eqn{H} is not supplied, the bootstrapped \eqn{W} matrix might have different prototype orders from the outputs of the function ssmf( ).
#' @param mtimes Integer, number of bootstrap samples. Default number is 50.
#' @param ncore The number of cores to use for parallel execution.
#'
#' @details
#' Create bootstrap samples of size \eqn{n} by sampling from the data set with replacement and repeat the steps \eqn{M} times.
#' The \eqn{m^{th}} bootstrap sample is denoted as
#' \deqn{X^{{\ast}(m)}=(x_1^{{\ast}(m)}, x_2^{{\ast}(m)},\ldots,x_n^{{\ast}(m)}),}
#'
#' where each \eqn{x_i^{{\ast}(m)}} is a random sample (with replacement) from the data set.
#'
#' Then, apply the SSMF algorithm to each bootstrap sample and calculate the \eqn{m^{th}} bootstrap replicate of the prototypes matrix,
#' which is denoted as \eqn{W^{{\ast}(m)}}.
#'
#' The estimate standard deviation of \eqn{M} bootstrap replicates can be calculated by
#'
#' \deqn{sd(W^{\ast}) =\sqrt {\frac{1}{M-1} \sum_{m=1}^{M} [W^{{\ast}(m)}-\overline{W}^{\ast}]^2 },}
#'
#' where \eqn{\overline{W}^{\ast}=\frac{1}{M} \sum_{m=1}^{M} W^{{\ast}(m)}}. Therefore, the 95\% CIs for the prototypes can be calculated by
#'
#' \deqn{(\overline{W}^{\ast}-t_{(0.025, M-1)} \cdot sd(W^{\ast}),\ \overline{W}^{\ast}+t_{(0.975, M-1)} \cdot sd(W^)),}
#
#' where \eqn{t_{(0.025, n-1)}} and \eqn{t_{(0.975, n-1)}} is the quantiles of student \eqn{t} distribution with 95\% significance and \eqn{(M-1)} degrees of freedom.
#'
#' @import iterators
#' @import foreach
#' @import doParallel
#'
#' @return \code{W.est} The \eqn{W} matrix estimated by bootstrap.
#' @return \code{lower} Lower bound of confidence intervals.
#' @return \code{upper} Upper bound of confidence intervals.
#' @references Stine, R. (1989). An Introduction to Bootstrap Methods: Examples and Ideas. Sociological Methods & Research, 18(2-3), 243-291. <doi:10.1177/0049124189018002003>
#' @author Wenxuan Liu
#' @export
#' @examples
#' # example code
#' \donttest{
#' data <- SimulatedDataset
#'
#' k <- 4
#'
#' fit <- ssmf(data = data, k = k)
#'
#' bootstrap(data = data , k = k, H = fit$H)
#' }

bootstrap <- function(data, k, H, mtimes = 50, ncore = 2) {
  result <- list()
  n <- dim(data)[1]
  registerDoParallel(cores = ncore)
  out <- foreach(icount(mtimes)) %dopar% {
    ind <- sample(1:n, replace = TRUE, size = n)
    data_boot <- data[ind, ]
    H_boot <- H[ind,]
    res_boot <- ssmf(as.matrix(data_boot), k = k, H = H_boot, lr = 0.001)
    res_boot
  }
  W_list <- lapply(out, function(x) x$W)

  W_m <- matrix(0, k, dim(data)[2])
  for (i in 1:mtimes) {
    W_m <- W_m + W_list[[i]]
  }
  W_mMean <- W_m / mtimes


  W_mVar <- matrix(0, k, dim(data)[2])
  for (m in 1:mtimes) {
    W_mVar <- W_mVar + (W_list[[m]] - W_mMean)^2
  }
  W_mSE <- sqrt((W_mVar) / (mtimes-1))

  upper <- W_mMean + 1.96 * W_mSE
  lower <- W_mMean - 1.96 * W_mSE

  result$W.est <- W_mMean
  result$lower <- lower
  result$upper <- upper

  return(result)
}


#Soft adjusted Rand index
#' @title Soft adjusted Rand index.
#' @description Soft adjusted Rand index, a soft agreement measure for class partitions incorporating assignment probabilities.
#' @param partition1 Numeric matrix/data frame of the probabilities of assignment of observations in partition 1 (membership matrix).
#' @param partition2 Numeric matrix/data frame of the probabilities of assignment of observations in partition 2 (membership matrix).
#'
#'
#' @return Soft adjusted Rand index.
#' @references Flynt, A., Dean, N. & Nugent, R. (2019) sARI: a soft agreement measure for class partitions incorporating assignment probabilities. Adv Data Anal Classif 13, 303–323 (2019). <doi:10.1007/s11634-018-0346-x>
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


#SimulatedDataset
#' A simulated metabolomic dataset.
#'
#' A simulated metabolomic data set containing 138 variables for 177 individuals.
#'
#' @format A data frame with 177 rows and 138 columns.
#'
#' @usage data(SimulatedDataset)
"SimulatedDataset"

#SimulatedMemberships
#' A simulated membership matrix.
#'
#' A simulated membership matrix containing 4 cluster memberships for 177 individuals.
#'
#' @format A data frame with 177 rows and 4 columns.
#'
#' @usage data(SimulatedMemberships)
"SimulatedMemberships"

#SimulatedPrototypes
#' A simulated prototype matrix.
#'
#' A simulated prototype matrix containing 4 cluster prototypes.
#'
#' @format A data frame with 4 rows and 138 columns.
#'
#' @usage data(SimulatedPrototypes)
"SimulatedPrototypes"

#fit_SSMF
#' Example results of SSMF.
#'
#' A list of the results for SSMF example for \eqn{k=1, 2, ..., 10}.
#'
#' @format A list with 10 items, each item is a results of SSMF,
#' containing the values of the estimated prototype matrix (\eqn{W}) and
#' the estimated membership matrix (\eqn{H}) matrix and the value of
#' the residuals sum of square (SSE).
#'
#' @usage fit_SSMF
"fit_SSMF"

#fit_gap
#' Example results of gap statistic.
#'
#' A list of the results for gap statistic example for \eqn{k=1, 2, ..., 10}.
#'
#' @format A list of gap statistic result, including the gap value vector,
#' the optimal number of prototypes/clusters and the Standard error vector.
#'
#' @usage fit_gap
"fit_gap"

#fit_boot
#' Example results of bootstrap.
#'
#' A list of the results for bootstrap example.
#'
#' @format A list of bootstrap result, including the values of estimated prototype matrix (\eqn{W}),
#' the lower bound of confidence intervals and the upper bound of confidence intervals.
#'
#' @usage fit_boot
"fit_boot"
