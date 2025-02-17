## ----setup, include = FALSE----------------------------------------------------------------
knitr::opts_chunk$set(
   collapse = TRUE,
   comment = "#>"
 )
knitr::opts_chunk$set(fig.width=7, fig.height=6)


## ----LoadPackage, include = FALSE----------------------------------------------------------
library(MetabolSSMF)
set.seed(12345)


## ----eval=FALSE----------------------------------------------------------------------------
## # Install the package
## devtools::install_github("WenxuanLiu1996/MetabolSSMF")
## 
## # Load the package
## library(MetabolSSMF)


## ----Data----------------------------------------------------------------------------------
# Simulated data set X
X <- as.matrix(SimulatedDataset)

# Simulated W matrix (prototype matrix)
W <-  as.matrix(SimulatedPrototypes)

# Simulated H matrix (membership matrix)
H <- as.matrix(SimulatedMemberships)


## ----Visulisation1, warning=FALSE----------------------------------------------------------
# Define 4 different colors
color <- c(ggsci::pal_jco("default", alpha=0.8)(7)[2], ggsci::pal_jama("default", alpha=0.8)(7)[c(3:5)])

par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
matplot(t(W), type='l', lty = 1, ylab='Prototypes', col = color, lwd=2, ylim = c(0,1))
legend('topright', inset=c(-0.15,0), legend = 1:4, fill=color[1:4])


## ----Visulisation2, warning=FALSE----------------------------------------------------------
pca_X <- prcomp(X, center = F, scale. = F) 
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
caroline::pies(lapply(apply(H, 1, list), unlist), color.table = caroline::nv(color[1:4], paste0('Pty', 1:4)), x0=pca_X$x[,1], y0=pca_X$x[,2], radii = 2, ylab='PC2', xlab='PC1')
legend('topright', inset=c(-0.15,0), legend = 1:4, fill=color[1:4])


## ----SSMF, eval=FALSE----------------------------------------------------------------------
## # Set the maximum k that is used in loop.
## K <- 10
## 
## # Run the SSMF algorithm with various k and save the results as a list
## fit_SSMF <- list()
## for(k in 1:K){
##   fit <- ssmf(X, k = k, lr = 0.001, meth='kmeans')
##   fit_SSMF[[k]] <- fit
## }


## ----include=FALSE-------------------------------------------------------------------------
K <- 10
load("~/Desktop/Projects/PhD Project/Wenxuan_Works/Github/MetabolSSMF/vignettes/results_ssmf.RData")


## ----RSS-----------------------------------------------------------------------------------
# Extract the RSS and save as a vector
rss_SSMF <- unlist(lapply(fit_SSMF, function(x) x$SSE))

# Plot RSS
plot(1:K, rss_SSMF, type="b", xlab = "K", ylab = "RSS", main='Elbow Criterion')


## ----Gap, eval=FALSE-----------------------------------------------------------------------
## # Apply the gap statistic to the simulated data
## fit_gap <- gap(X, rss = rss_SSMF)
## 
## # Visualise the results
## par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
## plot(x = c(1:K), y = fit_gap$gap, type='b', ylim=range((fit_gap$gap), (fit_gap$gap - fit_gap$standard.error)[-1]), xlab='K', ylab='GAP', main='Gap statistic')
## lines(x = c(1:(K-1)), y = (fit_gap$gap - fit_gap$standard.error)[-1], col='blue', lty=2, type='b')
## lines(x = rep(fit_gap$optimal.k, 2), y = range(fit_gap$gap), lty = 2, col='red')
## legend('topright', inset=c(-0.35,0), col=c('black', 'blue', 'red'), lty=c(1, 2, 2), legend = c('Gap', '-1 SE', 'Optimal K'), cex=0.8)


## ----Visulisation3, include=FALSE----------------------------------------------------------
load("~/Desktop/Projects/PhD Project/Wenxuan_Works/Github/MetabolSSMF/vignettes/results_gap.RData")

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x = c(1:K), y = fit_gap$gap, type='b', ylim=range((fit_gap$gap), (fit_gap$gap - fit_gap$standard.error)[-1]), xlab='K', ylab='GAP', main='Gap statistic')
lines(x = c(1:(K-1)), y = (fit_gap$gap - fit_gap$standard.error)[-1], col='blue', lty=2, type='b')
lines(x = rep(fit_gap$optimal.k, 2), y = range(fit_gap$gap), lty = 2, col='red')
legend('topright', inset=c(-0.35,0), col=c('black', 'blue', 'red'), lty=c(1, 2, 2), legend = c('Gap', '-1 SE', 'Optimal K'), cex=0.8)


## ----Bootstrap, eval=FALSE-----------------------------------------------------------------
## # Set the number of times to bootstrap
## M <- 1000
## 
## # Bootstrap
## # Initial H (membership) matrix to start the algorithm
## # Based on the results above, k=4 here
## initialisation <- init(data = X, k = 4, method = 'kmeans')
## fit_boot <- bootstrap(data = X, k = 4, H = initialisation$H, mtimes=M)


## ----include=FALSE-------------------------------------------------------------------------
load("~/Desktop/Projects/PhD Project/Wenxuan_Works/Github/MetabolSSMF/vignettes/results_bootstrap.RData")


## ----Visulisation4-------------------------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
matplot(t(fit_boot$W.est[c(2,4,1,3),]), type='n', xaxt='n', ylab='Features', ylim = range(fit_boot$lower,fit_boot$upper))
for(i in 1:4){
  alpha <- 0.2
  plat <- c(ggsci::pal_jco("default",alpha=alpha)(7)[2], ggsci::pal_jama("default", alpha=alpha)(7)[c(3:5)])
  col <- switch(i, plat[1], plat[2], plat[3], plat[4])
  polygon(c(1:138, 138:1), c(fit_boot$upper[c(2,4,1,3),][i,], rev(fit_boot$lower[c(2,4,1,3),][i,])), col = col, border= col, lty = 1)
}
matplot(t(fit_boot$W.est[c(2,4,1,3),]), type='l', lty = 1, ylab='Features', add=T, col=color, lwd=2)
legend('topright', inset=c(-0.15,0), legend = 1:4, fill=color[1:4])
axis(1, at=c(1, seq(10, 130, 10), 138), labels = c(1, seq(10, 130, 10), 138), cex=0.8)


## ----Visulisation5, warning=FALSE----------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
caroline::pies(lapply(apply(fit_SSMF[[4]]$H, 1, list), unlist), color.table = caroline::nv(color[1:4], c(4,3,2,1)), x0=pca_X$x[,1], y0=pca_X$x[,2], radii = 2, ylab='PC2', xlab='PC1') 
legend('topright', inset=c(-0.15,0), legend = 1:4, fill=color[1:4])


## ----sARI----------------------------------------------------------------------------------
# Compare the true and estimated soft membership matrix
sARI(fit_SSMF[[4]]$H, H)

# Self comparison of the true soft membership matrix
sARI(H, H)


## ----Shannon-------------------------------------------------------------------------------
E <- rep(NA, nrow(fit_SSMF[[4]]$H))
for(i in 1:nrow(fit_SSMF[[4]]$H)){
  E[i] <- diversity(fit_SSMF[[4]]$H[i,], two.power=T)
}
round(mean(E), 1)

