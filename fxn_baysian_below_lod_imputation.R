# R script/function to perform Bayesian imputation for truncated missing due to values being below detection limit
library(rjags)

impute.X.BDL <- function(X.obs=NULL, LOD=NULL) {
  X.obs <- as.vector(X.obs)
  N = length(X.obs)
  X.obs <- ifelse(X.obs>LOD, yes = X.obs, no = NA)
  R <- ifelse(is.na(X.obs),yes =  1,no = 0)
  imputeX <- 
    "model {
      for(i in 1:N) {
        X.obs[i] ~ dnorm(X.true[i],prec) 
        X.true[i] <- X.notmiss[i]*(1-R[i]) + X.miss[i]*R[i]
        X.notmiss[i] ~ dnorm(mu, tau)
        X.miss[i] ~ dnorm(mu, tau)T( , LOD)
      }
      tau <- 1/(sigma*sigma)
      sigma ~ dunif(0,5)
      mu ~ dnorm(0, 1.0E-06)
      prec <- 10000
      }"
  
  jdata <- list(N=length(X.obs), X.obs=X.obs, R=R, LOD=LOD)
  var.s <- c("X.true","mu", "sigma")
  model.fit <- jags.model(file=textConnection(imputeX), data=jdata, n.chains=2, n.adapt=1000, quiet=T)
  update(model.fit, n.iter=1000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
  
  r <- summary(model.fit)
  X.est <- r$quantiles[1:N, "50%"]
  return(X.est)
}
