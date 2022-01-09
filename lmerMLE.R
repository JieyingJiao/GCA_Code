library(lme4)
library(splines2)
library(mvtnorm)
library(sparseMVN)

source("lmerIter.R") # simulate data


##' fix model with MLE
lmerMLE <- function(formula, data = NULL, df, degree) {
  lmod <- lFormula(formula, data, REML = FALSE)
  mixmod.0 <- lmer(formula, data, REML = FALSE)
  X <- lmod$X
  n.random <- ncol(lmod$fr) - ncol(X) 
  lower <- min(X %*% fixef(mixmod.0))
  upper <- max(X %*% fixef(mixmod.0))
  lower.ad <- lower - 0.05 * (upper-lower)
  upper.ad <- upper + 0.05 * (upper-lower)
  isMat <- iSpline(c(lower.ad, upper.ad), df = df, degree = degree)
  sfix <- fixef(mixmod.0)
  svarR <- data.frame(VarCorr(mixmod.0))$vcov[1:n.random]
  sis <- rep(1, df+1)
  startpara <- c(sfix, svarR, sis)
  index <- cumsum(c(length(sfix), length(svarR), length(sis)))
  fit.mle <- constrOptim(startpara, f = neloglh, formula = formula, data = data, 
                         isMat = isMat, grad = NULL, 
                         ui = cbind(matrix(0, ncol = ncol(X), nrow=df+1+n.random), 
                                    diag(1, df+1+n.random)), 
                         ci = rep(0, df+1+n.random), control = list(maxit = 5000))
  return(list(fit = fit,
              isMat = isMat))
}


##' negative loglikelihood function
##' parameter = c(fix, varR, iSpline)
neloglh <- function(parameter, formula, data = NULL, isMat) {
  lmod <- lFormula(formula, data, REML = FALSE)
  y <- model.response(lmod$fr)
  X <- lmod$X
  Zt <- lmod$reTrms$Zt
  para.fix <- parameter[1:ncol(X)]
  para.varR <- parameter[(ncol(X)+1):ncol(lmod$fr)]
  para.is <- parameter[(ncol(lmod$fr)+1):length(parameter)]
  covR <- para.varR * crossprod(Zt) # FIX when more than one random effects
  errVar <- cbind(rep(1, length(y)), predict(isMat, X %*% para.fix))%*%para.is
  covE <- diag(as.vector(errVar))
  ch <- Cholesky(covR + covE)
  loglh <- dmvn.sparse(y, mu = X %*% para.fix, CH = ch, prec = FALSE, 
                       log = TRUE)
  return(-loglh)
}


