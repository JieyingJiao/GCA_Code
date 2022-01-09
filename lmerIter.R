library(lme4)
library(nlme)
library(splines2)
library(matrixStats)
library(mvtnorm)
library(sparseMVN)

##' MLE and ispline to fit the variance structure of y
##' on sigma scale
isplineMLE <- function(y, x, df, degree, isMat) {
  basis <- predict(isMat, x)
  basis <- cbind(1, basis)
  neloglh <- function(beta, y, basis) {
    errvar <- (basis %*% beta)^2
    sum(log(2*pi*errvar)/2 + y^2/2/errvar)
  }
  startbeta <- c(sqrt(mean(y^2)), rep(0.001, df))
  fit <- constrOptim(theta = startbeta, f = neloglh, y = y, basis = basis,
                     grad = NULL, ui =diag(1, df+1), ci = rep(0, df+1),
                     control = list(maxit = 5000))
  return(fit)
}

##' compute AIC and BIC
aic <- function(loglh, k) 2*k-2*loglh
bic <- function(loglh, k, n) k*log(n) - 2*loglh

##' self-defined sampling with cluster function
sampleCluster <- function(data, clustername, size, replace = TRUE) {
  cluster <- data[, which(names(data) == clustername)]
  ind <- sample(unique(cluster), size = size, replace = replace)
  subset.cluster <- function(x) subset(data, cluster == x)
  return(do.call(rbind, lapply(ind, subset.cluster)))
}


##' iteration method to fit a model with specified ispline df
##' fix the boundries of splines basis
iter <- function(formula, data = NULL, df, degree, tolerance = 1e-4, 
                 random.include = FALSE) {
  mixmod.0 <- lmer(formula=formula, data=data, REML=FALSE)
  X <- getME(mixmod.0, name = "X")
  fiterr <- resid(mixmod.0)
  if (!random.include) fitmean <- X %*% fixef(mixmod.0) else 
    fitmean <- fitted(mixmod.0)
  response <- fiterr + fitted(mixmod.0)
  dist.fix <- rep(10, 10)
  para.fix <- fixef(mixmod.0)
  
  ## extend the range of mean as 1.5 times of the naive fit mean range
  lower <- min(fitmean) - 0.1 * diff(range(fitmean))
  upper <- max(fitmean) + 0.1 * diff(range(fitmean))
  isMat <- iSpline(c(lower, upper), df = df, degree = degree, intercept = TRUE)
  
  for (i in 1:10) {
    sfit <- isplineMLE(fiterr, fitmean, df = df, degree = degree, isMat = isMat)
    basis <- predict(isMat, fitmean)
    var.fit <- (cbind(1, basis) %*% sfit$par)^2
    w <- 1/as.vector(var.fit)
    w <- pmin(w, rep(1e8, length(w)))
    mixmod.1 <- eval(substitute(lmer(formula=formula, data=data, REML=FALSE, 
                                     weights=dummy),list(dummy=w)))
    para.fixnew <- fixef(mixmod.1)
    dist.fix[i] <- sqrt(sum((para.fixnew - para.fix)^2))
    if (dist.fix[i] <= tolerance) break
    para.fix <- para.fixnew
    if (!random.include) fitmean <- X %*% fixef(mixmod.1) else 
      fitmean <- fitted(mixmod.1)
    fitmean <- pmax(fitmean, lower)
    fitmean <- pmin(fitmean, upper)
    if (!random.include) fiterr <- response -(fitmean + fitted(mixmod.1) - X %*% fixef(mixmod.1)) else
      fiterr <- response - fitmean
  }
  
  
  return(list(mixmod.0 = mixmod.0,
              mixmod.1 = mixmod.1,
              dist.fix = dist.fix[1:i],
              sfit = sfit, 
              bounds = c(lower, upper)))
}

##' get fitted error variance
##' input object from iter function
errVar <- function(object, degree, random.include = FALSE){
  X <- getME(object$mixmod.1, name = "X")
  df <- length(object$sfit$par) - 1
  if(!random.include) fitmean <- X %*% fixef(object$mixmod.1) else
    fitmean <- fitted(object$mixmod.1)
  fitmean <- pmax(fitmean, object$bounds[1])
  fitmean <- pmin(fitmean, object$bounds[2])
  isMat <- iSpline(object$bounds, df = df, degree = degree, intercept = TRUE)
  basis <- predict(isMat, fitmean)
  var.fit <- data.frame(VarCorr(object$mixmod.1))$vcov[2] * (cbind(1, basis) %*% 
                                                               object$sfit$par)^2
  return(var.fit)
}

##' select df
fitSelect <- function(formula, data = NULL, df, degree, tolerance = 1e-4,
                      random.include = FALSE) {
  df.n <- length(df)
  aic.n <- rep(Inf, df.n)
  bic.n <- rep(Inf, df.n)
  df.aic <- rep(0, df.n)
  df.bic <- rep(0, df.n)
  loglh <- rep(0, df.n)
  fit <- vector("list", df.n)
  lmod <- lFormula(formula=formula, data=data, REML=FALSE)
  y <- model.response(lmod$fr)
  X <- lmod$X
  Zt <- lmod$reTrms$Zt
  count <- 0
  for (i in 1:df.n) {
    t <- try(iter(formula=formula, data=data, df=df[i], degree=degree,
                  tolerance = tolerance, random.include = random.include), TRUE)
    if (is.list(t)) {
      count <- count + 1
      n.random <- getME(t$mixmod.1, name = "k")
      fit[[count]] <- t
      df.aic[count] <- df[i]
      df.bic[count] <- df[i]
      parameter <- list(para.fix = fixef(t$mixmod.1), 
                        para.sdR = data.frame(VarCorr(t$mixmod.1))$
                          sdcor[1:n.random],
                        para.sdE = data.frame(VarCorr(t$mixmod.1))$
                          sdcor[n.random+1],
                        para.is = t$sfit$par)
      covR <- parameter$para.sdR^2 * crossprod(Zt)
      covE <- diag(as.vector(errVar(t, degree = degree, random.include = random.include)))
      ch <- Cholesky(covR + covE)
      loglh[count] <- dmvn.sparse(y, mu=X %*% parameter$para.fix, CH=ch, 
                                  prec=FALSE)
      aic.n[count] <- aic(loglh[count], sum(lengths(parameter)))
      bic.n[count] <- bic(loglh[count], sum(lengths(parameter)), length(y))
    }
  }
  if (count == 0) {
    print("unable to fit")
    return("unable to fit")
  } else {
    return(list(df.aic = df.aic[order(aic.n[1:count])[1]],
                df.bic = df.bic[order(bic.n[1:count])[1]],
                loglh.aic = loglh[order(aic.n[1:count])[1]],
                loglh.bic = loglh[order(bic.n[1:count])[1]],
                model.aic = fit[[order(aic.n[1:count])[1]]],
                model.bic = fit[[order(bic.n[1:count])[1]]]))
  }
}

##' non-para bootstrap method to compute sd
## bootsd <- function(formula, data = NULL, df, degree, nsim, iter.n = 5,
##                    random.include = FALSE) {
##   lmod <- lFormula(formula = formula, data = data, REML = FALSE)
##   cluster <- lmod$fr[, ncol(lmod$fr)]
##   cluster.name <- names(lmod$fr)[ncol(lmod$fr)]
##   beta <- matrix(0, nrow = ncol(lmod$X), ncol = nsim)
##   sdR <- rep(0, nsim)
##   sdE <- rep(0, nsim)
##   theta <- matrix(0, nrow = df+1, ncol = nsim)
##   count <- 0

##   if (nsim > 1) {
##     for (i in 1:nsim) {
##       data.boot <- sampleCluster(data = data, cluster.name, 
##                                  length(unique(cluster)))
##       t <- try(iter(formula = formula, data = data.boot, df = df, 
##                     degree = degree, iter.n = iter.n, 
##                     random.include = random.include), TRUE)
##       if (is.list(t)) {
##         beta[, count+1] <- fixef(t$mixmod.1)
##         sdR[count+1] <- data.frame(VarCorr(t$mixmod.1))$sdcor[1]
##         sdE[count+1] <- data.frame(VarCorr(t$mixmod.1))$sdcor[2]
##         theta[, count+1] <- t$sfit$par
##         count <- count + 1
##       }
##     }
##     bootsd.beta <- rowSds(beta[, 1:count])
##     bootsd.sdR <- sd(sdR[1:count])
##     bootsd.sdE <- sd(sdE[1:count])
##     bootsd.theta <- rowSds(theta[, 1:count])
##     return(list(bootsd.beta = bootsd.beta,
##                 bootsd.sdR = bootsd.sdR,
##                 bootsd.sdE = bootsd.sdE,
##                 bootsd.theta = bootsd.theta,
##                 rate = count/nsim))
##   } else {
##     print("bootstrap replication is smaller than 2")
##     return(NA)
##   }
## }


##' parametric bootstrap method to compute sd
##' input iter object
bootsd.para <- function(formula, object, data, response, ID, degree, 
                        tolerance = 1e-4, random.include = FALSE, nsim) {
  X <- getME(object$mixmod.1, name = "X")
  beta <- matrix(0, nrow = ncol(X), ncol = nsim)
  sdR <- rep(0, nsim)
  sdE <- rep(0, nsim)
  df <- length(object$sfit$par) - 1
  theta <- matrix(0, nrow = df + 1, ncol = nsim)
  count <- 0
  
  varR <- data.frame(VarCorr(object$mixmod.1))$vcov[1]
  effect.fix <- X %*% fixef(object$mixmod.1)
  isMat <- iSpline(object$bounds, df = df, degree = degree, intercept = TRUE)
  
  if (nsim > 1) {
    for (i in 1:nsim) {
      effect.random <- rnorm(length(unique(data[, ID])), mean = 0, sd = sqrt(varR))
      effect.random <- rep(effect.random, data.frame(table(data[, ID]))$Freq)
      if (!random.include) fitmean <- X %*% fixef(object$mixmod.1) else
        fitmean <- X %*% fixef(object$mixmod.1) + effect.random
      fitmean <- pmax(fitmean, object$bounds[1])
      fitmean <- pmin(fitmean, object$bounds[2])
      basis <- predict(isMat, fitmean)
      varE <- data.frame(VarCorr(object$mixmod.1))$vcov[2] * (cbind(1, basis) %*% 
                                                                object$sfit$par)^2
      effect.error <- rnorm(length(varE), mean = 0, sd = sqrt(varE))
      response.boot <- effect.fix + effect.random + effect.error
      data.boot <- data
      data.boot[, response] <- response.boot
      t <- try(iter(formula = formula, data = data.boot, df = df, 
                    degree = degree, tolerance = tolerance, 
                    random.include = random.include), TRUE)
      if (is.list(t)) {
        beta[, count+1] <- fixef(t$mixmod.1)
        sdR[count+1] <- data.frame(VarCorr(t$mixmod.1))$sdcor[1]
        sdE[count+1] <- data.frame(VarCorr(t$mixmod.1))$sdcor[2]
        theta[, count+1] <- t$sfit$par
        count <- count + 1
      }
    }
    bootsd.beta <- rowSds(beta[, 1:count])
    bootsd.sdR <- sd(sdR[1:count])
    bootsd.sdE <- sd(sdE[1:count])
    bootsd.theta <- rowSds(theta[, 1:count])
    return(list(bootsd.beta = bootsd.beta,
                bootsd.sdR = bootsd.sdR,
                bootsd.sdE = bootsd.sdE,
                bootsd.theta = bootsd.theta,
                rate = count/nsim))
  } else {
    print("bootstrap replication is smaller than 2")
    return(NA)
  }
}

##' model fitting with df choice and bootstrap
lmerIter <- function(formula, response, ID, data = NULL, df, degree, nsim, 
                     tolerance = 1e-4, random.include = FALSE) {
  fit <- fitSelect(formula=formula, data=data, df=df, degree=degree, 
                   tolerance = tolerance, random.include = random.include)
  if (is.list(fit)) {
    if (fit$df.aic == fit$df.bic) {
      sd.paraboot <- bootsd.para(formula = formula, object = fit$model.aic, 
                                 data = data, response = response, ID, 
                                 degree = degree, tolerance = tolerance, 
                                 random.include = random.include, 
                                 nsim = nsim)
      
      AICsd.paraboot <- sd.paraboot
      BICsd.paraboot <- sd.paraboot
    } else {
      AICsd.paraboot <- bootsd.para(formula = formula, object = fit$model.aic, 
                                    data = data, response = response, ID,
                                    degree = degree, tolerance = tolerance, 
                                    random.include = random.include, nsim = nsim)
      BICsd.paraboot <- bootsd.para(formula = formula, object = fit$model.bic, 
                                    data = data, response = response, ID,
                                    degree = degree, tolerance = tolerance,
                                    random.include = random.include, 
                                    nsim = nsim)
      
    }
    return(list(fit = fit,
                AICsd.paraboot = AICsd.paraboot,
                BICsd.paraboot = BICsd.paraboot))
  } else return(fit)
}

##' simulate data, random error variance only depends on fixed effects
##' X1 ~ Bernoulli(p), X2 ~ U(0, u)
##' power function: var(error) = sdE^2 * (fixed)^(2*power)
##' self-defined function: var(error) = sdE^2 * f(fixed)^2
simu <- function(beta, sdE, sdR, n, k, l, u, p, random.include = FALSE, 
                 f = NULL, ...) {
  ## n: #object, k: #fixed effects level, l: #replication
  X1 <- rep(rbinom(n, size = 1, prob = p), each = k*l)
  X2 <- rep(runif(n * k, min = 0, max = u), each = l)
  ID <- rep(1:n, each = k*l)
  alpha <- rep(rnorm(n, mean = 0, sd = sdR), each = k*l)
  X <- cbind(rep(1, length(X1)), X1, X2)
  fixed <- X %*% beta
  eta <- fixed + alpha
  
  if (!is.null(f)) {
    if (!random.include){
      errsd <- sdE * f(fixed, ...)
    } else errsd <- sdE * f(eta, ...)
  } else errsd <- sdE
  
  err <- rnorm(length(errsd), mean = 0, sd = pmax(0, errsd))
  Y <- eta + err
  simudata <- data.frame(Y, ID, X1, X2)
  
  return(simudata)
}

##' power function
powerfunc1 <- function(x, power, s, c) {
  y <- s * x^power + c
  y <- pmax(y, 0)
}

##' logistic function
f5 <- function(x, a, b, location, scale) {
  y <- a * plogis(x, location = location, scale = scale) + b
}

##' normal cdf function
f6 <- function(x, a, b, mean, sd) {
  y <- a * pnorm(x, mean = mean, sd = sd) + b
}

##' Iterative method to fit lm
##' fix the boundries of splines basis
iter.lm <- function(formula, data = NULL, df, degree, tolerance = 1e-8) {
  mod.0 <- lm(formula = formula, data = data)
  fitmean <- fitted(mod.0)
  fiterr <- resid(mod.0)
  response <- fitmean + fiterr
  para <- coef(mod.0)
  dist.fix <- rep(10, 10)
  ## extend the range of mean as 1.5 times of the naive fit mean range
  lower <- min(fitmean) - 0.1 * diff(range(fitmean))
  upper <- max(fitmean) + 0.1 * diff(range(fitmean))
  isMat <- iSpline(c(lower, upper), df = df, degree = degree, intercept = TRUE)
  
  for (i in 1:10) {
    sfit <- isplineMLE(fiterr, fitmean, df = df, degree = degree, isMat = isMat)
    basis <- predict(isMat, fitmean)
    var.fit <- (cbind(1, basis) %*% sfit$par)^2
    w <- 1/as.vector(var.fit)
    w <- pmin(w, rep(1e8, length(w)))
    mod.1 <- eval(substitute(lm(formula=formula, data=data, 
                                weights=dummy),list(dummy=w)))
    para.new <- coef(mod.1)
    dist.fix[i] <- sqrt(sum((para - para.new)^2))
    if (dist.fix[i] <= tolerance) break
    para <- para.new
    fiterr <- resid(mod.1)
    fitmean <- fitted(mod.1)
    fitmean <- pmax(fitmean, lower)
    fitmean <- pmin(fitmean, upper)
    fiterr <- response - fitmean
  }
  
  return(list(mod.0 = mod.0,
              mod.1 = mod.1,
              dist.fix = dist.fix[1:i],
              sfit = sfit, 
              bounds = c(lower, upper)))
}

##' parametric bootstrap for independent data
##' input formula, iter.lm object
bootlm.para <- function(formula, object, degree, tolerance = 1e-8, nsim) {
  df <- length(object$sfit$par) - 1
  XMat <- model.matrix(object$mod.1)
  mean <- XMat %*% coef(object$mod.1)
  mean <- pmax(mean, object$bounds[1])
  mean <- pmin(mean, object$bounds[2])
  isMat <- iSpline(object$bounds, df = df, degree = degree, intercept = TRUE)
  basis <- predict(isMat, mean)
  errsd <- sigma(object$mod.1) * cbind(1, basis) %*% object$sfit$par
  
  beta <- matrix(0, nrow = ncol(XMat), ncol = nsim)
  theta <- matrix(0, nrow = df + 1, ncol = nsim)
  sdE <- rep(0, nsim)
  count <- 0
  
  if (nsim > 1) {
    for (i in 1:nsim) {
      error <- rnorm(n = length(errsd), mean = 0, sd = errsd)
      Y <- mean + error
      data.boot <- data.frame(Y, XMat)
      colnames(data.boot)[1] <- all.vars(formula)[1]
      t <- try(iter.lm(formula = formula, data = data.boot, df = df, 
                       degree = degree, tolerance = tolerance), TRUE)
      if (is.list(t)) {
        beta[, count + 1] <- coef(t$mod.1)
        theta[, count + 1] <- t$sfit$par
        sdE[count + 1] <- sigma(t$mod.1)
        count <- count + 1
      }
    }
    bootsd.beta <- rowSds(beta[, 1:count])
    bootsd.sdE <- sd(sdE[1:count])
    bootsd.theta <- rowSds(theta[, 1:count])
    return(list(bootsd.beta = bootsd.beta,
                bootsd.sdE = bootsd.sdE,
                bootsd.theta = bootsd.theta,
                rate = count/nsim))
  }else {
    print("bootstrap replication is smaller than 2")
    return(NA)
  }
}

lmIter <- function(formula, data = NULL, df, degree, tolerance = 1e-8, nsim) {
  fit <- try(iter.lm(formula = formula, data = data, df = df, 
                     degree = degree, tolerance = tolerance), TRUE)
  if (is.list(fit)) {
    bootsd <- bootlm.para(formula = formula, object = fit, degree = degree, 
                          tolerance = tolerance, nsim = nsim)
    return(list(fit = fit, bootsd = bootsd))
  } else {
    print("unable to fit")
    return(NA)
  }
}

##' simulate independent data
simu.ind3 <- function(beta, sdE, n, u, p, f = NULL, ...) {
  X1 <- rbinom(n, size = 1, p = p)
  X2 <- runif(n, min = 0, max = u)
  XMat <- cbind(1, X1, X2)
  meanlevel <- XMat %*% beta
  errsd <- sdE * f(meanlevel, ...)
  
  err <- rnorm(n, mean = 0, sd = pmax(errsd, 0))
  Y <- meanlevel + err
  
  simudata <- data.frame(Y, X1, X2)
}