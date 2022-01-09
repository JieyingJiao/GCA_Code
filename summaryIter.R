## functions to get estimation from iteration method from output of lmerIter

library(lme4)
library(splines2)

##' object from lmerIter function
fitIter <- function(object, degree, random.include = FALSE) {
 if (is.list(object)) {
    ## get number of random effects terms
    n.random <- getME(object$fit$model.aic$mixmod.1, name = "k")
    ## get the design matrix for fixed effects
    X <- getME(object$fit$model.aic$mixmod.1, name = "X")
    
    vc1 <- data.frame(VarCorr(object$fit$model.aic$mixmod.1))
    vc2 <- data.frame(VarCorr(object$fit$model.bic$mixmod.1))
    para.fix1 <- fixef(object$fit$model.aic$mixmod.1)
    para.fix2 <- fixef(object$fit$model.bic$mixmod.1)
    para.sdR1 <- vc1$sdcor[1:n.random]
    para.sdR2 <- vc2$sdcor[1:n.random]
    para.sdE1 <- vc1$sdcor[n.random+1]
    para.sdE2 <- vc2$sdcor[n.random+1]
    para.is1 <- object$fit$model.aic$sfit$par
    para.is2 <- object$fit$model.bic$sfit$par
    df1 <- object$fit$df.aic
    df2 <- object$fit$df.bic
    loglh1 <- object$fit$loglh.aic
    loglh2 <- object$fit$loglh.bic
    if (is.list(object$AICsd.paraboot)) {
      bootsd.fix1 <- object$AICsd.paraboot$bootsd.beta
      bootsd.fix2 <- object$BICsd.paraboot$bootsd.beta
      bootsd.sdR1 <- object$AICsd.paraboot$bootsd.sdR
      bootsd.sdR2 <- object$BICsd.paraboot$bootsd.sdR
      bootsd.sdE1 <- object$AICsd.paraboot$bootsd.sdE
      bootsd.sdE2 <- object$BICsd.paraboot$bootsd.sdE
      bootsd.is1 <- object$AICsd.paraboot$bootsd.theta
      bootsd.is2 <- object$BICsd.paraboot$bootsd.theta
      bootsd.fix = data.frame(AIC = bootsd.fix1, BIC = bootsd.fix2)
      bootsd.sdR = data.frame(AIC = bootsd.sdR1, BIC = bootsd.sdR2)
      bootsd.sdE = data.frame(AIC = bootsd.sdE1, BIC = bootsd.sdR2)
      bootsd.is = list(AIC = bootsd.is1, BIC = bootsd.is2)
    } else {
      bootsd.fix <- NA
      bootsd.sdR <- NA
      bootsd.sdE <- NA
      bootsd.is <- NA
    }
    
    if (!random.include) {
      fitmean1 <- X %*% para.fix1
      fitmean2 <- X %*% para.fix2
    } else {
      fitmean1 <- fitted(object$fit$model.aic$mixmod.1)
      fitmean2 <- fitted(object$fit$model.bic$mixmod.1)
    }
    
    isMat1 <- iSpline(object$fit$model.aic$bounds, df = df1, degree = degree, intercept = TRUE)
    isMat2 <- iSpline(object$fit$model.bic$bounds, df = df2, degree = degree, intercept = TRUE)
    
    isMat1 <- predict(isMat1, fitmean1)
    isMat2 <- predict(isMat2, fitmean2)
    
    return(list(para.fix = data.frame(AIC = para.fix1, BIC = para.fix2),
                para.sdR = data.frame(AIC = para.sdR1, BIC = para.sdR2),
                para.sdE = data.frame(AIC = para.sdE1, BIC = para.sdE2),
                para.is = list(AIC = para.is1, BIC = para.is2),
                df = data.frame(AIC = df1, BIC = df2),
                loglh = data.frame(AIC = loglh1, BIC = loglh2),
                bootsd.fix = bootsd.fix,
                bootsd.sdR = bootsd.sdR,
                bootsd.sdE = bootsd.sdE,
                bootsd.is = bootsd.is,
                isMat = list(AIC = isMat1, BIC = isMat2)))
  } else return(object)
}

##' get confidence interval for fix effects estimation and random variance
CIIter <- function(object, degree, random.include, alpha = 0.05) {
  if (is.list(object)) {
    z <- qnorm(alpha/2, lower.tail = FALSE)
    fit <- fitIter(object, degree = degree, random.include = random.include)
    para.fix <- fit$para.fix
    bootsd.fix <- fit$bootsd.fix
    para.sdR <- fit$para.sdR
    bootsd.sdR <- fit$bootsd.sdR
    lower.fix <- para.fix - z * bootsd.fix
    upper.fix <- para.fix + z * bootsd.fix
    lower.sdR <- para.sdR - z * bootsd.sdR
    upper.sdR <- para.sdR + z * bootsd.sdR
    rnames <- c(rownames(fit$para.fix), rownames(fit$para.sdR))
    return(list(AIC = data.frame(lower = c(lower.fix$AIC, lower.sdR$AIC), 
                                   upper = c(upper.fix$AIC, upper.sdR$AIC),
                                   row.names = rnames),
                BIC = data.frame(lower = c(lower.fix$BIC, lower.sdR$BIC),
                                   upper = c(upper.fix$BIC, upper.sdR$BIC),
                                   row.names = rnames)))
    } else return(object)
}

##' error sd estimation
##' input lmerIter object
eSd <- function(object, degree, random.include = FALSE) {
  X <- getME(object$fit$model.aic$mixmod.1, name = "X")
  fitsum <- fitIter(object = object, degree = degree, 
                 random.include = random.include)
  errsd1 <- cbind(1, fitsum$isMat$AIC) %*% fitsum$para.is$AIC * 
    data.frame(VarCorr(object$fit$model.aic$mixmod.1))$sdcor[2]
  errsd2 <- cbind(1, fitsum$isMat$BIC) %*% fitsum$para.is$BIC * 
    data.frame(VarCorr(object$fit$model.bic$mixmod.1))$sdcor[2]
  if (!random.include) {
    fitmean1 <- X %*% fitsum$para.fix$AIC
    fitmean2 <- X %*% fitsum$para.fix$BIC
  } else {
    fitmean1 <- fitted(object$fit$model.aic$mixmod.1)
    fitmean2 <- fitted(object$fit$model.bic$mixmod.1)
  }
  AIC <- list(fitmean = fitmean1, esd = errsd1)
  BIC <- list(fitmean = fitmean2, esd = errsd2)
  return(list(AIC = AIC, BIC = BIC))
}



##' get residual
residIter <- function(object) {
  AIC <- residuals(object$fit$model.aic$mixmod.1)
  BIC <- residuals(object$fit$model.bic$mixmod.1)
  return(data.frame(AIC, BIC))
}

