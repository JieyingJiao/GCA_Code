## likelihood of mixed-effects model when error variance change with conditional mean
## use numerical integration

## input: object from function iter()
LL_func <- function(t, formula, data, degree) {
  n.random <- getME(t$mixmod.1, name = "k")
  parameter <- list(para.fix = fixef(t$mixmod.1), 
                    para.sdR = data.frame(VarCorr(t$mixmod.1))$
                      sdcor[1:n.random],
                    para.sdE = data.frame(VarCorr(t$mixmod.1))$
                      sdcor[n.random+1],
                    para.is = t$sfit$par)
  lmod <- lFormula(formula=formula, data=data, REML=FALSE)
  y <- model.response(lmod$fr)
  X <- lmod$X
  Zt <- lmod$reTrms$Zt
  df <- length(t$sfit$par) - 1
  bounds <- t$bounds
  # covR <- diag(parameter$para.sdR^2, nrow = dim(Zt)[1])
  isMat <- iSpline(bounds, df = df, degree = degree, intercept = TRUE)
  
  ## use the product of numerical integration, since random effects are independent of each other
  ll_uni <- rep(0, dim(Zt)[1])
  j <- 1
  ## record dimension of X_uni
  for(i in levels(data$Chick)) {
    ## data for i_th chick
    y_uni <- data$weight[data$Chick == strtoi(i)]
    Z_uni <- matrix(data$Time[data$Chick == strtoi(i)], ncol = 1)
    X_uni <- X[data$Chick == strtoi(i), ]
    ll_uni[j] <- integrate(ll_int_uni, lower = -Inf, upper = Inf, y_uni = y_uni, 
                           X_uni = X_uni, Z_uni = Z_uni, sdR = parameter$para.sdR, 
                           para.fix = parameter$para.fix, isMat = isMat, 
                           para.sdE = parameter$para.sdE, para.is = parameter$para.is, 
                           bounds = bounds)$value
    j <- j+1
  }
  
  ll <- sum(log(ll_uni))
  
  return(ll)
}


## define the integrand function for univariate situation
## b could be a vector
ll_int_uni <- function(b, y_uni, X_uni, Z_uni, sdR, para.fix, isMat, para.sdE, para.is, bounds) {
  ll_int <- rep(0, length(b))
  for(i in 1:length(b)) {
    fitmean <- X_uni %*% para.fix + b[i] * Z_uni
    fitmean <- pmax(fitmean, bounds[1])
    fitmean <- pmin(fitmean, bounds[2])
    
    basis <- predict(isMat, fitmean[, 1])
    covE <- diag(as.vector(para.sdE^2 * (cbind(1, basis) %*% para.is)^2))
    ll_int[i] <- dmvnorm(y_uni, fitmean[, 1], covE, log=FALSE)*dnorm(b[i], 0, sdR, log=FALSE)
  }
  
  
  
  return(ll_int)
}
