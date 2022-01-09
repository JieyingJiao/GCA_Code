library(lme4)
library(nlme)
library(splines2)
library(matrixStats)
library(mvtnorm)
library(sparseMVN)
library(ggplot2)

##' normal cdf function
f6 <- function(x, a, b, mean, sd) {
  y <- a * pnorm(x, mean = mean, sd = sd) + b
}

powerfunc1 <- function(x, power, s, c) {
  y <- s * x^power + c
  y <- pmax(y, 0)
}


## independeent data results
## output111, 110
LMsimusummary <- function(datafile, alpha = 0.05, beta = c(1, 1, 1), order) {
  coef.naive <- matrix(0, nrow = length(beta), ncol = 1000)
  coef.iter <- coef.naive
  sehat.naive <- coef.naive
  sehat.paraboot <- coef.naive
  cover.naive <- coef.naive
  cover.paraboot <- coef.naive
  count <- 0
  for (i in 1:1000) {
    path <- paste0("/Users/jieyingjiao/Desktop/PIproject/simulation/output/", 
                   datafile, "/lmm_", i, ".rdata")
    t <- try(load(path), TRUE)
    if (t == "output") {
      output <- output[[order]]
      if (class(output) == "list") {
        count <- count+1
        coef.naive[, count] <- coef(output$fit$mod.0)
        coef.iter[, count] <- coef(output$fit$mod.1)
        sehat.naive[, count] <- summary(output$fit$mod.0)$coef[, 2]
        sehat.paraboot[, count] <- output$bootsd$bootsd.beta
        lower.naive <- coef.naive[, count]-qnorm(alpha/2, lower.tail = FALSE)*sehat.naive[, count]
        upper.naive <- coef.naive[, count]+qnorm(alpha/2, lower.tail = FALSE)*sehat.naive[, count]
        lower.paraboot <- coef.iter[, count]-qnorm(alpha/2, lower.tail = FALSE)*sehat.paraboot[, count]
        upper.paraboot <- coef.iter[, count]+qnorm(alpha/2, lower.tail = FALSE)*sehat.paraboot[, count]
        for (j in 1:length(beta)) {
          if (lower.naive[j] <= beta[j] && upper.naive[j] >= beta[j]) cover.naive[j, count] <- 1
          if (lower.paraboot[j] <= beta[j] && upper.paraboot[j] >= beta[j]) cover.paraboot[j, count] <- 1
        }
      }
    }
  }
  coefnaive <- rowMeans(coef.naive[, 1:count])
  biasnaive <- coefnaive - 1
  coefiter <- rowMeans(coef.iter[, 1:count])
  biasiter <- coefiter - 1
  se.naive <- rowSds(coef.naive[, 1:count])
  se.iter <- rowSds(coef.iter[, 1:count])
  hat.naive <- rowMeans(sehat.naive[, 1:count])
  hat.paraboot <- rowMeans(sehat.paraboot[, 1:count])
  cover.naive <- rowMeans(cover.naive[, 1:count])
  cover.paraboot <- rowMeans(cover.paraboot[, 1:count])
  naive <- data.frame(coefnaive, biasnaive, se.naive, hat.naive, cover.naive)
  parabootout <- data.frame(coefiter, biasiter, se.iter, hat.paraboot, cover.paraboot)
  return(list(naive = naive, parabootout = parabootout,
              count = count))
}

LMerrsdSum <- function(x, datafile, degree, order, p = 0.05) {
  y <- matrix(0, nrow = length(x), ncol = 1000)
  count <- 0
  for (i in 1:1000) {
    path <- paste0("/Users/jieyingjiao/Desktop/PIproject/simulation/output/", 
                   datafile, "/lmm_", i, ".rdata")
    t <- try(load(path), TRUE)
    if (t == "output") {
      count <- count + 1
      output <- output[[order]]
      df <- length(output$fit$sfit$par) - 1
      isMat <- iSpline(output$fit$bounds, df = df, degree = degree, 
                       intercept = TRUE)
      basis <- predict(isMat, x)
      y[, count] <- cbind(1, basis) %*% output$fit$sfit$par * 
        sigma(output$fit$mod.1)
    }
  }
  errsd.fit <- rowMeans(y[, 1:count])
  errsd.quantile <- apply(y[, 1:count], MARGIN = 1, FUN = quantile, prob = c(p/2, 1-p/2))
  return(list(Estimate = errsd.fit, Quantile = errsd.quantile, y = y[, 1:count]))
}

summary_stat <- LMsimusummary(datafile = 'output110', order = 1)

x <- seq(1, 4, by = 0.01)
err_est1 <- LMerrsdSum(x  = x, datafile = 'output111', degree = 1, order = 1)
err_est2 <- LMerrsdSum(x  = x, datafile = 'output111', degree = 1, order = 2)
err_est3 <- LMerrsdSum(x  = x, datafile = 'output111', degree = 1, order = 3)
err_est4 <- LMerrsdSum(x  = x, datafile = 'output111', degree = 1, order = 4)
err_est5 <- LMerrsdSum(x  = x, datafile = 'output110', degree = 1, order = 1)
err_est6 <- LMerrsdSum(x  = x, datafile = 'output110', degree = 1, order = 2)

y_true1 <- 0.25 * powerfunc1(x = x, power = 1, s = 1, c = -0.9)
y_true2 <- 0.02 * powerfunc1(x = x, power = 3, s = 1, c = 1.2)
y_true3 <- 0.1 * f6(x = x, a = 5, b = 1, mean = 2, sd = 0.3)

fig1_data <- data.frame(X = x, True = y_true1, Estimate = err_est1$Estimate,
                        Lower = err_est1$Quantile[1, ], Upper = err_est1$Quantile[2, ],
                        setting = 'n=100, g=g1')
fig1_data <- rbind(fig1_data,
                   data.frame(X = x, True = y_true1, Estimate = err_est2$Estimate,
                              Lower = err_est2$Quantile[1, ], Upper = err_est2$Quantile[2, ],
                              setting = 'n=200, g=g1'))
fig1_data <- rbind(fig1_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est3$Estimate,
                              Lower = err_est3$Quantile[1, ], Upper = err_est3$Quantile[2, ],
                              setting = 'n=100, g=g2'))
fig1_data <- rbind(fig1_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est4$Estimate,
                              Lower = err_est4$Quantile[1, ], Upper = err_est4$Quantile[2, ],
                              setting = 'n=200, g=g2'))
fig1_data <- rbind(fig1_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est5$Estimate,
                              Lower = err_est5$Quantile[1, ], Upper = err_est5$Quantile[2, ],
                              setting = 'n=100, g=g3'))
fig1_data <- rbind(fig1_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est6$Estimate,
                              Lower = err_est6$Quantile[1, ], Upper = err_est6$Quantile[2, ],
                              setting = 'n=200, g=g3'))
save(fig1_data, file = '/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig1Data.RData')

load('/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig1Data.RData')
p <- ggplot(data = fig1_data, aes(x = X, y = True)) + 
  geom_line(linetype = "dashed", show.legend = FALSE) + geom_line(aes(y = Estimate)) + 
  facet_wrap(facets = ~ setting, labeller = "label_both", scales = "free_y") + 
  geom_ribbon(aes(min = Lower, max = Upper), fill = "grey", alpha = 0.5) + 
  ylab(expression(sigma(epsilon))) +
  xlab('')
  theme(text = element_text(size = 15), axis.title = element_text(size = 15), axis.text.x=element_blank()) + 
  theme_bw()
ggsave('../manuscript/plots/errplotlm.pdf', plot = p, width = 10, height = 5)



## clustered data
LMMerrsdSum <- function(x, datafile, degree, order, p = 0.05) {
  y <- matrix(0, nrow = length(x), ncol = 1000)
  count <- 0
  for (i in 1:1000) {
    path <- paste0("/Users/jieyingjiao/Desktop/PIproject/simulation/output/", 
                   datafile, "/lmm_", i, ".rdata")
    t <- try(load(path), TRUE)
    if (t == "output") {
      count <- count + 1
      output <- output[[order]]
      df <- length(output$fit$model.bic$sfit$par) - 1
      isMat <- iSpline(output$fit$model.bic$bounds, df = df, degree = degree, 
                       intercept = TRUE)
      basis <- predict(isMat, x)
      y[, count] <- cbind(1, basis) %*% output$fit$model.bic$sfit$par * 
        sigma(output$fit$model.bic$mixmod.1)
    }
  }
  errsd.fit <- rowMeans(y[, 1:count])
  errsd.quantile <- apply(y[, 1:count], MARGIN = 1, FUN = quantile, prob = c(p/2, 1-p/2))
  return(list(Estimate = errsd.fit, Quantile = errsd.quantile, y = y[, 1:count]))
}

sdRsimu <- function(datafile, order, sdR, alpha =  0.05) {
  random.naive <- rep(0, 1000)
  random.iter <- random.naive
  random.lme <- random.naive
  sehat.iter <- random.naive
  cover.iter <- rep(0, 1000)
  count <- 0
  count.lme <- 0
  for (i in 1:1000) {
    path <- paste0("/Users/jieyingjiao/Desktop/PIproject/simulation/output/", 
                   datafile, "/lmm_", i, ".rdata")
    t <- try(load(path), TRUE)
    if (t == "output") {
      output1 <- output[[paste0(order, ".1")]]
      count <- count+1
      random.naive[count] <- data.frame(VarCorr(output1$fit$model.bic$mixmod.0))$sdcor[1]
      random.iter[count] <- data.frame(VarCorr(output1$fit$model.bic$mixmod.1))$sdcor[1]
      sehat.iter[count] <- output1$BICsd.paraboot$bootsd.sdR
      lower.iter <- sdR - qnorm(alpha/2, lower.tail = FALSE) * sehat.iter[count]
      upper.iter <- sdR + qnorm(alpha/2, lower.tail = FALSE) * sehat.iter[count]
      cover.iter[count] <- as.numeric(random.iter[count] >= lower.iter & 
                                        random.iter[count] <= upper.iter)
      output2 <- output[[paste0(order, ".2")]]
      if (class(output2) != "try-error") {
        count.lme <- count.lme + 1
        random.lme[count.lme] <- as.numeric(VarCorr(output2)[1, 2])
      }
    }
  }
  fit.naive <- mean(random.naive[1:count])
  fit.iter <- mean(random.iter[1:count])
  fit.lme <- mean(random.lme[1:count.lme])
  se.naive <- sd(random.naive[1:count])
  se.iter <- sd(random.iter[1:count])
  se.lme <- sd(random.lme[1:count.lme])
  sehat.iter <- mean(sehat.iter[1:count])
  cover.iter <- mean(cover.iter[1:count])
  naive <- data.frame(fit.naive, se.naive)
  iter <- data.frame(fit.iter, se.iter, sehat.iter, cover.iter)
  lme <- data.frame(fit.lme, se.lme)
  return(list(naive = naive, iter = iter, lme = lme, count = count, count.lme = count.lme))
}

## random effects not includede (marginal mean)
sdR_results1 <- sdRsimu(datafile = 'output118', order = 'output1', sdR = 0.1, alpha = 0.05)
sdR_results2 <- sdRsimu(datafile = 'output118', order = 'output2', sdR = 0.1, alpha = 0.05)
sdR_results3 <- sdRsimu(datafile = 'output118', order = 'output3', sdR = 0.1, alpha = 0.05)
sdR_results4 <- sdRsimu(datafile = 'output118', order = 'output4', sdR = 0.1, alpha = 0.05)
sdR_results5 <- sdRsimu(datafile = 'output119', order = 'output5', sdR = 0.1, alpha = 0.05)
sdR_results6 <- sdRsimu(datafile = 'output119', order = 'output6', sdR = 0.1, alpha = 0.05)

x <- seq(1, 7, by = 0.01)
err_est1 <- LMMerrsdSum(x  = x, datafile = 'output118', degree = 1, order = 1)
err_est2 <- LMMerrsdSum(x  = x, datafile = 'output118', degree = 1, order = 3)
err_est3 <- LMMerrsdSum(x  = x, datafile = 'output118', degree = 1, order = 5)
err_est4 <- LMMerrsdSum(x  = x, datafile = 'output118', degree = 1, order = 7)
err_est5 <- LMMerrsdSum(x  = x, datafile = 'output119', degree = 1, order = 1)
err_est6 <- LMMerrsdSum(x  = x, datafile = 'output119', degree = 1, order = 3)

y_true1 <- 0.2 * powerfunc1(x = x, power = 1, s = 1, c = 1)
y_true2 <- 0.1 * powerfunc1(x = x, power = 3, s = 0.08, c = 2)
y_true3 <- 0.2 * f6(x = x, a = 5, b = 2, mean = 4, sd = 0.6)

fig2_data <- data.frame(X = x, True = y_true1, Estimate = err_est1$Estimate,
                        Lower = err_est1$Quantile[1, ], Upper = err_est1$Quantile[2, ],
                        setting = 'n=100, g=g1')
fig2_data <- rbind(fig2_data,
                   data.frame(X = x, True = y_true1, Estimate = err_est2$Estimate,
                              Lower = err_est2$Quantile[1, ], Upper = err_est2$Quantile[2, ],
                              setting = 'n=200, g=g1'))
fig2_data <- rbind(fig2_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est3$Estimate,
                              Lower = err_est3$Quantile[1, ], Upper = err_est3$Quantile[2, ],
                              setting = 'n=100, g=g2'))
fig2_data <- rbind(fig2_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est4$Estimate,
                              Lower = err_est4$Quantile[1, ], Upper = err_est4$Quantile[2, ],
                              setting = 'n=200, g=g2'))
fig2_data <- rbind(fig2_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est5$Estimate,
                              Lower = err_est5$Quantile[1, ], Upper = err_est5$Quantile[2, ],
                              setting = 'n=100, g=g3'))
fig2_data <- rbind(fig2_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est6$Estimate,
                              Lower = err_est6$Quantile[1, ], Upper = err_est6$Quantile[2, ],
                              setting = 'n=200, g=g3'))
save(fig2_data, file = '/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig2Data.RData')

load('/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig2Data.RData')
p <- ggplot(data = fig2_data, aes(x = X, y = True)) + 
  geom_line(linetype = 'dashed', show.legend = FALSE) + geom_line(aes(y = Estimate)) + 
  facet_wrap(facets = ~ setting, labeller = "label_both", scales = "free_y") + 
  geom_ribbon(aes(min = Lower, max = Upper), fill = "grey", alpha = 0.5) + 
  ylab(expression(sigma(epsilon))) + 
  xlab('') + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15)) + 
  theme_bw()
ggsave('../manuscript/plots/errplot1.pdf', plot = p, width = 10, height = 5)



## random effects included (conditional mean)
sdR_results1 <- sdRsimu(datafile = 'output115', order = 'output7', sdR = 0.1, alpha = 0.05)
sdR_results2 <- sdRsimu(datafile = 'output115', order = 'output8', sdR = 0.1, alpha = 0.05)
sdR_results3 <- sdRsimu(datafile = 'output116', order = 'output9', sdR = 0.1, alpha = 0.05)
sdR_results4 <- sdRsimu(datafile = 'output116', order = 'output10', sdR = 0.1, alpha = 0.05)
sdR_results5 <- sdRsimu(datafile = 'output119', order = 'output11', sdR = 0.1, alpha = 0.05)
sdR_results6 <- sdRsimu(datafile = 'output119', order = 'output12', sdR = 0.1, alpha = 0.05)

x <- seq(1, 7, by = 0.01)
err_est1 <- LMMerrsdSum(x  = x, datafile = 'output115', degree = 1, order = 1)
err_est2 <- LMMerrsdSum(x  = x, datafile = 'output115', degree = 1, order = 3)
err_est3 <- LMMerrsdSum(x  = x, datafile = 'output116', degree = 1, order = 1)
err_est4 <- LMMerrsdSum(x  = x, datafile = 'output116', degree = 1, order = 3)
err_est5 <- LMMerrsdSum(x  = x, datafile = 'output119', degree = 1, order = 5)
err_est6 <- LMMerrsdSum(x  = x, datafile = 'output119', degree = 1, order = 7)

y_true1 <- 0.2 * powerfunc1(x = x, power = 1, s = 1, c = 1)
y_true2 <- 0.1 * powerfunc1(x = x, power = 3, s = 0.08, c = 2)
y_true3 <- 0.2 * f6(x = x, a = 5, b = 2, mean = 4, sd = 0.6)

fig3_data <- data.frame(X = x, True = y_true1, Estimate = err_est1$Estimate,
                        Lower = err_est1$Quantile[1, ], Upper = err_est1$Quantile[2, ],
                        setting = 'n=100, g=g1')
fig3_data <- rbind(fig3_data,
                   data.frame(X = x, True = y_true1, Estimate = err_est2$Estimate,
                              Lower = err_est2$Quantile[1, ], Upper = err_est2$Quantile[2, ],
                              setting = 'n=200, g=g1'))
fig3_data <- rbind(fig3_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est3$Estimate,
                              Lower = err_est3$Quantile[1, ], Upper = err_est3$Quantile[2, ],
                              setting = 'n=100, g=g2'))
fig3_data <- rbind(fig3_data,
                   data.frame(X = x, True = y_true2, Estimate = err_est4$Estimate,
                              Lower = err_est4$Quantile[1, ], Upper = err_est4$Quantile[2, ],
                              setting = 'n=200, g=g2'))
fig3_data <- rbind(fig3_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est5$Estimate,
                              Lower = err_est5$Quantile[1, ], Upper = err_est5$Quantile[2, ],
                              setting = 'n=100, g=g3'))
fig3_data <- rbind(fig3_data,
                   data.frame(X = x, True = y_true3, Estimate = err_est6$Estimate,
                              Lower = err_est6$Quantile[1, ], Upper = err_est6$Quantile[2, ],
                              setting = 'n=200, g=g3'))
save(fig3_data, file = '/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig3Data.RData')

load('/Users/jieyingjiao/Desktop/PIproject/simulation/output/fig3Data.RData')
p <- ggplot(data = fig3_data, aes(x = X, y = True)) + 
  geom_line(linetype = 'dashed', show.legend = FALSE) + geom_line(aes(y = Estimate)) + 
  facet_wrap(facets = ~ setting, labeller = "label_both", scales = "free_y") + 
  geom_ribbon(aes(min = Lower, max = Upper), fill = "grey", alpha = 0.5) + 
  ylab(expression(sigma(epsilon))) + 
  xlab('') +
  theme(text = element_text(size = 15), axis.title = element_text(size = 15)) + 
  theme_bw()
ggsave('../manuscript/plots/errplot2.pdf', plot = p, width = 10, height = 5)



