data("ChickWeight")

source("summaryIter.R")
source("lmerIter.R")
source('utils.R')
source('LL_mixed_conditional_mean.R')

library(car) # qqPlot function
library(ggplot2)
library(gridExtra) # grid.arrange function
library(MASS) # boxcox function
library(tidyr) # gather function
library(tidyverse)
library(sas7bdat)
library(matrixStats)

## remove outliers
ChickWeight <- ChickWeight[ChickWeight$Chick != 24, ]
# ChickWeight$Chick[ChickWeight$Chick == 50] <- 24
# ChickWeight$Chick <- droplevels(ChickWeight$Chick)

(ggplot(ChickWeight, aes(x = Time, y = weight, group = Chick))+geom_line() + 
    geom_point()+theme_bw() + facet_wrap( ~ Diet, nrow = 1, labeller = label_both) +
    ylab("Weight (g)") + xlab("Time (days)")) %>% ggsave(plot = .,
                                                         filename = "../manuscript/plots/ChickWeight.pdf",
                                                         width = 9, height= 3)

##============================================================================

## fit 4 diet groups together

##============================================================================


isMat <- iSpline(ChickWeight$Time, df = 3, degree = 0, intercept = TRUE)
## create dummy variable for Diet. Use Diet 1 as baseline
ChickWeight$D2 <- 0
ChickWeight$D3 <- 0
ChickWeight$D4 <- 0
for (i in 1:nrow(ChickWeight)) {
  if (ChickWeight$Diet[i] == 2) {
    ChickWeight$D2[i] <- 1
  }
  else if (ChickWeight$Diet[i] == 3) {
    ChickWeight$D3[i] <- 1
  }
  else if (ChickWeight$Diet[i] == 4) {
    ChickWeight$D4[i] <- 1
  }
}

## create ispline basis
isMat <- iSpline(ChickWeight$Time, df = 3, degree = 0, intercept = TRUE)
ChickWeight <- data.frame(ChickWeight, isMat)

ChickWeight$Log_weight <- log(ChickWeight$weight)

## use the model that different diet groups still share same intercept, but has different slope
chickfit0 <- lmer(weight ~ X1 + X2 + X3 + X1*D2 + X1*D3 + X1*D4 + X2*D2 + X2*D3 + X2*D4 + 
                    X3*D2 + X3*D3 + X3*D4 + (Time - 1 |Chick),
                  data = ChickWeight, REML = FALSE)
p1 <- gg_qq(resid(chickfit0), labels = NA) + theme_bw()
p1

## get naive model likelihood
chickfit0
## -2315.220


## use marginal mean
naivefit1 <- data.frame(x = getME(chickfit0, name = "X") %*% fixef(chickfit0), 
                        y = resid(chickfit0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit1 <- lowess(naivefit1$x, naivefit1$y)
p2 <- ggplot(data = naivefit1, aes(x = x, y = y)) + xlab("Fitted Mean") + 
  ylab("Residual") + geom_point() + 
  # geom_line(aes(x = lowessfit1$x, y = lowessfit1$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) + theme_bw()
p2

## use conditional mean
naivefit1.2 <- data.frame(x = fitted(chickfit0), 
                          y = resid(chickfit0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit1.2 <- lowess(naivefit1.2$x, naivefit1.2$y)
p2.2 <- ggplot(data = naivefit1.2, aes(x = x, y = y)) + xlab("Fitted Mean") + 
  ylab("Residual") + geom_point() + 
  # geom_line(aes(x = lowessfit1.2$x, y = lowessfit1.2$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15),  
        axis.title.x = element_blank())+theme_bw()
p2.2


## use the iteration method

## marginal mean model
chickfit.1 <- lmerIter(weight ~ X1 + X2 + X3 + X1*D2 + X1*D3 + X1*D4 + X2*D2 + X2*D3 + X2*D4 + 
                         X3*D2 + X3*D3 + X3*D4 + (Time - 1 | Chick), 
                       response = "weight", ID = "Chick", data = ChickWeight, 
                       df = 9, degree = 1, nsim = 2, tolerance = 1e-4, 
                       random.include = FALSE)

## df: 12:20, degree: 3
chickerr <- eSd(chickfit.1, degree = 1, random.include = FALSE)
chickresid <- data.frame(x = chickerr$BIC$fitmean, 
                         y = residIter(chickfit.1)$BIC / chickerr$BIC$esd, 
                         x1 = fitted(chickfit.1$fit$model.bic$mixmod.1))


p3 <- gg_qq2(residIter(chickfit.1)$BIC / chickerr$BIC$esd)+theme_bw()
p3

fitmean <- chickerr$BIC$fitmean

alpha <- 0.05
lower <- chickerr$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr <- data.frame(x = fitmean, y = resid(chickfit.1$fit$model.bic$mixmod.1), 
                            lower = lower, upper = upper)

p4 <- ggplot(chickresiderr, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  xlab("Fitted Marginal Mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 20), 
        axis.title = element_text(size = 20))+theme_bw()
p4

# p <- grid.arrange(p1, p2, p3, p4, nrow = 2)
# ggsave('../manuscript/plots/chickfit.pdf', p, width = 9, height = 7)

## get log likelihood, marginal situation
print(chickfit.1$fit$loglh.bic)
## -2199.183
chickfit.1$fit$df.bic
## df.bic = 9


## fit conditional mean model
chickfit.2 <- lmerIter(weight ~ X1 + X2 + X3 + X1*D2 + X1*D3 + X1*D4 + X2*D2 + X2*D3 + X2*D4 + 
                         X3*D2 + X3*D3 + X3*D4 + (Time - 1 | Chick), 
                       response = "weight", ID = "Chick", data = ChickWeight, 
                       df = 9, degree = 1, nsim = 2, tolerance = 1e-4, 
                       random.include = TRUE)

# save(chickfit.2, file = '~/Desktop/chickfit2.rdata')

chickerr.2 <- eSd(chickfit.2, degree = 1, random.include = TRUE)
chickresid.2 <- data.frame(x = chickerr.2$BIC$fitmean, 
                           y = residIter(chickfit.2)$BIC / chickerr.2$BIC$esd, 
                           x1 = fitted(chickfit.2$fit$model.bic$mixmod.1))


p3.2 <- gg_qq2(residIter(chickfit.2)$BIC / chickerr.2$BIC$esd)+theme_bw()
p3.2

fitmean.2 <- chickerr.2$BIC$fitmean

alpha <- 0.05
lower <- chickerr.2$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr.2 <- data.frame(x = fitmean.2, y = resid(chickfit.2$fit$model.bic$mixmod.1), 
                              lower = lower, upper = upper)

p4.2 <- ggplot(chickresiderr.2, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  xlab("Fitted Conditional Mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 20), 
        axis.title = element_text(size = 20))+theme_bw()
p4.2

# p <- grid.arrange(p1, p2.2, p3.2, p4.2, nrow = 2)
# ggsave('../manuscript/plots/chickfit_randominclude.pdf', p, width = 9, height = 7)

p <- grid.arrange(p1, p3, p3.2, p2, p4, p4.2, nrow = 2)
ggsave('../manuscript/plots/ChickWeight_diagnosis.pdf', p, width = 18, height = 12)


ll2 <- LL_func(chickfit.2$fit$model.bic, formula = weight ~ X1 + X2 + X3 + X1*D2 + X1*D3 + X1*D4 + X2*D2 + X2*D3 + X2*D4 + 
                 X3*D2 + X3*D3 + X3*D4 + (Time - 1 | Chick), 
               data = ChickWeight, degree = 1)
ll2
## -2164.899
chickfit.2$fit$df.bic
## df.bic = 9

## calculate the response variance, to construct fitted growth curve quantiles.

## use bootstrap sample
## input iter object and time (a full list of time points) to be evaluated
## n_object, each evaluated at timepoints of "time"
response_boot_sample <- function(object, time, n_object, degree, Diet, random.include = TRUE) {
  df <- length(object$sfit$par) - 1
  n.random <- getME(object$mixmod.1, name = "k")
  parameter <- list(para.fix = fixef(object$mixmod.1), 
                    para.sdR = data.frame(VarCorr(object$mixmod.1))$
                      sdcor[1:n.random],
                    para.sdE = data.frame(VarCorr(object$mixmod.1))$
                      sdcor[n.random+1],
                    para.is = object$sfit$par)
  
  response_sample <- matrix(0, nrow = n_object, ncol = length(time))
  
  isMat_mean <- iSpline(ChickWeight$Time, df = 3, degree = 0, intercept = TRUE)
  isMat <- iSpline(object$bounds, df = df, degree = degree, intercept = TRUE)
  varR <- data.frame(VarCorr(object$mixmod.1))$vcov[1]
  for (i in 1:n_object) {
    effect.random <- rnorm(n = 1, mean = 0, sd = sqrt(varR))
    for (j in 1:length(time)) {
      X_eval <- predict(isMat_mean, time[j])
      X_fitted <- cbind(1, X_eval[1], X_eval[2], X_eval[3], Diet[1], Diet[2], Diet[3], 
                        X_eval[1]*Diet[1], X_eval[1]*Diet[2], X_eval[1]*Diet[3],
                        X_eval[2]*Diet[1], X_eval[2]*Diet[2], X_eval[2]*Diet[3],
                        X_eval[3]*Diet[1], X_eval[3]*Diet[2], X_eval[3]*Diet[3])
      fixed_eval <- X_fitted %*% parameter$para.fix
      fixed_eval <- fixed_eval[1, 1]
      if (!random.include) fitmean <- fixed_eval else
        fitmean <- fixed_eval + effect.random * time[j]
      fitmean <- pmax(fitmean, object$bounds[1])
      fitmean <- pmin(fitmean, object$bounds[2])
      basis <- predict(isMat, fitmean)
      varE <- data.frame(VarCorr(object$mixmod.1))$vcov[2] * (cbind(1, basis) %*% 
                                                                parameter$para.is)^2
      varE <- varE[1, 1]
      effect.error <- rnorm(1, mean = 0, sd = sqrt(varE))
      response_sample[i, j] <- fixed_eval + effect.random*time[j] + effect.error
      
    }
  }
  return(response_sample)
}

time <- seq(0, 21, by = 1)
response_boot_sample_1 <- response_boot_sample(chickfit.2$fit$model.bic, time = time, 
                                               n_object = 10000, degree = 1, Diet = c(0, 0, 0))
response_boot_sample_2 <- response_boot_sample(chickfit.2$fit$model.bic, time = time, 
                                               n_object = 10000, degree = 1, Diet = c(1, 0, 0))
response_boot_sample_3 <- response_boot_sample(chickfit.2$fit$model.bic, time = time, 
                                               n_object = 10000, degree = 1, Diet = c(0, 1, 0))
response_boot_sample_4 <- response_boot_sample(chickfit.2$fit$model.bic, time = time, 
                                               n_object = 10000, degree = 1, Diet = c(0, 0, 1))



quantile_plot_1 <- data.frame(x = time, 
                              lower1 = colQuantiles(response_boot_sample_1, probs = 0.025),
                              upper1 = colQuantiles(response_boot_sample_1, probs = 0.975),
                              lower2 = colQuantiles(response_boot_sample_1, probs = 0.05),
                              upper2 = colQuantiles(response_boot_sample_1, probs = 0.95), 
                              Diet = 1)
quantile_plot_2 <- data.frame(x = time, 
                              lower1 = colQuantiles(response_boot_sample_2, probs = 0.025),
                              upper1 = colQuantiles(response_boot_sample_2, probs = 0.975),
                              lower2 = colQuantiles(response_boot_sample_2, probs = 0.05),
                              upper2 = colQuantiles(response_boot_sample_2, probs = 0.95), 
                              Diet = 2)
quantile_plot_3 <- data.frame(x = time, 
                              lower1 = colQuantiles(response_boot_sample_3, probs = 0.025),
                              upper1 = colQuantiles(response_boot_sample_3, probs = 0.975),
                              lower2 = colQuantiles(response_boot_sample_3, probs = 0.05),
                              upper2 = colQuantiles(response_boot_sample_3, probs = 0.95), 
                              Diet = 3)
quantile_plot_4 <- data.frame(x = time, 
                              lower1 = colQuantiles(response_boot_sample_4, probs = 0.025),
                              upper1 = colQuantiles(response_boot_sample_4, probs = 0.975),
                              lower2 = colQuantiles(response_boot_sample_4, probs = 0.05),
                              upper2 = colQuantiles(response_boot_sample_4, probs = 0.95), 
                              Diet = 4)
ChickWeight$fitted <- fitted(chickfit.2$fit$model.bic$mixmod.1)
ChickWeight$fixed_fit <- getME(chickfit.2$fit$model.bic$mixmod.1, name = "X") %*% fixef(chickfit.2$fit$model.bic$mixmod.1)

temp <- ChickWeight[ChickWeight$Diet == 1, ]
p1_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p1_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)


temp <- ChickWeight[ChickWeight$Diet == 2, ]
p2_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p2_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

temp <- ChickWeight[ChickWeight$Diet == 3, ]
p3_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p3_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

temp <- ChickWeight[ChickWeight$Diet == 4, ]
p4_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p4_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)


p_quantile <- grid.arrange(p1_quantile, p2_quantile, p3_quantile, p4_quantile, nrow = 1)
ggsave('../manuscript/plots/ChickWeight_quantile_bootsample.pdf', p_quantile, width = 9, height = 3)
p_quantile_fitted <- grid.arrange(p1_quantile_fitted, p2_quantile_fitted, p3_quantile_fitted, p4_quantile_fitted, nrow = 1)
ggsave('../manuscript/plots/ChickWeight_fitted_quantile_bootsample.pdf', p_quantile_fitted, width = 9, height = 3)




time <- seq(0, 21, by = 1)
varY_1 <- rep(0, length(time))
varY_2 <- rep(0, length(time))
varY_3 <- rep(0, length(time))
varY_4 <- rep(0, length(time))
fixed_1 <- rep(0, length(time))
fixed_2 <- rep(0, length(time))
fixed_3 <- rep(0, length(time))
fixed_4 <- rep(0, length(time))
for (i in 1:length(time)) {
  temp <- var_response(chickfit.2$fit$model.bic, Diet = c(0, 0, 0), time = time[i], degree = 1)
  varY_1[i] <- temp$varY
  fixed_1[i]<- temp$fixed_eval
  temp <- var_response(chickfit.2$fit$model.bic, Diet = c(1, 0, 0), time = time[i], degree = 1)
  varY_2[i] <- temp$varY
  fixed_2[i]<- temp$fixed_eval
  temp <- var_response(chickfit.2$fit$model.bic, Diet = c(0, 1, 0), time = time[i], degree = 1)
  varY_3[i] <- temp$varY
  fixed_3[i]<- temp$fixed_eval
  temp <- var_response(chickfit.2$fit$model.bic, Diet = c(0, 0, 1), time = time[i], degree = 1)
  varY_4[i] <- temp$varY
  fixed_4[i]<- temp$fixed_eval
}

q1 <- qnorm(0.05/2, lower.tail = FALSE)
q2 <- qnorm(0.1/2, lower.tail = FALSE)
quantile_plot_1 <- data.frame(x = time, 
                              lower1 = fixed_1 - q1 * sqrt(varY_1),
                              upper1 = fixed_1 + q1 * sqrt(varY_1),
                              lower2 = fixed_1 - q2 * sqrt(varY_1),
                              upper2 = fixed_1 + q2 * sqrt(varY_1), 
                              Diet = 1)
quantile_plot_2 <- data.frame(x = time, 
                              lower1 = fixed_2 - q1 * sqrt(varY_2),
                              upper1 = fixed_2 + q1 * sqrt(varY_2),
                              lower2 = fixed_2 - q2 * sqrt(varY_2),
                              upper2 = fixed_2 + q2 * sqrt(varY_2), 
                              Diet = 2)
quantile_plot_3 <- data.frame(x = time, 
                              lower1 = fixed_3 - q1 * sqrt(varY_3),
                              upper1 = fixed_3 + q1 * sqrt(varY_3),
                              lower2 = fixed_3 - q2 * sqrt(varY_3),
                              upper2 = fixed_3 + q2 * sqrt(varY_3), 
                              Diet = 3)
quantile_plot_4 <- data.frame(x = time, 
                              lower1 = fixed_4 - q1 * sqrt(varY_4),
                              upper1 = fixed_4 + q1 * sqrt(varY_4),
                              lower2 = fixed_4 - q2 * sqrt(varY_4),
                              upper2 = fixed_4 + q2 * sqrt(varY_4), 
                              Diet = 4)

quantile_plot <- rbind(quantile_plot_1, quantile_plot_2, quantile_plot_3, quantile_plot_4)

## input object: object from function iter()
## input Diet dummy variable, default is for Diet 1
## input time stamp, an integer from 0 to 21
var_response <- function(object, Diet = c(0, 0, 0), time, degree) {
  n.random <- getME(object$mixmod.1, name = "k")
  parameter <- list(para.fix = fixef(object$mixmod.1), 
                    para.sdR = data.frame(VarCorr(object$mixmod.1))$
                      sdcor[1:n.random],
                    para.sdE = data.frame(VarCorr(object$mixmod.1))$
                      sdcor[n.random+1],
                    para.is = object$sfit$par)
  bounds <- object$bounds
  df <- length(object$sfit$par) - 1
  
  isMat <- iSpline(bounds, df = df, degree = degree, intercept = TRUE)
  isMat_mean <- iSpline(ChickWeight$Time, df = 3, degree = 0, intercept = TRUE)
  
  X_eval <- predict(isMat_mean, time)
  X_fitted <- cbind(1, X_eval[1], X_eval[2], X_eval[3], Diet[1], Diet[2], Diet[3], 
                    X_eval[1]*Diet[1], X_eval[1]*Diet[2], X_eval[1]*Diet[3],
                    X_eval[2]*Diet[1], X_eval[2]*Diet[2], X_eval[2]*Diet[3],
                    X_eval[3]*Diet[1], X_eval[3]*Diet[2], X_eval[3]*Diet[3])
  fixed_eval <- X_fitted %*% parameter$para.fix
  fixed_eval <- fixed_eval[1, 1]
  
  varY <- integrate(err_var_int, lower = -Inf, upper = Inf, 
                    fixed_eval = fixed_eval, parameter = parameter, 
                    isMat = isMat, time = time, bounds = bounds)$value
  varY <- varY + time^2*parameter$para.sdR^2
  
  return (list(varY = varY, fixed_eval = fixed_eval))
  
}

## define error variance as a function of random effects and time, as integrand
err_var_int <- function(b, fixed_eval, parameter, isMat, time, bounds) {
  err_var_int_val <- rep(0, length(b))
  for (i in 1:length(b)) {
    fitmean <- fixed_eval + time*b[i]
    fitmean <- max(fitmean, bounds[1])
    fitmean <- min(fitmean, bounds[2])
    
    basis <- predict(isMat, fitmean)
    varE <- parameter$para.sdE^2 * (cbind(1, basis) %*% parameter$para.is)^2
    varE <- varE[1, 1]
    err_var_int_val[i] <- varE*dnorm(b[i], 0, parameter$para.sdR, log=FALSE)
  }
  
  return(err_var_int_val)
}

## get fitted mean curve for each diet group
## over-lay quantile band
## include random effects
ChickWeight$fitted <- fitted(chickfit.2$fit$model.bic$mixmod.1)
ChickWeight$fixed_fit <- getME(chickfit.2$fit$model.bic$mixmod.1, name = "X") %*% fixef(chickfit.2$fit$model.bic$mixmod.1)

(ggplot(ChickWeight, aes(x = Time, y = fitted, group = Chick))+geom_line(linetype = 'dashed') + 
    theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + facet_wrap( ~ Diet, nrow = 1, labeller = label_both) +
    ylab("Weight")) %>% ggsave(plot = .,
                               filename = "../manuscript/plots/ChickWeight_fit.pdf",
                               width = 9, height= 3)

temp <- ChickWeight[ChickWeight$Diet == 1, ]
p1_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p1_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_1, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)


temp <- ChickWeight[ChickWeight$Diet == 2, ]
p2_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p2_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_2, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

temp <- ChickWeight[ChickWeight$Diet == 3, ]
p3_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p3_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_3, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

temp <- ChickWeight[ChickWeight$Diet == 4, ]
p4_quantile <- ggplot(temp)+geom_line(aes(x = Time, y = weight, group = Chick)) + 
  geom_point(aes(x = Time, y = weight)) + theme_bw() + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)

## over-lay the fitted curve
p4_quantile_fitted <- ggplot(temp)+geom_line(linetype = 'dashed', aes(x = Time, y = fitted, group = Chick)) + 
  theme_bw() + geom_line(aes(x = Time, y = fixed_fit)) + ylab("Weight") + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower1, ymax = upper1), fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = quantile_plot_4, aes(x = time, ymin = lower2, ymax = upper2), fill = "grey", alpha = 0.8)


p_quantile <- grid.arrange(p1_quantile, p2_quantile, p3_quantile, p4_quantile, nrow = 1)
ggsave('../manuscript/plots/ChickWeight_quantile.pdf', p_quantile, width = 9, height = 3)
p_quantile_fitted <- grid.arrange(p1_quantile_fitted, p2_quantile_fitted, p3_quantile_fitted, p4_quantile_fitted, nrow = 1)
ggsave('../manuscript/plots/ChickWeight_fitted_quantile.pdf', p_quantile_fitted, width = 9, height = 3)
