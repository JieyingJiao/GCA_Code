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


(ggplot(ChickWeight, aes(x = Time, y = weight, group = Chick))+geom_line() + 
    geom_point()+theme_bw() + facet_wrap( ~ Diet, nrow = 1, labeller = label_both) +
    ylab("Weight")) %>% ggsave(plot = .,
                               filename = "../manuscript/plots/ChickWeight.pdf",
                               width = 9, height= 3)
chickdata.1 <- subset(ChickWeight, Diet == 1)
isMat <- iSpline(chickdata.1$Time, df = 3, degree = 0, intercept = TRUE)
chickdata.1 <- data.frame(chickdata.1, isMat)
chickfit1.0 <- lmer(weight ~ X1 + X2 + X3 + (Time - 1 |Chick),
                    data = chickdata.1, REML = FALSE)

p1 <- gg_qq(resid(chickfit1.0), labels = NA)
naivefit1 <- data.frame(x = getME(chickfit1.0, name = "X") %*% fixef(chickfit1.0), 
                        y = resid(chickfit1.0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit1.0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit1 <- lowess(naivefit1$x, naivefit1$y)
p2 <- ggplot(data = naivefit1, aes(x = x, y = y)) + xlab("Fitted mean") + 
  ylab("Residual") + geom_point() + 
  geom_line(aes(x = lowessfit1$x, y = lowessfit1$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15),  
        axis.title.x = element_blank())

naivefit1.2 <- data.frame(x = fitted(chickfit1.0), 
                        y = resid(chickfit1.0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit1.0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit1.2 <- lowess(naivefit1.2$x, naivefit1.2$y)
p2.2 <- ggplot(data = naivefit1.2, aes(x = x, y = y)) + xlab("Fitted mean") + 
  ylab("Residual") + geom_point() + 
  geom_line(aes(x = lowessfit1.2$x, y = lowessfit1.2$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15),  
        axis.title.x = element_blank())


## use the iteration method
chickfit1.1 <- lmerIter(weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
                        response = "weight", ID = "Chick", data = chickdata.1, 
                        df = 3:10, degree = 2, nsim = 2, tolerance = 1e-4, 
                        random.include = FALSE)
chickerr1 <- eSd(chickfit1.1, degree = 2, random.include = FALSE)
chickresid1 <- data.frame(x = chickerr1$BIC$fitmean, 
                          y = residIter(chickfit1.1)$BIC / chickerr1$BIC$esd, 
                          x1 = fitted(chickfit1.1$fit$model.bic$mixmod.1))


p3 <- gg_qq2(residIter(chickfit1.1)$BIC / chickerr1$BIC$esd)

fitmean1 <- chickerr1$BIC$fitmean

alpha <- 0.05
lower <- chickerr1$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr1 <- data.frame(x = fitmean1, y = resid(chickfit1.1$fit$model.bic$mixmod.1), 
                             lower = lower, upper = upper)

p4 <- ggplot(chickresiderr1, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1) + 
  xlab("Fitted mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 15), 
        axis.title = element_text(size = 15))
p <- grid.arrange(p1, p2, p3, p4, nrow = 2)
ggsave('../manuscript/plots/chickfit1.pdf', p, width = 9, height = 7)

## get log likelihood, marginal situation
print(chickfit1.1$fit$loglh.bic)
## -860.783


## fit conditional mean situation
chickfit1.2 <- lmerIter(weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
                        response = "weight", ID = "Chick", data = chickdata.1, 
                        df = 3:10, degree = 2, nsim = 2, tolerance = 1e-4, 
                        random.include = TRUE)
chickerr1.2 <- eSd(chickfit1.2, degree = 2, random.include = TRUE)
chickresid1.2 <- data.frame(x = chickerr1.2$BIC$fitmean, 
                          y = residIter(chickfit1.2)$BIC / chickerr1.2$BIC$esd, 
                          x1 = fitted(chickfit1.2$fit$model.bic$mixmod.1))


p3.2 <- gg_qq2(residIter(chickfit1.2)$BIC / chickerr1.2$BIC$esd)

fitmean1.2 <- chickerr1.2$BIC$fitmean

alpha <- 0.05
lower <- chickerr1.2$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr1.2 <- data.frame(x = fitmean1.2, y = resid(chickfit1.2$fit$model.bic$mixmod.1), 
                             lower = lower, upper = upper)

p4.2 <- ggplot(chickresiderr1.2, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1) + 
  xlab("Fitted mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 15), 
        axis.title = element_text(size = 15))
p <- grid.arrange(p1, p2.2, p3.2, p4.2, nrow = 2)
ggsave('../manuscript/plots/chickfit1_randominclude.pdf', p, width = 9, height = 7)

ll2 <- LL_func(chickfit1.2$fit$model.bic, formula = weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
               data = chickdata.1, degree = 2)
ll2
## -850.8608



##likelihood ratio test
# pvalue.1 <- pchisq(2 * (chickfit1.1$fit$loglh.bic - logLik(chickfit1.0)), 
#                    df = 3, lower.tail = FALSE)
# pvalue.1


## group 4
chickdata.4 <- subset(ChickWeight, Diet == 4)
isMat <- iSpline(chickdata.4$Time, df = 3, degree = 0, intercept = TRUE)
chickdata.4 <- data.frame(chickdata.4, isMat)
chickfit4.0 <- lmer(weight ~ X1 + X2 + X3 + (Time - 1 |Chick),
                    data = chickdata.4, REML = FALSE)
p5 <- gg_qq(resid(chickfit4.0), labels = NA)
naivefit4 <- data.frame(x = getME(chickfit4.0, name = "X") %*% fixef(chickfit4.0), 
                        y = resid(chickfit4.0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit4.0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit4 <- lowess(naivefit4$x, naivefit4$y)
p6 <- ggplot(data = naivefit4, aes(x = x, y = y)) + xlab("Fitted mean") + 
  ylab("Residual") + geom_point() + 
  geom_line(aes(x = lowessfit4$x, y = lowessfit4$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15),  
        axis.title.x = element_blank())

naivefit4.2 <- data.frame(x = fitted(chickfit4.0), 
                          y = resid(chickfit4.0))
cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
cutoff2 <- qnorm(0.05/2, lower.tail = FALSE) * data.frame(VarCorr(chickfit4.0))$sdcor[2]
cutoff3 <- -cutoff2
cutoff2 <- data.frame(yintercept=cutoff2, cutoff=factor("bound"))
cutoff3 <- data.frame(yintercept=cutoff3, cutoff=factor("bound"))
lowessfit4.2 <- lowess(naivefit4.2$x, naivefit4.2$y)
p6.2 <- ggplot(data = naivefit4.2, aes(x = x, y = y)) + xlab("Fitted mean") + 
  ylab("Residual") + geom_point() + 
  geom_line(aes(x = lowessfit4.2$x, y = lowessfit4.2$y), col = "red", show.legend = FALSE) +
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff2, show.legend = FALSE) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff3, show.legend = FALSE) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15),  
        axis.title.x = element_blank())


chickfit4.1 <- lmerIter(weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
                        response = "weight", ID = "Chick", data = chickdata.4, 
                        df = 3:10, degree = 2, nsim = 2, tolerance = 1e-4, 
                        random.include = FALSE)
chickerr4 <- eSd(chickfit4.1, degree = 2, random.include = FALSE)
chickresid4 <- data.frame(x = chickerr4$BIC$fitmean, 
                          y = residIter(chickfit4.1)$BIC / chickerr4$BIC$esd, 
                          x1 = fitted(chickfit4.1$fit$model.bic$mixmod.1))


p7 <- gg_qq2(residIter(chickfit4.1)$BIC / chickerr4$BIC$esd)

fitmean4 <- chickerr4$BIC$fitmean

alpha <- 0.05
lower <- chickerr4$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr4 <- data.frame(x = fitmean4, y = resid(chickfit4.1$fit$model.bic$mixmod.1), 
                             lower = lower, upper = upper)

p8 <- ggplot(chickresiderr4, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1) + 
  xlab("Fitted mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 15), 
        axis.title = element_text(size = 15))

p <- grid.arrange(p5, p6, p7, p8, nrow = 2)
ggsave('../manuscript/plots/chickfit4.pdf', p, width = 9, height = 7)

## get log likelihood, marginal situation
print(chickfit4.1$fit$loglh.bic)
## -429.6968

## fit conditional mean situation
chickfit4.2 <- lmerIter(weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
                        response = "weight", ID = "Chick", data = chickdata.4, 
                        df = 3:10, degree = 2, nsim = 2, tolerance = 1e-4, 
                        random.include = TRUE)
chickerr4.2 <- eSd(chickfit4.2, degree = 2, random.include = FALSE)
chickresid4.2 <- data.frame(x = chickerr4.2$BIC$fitmean, 
                          y = residIter(chickfit4.2)$BIC / chickerr4.2$BIC$esd, 
                          x1 = fitted(chickfit4.2$fit$model.bic$mixmod.1))


p7.2 <- gg_qq2(residIter(chickfit4.2)$BIC / chickerr4.2$BIC$esd)

fitmean4.2 <- chickerr4.2$BIC$fitmean

alpha <- 0.05
lower <- chickerr4.2$BIC$esd * qnorm(alpha/2, lower.tail = FALSE)
upper <- -lower

chickresiderr4.2 <- data.frame(x = fitmean4.2, y = resid(chickfit4.2$fit$model.bic$mixmod.1), 
                             lower = lower, upper = upper)

p8.2 <- ggplot(chickresiderr4.2, aes(x = x, y = y)) + 
  geom_point(aes(y = y), col = "black") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1) + 
  xlab("Fitted mean") + ylab("Residual") + 
  theme(legend.position="none",text = element_text(size = 15), 
        axis.title = element_text(size = 15))


p <- grid.arrange(p5, p6.2, p7.2, p8.2, nrow = 2)
ggsave('../manuscript/plots/chickfit4_randominclude.pdf', p, width = 9, height = 7)

ll2 <- LL_func(chickfit4.2$fit$model.bic, formula = weight ~ X1 + X2 + X3 + (Time - 1 | Chick), 
               data = chickdata.4, degree = 2)
ll2
## -428.9262




## likelihood ratio test
# pvalue.4 <- pchisq(2 * (chickfit4.1$fit$loglh.bic - logLik(chickfit4.0)), 
#                    df = 3, lower.tail = FALSE)
# pvalue.4


## overlaid growth curve plot for 4 diet groups
Time <- unique(ChickWeight$Time)
isMat <- iSpline(Time, df = 3, degree = 0, intercept = TRUE)
meancurve1 <- cbind(1, isMat) %*% fixef(chickfit1.1$fit$model.bic$mixmod.1)
meancurve2 <- cbind(1, isMat) %*% fixef(chickfit4.1$fit$model.bic$mixmod.1)
fitsum1 <- fitIter(chickfit1.1, degree = 2, random.include = FALSE)
isMaterr1 <- predict(fitsum1$isMat$BIC, meancurve1)
errsd1 <- cbind(1, isMaterr1) %*% fitsum1$para.is$BIC * 
  data.frame(VarCorr(chickfit1.1$fit$model.bic$mixmod.1))$sdcor[2]
sdcurve1 <- sqrt(Time^2 * data.frame(VarCorr(chickfit1.1$fit$model.bic$mixmod.1))$vcov[1] + 
                   errsd1^2)

fitsum4 <- fitIter(chickfit4.1, degree = 2, random.include = FALSE)
isMaterr