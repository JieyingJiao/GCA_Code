
pancrea  <- read.csv("../cleanData/pancrea.csv", header = TRUE)

## plot function in utilis
source('utils.R')

library(car) # qqPlot function
library(ggplot2)
library(gridExtra) # grid.arrange function
library(MASS) # boxcox function
library(tidyr) # gather function
library(tidyverse)
library(sas7bdat)
library(matrixStats)
library(splines2)
library(quadprog)

time.week <- as.numeric(substr(pancrea$Duration, 1, 2))
time.day <- as.numeric(substr(pancrea$Duration, 4, 4))
time <- 7 * time.week + time.day
pancrea <- data.frame(pancrea, time)
pancrea.long <- gather(pancrea, key = "Type", value = "measurement", Length, 
                       Top, Body, Tail, factor_key = TRUE)
ggplot(pancrea.long, aes(x = time, y = measurement)) + geom_point() + 
  theme_bw() + facet_wrap(~Type, nrow = 1, scale = "free_y")





##' MLE and ispline to fit mean and variance structure of y, normal dist
##' depend on x1 and x2 respectively
##' on variance scale
##' output ispline basis coefficient
##' degree1 and degree2 for mean and variance, respectively
csplineMLE.1 <- function(y, x1, x2, df1, df2, degree1, degree2) {
  N <- length(y)
  csMat1 <- cSpline(x1, df = df1, degree = degree1, intercept = TRUE)
  isMat2 <- iSpline(x2, df = df2, degree = degree2, intercept = TRUE)
  
  basis1 <- cbind(1, x1, csMat1)
  basis2 <- cbind(1, isMat2)
  
  basis <- list(basis1, basis2)
  neloglh <- function(beta, y, basis, df1, df2) {
    mean <- basis[[1]] %*% beta[1:(df1+2)]
    var <- (basis[[2]] %*% beta[(df1+3):(df1+df2+3)])^2
    sum(log(2*pi*var)/2 + (y-mean)^2/2/var)
  }
  startbeta <- c(mean(y), 0.001, rep(-0.001, df1), sd(y), rep(0.001, df2))
  
  ## constraints: mean increasing and concave, error standard deviation increasing
  ui1 <- rbind(c(0, 1, deriv(predict(csMat1, max(x1)))), cbind(matrix(0, nrow = df1, ncol = 2), diag(-1, df1))) 
  ui2 <- diag(1, df2+1)
  ui <- rbind(cbind(ui1, matrix(0, nrow = nrow(ui1), ncol = ncol(ui2))), cbind(matrix(0, nrow = nrow(ui2), ncol = ncol(ui1)), ui2))
  fit <- constrOptim(theta = startbeta, f = neloglh, y = y, basis = basis,
                     df1 = df1, df2 = df2,grad = NULL, ui = ui,
                     ci = c(rep(0, df1+1), 1, rep(0, df2)),
                     control = list(maxit = 5000))
  return(list(mean.par = fit$par[1:(df1+2)], 
              sd.par = fit$par[(df1+3):(df1+df2+3)],
              loglh = - fit$value,
              csMat.mean = csMat1,
              csMat.sd = isMat2, 
              mean_x = x1,
              sd_x = x2))
}

##' calculate fitted mean and variance at given x value:x.new
##' input object from isplineMLE.1
fitsum <- function(x_mean.new, x_sd.new, object) {
  mean <- cbind(1, x_mean.new, predict(object$csMat.mean, x_mean.new)) %*% object$mean.par
  sd <- cbind(1, predict(object$csMat.sd, x_sd.new)) %*% object$sd.par
  return(list(mean.fit = mean, sd.fit = sd))
}

##' compute AIC and BIC
aic <- function(loglh, k) 2*k-2*loglh
bic <- function(loglh, k, n) k*log(n) - 2*loglh

##' df selection using BIC criteria. df1 and df2 are vectors
fitSelect.2 <- function(y, x1, x2, df1, df2, degree1, degree2) {
  n <- length(y)
  bicMat <- matrix(0, nrow = length(df1), ncol = length(df2))
  
  for (i in 1:length(df1)) {
    for (j in 1:length(df2)) {
      fit.mle <- csplineMLE.1(y = y, x1 = x1, x2 = x2, df1 = df1[i], df2 = df2[j], 
                              degree1 = degree1, degree2 = degree2)
      bicMat[i, j] <- bic(fit.mle$loglh, k = df1[i] + df2[j] + 2, n = n)
    }
  }
  rownames(bicMat) <- df1
  colnames(bicMat) <- df2
  return(bicMat)
}

pancrea.length <- na.omit(pancrea[, c("time", "Length")])
pancrea.length <- data.frame(pancrea.length, type = "length")
## use standard time?
pancrea.length$time <- pancrea.length$time - min(pancrea.length$time)

degree1 <- 0
degree2 <- 1
bicMat <- fitSelect.2(pancrea.length$Length, pancrea.length$time, 
                      pancrea.length$time, df1 = 4:5, df2 = 4:5, 
                      degree1 = degree1, degree2 = degree2)

## '+3' since the df selection starting from 4
df.bic <- which(bicMat == min(bicMat), arr.ind = TRUE) + 3
df.bic
# both selected df as 4
# minimal bic is 308.9445

lengthfit <- csplineMLE.1(pancrea.length$Length, pancrea.length$time, 
                          pancrea.length$time, df1 = df.bic[1], df2 = df.bic[2], degree1 = degree1, 
                          degree2 = degree2)
length.naive <- lm(Length ~ time, data = pancrea.length)

## p-value for likelihood ratio test
pchisq(2*(lengthfit$loglh - logLik(length.naive)), df = 3, lower.tail = FALSE)

x.length <- seq(from = min(pancrea.length$time), to = max(pancrea.length$time), 
                by = 0.01)
length.fit <- fitsum(x.length, x.length, lengthfit)
## loglikelihood: -135.5513

## pancrea.top <- na.omit(pancrea[, c("time", "Top")])
## pancrea.top <- data.frame(pancrea.top, type = "top")
## topfit <- isplineMLE.1(pancrea.top$Top, pancrea.top$time, pancrea.top$time,
##                        df1 = 3, df2 = 3, degree1 = 2, degree2 = 2)
## x.top <- seq(from = min(pancrea.top$time), to = max(pancrea.top$time), by = 0.01)
## top.fit <- fitsum(x.top, topfit)

## pancrea.body <- na.omit(pancrea[, c("time", "Body")])
## pancrea.body <- data.frame(pancrea.body, type = "body")
## bodyfit <- isplineMLE.1(pancrea.body$Body, pancrea.body$time, pancrea.body$time,
##                        df1 = 3, df2 = 3, degree1 = 2, degree2 = 2)
## x.body <- seq(from = min(pancrea.body$time), to = max(pancrea.body$time), by = 0.01)
## body.fit <- fitsum(x.body, bodyfit)

## pancrea.tail <- na.omit(pancrea[, c("time", "Tail")])
## pancrea.tail <- data.frame(pancrea.tail, type = "tail")
## tailfit <- isplineMLE.1(pancrea.tail$Tail, pancrea.tail$time, pancrea.tail$time,
##                        df1 = 3, df2 = 3, degree1 = 2, degree2 = 2)
## x.tail <- seq(from = min(pancrea.tail$time), to = max(pancrea.tail$time), by = 0.01)
## tail.fit <- fitsum(x.tail, tailfit)

q1 <- qnorm(0.05/2, lower.tail = FALSE)
q2 <- qnorm(0.1/2, lower.tail = FALSE)
length.plot <- data.frame(time = x.length, mean = length.fit$mean.fit, 
                          lower1 = length.fit$mean.fit - q1 * length.fit$sd.fit,
                          upper1 = length.fit$mean.fit + q1 * length.fit$sd.fit,
                          lower2 = length.fit$mean.fit - q2 * length.fit$sd.fit,
                          upper2 = length.fit$mean.fit + q2 * length.fit$sd.fit,
                          type = "length")

length.resid <- fitsum(pancrea.length$time, pancrea.length$time, lengthfit)
resid.plot <- data.frame(fitted = length.resid$mean.fit, 
                         errsd = length.resid$sd.fit, 
                         resid = pancrea.length$Length - length.resid$mean.fit, 
                         lower = - q1 * length.resid$sd.fit, 
                         upper = q1 * length.resid$sd.fit)

p1 <- ggplot(data = length.plot, aes(x = time, y = mean)) + geom_line() + 
  ylab("measure") + xlab("Duration") + 
  geom_ribbon(data = length.plot, aes(ymin = lower1, ymax = upper1),
              fill = "grey", alpha = 0.5) + 
  geom_ribbon(data = length.plot, aes(ymin = lower2, ymax = upper2),
              fill = "grey", alpha = 0.8) +
  geom_point(data = pancrea.length, aes(y = Length)) + 
  theme(text = element_text(size = 15)) + 
  ggtitle("(I)") + ylab("Length (mm)") + xlab('Prenatal period (days)') + theme_bw()
p1


p2 <- gg_qq2(resid.plot$resid / resid.plot$errsd,
             labels = "Standardized Residual Sample Quantile") + 
  theme_bw() + ggtitle("(III)")
p2

cutoff1 <- data.frame(yintercept=0, cutoff=factor(0))
p3 <- ggplot(data = resid.plot, aes(x = fitted, y = resid)) + 
  geom_point() + theme(text = element_text(size = 15), axis.title = element_text(size = 15)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) + 
  geom_point() + 
  geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff1, show.legend = FALSE) + 
  xlab("Fitted Value") + ylab("Residual") + 
  theme_bw() + ggtitle("(II)")
p3

(grid.arrange(p1, p3, p2, nrow = 1)) %>% 
  ggsave(plot = .,filename = "../manuscript/plots/pancrealength_cspline.pdf",
         width = 8, height = 3)


## top.fit <- data.frame(time = x.top, mean = top.fit$mean.fit,
##                       lower1 = top.fit$mean.fit - q1 * sqrt(top.fit$var.fit),
##                       upper1 = top.fit$mean.fit + q1 * sqrt(top.fit$var.fit),
##                       lower2 = top.fit$mean.fit - q2 * sqrt(top.fit$var.fit),
##                       upper2 = top.fit$mean.fit + q2 * sqrt(top.fit$var.fit),
##                       type = "top")
## body.fit <- data.frame(time = x.body, mean = body.fit$mean.fit,
##                        lower1 = body.fit$mean.fit - q1 * sqrt(body.fit$var.fit),
##                        upper1 = body.fit$mean.fit + q1 * sqrt(body.fit$var.fit),
##                        lower2 = body.fit$mean.fit - q2 * sqrt(body.fit$var.fit),
##                        upper2 = body.fit$mean.fit + q2 * sqrt(body.fit$var.fit),
##                        type = "body")
## tail.fit <- data.frame(time = x.tail, mean = tail.fit$mean.fit,
##                        lower1 = tail.fit$mean.fit - q1 * sqrt(tail.fit$var.fit),
##                        upper1 = tail.fit$mean.fit + q1 * sqrt(tail.fit$var.fit),
##                        lower2 = tail.fit$mean.fit - q2 * sqrt(tail.fit$var.fit),
##                        upper2 = tail.fit$mean.fit + q2 * sqrt(tail.fit$var.fit),
##                        type = "tail")

## pancrea.fit <- rbind(length.fit, top.fit, body.fit, tail.fit)
## ggplot(data = pancrea.fit, aes(x = time, y = mean)) + 
##   facet_wrap(~type, nrow = 1, scales = "free_y") + 
##   geom_line() + ylab("measure") + xlab("Duration") + 
##   geom_ribbon(data = length.fit, aes(ymin = lower1, ymax = upper1),
##               fill = "grey", alpha = 0.5) + 
##   geom_ribbon(data = length.fit, aes(ymin = lower2, ymax = upper2),
##               fill = "grey", alpha = 0.8) +
##   geom_ribbon(data = top.fit, aes(ymin = lower1, ymax = upper1),
##               fill = "grey", alpha = 0.5) + 
##   geom_ribbon(data = top.fit, aes(ymin = lower2, ymax = upper2),
##               fill = "grey", alpha = 0.8) + 
##   geom_ribbon(data = body.fit, aes(ymin = lower1, ymax = upper1),
##               fill = "grey", alpha = 0.5) + 
##   geom_ribbon(data = body.fit, aes(ymin = lower2, ymax = upper2),
##               fill = "grey", alpha = 0.8) +
##   geom_ribbon(data = tail.fit, aes(ymin = lower1, ymax = upper1),
##               fill = "grey", alpha = 0.5) + 
##   geom_ribbon(data = tail.fit, aes(ymin = lower2, ymax = upper2),
##               fill = "grey", alpha = 0.8) + 
##   geom_point(data = pancrea.length, aes(y = Length)) + 
##   geom_point(data = pancrea.top, aes(y = Top)) + 
##   geom_point(data = pancrea.body, aes(y = Body)) + 
##   geom_point(data = pancrea.tail, aes(y = Tail))