library(tsutils)
library(forecast)
library(Metrics)
library(TTR)
library(RColorBrewer)
library(lubridate)
library(tidyverse)
library(ggplot2)
library(smooth)
source("MSci381_exploration_functions.r")

data <- read.csv("Dataset3 (4).csv")
date_daily = seq(from = as.Date("2019-11-24") - 99, to = as.Date("2019-11-24"), by = "day")

startW <- as.numeric(strftime(head(date_daily, 1), format = "%U"))
startD <- as.numeric(strftime(head(date_daily, 1), format = "%w")) 

endW <- as.numeric(strftime(tail(date_daily, 1), format = "%U"))
endD <- as.numeric(strftime(tail(date_daily, 1), format = "%w")) 

data_p1 <- ts(data$ProductP1, frequency = 52, end = decimal_date(ymd("2019-11-10")))
data_p2 <- ts(data$ProductP2, frequency = 52, end = decimal_date(ymd("2019-11-10")))
data_p3 <- ts(data$ProductP3, frequency = 7, start = c(startW, startD))
data_p4 <- ts(data$ProductP4, frequency = 7, start = c(startW, startD))
data_p5 <- ts(data$ProductP5, frequency = 52, end = decimal_date(ymd("2019-08-30")))

wday(date_daily[1:7], label = TRUE)

summary(data_p1)
summary(data_p2)
summary(data_p3)
summary(data_p4)
summary(data_p5)

pacf(data_p3)
pacf(data_p4)

plot(data$price, type = "l")

hist(data_p1)
hist(data_p2)
hist(data_p3)
hist(data_p4)
hist(data_p5)

plot(data_p1, main = "P1 Demand")
plot(data_p2, main = "P2 Demand")
plot(data_p3, main = "P3 Demand", xlab = "Week")
plot(data_p4, main = "P4 Demand")
plot(data_p5, main = "P5 Demand")


cma_p1 <- cmav(data_p1, outplot = 1)
cma_p2 <- cmav(data_p2, outplot = 1)
cma_p3 <- cmav(data_p3, outplot = 1)
cma_p4 <- cmav(data_p4, outplot = 1)
cma_p5 <- cmav(data_p5, outplot = 1)

#seasonality

seasplot(data_p1, outplot = 4)
seasplot(data_p2, outplot = 1)
seasplot(data_p3, outplot = 2)
seasplot(data_p4, outplot = 1)
seasplot(data_p5, outplot = 1)


#Decomposition

decomp_p1 <- decomp(data_p1, decomposition = "additive", outplot = 1) #returns error
decomp_p2 <- decomp(data_p2, decomposition = "additive", outplot = 1) #returns error
decomp_p3 <- decomp(data_p3, decomposition = "additive", outplot = 1)
decomp_p4 <- decomp(data_p4, decomposition = "additive", outplot = 1)
decomp_p5 <- decomp(data_p5, decomposition = "additive", outplot = 1)

#forecast P1

p1.trn <- head(data_p1, 80)
p1.tst <- tail(data_p1, 20)


p1_fit1 <- ets(p1.trn, model = "ANN")
plot(p1.trn)
lines(p1_fit1$fitted, col="red")
p1_fit1

p1_fit2 <- ets(p1.trn, model = "ANN", alpha = 0.1)
plot(p1.trn)
lines(p1_fit2$fitted, col="red")
p1_fit2

p1_fit1$mse
p1_fit2$mse

p1_frc1 <- forecast(p1_fit1, h=10)
p1_frc2 <- forecast(p1_fit2, h=10)
plot(p1_frc1)
lines(p1_fit1$fitted,col="blue")
lines(p1_frc2$mean, col="red")
lines(p1_fit2$fitted,col="red")
lines(p1.tst, lty = 2)
lines(p1_frc_lag, col = "magenta")
lines(p1_frc_lag_1, col = "magenta")
legend("bottomleft", c("alpha = 0.3646", "alpha = 0.1", "Autoregressive"), col = c("blue", "red", "magenta"), lty = 1)


crit_p1 <- array(NA,c(2,3),dimnames=list(c("Forecast 1", "Forecast 2"),
                                      c("AIC","AICc","BIC")))

crit_p1[1,1] <- p1_fit1$aic
crit_p1[1,2] <- p1_fit1$aicc
crit_p1[1,3] <- p1_fit1$bic
crit_p1[2,1] <- p1_fit2$aic
crit_p1[2,2] <- p1_fit2$aicc
crit_p1[2,3] <- p1_fit2$bic

crit_p1

MAE_p1.1 <- mean(abs(p1.tst - p1_frc1$mean))
MAE_p1.2 <- mean(abs(p1.tst - p1_frc2$mean))
MAE_p1 <- c(MAE_p1.1, MAE_p1.2)
names(MAE_p1) <- paste0("Forecast ",1:2) 
round(MAE_p1,3)

MAE_p1_1step <- c(mean(abs(p1.tst[1] - p1_frc1$mean[1])), mean(abs(p1.tst[1] - p1_frc2$mean[1])))
MAE_p1_1step

MAE_p1_insample <- c(mean(abs(p1.trn - p1_frc1$fitted)), mean(abs(p1.trn - p1_frc2$fitted)))
MAE_p1_insample

MSE_p1.1 <- mean((p1.tst - p1_frc1$mean)^2)
MSE_p1.2 <- mean((p1.tst - p1_frc2$mean)^2)

MSE_p1_1step <- c(mean((p1.tst[1] - p1_frc1$mean[1])^2), mean((p1.tst[1] - p1_frc2$mean[1])^2))
sqrt(MSE_p1_1step)

#SMA for p1

p1_sma <- sma(p1.trn, h = 10, silent = FALSE)
lines(p1.tst, lty = 2)
summary(p1_sma)
names(p1_sma)

MAE_p1.3 <- mean(abs(p1.tst - p1_sma$forecast))
MAE_p1.3
MAE_p1.3.insample <- mean(abs(p1.trn - p1_sma$fitted))
MAE_p1.3.insample
MAE_p1.3.1step <- mean(abs(p1.tst[1] - p1_sma$forecast[1]))
MAE_p1.3.1step

RMSE_p1.3 <- sqrt(mean((p1.tst - p1_sma$forecast)^2))
RMSE_p1.3.insample <- sqrt(mean((p1.trn - p1_sma$fitted)^2))
RMSE_p1.3.1step <- sqrt(mean((p1.tst[1] - p1_sma$forecast[1])^2))

RMSE_p1.3.1step
RMSE_p1.3
RMSE_p1.3.insample


#lag model for P1

pacf(data_p1)

P1 <- array(NA, c(length(p1.trn), 2))

for (i in 1:2) {
  P1[i:length(p1.trn), i] <- p1.trn[1:(length(p1.trn) - i + 1)]
}
colnames(P1) <- c("p1", "lag 1")
P1 <- as.data.frame(P1)
P1
plot(P1)

p1_fit_lag <- step(lm(p1~., P1))
summary(p1_fit_lag)
par(mfrow=c(2,2)) # Split into 2 by 2
plot(p1_fit_lag)
par(mfrow=c(1,1))


plot(P1$p1, type = "l")
p1_frc_lag <- predict(p1_fit_lag, P1)
lines(p1_frc_lag, col = "red")

prod1new <- array(tail(p1.trn, 1), c(1,1))
colnames(prod1new) <- paste0("lag ", 1:1)
prod1new <- as.data.frame(prod1new)
prod1new

predict(p1_fit_lag, prod1new, h = 10)

p1_frc_lag_1 <- array(NA, c(10,1))

for (i in 1:10) {
  prod1new <- tail(p1.trn, 1)
  prod1new <- c(prod1new, p1_frc_lag_1)
  prod1new <- prod1new[i:(0+i)]
  prod1new <- prod1new[1:1]
  prod1new <- array(prod1new, c(1,1))
  colnames(prod1new) <- paste0("lag ", 1:1)
  prod1new <- as.data.frame(prod1new)
  p1_frc_lag_1[i] <- predict(p1_fit_lag, prod1new)
}

p1_frc_lag_1

p1_frc_lag_1 <- ts(p1_frc_lag_1, frequency = frequency(p1.tst), start = start(p1.tst))
p1_frc_lag <- ts(p1_frc_lag, frequency = frequency(p1.trn), start = start(p1.trn))
ts.plot(p1.trn, p1.tst, p1_frc_lag,p1_frc_lag_1, col = c("black", "black", "red", "red"))

MSE_p1lag <- mean((p1.tst-p1_frc_lag_1)^2)
MSE_p1 <- c(MSE_p1.1, MSE_p1.2, MSE_p1lag) 
RMSE_p1 <- sqrt(MSE_p1) 
names(RMSE_p1) <- paste0("Forecast ",1:3) 
round(RMSE_p1,3) 

MSE_p1lag_1step <- mean((p1.tst[1] - p1_frc_lag_1[1])^2)
sqrt(MSE_p1lag_1step)

MSE_p1_insample <- c(mean((p1.trn - p1_fit1$fitted)^2), mean((p1.trn - p1_fit2$fitted)^2), mean((p1.trn - p1_fit_lag$fitted.values)^2))
sqrt(MSE_p1_insample)
MSE_p1_insample

MAE_p1lag <- mean(abs(p1.tst - p1_frc_lag_1))
MAE_p1lag
MAE_p1lag1step <- mean(abs(p1.tst[1] - p1_frc_lag_1[1]))
MAE_p1lag1step
MAE_p1laginsample <- mean(abs(tail(p1.trn, 79) - p1_fit_lag$fitted.values))
MAE_p1laginsample


p1_frc1 <- forecast(p1_fit1, h=10)
p1_frc2 <- forecast(p1_fit2, h=10)
plot(p1_frc1)
lines(p1_fit1$fitted,col="blue")
lines(p1_frc2$mean, col="red")
lines(p1_fit2$fitted,col="red")
lines(p1_frc2$lower[,2],col="red") 
lines(p1_frc2$upper[,2],col="red")
lines(p1.tst, lty = 2)
lines(p1_frc_lag, col = "magenta")
lines(p1_frc_lag_1, col = "magenta")
legend("bottomleft", c("alpha = 0.3646", "alpha = 0.1", "Autoregressive"), col = c("blue", "red", "magenta"), lty = 1)



#forecast P2

p2.trn <- head(data_p2, 80)
p2.tst <- tail(data_p2, 20)
pacf(p2.trn)


p2_fit1 <- ets(p2.trn, model = "MAN")
p2_fit1
plot(p2.trn)
lines(p2_fit1$fitted, col = "red")

p2_fit2 <- ets(p2.trn, model = "MAN", alpha = 0.005, beta = 0.005)
p2_fit2
plot(p2.trn)
lines(p2_fit2$fitted, col = "red")


sqrt(p2_fit1$mse)
sqrt(p2_fit2$mse)

p2_frc1 <- forecast(p2_fit1, h = 10)
p2_frc2 <- forecast(p2_fit2, h = 10)
plot(p2_frc1)
lines(p2_fit1$fitted,col="blue")
lines(p2_frc2$mean, col = "red")
lines(p2_fit2$fitted,col="red")
lines(p2_frc2$lower[,2],col="red") 
lines(p2_frc2$upper[,2],col="red")
lines(p2.tst, lty = 2)
legend("bottomright", c("Optimal Values EXSM", "More reactive EXSM"), col = c("blue", "red"), lty = 1)

MAE_p2.1 <- mean(abs(p2.tst - p2_frc1$mean))
MAE_p2.2 <- mean(abs(p2.tst - p2_frc2$mean))
MAE_p2 <- c(MAE_p2.1, MAE_p2.2)
names(MAE_p2) <- paste0("Forecast ",1:2) 
round(MAE_p2,3)

MSE_p2.1 <- mean((head(p2.tst) - p2_frc1$mean)^2)
MSE_p2.2 <- mean((head(p2.tst) - p2_frc2$mean)^2)
MSE_p2 <- c(MSE_p2.1, MSE_p2.2) 
RMSE_p2 <- sqrt(MSE_p2) 
names(RMSE_p2) <- paste0("Forecast ",1:2) 
round(RMSE_p2,3)

MSE_p2_1step <- c(mean((p2.tst[1] - p2_frc1$mean[1])^2), mean((p2.tst[1] - p2_frc2$mean[1])^2))
sqrt(MSE_p2_1step)

#naive trend

p2_naive <- rwf(p2.trn, h = 10, drift = TRUE)
p2_naive$mean

MSE_p2naive <- c(mean((tail(p2.trn, 79) - p2_naive$fitted)^2), mean((p2.tst - p2_naive$mean)^2), mean((p2.tst[1] - p2_naive$mean[1])^2))
RMSE_p2naive <- sqrt(MSE_p2naive)
RMSE_p2naive

plot(p2_frc1)
lines(p2_fit1$fitted,col="blue")
lines(p2_frc2$mean, col = "red")
lines(p2_fit2$fitted,col="red")
lines(p2_frc2$lower[,2],col="red") 
lines(p2_frc2$upper[,2],col="red")
lines(p2_naive$mean, col = "green")
lines(p2.tst, lty = 2)
legend("bottomright", c("Optimal Values EXSM", "More reactive EXSM", "Naive Trend"), col = c("blue", "red", "green"), lty = 1)


#lag model for P2

P2 <- array(NA, c(length(p2.trn), 4))

for (i in 1:4) {
  P2[i:length(p2.trn), i] <- p2.trn[1:(length(p2.trn) - i + 1)]
}
colnames(P2) <- c("p2", paste0("lag", 1:3))
P2 <- as.data.frame(P2)
P2
plot(P2)


p2_fit_lag <- step(lm(p2~., P2))
summary(p2_fit_lag)
par(mfrow=c(2,2)) # Split into 2 by 2
plot(p2_fit_lag)
par(mfrow=c(1,1))

plot(P2$p2, type = "l")
p2_frc_lag <- predict(p2_fit_lag, P2)
lines(p2_frc_lag, col = "red")

prod2new <- array(tail(p2.trn, 1), c(1,3))
colnames(prod2new) <- paste0("lag", 1:3)
prod2new <- as.data.frame(prod2new)
prod2new

predict(p2_fit_lag, prod2new, h = 10)

p2_frc_lag_1 <- array(NA, c(10,1))

for (i in 1:10) {
  prod2new <- tail(p2.trn, 3)
  prod2new <- c(prod2new, p2_frc_lag_1)
  prod2new <- prod2new[i:(2+i)]
  prod2new <- prod2new[3:1]
  prod2new <- array(prod2new, c(1,3))
  colnames(prod2new) <- paste0("lag", 1:3)
  prod2new <- as.data.frame(prod2new)
  p2_frc_lag_1[i] <- predict(p2_fit_lag, prod2new)
}

p2_frc_lag_1

p2_frc_lag_1 <- ts(p2_frc_lag_1, frequency = frequency(p2.tst), start = start(p2.tst))
p2_frc_lag <- ts(p2_frc_lag, frequency = frequency(p2.trn), start = start(p2.trn))
ts.plot(p2.trn, p2.tst, p2_frc_lag,p2_frc_lag_1, col = c("black", "black", "red", "red"))

MSE_p2lag <- mean((p2.tst-p2_frc_lag_1)^2)
RMSE_p2lag <- sqrt(MSE_p2lag)
RMSE_p2lag

sqrt(mean((p2.trn[4:80] - p2_fit_lag$fitted.values)^2))

MSE_p2_final <- c(MSE_p2.1, MSE_p2.2, MSE_p2lag) 
RMSE_p2_final <- sqrt(MSE_p2_final) 
names(RMSE_p2_final) <- paste0("Forecast ",1:3) 
round(RMSE_p2_final,3) 

MSE_p2lag_1step <- mean((p2.tst[1] - p2_frc_lag_1[1])^2)
sqrt(MSE_p2lag_1step)

#lagged model with handled trend P2

P2_diff <- P2
for (i in 1:ncol(P2_diff)){
  P2_diff[,i] <- c(NA, diff(P2_diff[,i]))
}
P2_diff


p2_fit_lag_diff <- step(lm(p2~., P2_diff))
summary(p2_fit_lag_diff)

par(mfrow=c(2,2)) # Split into 2 by 2
plot(p2_fit_lag_diff)
par(mfrow=c(1,1))


p2_frc_lag_diff <- array(NA,c(10,1))

for (i in 1:10){
  p2.diff <- diff(p2.trn)
  P2_new <- tail(p2.diff, 3)
  P2_new <- c(P2_new, p2_frc_lag_diff)
  P2_new <- P2_new[i:(2+i)]
  P2_new <- P2_new[3:1]
  P2_new <- array(P2_new, c(1,3))
  colnames(P2_new) <- paste0("lag", 1:3)
  P2_new <- as.data.frame(P2_new)
  p2_frc_lag_diff[i] <- predict(p2_fit_lag_diff, P2_new)
}


p2_frc_lag_diff

p2_frc_lag_ud <- cumsum(c(tail(p2.trn, 1), p2_frc_lag_diff))
p2_frc_lag_ud <- p2_frc_lag_ud[-1]

p2_frc_lag_ud <- ts(p2_frc_lag_ud, frequency = frequency(p2.tst), start = start(p2.tst))
ts.plot(p2.trn, p2.tst, p2_frc_lag_ud, col = c("black","black","red"))

MSE_p2lag_1 <- mean((p2.tst-p2_frc_lag_ud)^2)
MSE_p2_final <- c(MSE_p2.1, MSE_p2.2, MSE_p2lag_1) 
RMSE_p2_final <- sqrt(MSE_p2_final) 
names(RMSE_p2_final) <- paste0("Forecast ",1:3) 
round(RMSE_p2_final,3) 

MSE_p2lag_1step <- mean((p2.tst[1] - p2_frc_lag_ud[1])^2)
sqrt(MSE_p2lag_1step)

MSE_p2_insample <- c(mean((p2.trn - p2_fit1$fitted)^2), mean((p2.trn - p2_fit2$fitted)^2), mean((p2.trn - p2_fit_lag_diff$fitted.values)^2))
sqrt(MSE_p2_insample)
MSE_p2_insample

p2_fit_lag_diff_ts <- cumsum(c(p2.trn[4], p2_fit_lag_diff$fitted.values))
p2_fit_lag_diff_ts <- p2_fit_lag_diff_ts[-1]
p2_fit_lag_diff_ts <- ts(p2_fit_lag_diff_ts, frequency = frequency(p2.trn), start = start(tail(p2.trn, 75)))


p2_frc1 <- forecast(p2_fit1, h = 10)
p2_frc2 <- forecast(p2_fit2, h = 10)
plot(p2_frc1)
lines(p2_fit1$fitted,col="blue")
lines(p2_frc2$mean, col = "red")
lines(p2_fit2$fitted,col="red")
lines(p2_frc2$lower[,2],col="red") 
lines(p2_frc2$upper[,2],col="red")
lines(p2.tst, lty = 2)
lines(p2_fit_lag_diff_ts, col = "magenta")
lines(p2_frc_lag_ud, col = "magenta")
legend("bottomright", c("Optimal Values EXSM", "More reactive EXSM", "Differences"), col = c("blue", "red", "magenta"), lty = 1)


p2_frc_lag_x <- ts(p2_frc_lag, frequency = frequency(p2.trn), start = start(p2.trn))

p2_frc1 <- forecast(p2_fit1, h = 10)
p2_frc2 <- forecast(p2_fit2, h = 10)
plot(p2_frc1)
lines(p2_fit1$fitted,col="blue")
lines(p2_frc2$mean, col = "red")
lines(p2_fit2$fitted,col="red")
lines(p2_frc2$lower[,2],col="red") 
lines(p2_frc2$upper[,2],col="red")
lines(p2.tst, lty = 2)
lines(p2_frc_lag_x, col = "magenta")
lines(p2_frc_lag_1, col = "magenta")
legend("bottomright", c("Optimal Values EXSM", "More reactive EXSM", "Autoregressive"), col = c("blue", "red", "magenta"), lty = 1)



#forecast for P3

p3.trn <- head(data_p3, 77)
p3.tst <- tail(data_p3, 23)



p3_fit1 <- ets(p3.trn, model = "AAA")
p3_fit1
plot(p3.trn)
lines(p3_fit1$fitted, col = "red")

p3_fit2 <- ets(p3.trn, model = "AAA", alpha = 0.01, beta = 0.005, gamma = 0.001)
p3_fit2
plot(p3.trn)
lines(p3_fit2$fitted, col = "red")
lines(p3_fit1$fitted, col = "blue")
legend("bottomright", c("Optimal values EXSM", "More reactive EXSM"), col = c("blue" , "red"), lty = 1)

p3_frc1 <- forecast(p3_fit1, h = 10)
p3_frc2 <- forecast(p3_fit2, h = 10)

plot(p3_frc2)
lines(p3_fit2$fitted, col="blue")
lines(p3_frc3, col = "red")
lines(p3_frc4, col = "red")
lines(p3.tst, lty = 2)
legend("bottomright", c("More reactive EXSM", "Autoregressive"), col = c("blue", "red"), lty = 1)

MAE_p3.1 <- mean(abs(p3.tst - p3_frc1$mean))
MAE_p3.2 <- mean(abs(p3.tst - p3_frc2$mean))
MAE_p3 <- c(MAE_p3.1, MAE_p3.2)
names(MAE_p3) <- paste0("Forecast ",1:2) 
round(MAE_p3,3)

MSE_p3.1 <- mean((p3.tst - p3_frc1$mean)^2)
MSE_p3.2 <- mean((p3.tst - p3_frc2$mean)^2)
MSE_p3 <- c(MSE_p3.1, MSE_p3.2) 
RMSE_p3 <- sqrt(MSE_p3) 
names(RMSE_p3) <- paste0("Forecast ",1:2) 
round(RMSE_p3,3)

MSE_p3_1step <- c(mean((p3.tst[1] - p3_frc1$mean[1])^2), mean((p3.tst[1] - p3_frc2$mean[1])^2))
sqrt(MSE_p3_1step)

#linear regression for P3

x <- data[, c(4,5)] #dataframe for p3 and p4
plot(x)
cor(data_p3, data_p4)

x1 <- data[, c(4, 8)]
plot(x1)
cor(x1)

yy <- cbind(diff(data$ProductP3), diff(data$ProductP4))
cor(yy)
plot(yy)

yy1 <- cbind(diff(data$ProductP3), diff(data$temperature))
plot(yy1)
cor(yy1)

ccf(data$ProductP3, data$ProductP4)
ccf(diff(data$ProductP3), diff(data$ProductP4))



#lag model for P3

P3 <- array(NA, c(length(p3.trn), 9))

for (i in 1:9) {
  P3[i:length(p3.trn), i] <- p3.trn[1:(length(p3.trn) - i + 1)]
}
colnames(P3) <- c("P3", paste0("lag", 1:8))
P3 <- as.data.frame(P3)
plot(P3)



pacf(p3.trn)

n <- length(p3.trn)
n

prod3 <- array(NA, c(n,9))
for(i in 1:9){
  prod3[i:n,i] <- p3.trn[1:(n-i+1)]
  
}
colnames(prod3) <- c("p3",paste0("lag",1:8))
prod3[1:10,]
prod3[(n-9):n,]

prod3 <- as.data.frame(prod3)
plot(prod3)

p3_fit4 <- lm(p3~.,data = prod3)
summary(p3_fit4)

p3_fit5 <- step(p3_fit4)
summary(p3_fit5)
par(mfrow=c(2,2)) # Split into 2 by 2
plot(p3_fit5)
par(mfrow=c(1,1))

c(AIC(p3_fit1),AIC(p3_fit2),AIC(p3_fit3),AIC(p3_fit4),AIC(p3_fit5))

plot(prod3$p3, type = "l")
p3_frc3 <- predict(p3_fit5, prod3)
lines(p3_frc3, col = "red")

prod3new <- array(tail(p3.trn, 8), c(1,8))
colnames(prod3new) <- paste0("lag", 8:1)
prod3new <- as.data.frame(prod3new)
prod3new

p3_frc_lag <- predict(p3_fit5, prod3)

p3_frc4 <- array(NA,c(10,1))

for (i in 1:10) {
  prod3new <- tail(p3.trn, 8)
  prod3new <- c(prod3new, p3_frc4)
  prod3new <- prod3new[i:(7+i)]
  prod3new <- prod3new[8:1]
  prod3new <- array(prod3new, c(1,8))
  colnames(prod3new) <- paste0("lag", 1:8)
  prod3new <- as.data.frame(prod3new)
  p3_frc4[i] <- predict(p3_fit5, prod3new)
}

p3_frc3 <- ts(p3_frc3, frequency = frequency(p3.trn), start = start(p3.trn))
p3_frc4 <- ts(p3_frc4, frequency = frequency(p3.tst), start = start(p3.tst))
ts.plot(p3.trn,p3.tst,p3_frc4,p3_frc3, col = c("black", "black", "red","red"))

MSE_p3_lag <- mean((p3.tst - p3_frc4)^2)
sqrt(MSE_p3_lag)

MSE_p3_lag_1step <- mean((p3.tst[1] - p3_frc4[1])^2)
sqrt(MSE_p3_lag_1step)

#lag model with trend handling for P3

prod3diff <- prod3
for (i in 1:ncol(prod3diff)) {
  prod3diff[,i] <- c(NA, diff(prod3diff[,i]))
}
prod3diff

p3_fit7 <- step(lm(p3~., prod3diff[-(1:9),]))
summary(p3_fit7)
test <- cumsum(c(p3.trn[9], p3_fit7$fitted.values))
test <- test[-1]
length(test)

test_mse <- mean((tail(p3.trn,68) - test)^2)
sqrt(test_mse)
sqrt(p3_fit2$mse)

p3_frc6 <- array(NA,c(10,1))
for (i in 1:10){
  
  p3.diff1 <- diff(p3.trn)
  prod3new3 <- tail(p3.diff1,8)
  prod3new3 <- c(prod3new3,p3_frc6)
  prod3new3 <- prod3new3[i:(7+i)]
  prod3new3 <- prod3new3[8:1]
  prod3new3 <- array(prod3new3, c(1,8))
  colnames(prod3new3) <- paste0("lag",1:8)
  prod3new3 <- as.data.frame(prod3new3)
  p3_frc6[i] <- predict(p3_fit7,prod3new3)
}

p3_frc6ud <- cumsum(p3_frc6) + as.vector(tail(p3.trn, 1))

p3_frc6ud <- ts(p3_frc6ud, frequency = frequency(p3.tst), start = start(p3.tst))
ts.plot(p3.trn, p3.tst, p3_frc4, p3_frc6ud, col = c("black", "black", "red", "blue"))
MSE_p3_diff <- mean((p3.tst - p3_frc6ud)^2)
sqrt(MSE_p3_diff)

MSE_p3_diff_1step <- mean((p3.tst[1] - p3_frc6ud[1])^2)

#dummies for seasonality P3

startD

D_p3 <- as.factor(rep(c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday"),11))
D_p3

p3_dummies <- cbind(prod3, D_p3)
colnames(p3_dummies) <- c(colnames(p3_dummies)[1:9],"D")
p3_dummies

p3_fit_d <- lm(p3~., p3_dummies)
summary(p3_fit_d)

idx <- is.na(p3_dummies)
idx <- rowSums(idx)
idx <- idx == 0 
idx

p3_dummies[idx,]

p3_fit_d_x <- lm(p3~., data = p3_dummies[idx,])
p3_fit_d_1 <- step(p3_fit_d_x)
summary(p3_fit_d_1)

c(AIC(p3_fit7), AIC(p3_fit_d_1), AIC(p3_fit5))

p3_frc_dum <- predict(p3_fit_d_1, p3_dummies)
ts.plot(p3.trn, p3_frc_dum, p3_frc_lag, p3_fit2$fitted,col = c("black", "green", "red", "blue"))
legend("bottomright", c("More reactive EXSM", "Autoregressive", "Dummies"), col = c("blue", "red", "green"), lty = 1)


p3_frc_dum_1 <- array(NA,c(10,1))
for (i in 1:10){
  
  p3_dum <- tail(p3.trn,8)
  p3_dum <- c(p3_dum, p3_frc_dum_1)
  p3_dum <- p3_dum[i:(7+i)]
  p3_dum <- p3_dum[8:1]
  p3_dum <- array(p3_dum, c(1,8))
  colnames(p3_dum) <- paste0("lag",1:8)
  p3_dum <- as.data.frame(p3_dum)
  
  D <- as.factor(rep(c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday"),2)[i])
  p3_dum <- cbind(p3_dum, D)
  p3_frc_dum_1[i] <- predict(p3_fit_d_1, p3_dum)
}

p3_frc_dum_1 <- ts(p3_frc_dum_1, frequency = frequency(p3.tst), start = start(p3.tst))
ts.plot(p3.trn, head(p3.tst,10), tail(p3.tst, 13), p3_frc_dum_1, p3_frc4, p3_frc6ud, col = c("black", "black","black", "red", "blue", "magenta"))

MSE_p3_dum <- mean((p3.tst - p3_frc_dum_1)^2)
sqrt(MSE_p3_dum)

MSE_p3_all <- c(MSE_p3.1, MSE_p3.2, MSE_p3_lag, MSE_p3_diff, MSE_p3_dum) 
sqrt(MSE_p3_all)

MSE_p3_dum_1step <- mean((p3.tst[1] - p3_frc_dum_1[1])^2)
MSE_p3_all_1step <- c(MSE_p3_1step, MSE_p3_lag_1step, MSE_p3_diff_1step, MSE_p3_dum_1step)
sqrt(MSE_p3_all_1step)


plot(p3_frc2, xlim = c(42,45))
lines(p3_frc4,col="red")
lines(p3_frc_dum_1, col = "green")
lines(p3_frc6ud, col="magenta")
lines(p3.tst, lty = 2)
legend("bottomright", c("More reactive EXSM", "Autoregressive", "Dummies", "Differences"), col = c("blue", "red", "green","magenta"), lty = 1)


#forecast for P4

p4.trn <- head(data_p4, 77)
p4.tst <- tail(data_p4, 23)


p4_fit7 <- ets(p4.trn, model = "ANA")
p4_fit7
plot(p4.trn)
lines(p4_fit7$fitted, col = "red")


p4_fit8 <- ets(p4.trn, model = "ANA", alpha = 0.1, gamma = 0.05)
p4_fit8
plot(p4.trn)
lines(p4_fit8$fitted, col = "red")
lines(p4_fit7$fitted, col = "blue")

sqrt(p4_fit7$mse)
sqrt(p4_fit8$mse)

p4_frc7 <- forecast(p4_fit7, h = 10)
p4_frc8 <- forecast(p4_fit8, h = 10)



plot(p4_frc8)
lines(p4_fit8$fitted,col="blue")
lines(p4_frc3, col = "red")
lines(p4_frc_lag, col = "red")
lines(p4.tst, lty = 2)

MAE_p4.1 <- mean(abs(p4.tst - p4_frc7$mean))
MAE_p4.2 <- mean(abs(p4.tst - p4_frc8$mean))
MAE_p4 <- c(MAE_p4.1, MAE_p4.2)
names(MAE_p4) <- paste0("Forecast ",1:2) 
round(MAE_p4,3)

MSE_p4.1 <- mean((p4.tst - p4_frc7$mean)^2)
MSE_p4.2 <- mean((p4.tst - p4_frc8$mean)^2)
MSE_p4 <- c(MSE_p4.1, MSE_p4.2) 
RMSE_p4 <- sqrt(MSE_p4) 
names(RMSE_p4) <- paste0("Forecast ",1:2) 
round(RMSE_p4,3)

#correlation for P4

y <- data[, c(5, 7, 8)]
plot(y)
cor(y)

z <- cbind(diff(data$ProductP4), diff(data$price), diff(data$temperature))
plot(z)
cor(z)

#linear regression for P4

p4_fit1 <- lm(ProductP4 ~ price, data=data)
print(p4_fit1)
sum_p4 <- summary(p4_fit1)
sum_p4
pval <- sum_p4$coefficients[,4]
pval <= 0.05

par(mfrow=c(2,2)) # Split into 2 by 2
plot(p4_fit1)
par(mfrow=c(1,1))

newprice <- as.data.frame(data$price[78:87])
newprice

p4_frc1 <- predict(p4_fit1, data)
p4_frc1 <- head(p4_frc1, 87)
p4_frc1 <- ts(p4_frc1, frequency = frequency(p4.trn), start = start(p4.trn))

plot(p4_frc8)
lines(p4_fit7$fitted,col="blue")
lines(p4_frc1, col = "red")
lines(p4.tst, lty = 2)


c(AIC(p4_fit1), AIC(p4_fit2))

#lag model for P4


pacf(p4.trn)

n <- length(p4.trn)
n

prod4 <- array(NA, c(n,8))
for(i in 1:8){
  prod4[i:n,i] <- p4.trn[1:(n-i+1)]
  
}
colnames(prod4) <- c("p4",paste0("lag",1:7))
prod4[1:10,]
prod4[(n-8):n,]

prod4 <- as.data.frame(prod4)
plot(prod4)


p4_fit3 <- lm(p4~.,data = prod4)
summary(p4_fit3)

p4_fit4 <- step(p4_fit3)
summary(p4_fit4)

c(AIC(p4_fit1),AIC(p4_fit2),AIC(p4_fit3),AIC(p4_fit4))

plot(prod4$p4, type = "l")
p4_frc3 <- predict(p4_fit4, prod4)
lines(p4_frc3, col = "red")

prod4new <- array(tail(p4.trn, 7), c(1,7))
colnames(prod4new) <- paste0("lag", 7:1)
prod4new <- as.data.frame(prod4new)
prod4new

predict(p4_fit4, prod4new)

p4_frc_lag <- array(NA,c(10,1))

for (i in 1:10) {
  prod4new <- tail(p4.trn, 7)
  prod4new <- c(prod4new, p4_frc_lag)
  prod4new <- prod4new[i:(6+i)]
  prod4new <- prod4new[7:1]
  prod4new <- array(prod4new, c(1,7))
  colnames(prod4new) <- paste0("lag", 1:7)
  prod4new <- as.data.frame(prod4new)
  p4_frc_lag[i] <- predict(p4_fit4, prod4new)
}

p4_frc_lag <- ts(p4_frc_lag, frequency = frequency(p4.tst), start = start(p4.tst))
p4_frc3 <- ts(p4_frc3, frequency = frequency(p4.trn), start = start(p4.trn))
ts.plot(p4.trn,p4.tst,p4_frc_lag, p4_frc3, col = c("black", "black", "red", "red"))

ts.plot(p4.trn, p4_frc3,p4_fit8$fitted, p4_fit_d_1$fitted.values, col = c("black", "red", "blue", "green"))
legend("topleft", c("More reactive EXSM", "Autoregressive"), col = c("blue", "red"), lty = 1)

MSE_p4_lag_insample <- mean((tail(p4.trn, 70) - p4_frc3)^2)
sqrt(MSE_p4_lag_insample)
sqrt(p4_fit8$mse)

MSE_p4_lag <- mean((p4.tst - p4_frc_lag)^2)
sqrt(MSE_p4_lag)

MSE_p4_lag_1step <- mean((p4.tst[1] - p4_frc_lag[1])^2)
sqrt(MSE_p4_lag_1step)




#lags with price P4

plot(as.vector(data$price), as.vector(data_p4), ylab = "P4 Sales" , xlab = "price")
abline(lm(data_p4~data$price), col = "red")

price <- c(data$price[1:(length(p4.trn))])

prod4_2 <- cbind(prod4, price)
p4_fit5 <- step(lm(p4~., prod4_2))

summary(p4_fit5)

p4_frc_prc <- predict(p4_fit5, prod4_2)

p4_frc_price <- array(NA,c(10,1))
for (i in 1:10){
  
  prod4new2 <- tail(p4.trn,7)
  prod4new2 <- c(prod4new2,p4_frc_price)
  prod4new2 <- prod4new2[i:(6+i)]
  prod4new2 <- prod4new2[7:1]
  prod4price <- tail(price,24)
  prod4price <- prod4price[i]
  prod4new2 <- c(prod4new2,prod4price)
  prod4new2 <- array(prod4new2, c(1,8))
  colnames(prod4new2) <- c(paste0("lag",1:7),"price")
  prod4new2 <- as.data.frame(prod4new2)
  p4_frc_price[i] <- predict(p4_fit5,prod4new2)
}


p4_frc_price <- ts(p4_frc_price, frequency = frequency(p4.tst), start = start(p4.tst))
ts.plot(p4.trn, p4.tst, p4_frc_lag, p4_frc_price, col = c("black", "black", "red", "blue"))

MSE <- c(mean((p4.tst - p4_frc_lag)^2), mean((p4.tst - p4_frc_price)^2))
sqrt(MSE)

MSE_p4_price_insample <- mean((tail(p4.trn, 70) - tail(p4_frc_prc,70))^2)
sqrt(MSE_p4_price_insample)

#seasonality for P4


D_p4 <- as.factor(rep(c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday"),11))
D_p4

p4_dummies <- cbind(prod4, D_p4)
colnames(p4_dummies) <- c(colnames(p4_dummies)[1:8],"D")
p4_dummies

p4_fit_d <- lm(p4~., p4_dummies)
summary(p4_fit_d)

idx1 <- is.na(p4_dummies)
idx1 <- rowSums(idx1)
idx1 <- idx1 == 0 
idx1

p4_dummies[idx1,]

p4_fit_d_1 <- step(lm(p4~., data = p4_dummies[idx1,]))
summary(p4_fit_d_1)
p4_frc_dum <- predict(p4_fit_d_1, p4_dummies)


p4_frc_dum_1 <- array(NA,c(10,1))
for (i in 1:10){
  
  p4_dum <- tail(p4.trn,7)
  p4_dum <- c(p4_dum, p4_frc_dum_1)
  p4_dum <- p4_dum[i:(6+i)]
  p4_dum <- p4_dum[7:1]
  p4_dum <- array(p4_dum, c(1,7))
  colnames(p4_dum) <- paste0("lag",1:7)
  p4_dum <- as.data.frame(p4_dum)
  
  D <- as.factor(rep(c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday"),2)[i])
  p4_dum <- cbind(p4_dum, D)
  p4_frc_dum_1[i] <- predict(p4_fit_d_1, p4_dum)
}

p4_frc_dum_1 <- ts(p4_frc_dum_1, frequency = frequency(p4.tst), start = start(p4.tst))
ts.plot(p4.trn, head(p4.tst,10), tail(p4.tst, 13), p4_frc_dum_1, p4_frc_lag, p4_frc_price, col = c("black", "black","black", "red", "blue" , "green"))


ts.plot(p4.trn, p4_frc3,p4_fit8$fitted, p4_frc_dum, p4_frc_prc, col = c("black", "red", "blue", "green","magenta"))
legend("topleft", c("More reactive EXSM", "Autoregressive", "Dummies", "Price"), col = c("blue", "red", "green", "magenta"), lty = 1)

ts.plot(p4.trn, p4_frc3,p4_fit8$fitted, col = c("black", "red", "blue"))
legend("topleft", c("More reactive EXSM", "Autoregressive"), col = c("blue", "red"), lty = 1)


MSE_p4_dum_insample <- mean((p4.trn - p4_frc_dum)^2)
sqrt(MSE_p4_dum_insample)

c(AIC(p4_fit_d_1), AIC(p4_fit5), AIC(p4_fit4))

plot(p4_frc8, xlim = c(42, 45))
lines(p4.tst, lty = 2)
lines(p4_frc_lag, col = "red")
lines(p4_frc_dum_1, col = "green")
lines(p4_frc_price, col = "magenta")
legend("topleft", c("More reactive EXSM", "Autoregressive", "Dummies", "Price"), col = c("blue", "red", "green", "magenta"), lty = 1)

MSE_p4 <- c(mean((p4.tst - p4_frc8$mean)^2), mean((p4.tst - p4_frc_lag)^2), mean((p4.tst - p4_frc_dum_1)^2), mean((p4.tst - p4_frc_price)^2))
sqrt(MSE_p4)

MSE_p4_1step <- c(mean((p4.tst[1] - p4_frc8$mean[1])^2), mean((p4.tst[1] - p4_frc_lag[1])^2), mean((p4.tst[1] - p4_frc_dum_1[1])^2), mean((p4.tst[1] - p4_frc_price[1])^2))
sqrt(MSE_p4_1step)


#P5 forecast

library(tsintermittent)


p5.trn <- head(data_p5, 80)
p5.tst <- tail(data_p5, 20)

p5_fit1 <- crost(p5.trn, h = 10, type='croston', outplot=TRUE)


p5_cro_in <- ts(p5_fit1$frc.in, frequency = frequency(p5.trn), start = start(p5.trn))
p5_cro_out <- ts(p5_fit1$frc.out, frequency = frequency(p5.tst), start = start(p5.tst))


p5_fit2 <- crost(p5.trn, h = 10, type='sba', outplot=TRUE)


p5_sba_in <- ts(p5_fit2$frc.in, frequency = frequency(p5.trn), start = start(p5.trn))
p5_sba_out <- ts(p5_fit2$frc.out, frequency = frequency(p5.tst), start = start(p5.tst))

plot(data_p5)
lines(p5_cro_in, col = "red")
lines(p5_cro_out, col = "red")
lines(p5_sba_in, col = "blue")
lines(p5_sba_out, col = "blue")
legend("topleft", c("CRO","SBA"), col = c("red","blue"), lty = 1)



MSE_p5_1 <- mean((head(p5.tst,10) - p5_fit1$frc.out)^2)
MSE_p5_2 <- mean((head(p5.tst,10) - p5_fit2$frc.out)^2)
c(MSE_p5_1, MSE_p5_2)

MSE_p5_11step <- mean((p5.tst[1] - p5_fit1$frc.out[1])^2)
MSE_p5_21step <- mean((p5.tst[1] - p5_fit2$frc.out[1])^2)
c(MSE_p5_11step, MSE_p5_21step)

MSE_p5_1_insample <- mean(tail((p5.trn - p5_cro_in)^2), 77)
MSE_p5_2_insample <- mean(tail((p5.trn - p5_sba_in)^2), 77)
c(MSE_p5_1_insample, MSE_p5_2_insample)









