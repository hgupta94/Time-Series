##################### Load and clean data #####################

# Set working directory & load packages
#setwd("C:/Users/hirsh/OneDrive/Desktop/NC State/02 Courses/02 Fall/AA502-Analytics Methods & Aadflications I/Fall 2/Time Series 2/Data")

library(lubridate) #date/time manipulations
library(dplyr) # general data manipulation
library(forecast)
library(fma)
library(tseries)
library(expsmooth)
library(lmtest)
library(zoo)
library(ggplot2)
library(rucm)

# Load data
aq <- read.csv("PM_2_5_Raleigh2.csv")
co <- read.csv("CO_Raleigh.csv")
no <- read.csv("NO_Raleigh.csv")
so2 <- read.csv("SO2_Raleigh.csv")
weather <- read.csv("Weatherdata.csv")

# Convert dates to date format
aq$Date <- as.Date(aq$Date, format = "%m/%d/%Y")
co$Date <- as.Date(co$Date, format = "%m/%d/%Y")
no$Date <- as.Date(no$Date, format = "%m/%d/%Y")
so2$Date <- as.Date(so2$Date, format = "%m/%d/%Y")
weather$DATE <- as.Date(weather$DATE, format = "%m/%d/%Y")

# Roll up from daily to monthly
aq <- aq %>% 
  select(Date, Daily.Mean.PM2.5.Concentration) %>% 
  group_by(year=year(Date),month=month(Date)) %>% 
  summarise(avg.pm = mean(Daily.Mean.PM2.5.Concentration))

co <- co %>% 
  select(Date, Daily.Max.8.hour.CO.Concentration) %>% 
  group_by(year=year(Date), month=month(Date)) %>% 
  summarise(avg.co = mean(Daily.Max.8.hour.CO.Concentration))

no <- no %>% 
  select(Date, Daily.Max.1.hour.NO2.Concentration) %>% 
  group_by(year=year(Date), month=month(Date)) %>% 
  summarise(avg.no = mean(Daily.Max.1.hour.NO2.Concentration))

so2 <- so2 %>% 
  select(Date, Daily.Max.1.hour.SO2.Concentration) %>% 
  group_by(year=year(Date), month=month(Date)) %>% 
  summarise(avg.so2 = mean(Daily.Max.1.hour.SO2.Concentration))

# Convert to time series object
aq.ts <- ts(aq$avg.pm, frequency = 12)
co.ts <- ts(co$avg.co, frequency = 12)
no.ts <- ts(no$avg.no, frequency = 12)
so2.ts <- ts(so2$avg.so2, frequency = 12)

decomp <- stl(aq.ts, s.window = 7)
plot(decomp)

# Split into training and validation
train <- subset(aq.ts, end = length(aq.ts)-6)
test <- subset(aq.ts, start = length(aq.ts)-5)


##################### Fit exponential smoothing models #####################
# Single
aq.ses <- ses(train, initial = "optimal", h=6)
summary(aq.ses)
ses.pred <- forecast(aq.ses, h = 6)
ses.error <- test - ses.pred$mean
ses.mape <- mean(abs(ses.error)/abs(test))
ses.mape

plot(aq.ses)
lines(test, col="red")

# Trending (not damped)
aq.trend <- holt(train, initial = "optimal", h=6)
summary(aq.trend)
trend.pred <- forecast(aq.trend, h = 6)
trend.error <- test - trend.pred$mean
trend.mape <- mean(abs(trend.error)/abs(test))
trend.mape # wose than single

plot(aq.trend)
lines(test, col="red")

# Trending (damped)
aq.damped <- holt(train, initial = "optimal", h=6, damped = T)
summary(aq.damped)
damped.pred <- forecast(aq.damped, h = 6)
damped.error <- test - damped.pred$mean
damped.mape <- mean(abs(damped.error)/abs(test))
damped.mape # wose than single

plot(aq.damped)
lines(test, col="red")

# Holt-Winters (additive)
aq.hw_add <- hw(train, initial = "optimal", h=6, seasonal = "additive")
summary(aq.hw_add)
hw_add.pred <- forecast(aq.hw_add, h = 6)
hw_add.error <- test - hw_add.pred$mean
hw_add.mape <- mean(abs(hw_add.error)/abs(test))
hw_add.mape # wose than single

plot(aq.hw_add)
lines(test, col="red")

# Holt-Winters (multiplicative)
aq.hw_mult <- hw(train, initial = "optimal", h=6, seasonal = "multiplicative")
summary(aq.hw_mult)
hw_mult.pred <- forecast(aq.hw_mult, h = 6)
hw_mult.error <- test - hw_mult.pred$mean
hw_mult.mape <- mean(abs(hw_mult.error)/abs(test))
hw_mult.mape # wose than single

plot(aq.hw_mult)
lines(test, col="red")

# Combine all MAPE results
model <- c("Single", "Trend", "Trend (Damped)", "HW Additive", "HW Mult")
mape <- c(round(ses.mape,4), round(trend.mape,4), round(damped.mape,4), round(hw_add.mape,4), round(hw_mult.mape,4))
comb <- cbind(model, mape)
comb # Single ESM best model


##################### Fit ARMIA model #####################
ndiffs(train, test = "adf") # 1 difference
Acf(train) # 1 MA term
Pacf(train) # 1 AR term

mod1 <- Arima(train, order = c(1,0,1))
summary(mod1)

plot(aq.ts)
lines(stats::lag(mod1$fitted, 1), col="red")
lines(stats::lag(arima.pred$mean), col="blue")

white.lb <- rep(NA, 10)
for(i in 1:10){
  white.lb[i] <- Box.test(mod1$residuals, lag = i, type = "Lj", fitdf = 2)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # We have white noise

Box.test(mod1$residuals, lag = 5, type = "Lj", fitdf = 2)

plot(forecast(mod1, h = 6))

arima.pred <- forecast(mod1, h = 6)
arima.error <- test - stats::lag(arima.pred$mean)
arima.mape <- mean(abs(arima.error)/abs(test))
arima.mape

plot(aq.ts)
lines(stats::lag(mod1$fitted, 1), col="red")
lines(stats::lag(arima.pred$mean), col="blue")



##################### Fit Seasonal ARMIA model #####################
# Using trig functions
mod2 <- Arima(train, order = c(1,1,1), xreg = fourier(train, K=4))
summary(mod2) # better than regular ARIMA

s1.arima.pred <- forecast(mod2, h = 6, xreg = fourier(train, K=4))
s1.arima.error <- test - s1.arima.pred$mean
s1.arima.mape <- mean(abs(s1.arima.error)/abs(test))
s1.arima.mape # still worse than single esm

plot(aq.ts)
lines(mod2$fitted, col="red")
lines(s1.arima.pred$mean, col="blue")

white.lb <- rep(NA, 10)
for(i in 1:10){
  white.lb[i] <- Box.test(mod2$residuals, lag = i, type = "Lj", fitdf = 2)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # We have white noise


# Add seasonal terms
Acf(aq.ts) # looks like 1 seasonal MA term
Pacf(aq.ts) # looks like 1 seasonal AR term
nsdiffs(aq.ts) # no seasonal difference


mod3 <- Arima(train, order = c(1,1,1), seasonal = c(1,0,1), xreg = fourier(train, K=4))
summary(mod3) # better than just using trig

s2.arima.pred <- forecast(mod3, h = 6, xreg = fourier(train, K=4))
s2.arima.error <- test - s2.arima.pred$mean
s2.arima.mape <- mean(abs(s2.arima.error)/abs(test))
s2.arima.mape # still worse than single esm

white.lb <- rep(NA, 10)
for(i in 1:10){
  white.lb[i] <- Box.test(mod3$residuals, lag = i, type = "Lj", fitdf = 3)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # Not white noise

plot(aq.ts)
lines(mod3$fitted, col="red")
lines(s1.arima.pred$mean, col="blue")



mod4 <- Arima(train, order = c(1,1,1), seasonal = c(1,0,0), xreg = fourier(train, K=4)) # play with seasonal terms to see if we can get white noise
summary(mod4)

white.lb <- rep(NA, 10)
for(i in 1:10){
  white.lb[i] <- Box.test(mod4$residuals, lag = i, type = "Lj", fitdf = 3)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # White noise

s3.arima.pred <- forecast(mod4, h = 6, xreg = fourier(train, K=4))
s3.arima.error <- test - s3.arima.pred$mean
s3.arima.mape <- mean(abs(s3.arima.error)/abs(test))
s3.arima.mape # still worse than single esm

plot(aq.ts)
lines(mod4$fitted, col="red")
lines(s3.arima.pred$mean, col="blue")


##################### Fit ARMIAX model #####################
# Roll up weather data
weather <- weather[year(weather$DATE)<2019,-c(1,2)]

w.awnd <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.awnd=mean(AWND))
w.awnd.ts <- ts(w.awnd$avg.awnd, frequency = 12)

w.prcp <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(sum.prcp=sum(PRCP))
w.prcp.ts <- ts(w.prcp$sum.prcp, frequency = 12)

w.snow <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(sum.snow=sum(SNOW))
w.snow.ts <- ts(w.snow$sum.snow, frequency = 12)

w.snwd <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.snwd=mean(SNWD))
w.snwd.ts <- ts(w.snwd$avg.snwd, frequency = 12)

w.tavg <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.tavg=mean(TAVG))
w.tavg.ts <- ts(w.tavg$avg.tavg, frequency = 12)

w.tmax <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.tavg=mean(TMAX))
w.tmax.ts <- ts(w.tmax$avg.tavg, frequency = 12)

w.tmin <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.tavg=mean(TMIN))
w.tmin.ts <- ts(w.tmin$avg.tavg, frequency = 12)

w.wsf2 <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.tavg=mean(WSF2))
w.wsf2.ts <- ts(w.wsf2$avg.tavg, frequency = 12)

w.wsf5 <- weather %>% 
  group_by(year=year(DATE), month=month(DATE)) %>% 
  summarise(avg.tavg=mean(WSF5))
w.wsf5.ts <- ts(w.wsf5$avg.tavg, frequency = 12)
  
plot(w.tavg.ts)
# ADF test on all x variables
ndiffs(co.ts, test = "adf")
ndiffs(no.ts, test = "adf")
ndiffs(so2.ts, test = "adf")
ndiffs(w.awnd.ts, test = "adf")
ndiffs(w.prcp.ts, test = "adf")
ndiffs(w.snow.ts, test = "adf")
ndiffs(w.snwd.ts, test = "adf")
ndiffs(w.tavg.ts, test = "adf")
ndiffs(w.tmax.ts, test = "adf")
ndiffs(w.tmin.ts, test = "adf")
ndiffs(w.wsf2.ts, test = "adf")
ndiffs(w.wsf5.ts, test = "adf")
# no diffs


# Combine all variables into a dataframe
# removed snwd, tmin/tmax, wsf5 since they have high cor with snow/tavg; wsf5 had NAs
all.vars <- data.frame(aq.ts, co.ts, no.ts, so2.ts, w.awnd.ts, w.prcp.ts, w.snow.ts, w.tavg.ts, w.wsf2.ts)

# Split all.x into same train/test set as airquality variable

train.x <- all.vars[1:54,]
test.x <- all.vars[55:60,]


# Backwards regression on airquality to see which x vars to use
mod.x <- lm(aq.ts ~ ., data = train.x)
step(mod.x, direction = "backward")

final.x <- train.x[,c("co.ts", "so2.ts", "w.tavg.ts")]

xreg <- as.matrix(final.x)
mod5 <- Arima(train, order = c(1,1,1), xreg = xreg)
summary(mod5)

arimax.pred <- forecast(mod5, h = 6, xreg = xreg)
arimax.error <- test - arimax.pred$mean
arimax.mape <- mean(abs(arimax.error)/abs(test))
arimax.mape # still worse than single esm

white.lb <- rep(NA, 10)
for(i in 1:10){
  white.lb[i] <- Box.test(mod5$residuals, lag = i, type = "Lj", fitdf = 2)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # White noise?

plot(aq.ts)
lines(mod5$fitted, col="red")
lines(arimax.pred$mean, col="blue")


# Try to add seasonality?
mod6 <- Arima(train, order = c(1,1,1), seasonal = c(0,0,1), xreg = xreg)
summary(mod6)


arimax2.pred <- forecast(mod6, h = 6, xreg = xreg)
arimax2.error <- test - arimax2.pred$mean
arimax2.mape <- mean(abs(arimax2.error)/abs(test))
arimax2.mape # still worse than single esm

white.lb <- rep(NA, 10)
for(i in 6:15){
  white.lb[i] <- Box.test(mod6$residuals, lag = i, type = "Lj", fitdf = 3)$p.value
}
barplot(white.lb, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags") # White noise

plot(aq.ts)
lines(mod6$fitted, col="red")
lines(arimax2.pred$mean, col="blue")


##################### Fit Neural Network model #####################
nn.mod <- nnetar(train, p = 2, P = 1, size = 2)
nn.fcast <- forecast(nn.mod, h=6)
plot(nn.fcast)

xreg1 <- cbind(fourier(train, K=5), seq(1,length(train)))
colnames(xreg1) <- c('s1','c1','s2','c2','s3','c3','s4','c4','s5','c5','time')

mod7 <- Arima(diff(train), order = c(1,1,1), xreg = xreg1)
nn.mod2 <- nnetar(mod7$residuals, p=1,P=1,size = 2)
nn.fcast2 <- forecast(nn.mod2, h=6)
plot(nn.fcast2)
lines(nn.mod2$fitted, col="red")

