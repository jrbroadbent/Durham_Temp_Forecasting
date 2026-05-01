
# install.packages("forecast", dependencies = TRUE)
library(dplyr)
library(skimr)
library(lubridate)
library(tseries)
library(forecast)

#### Load Data ####
# NOTE: may need to change path to data
temp_data <- read.csv("~/Dev/EES/durhamtemp_1901_2019.csv", header=TRUE, sep = ",")
head(temp_data)
dim(temp_data)

# remove NAs 
skim(temp_data) # 762 missing
which(is.na(temp_data$Av.temp)) # 43465 onwards is NAs
temp_data$Date[43464] # => this is the last day of 2019 so can remove all NAs
temp_data <- temp_data |> 
  na.omit()
skim(temp_data) # confirm NAs removed 

# produce time series object
temp_data_ts <- ts(temp_data$Av.temp,  start=c(1901, 1), frequency=365)
is.ts(temp_data_ts)


# Figure 1 - plot of the time series 
temp_data$Date <- dmy(temp_data$Date)
plot(temp_data$Date, temp_data$Av.temp, type="l", ylab="Daily average temperature (degC)", xlab="Date (Year)")


# Augmented Dickey-Fuller test => check stationarity, unit root 
adf.test(temp_data$Av.temp) # p-value <0.01 => suggests stationarity and no unit root


# plot of recent years 2017-2019 
Idxs_2017 <- which(temp_data$Year == "2017")
Idxs_2018 <- which(temp_data$Year == "2018")
Idxs_2019 <- which(temp_data$Year == "2019")
Idxs <- c(Idxs_2017, Idxs_2018, Idxs_2019)
temp_recent <- temp_data$Av.temp[Idxs]
date_recent <- temp_data$Date[Idxs]
length(temp_recent) # check 3*365 days selected 
par(mfrow=c(1,2))
plot(date_recent, temp_recent, type="l", xlab="Date (year)", ylab="Daily average temperature (degC)")
plot(date_recent, TTR::SMA(temp_recent, n=30), type="l", xlab="Date (year)", ylab="Daily average temperature (degC)") 
# defintely a visible seasonal trend, moving average of 30 time points 


# seasonal decomposition of time series 
decomp <- decompose(temp_data_ts) # use stl instead for more sophisticated decomposition?
autoplot(decomp)  # autoplot from library forecast
# shows seasonal, trend and remainder components of time series 
# remainder looks to have constant variance 
skim(decomp$random) # decompose introduces NAs at start and end of decomps 


temp_ma <- TTR::SMA(temp_data_ts, n=365)
par(mfrow=c(2,1), mar=c(4,4,3,1), cex=1)
plot(temp_data$Date, temp_data_ts, type="l", xlab="Date (year)", ylab="Daily average temperature (degC)")
plot(temp_data$Date, temp_ma, type="l", xlab="Date (year)", ylab="Daily average temperature (degC)")


# abline(a=0, b=0, lty="dashed")

# decomp2 <- stl(temp_data_ts, s.window = 365)
# autoplot(decomp2)


#### Fourier analysis for seasonal component ####
temp_fourier <- fourier(x = temp_data_ts, K = 3) # K = maximum order of fourier terms


#### check for autocorrelations ####
par(mfrow=c(1,2))
acf(temp_data_ts, main="") # for determining if AR component needed; can see there is strong autocorrelation (ACF close to 1)
pacf(temp_data_ts, main="") # for determining order of AR component; partial ACF drops off sharply after 1, but appears to be a pattern 

# Visualising effect of subtracting fourier terms from time series
# temp <- temp_data_ts
# for(i in 1:6){ # 6 because 2*K fourier terms; K = 3 
#   temp <- temp - temp_fourier[,i]
# }
# plot(temp)
# lines(temp_data_ts, col=2)
# acf(temp)


#### Differencing #### 

# first difference 
temp_data_diff <- diff(temp_data_ts, lag=1)

L <- length(temp_data_diff)
plot(temp_data_diff[(L-(365*3)):L], type="l") # doesn't appear to have seasonal component after differencing to remove trend
plot(temp_data_ts[(L-(365*3)):L], type="l") # clear seasonal component here 

# seasonal diff
temp_seasonal <- diff(temp_data$Av.temp, lag=365)
L <- length(temp_seasonal)
plot(temp_seasonal[(L-(365*4)):(L-365)], type="l")

# do both? 
temp_both <- diff(temp_data_diff, lag=365)
L <- length(temp_both)
plot(temp_seasonal[(L-(365*4)):(L-365)], type="l")


# check ACF PACF after differencing with lag 1 
par(mfrow=c(1,2))
acf(temp_data_diff, main="", ylim=c(-0.15, 1))
pacf(temp_data_diff, main="", ylim=c(-1, 0.15))
# ACF now cuts off to zero after differencing 




#### ARIMA models ####

model1 <- arima(temp_data_ts, order=c(2,1,1))
model2 <- arima(temp_data_ts, order = c(2,1,2))

model3 <- Arima(temp_data_ts, order=c(2,1,1), xreg = fourier(temp_data_ts, K = 3))
fc <- forecast(model3, xreg = fourier(temp_data_ts, K = 3, h = 365))
plot(fc$mean)

model4 <- Arima(temp_data_ts, order=c(2,1,2), xreg = fourier(temp_data_ts, K = 3))
fc <- forecast(model4, xreg = fourier(temp_data_ts, K = 3, h = 365))
plot(fc$mean)
# plots indistinguishable by eye so prefer simpler model? 

model5 <- Arima(temp_data_ts, order=c(3,1,1), xreg = fourier(temp_data_ts, K = 3))
fc <- forecast(model5, xreg = fourier(temp_data_ts, K = 3, h = 365))
plot(fc$mean)

model6 <- auto.arima(temp_data_ts, xreg = fourier(temp_data_ts, K = 3), seasonal = FALSE)
model7 <- auto.arima(temp_data_ts, xreg = fourier(temp_data_ts, K = 6), seasonal = FALSE)
model8 <- auto.arima(temp_data_ts, xreg = fourier(temp_data_ts, K = 1), seasonal = FALSE)


#### ARIMA models: check residuals (ACF and white noise) #### 
# condider ACF of residuals => autocorrelations should be non-significant 
# residuals should look like white noise, normally distributed with mean of 0 
par(mfrow=c(1,2))
acf(model5$residuals, lag.max=365*5, main="") # autocorrelation very high at lag 1 => try increasing difference order
hist(model5$residuals, xlab="residuals", main="")
#plot(model3$residuals)
#plot(model3$fitted, model3$residuals) # fitted value vs residual


### AIC
model1$aic
model2$aic
model3$aic
model4$aic
model5$aic

#### Forecast #### 
# model3 selected for forecast 
model5 <- Arima(temp_data_ts, order=c(3,1,1), xreg = fourier(temp_data_ts, K = 3))
fc <- forecast(model5, xreg = fourier(temp_data_ts, K = 3, h = 365))
par(mfrow=c(1,1))
plot(fc$mean, ylim=c(min(fc$lower[,2]), max(fc$upper[,2])), ylab="Daily average temperature (degC)") # looks good 
lines(fc$upper[,2], col=2)
lines(fc$lower[,2], col=2) # 95% of data between lower and upper bounds
# lines(fc$upper[,1], col=3)
# lines(fc$lower[,1], col=3)









