rm(list = ls())
library(readr)
library(forecast)
library(tseries)
library(urca)
# install.packages('vars')
require(vars)
# install.packages('tsDyn')
library(tsDyn)
# install.packages('FCVAR')
require(FCVAR)

setwd("~/Brock University related/Master's/Thesis Research/Example Data")
mydata = read.csv('first_data.csv')
head(mydata)
mydata = mydata[,-c(8:53)]
mydata = mydata[-c(46:57),]
mydata
par(mfrow = c(2,2))
del_rates_ts = ts(mydata$del_rates)
# Compute r1

T = dim(mydata)[1]

y = mydata$del_rates[c(1:(T - 1))]
z = mydata$del_rates[c(2:T)]

cor(y, z)

### Correlation of 98.65%

# Check for k = 2
T = dim(mydata)[1]
k = 2
y = mydata$del_rates[c(1:(T - k))]
z = mydata$del_rates[c((1 + k):T)]
cor(y, z)

# strong autocorrelation at k = 2, 97.31%

acf(mydata$del_rates, lag = 50)
pacf(mydata$del_rates)


# Mean model
mean(mydata$del_rates)
mean_model = meanf(del_rates_ts, h = 10)
summary(mean_model)
plot(mean_model)

# AR1 model
T = dim(mydata)[1]
y = mydata$del_rates[c(2:T)]
z = mydata$del_rates[c(1:(T-1))]

ar1_model = lm(y ~ z)
summary(ar1_model)
eHat = residuals(ar1_model)
acf(eHat, lag = 50)


# Or use Arima function
ar_model = arima(del_rates_ts, order = c(1,0,0))
ar_forecast = forecast(ar_model, h = 50)
plot(ar_forecast)
lines(mydata$del_rates, col = 'blue')


# ARIMA model
arima_model = auto.arima(del_rates_ts)
summary(arima_model)
forecasted_arima = forecast(arima_model, h = 50)
plot(forecasted_arima)
lines(mydata$del_rates, col = 'red')

arima_model = arima(del_rates_ts, order = c(1,1,1))
summary(arima_model)
forecasted_arima = forecast(arima_model, h = 50)
plot(forecasted_arima)

# ---- Monday - Thursday ----
# Using lm for ARIMA(0,1,0)
z_t = diff(del_rates_ts)
model = lm(z_t ~ 1)
summary(model)
e = residuals(model)
acf(z_t, main = 'ACF of differenced default rates')
pacf(z_t, main = 'PACF of differenced default rates', ylim = c(-1, 1))
plot(model, which = 1:4)


# Training and Testing accuracy
set.seed(0)
n = length(del_rates_ts)
train_length = floor(0.75 * n)
training_y = del_rates_ts[1:train_length]
testing_y = del_rates_ts[(train_length + 1) : n]

diff_train = diff(training_y)
training_model = lm(diff_train ~ 1)
summary(training_model)

testing_pred = rep(tail(training_y, 1), length(testing_y))
forecasted = c(training_y, testing_pred)

training_resid = residuals(training_model)
tr_rmse = sqrt(mean(training_resid^2))
tr_rmse

test_rmse = sqrt(mean((testing_y - testing_pred)^2))
test_rmse


## Predicting next difference ----
next_dif = coef(training_model)[1]
next_dif
last_val = tail(training_y, 1)
next_val = last_val + next_dif

next_steps = 10
forecasted_steps = numeric(next_steps)
forecasted_steps[1] = next_val

for (i in 2:next_steps) {
  last_val = forecasted_steps[i - 1]
  forecasted_steps[i] = last_val + next_dif
}

future_mse = sqrt(mean((testing_y[1:next_steps] - forecasted_steps)^2))
future_mse # 0.08649
# ----- Is variance constant?? ----
mod = arima(del_rates_ts, order = c(0,1,0))
resids = residuals(mod)
plot.ts(resids)
mean(resids)
abline(h = -0.005, col = 'violetred')

fitted = del_rates_ts - resids
plot(fitted, resids, col = 'skyblue')
abline(h = 0)

plot.ts(resids^2)

FinTS::ArchTest(resids)
# Fail to reject null hypothesis, so no ARCH effects
# Plots also look mostly constant


# ---- TIME SERIES ----
# ---- VAR ----
# ---- Ignore this ----
delinq_ts = ts(mydata$del_rates, start = c(2012, 4), frequency = 4)
gdp_ts = ts(mydata$gdp_growth, start = c(2012, 3), frequency = 4)
infl_ts = ts(mydata$inflation_rate, start = c(2012, 3), frequency = 4)

# Plot to inspect stationarity
ts.plot(delinq_ts)
abline(reg = lm(delinq_ts ~ time(delinq_ts)), col = 'violetred')
## Plot looks like it's non-stationary because the mean decreases over time

# Augmented Dickey-Fuller Test
adf.test(delinq_ts)  # null hypothesis = non-stationary
# p-value is 0.5371, so we fail to reject null hypothesis

delinq_diff = diff(delinq_ts)
plot(delinq_diff)
adf.test(delinq_diff)

log_delinq = log(delinq_ts)
adf.test(diff(log_delinq))
# still non-stationary
# ---- Okay we back ----
# Assuming ts_data is a multivariate time series (e.g., 3 variables)
data = mydata[c(2:45),c(2,3,5,7)]
data$gdp_growth = as.numeric(data$gdp_growth)
data$inflation_rate = as.numeric(data$inflation_rate)

del_ts = ts(data$del_rates, start = c(2012, 4), frequency = 4)
unemploy_ts = ts(data$unemployment_rate, start = c(2012, 4), frequency = 4)
gdp_growth_ts = ts(data$gdp_growth, start = c(2012, 4), frequency = 4)
infl_rate_ts = ts(data$inflation_rate, start = c(2012, 4), frequency = 4)


# ---- Regular VAR model (Ignore this) ----
VAR_data = window(ts.union(del_ts, unemploy_ts, gdp_growth_ts, infl_rate_ts), start = c(2012, 4), end = c(2023, 3))
VAR_est = VAR(VAR_data)
VAR_est
###
summary(VAR_est$varresult$del_ts)$adj.r.squared # 0.9799
summary(VAR_est$varresult$unemploy_ts)$adj.r.squared # 0.3994
summary(VAR_est$varresult$gdp_growth_ts)$adj.r.squared # 0.2321
summary(VAR_est$varresult$infl_rate_ts)$adj.r.squared # 0.3875
###


# # But is model stationary? ----
# stat_model = lm(del_ts ~ unemploy_ts + gdp_growth_ts + infl_rate_ts)
# res_cont = residuals(stat_model)
# adf.test(res_cont)
# # p-value = 0.342, so we fail to reject H0, that data is nonstationary

VARselect(VAR_data)

# Check stationarity
ts.plot(VAR_data, col = 1:4)
legend("topleft", legend = colnames(VAR_data), col = 1:4, lty = 1)

d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2))
ts.plot(d_data, col = 5:8) # Obviously non-stationary due to spike around COVID


# K = 2
j_test = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(j_test)
# Failed to reject for r <=2 so this suggests r = 2
β = j_test@V[,c(1:2)]
α = j_test@W[,c(1:2)]
Π = α %*% t(β)
Π



#####
# Lag order of 2 (K), ecdet= intercept for cointegration, largest eigenvalue of 0.5085923
# 4 tests for r values, and critical values at 10, 5, and 1% levels of confidence

# First hypothesis, r = 0, is the presence of cointegration. Null hypothesis is no cointegration.
# Test statistic is 49.95,which is greater than all critical values at their various levels of confidence
# SO we reject H0: that there is no cointegration.

# Second: H0: r <= 1, vs Ha: r > 1, 20.11 is lesser than 23.52 (for 1% confidence interval), 
# so we fail to reject H0 at the 1% level, i.e. the relationship is at most 1 (since I rejected it be 0)

# Third: r <= 2 vs r > 2: test statistic is less than all critical values, so we fail to reject H0.
# Since second test showed r has to be at most less or equal than 1, and we failed to reject H0 (i.e r <= 2), 
# so estimate of matrix rank is 1??

# So I need linear combo of 1 time series to form a stationary series

# Linear combo has to be formed from the eigenvector components of the eigenvector associated with the largest eigenvalue.
# Largest eigenvector: 0.65, so eigenvector is (1.000, -1.457, -0.3570) 
#####
# s = 1.000*data$unemployment_rate -1.457*data$gdp_growth - 0.357*data$inflation_rate
# plot(s, type = 'l')
# adf.test(s) # p-value = 0.03062, so we reject H0 (that s is not stationary)
#######
vecm_mod = cajorls(j_test, r = 1)
summary(vecm_mod$rlm)
vecm = vec2var(j_test, r =  1)
forecast_vecm = predict(vecm, n.ahead = 10)
par(mfrow = c(1,1))
plot(forecast_vecm)

# Next stop, since cointegration exists: VECM!!
vecm_model2 = VECM(VAR_data, lag = 2, r = 1, include = "const", estim = "ML")
summary(vecm_model2)








