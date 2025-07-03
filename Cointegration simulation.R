rm(list = ls())
library(tseries)
library(vars)
library(urca)
library(forecast)
# Xt = ψ1X_t-1 + φ1Y_t-1 + ∈_1t
# Yt = ψ2X_t-1 + φ2Y_t-1 + ∈_2t
# Data Generation, T = 200
# ---- no mean ----
T = 200
set.seed(0)

X = numeric(T)
Y = numeric(T)

X[1] = rnorm(1,0,1)
Y[1] = rnorm(1,0,1)

ψ1 = 0.5
φ1 = 0.3

ψ2 = 0.2
φ2 = 0.8

e_1 = e_2 = numeric(T)

e_1[1] = rnorm(1, 0, 1)
e_2[1] = rnorm(1, 0, 1)

# set.seed(0)
for (t in 2:T) {
  e_1[t] = rnorm(1, 0, 1)
  e_2[t] = rnorm(1, 0, 1)
  
  X[t] = (ψ1 * X[t-1]) + (φ1 * Y[t-1]) + e_1[t]
  Y[t] = (ψ2 * X[t-1]) + (φ2 * Y[t-1]) + e_2[t]
}

A = matrix(c(ψ1, ψ2, φ1, φ2), nrow = 2)
eigen(A) # Eigenvalues are both < 1 so my model is stationary :P

# Okay, confirming that the coefficients of my A matrix 
# actually match with VAR model in R 
data = cbind(X,Y)
var_model = VAR(ts(data), p = 1, type = "none")
summary(var_model)
# Eigenvalues are 0.937 and 0.362, so model is stationary. 
# Training split, 80%
train = ts(data[1:160, ])
test  = ts(data[161:200, ])
var_fit = VAR(train, p = 1, type = "none")
preds = predict(var_fit, n.ahead = 100)

accuracy(preds$fcst$X[,1], test[,1])
accuracy(preds$fcst$Y[,1], test[,2])

# Now adding µ_t - this is a mean increasing in time to make model non-stationary
set.seed(0)

µ_t = 1:T # Random

# Regenerating my X and Y
X = Y = numeric(T)

e_1 = e_2 = numeric(T)

e_1[1] = rnorm(1, 0, 1)
e_2[1] = rnorm(1, 0, 1)

X[1] = rnorm(1, 0, 1)
Y[1] = rnorm(1, 0, 1)

ψ1 = 0.5
φ1 = 0.3

ψ2 = 0.2
φ2 = 0.8

# set.seed(0)
for (t in 2:T) {
  e_1[t] = rnorm(1, 0, 1)
  e_2[t] = rnorm(1, 0, 1)
  
  # Subtracting µ_t inside X_t-1 & Y_t-1 should take care of µ_t throughout??
  X[t] = (ψ1 * (X[t-1] - µ_t[t])) + (φ1 * (Y[t-1] - µ_t[t])) + µ_t[t] + e_1[t]
  Y[t] = (ψ2 * (X[t-1] - µ_t[t])) + (φ2 * (Y[t-1] - µ_t[t])) + µ_t[t] + e_2[t]
}
ts.plot(X)

data = cbind(X,Y)
var_model = VAR(ts(data), p = 1, type = "trend")
summary(var_model)

VARselect(data)
j_test = ca.jo(data, type = "eigen", ecdet = "trend", K = 2)
summary(j_test)
## FULL RANK: means no cointegrating vector

# Ignore
# Reject r = 0, fail to reject r <= 1 so rank = 1
#

# ---- If rank is not full ----
# Now pulling out α and β that decompose to form Π
α = j_test@W[1:2,] # alpha is the loading matrix
β = j_test@V[1:2,] # beta is eigenvectors matrix, both columns since r = 2
Π = α %*% t(β)
Π + diag(2) # diag(2) is just my identity matrix
# Reverse calculation of resulting alpha and beta from Johansen test

# OR
var_model = VAR(ts(data), p = 1, type = "trend")
summary(var_model)

# ---- Back to business ----

# plot(resid(var_model))
# acf(resid(var_model)[,1])
# roots(var_model)

plot(ts(X - Y))
abline(h = 4, col = 'violetred')

# fit = auto.arima((ts(X - Y)))
# summary(fit)
# # tsdisplay(residuals(fit))
# accuracy(fit)
# # forecast(fit, h = 100) |> autoplot()
# # Why does forecast plot look like this ??
adf.test(ts(X - Y))

# Plot looks like X - Y is stationary, my adf.test also 
# has p-value of less than 0.05 so confirm that X - Y is stationary


# Let's try real example - non full rank 
set.seed(0)

µ_t = cumsum(abs(rnorm(T))) # Random

# Regenerating my X and Y
X = Y = numeric(T)

e_1 = e_2 = numeric(T)

e_1[1] = rnorm(1, 0, 1)
e_2[1] = rnorm(1, 0, 1)

X[1] = rnorm(1, 0, 1)
Y[1] = rnorm(1, 0, 1)

ψ1 = 0.5
φ1 = 0.3

ψ2 = 0.2
φ2 = 0.8

# set.seed(0)
for (t in 2:T) {
  e_1[t] = rnorm(1, 0, 1)
  e_2[t] = rnorm(1, 0, 1)
  
  # Subtracting µ_t inside X_t-1 & Y_t-1 should take care of µ_t throughout??
  X[t] = (ψ1 * (X[t-1] - µ_t[t])) + (φ1 * (Y[t-1] - µ_t[t])) + µ_t[t] + e_1[t]
  Y[t] = (ψ2 * (X[t-1] - µ_t[t])) + (φ2 * (Y[t-1] - µ_t[t])) + µ_t[t] + e_2[t]
}
ts.plot(X)

data = cbind(X,Y)
VARselect(data) # Lags of 1/2
j_test = ca.jo(data, type = "eigen", ecdet = "trend", K = 2)
summary(j_test)

# Reject r = 0, fail to reject r <= 1 so rank = 1

# Now pulling out α and β that decompose to form Π
α = j_test@W[1:2,1] # alpha is the loading matrix
β = j_test@V[1:2,1] # beta is eigenvectors matrix, both columns since r = 2
Π = α %*% t(β)
Π + diag(2) # diag(2) is just my identity matrix
# Reverse calculation of resulting alpha and beta from Johansen test

# OR
var_model = VAR(ts(data), p = 1, type = "trend")
summary(var_model)

plot(ts(X - Y))
abline(h = 3, col = 'violetred')



# Try n = 3 (w/o mean term)
set.seed(0)
T = 200
X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = rnorm(1,0,1)
Y[1] = rnorm(1,0,1)
Z[1] = rnorm(1,0,1)

a1 = -0.5; b1 = 0.1; c1 = -0.3
a2 = 0.2;  b2 = -0.7; c2 = 0.4
a3 = 0.4;  b3 = 0.5;  c3 = 0.6

A = matrix(c(a1,a2,a3,b1,b2,b3,c1,c2,c3), nrow = 3,byrow = F)
eigen(A)

# Absolute values of eigenvalues are all < 1 so supposedly stationary ^_^

e_1 = e_2 = e_3 = numeric(T)

e_1[1] = rnorm(1, 0, 1)
e_2[1] = rnorm(1, 0, 1)
e_3[1] = rnorm(1, 0, 1)

# set.seed(0)
for (t in 2:T) {
  e_1[t] = rnorm(1, 0, 1)
  e_2[t] = rnorm(1, 0, 1)
  e_3[t] = rnorm(1, 0, 1)
  
  X[t] = (a1 * X[t-1]) + (b1 * Y[t-1]) + (c1 * Z[t-1]) + e_1[t]
  Y[t] = (a2 * X[t-1]) + (b2 * Y[t-1]) + (c2 * Z[t-1]) + e_2[t]
  Z[t] = (a3 * X[t-1]) + (b3 * Y[t-1]) + (c3 * Z[t-1]) + e_3[t]
}

data = cbind(X,Y,Z)
VARselect(data, type = 'trend')
var_model = VAR(ts(data), p = 1, type = "trend")
summary(var_model)
# Returned coefficients are also pretty close :D

# RESET
# Try n = 3 (w mean term)
# rm(list = ls())
set.seed(0)

X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = rnorm(1)
Y[1] = rnorm(1)
Z[1] = rnorm(1)

µ_t = 1:T

e_1 = numeric(T)
e_2 = numeric(T)
e_3 = numeric(T)

for (t in 2:T) {
  e_1[t] = rnorm(1)
  e_2[t] = rnorm(1)
  e_3[t] = rnorm(1)
  
  X[t] = a1 * (X[t-1] - µ_t[t]) + b1 * (Y[t-1] - µ_t[t]) + 
    c1 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_1[t]
  Y[t] = a2 * (X[t-1] - µ_t[t]) + b2 * (Y[t-1] - µ_t[t]) + 
    c2 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_2[t]
  Z[t] = a3 * (X[t-1] - µ_t[t]) + b3 * (Y[t-1] - µ_t[t]) + 
    c3 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_3[t]
}

VAR_data = ts(cbind(X, Y, Z))

VARselect(VAR_data, type = 'trend')

j_test1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(j_test1)

# Reject r <= 2, so rank is 3
# Full rank means 3 cointegrating vector

###
# IF NOT FULL RANK
# Now pulling out α and β that decompose to form Π
α = j_test1@W[,c(1:3)] # alpha is the loading matrix
β = j_test1@V[1:3,c(1:3)] # beta is eigenvectors matrix, both columns since r = 3
Π = α %*% t(β)
Π + diag(3)

var_model = VAR(ts(VAR_data), p = 1, type = "trend")
summary(var_model)
# Coefficients are close
####

u = X - Y
v = Y - Z
w = X - Z

plot(ts(u))
plot(ts(v))
plot(ts(w))

coin_data = cbind(u,v)
VARselect(coin_data)
j_test2 = ca.jo(coin_data, type = "trace", ecdet = "trend", K = 3)
summary(j_test2)
# Full rank
var_model1 = VAR(coin_data, p = 1, type = "trend")
summary(var_model1)

coin_data1 = cbind(v,w)
VARselect(coin_data1)
j_test3 = ca.jo(coin_data1, type = "trace", ecdet = "trend", K = 3)
summary(j_test3)
# Full rank
var_model2 = VAR(coin_data1, p = 1, type = "trend")
summary(var_model1)

coin_data2 = cbind(u,w)
j_test3 = ca.jo(coin_data2, type = "trace", ecdet = "const", K = 2)
summary(j_test3)
# Full rank
var_model3 = VAR(coin_data2, p = 1, type = "const")
summary(var_model3)


# non-full rank, r = 2
set.seed(0)
T = 100

X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = rnorm(1)
Y[1] = rnorm(1)
Z[1] = rnorm(1)

µ_t = (1:T)/100

e_1 = numeric(T)
e_2 = numeric(T)
e_3 = numeric(T)

a1 = 0.4; b1 = -0.3; c1 = 0.7
a2 = 0.2; b2 = 0.1; c2 = 0.3
a3 = 0.5; b3 = -0.2; c3 = 0.7

for (t in 2:T) {
  e_1[t] = rnorm(1)
  e_2[t] = rnorm(1)
  e_3[t] = rnorm(1)
  
  X[t] = a1 * (X[t-1] - µ_t[t]) + b1 * (Y[t-1] - µ_t[t]) + 
    c1 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_1[t]
  Y[t] = a2 * (X[t-1] - µ_t[t]) + b2 * (Y[t-1] - µ_t[t]) + 
    c2 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_2[t]
  Z[t] = a3 * (X[t-1] - µ_t[t]) + b3 * (Y[t-1] - µ_t[t]) + 
    c3 * (Z[t-1] - µ_t[t]) + µ_t[t] + e_3[t]
}
ts.plot(Z)
VAR_data = ts(cbind(X, Y, Z))

VARselect(VAR_data, type = 'trend') # Lag choice of 1's/2's across the board

j_test1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(j_test1)

# Fail to reject r <= 2, so rank is 2
α = j_test1@W[,c(1:2)] # alpha is the loading matrix
β = j_test1@V[1:3,c(1:2)] # beta is eigenvectors matrix, both columns since r = 3
Π = α %*% t(β)
Π + diag(3)
# Close enough lol

var_model = VAR(ts(VAR_data), p = 1, type = "trend")
summary(var_model)
# Coefficients are close
####

u = X - Y
v = Y - Z
w = X - Z

plot(ts(u))
plot(ts(v))
plot(ts(w))

coin_data = cbind(u,v)
VARselect(coin_data)
j_test2 = ca.jo(coin_data, type = "trace", ecdet = "trend", K = 2)
summary(j_test2)
# Full rank
var_model1 = VAR(coin_data, p = 1, type = "trend")
summary(var_model1)

coin_data1 = cbind(v,w)
VARselect(coin_data1)
j_test3 = ca.jo(coin_data1, type = "trace", ecdet = "trend", K = 2)
summary(j_test3)
# Full rank
var_model2 = VAR(coin_data1, p = 1, type = "trend")
summary(var_model1)

coin_data2 = cbind(u,w)
j_test3 = ca.jo(coin_data2, type = "trace", ecdet = "const", K = 2)
summary(j_test3)
# Full rank
var_model3 = VAR(coin_data2, p = 1, type = "const")
summary(var_model3)

# Non-full rank: Different means
# Rank 0
set.seed(0)
T = 200

rank = 0

trend1 = cumsum(rnorm(T, 0.1, 1))
trend2 = cumsum(rnorm(T, 0.05, 1))
trend3 = cumsum(rnorm(T, 0.02, 1))

if (rank == 0) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = trend3
} else if (rank == 1) {
  µ_t1 = trend1
  µ_t2 = trend1 + rnorm(T, 0, 0.1)
  µ_t3 = trend2
} else if (rank == 2) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = 0.5 * (trend1 + trend2)
}

a1 = 0.5; b1 = 0.2; c1 = 0.1
a2 = 0.1; b2 = 0.4; c2 = 0.2
a3 = 0.2; b3 = 0.1; c3 = 0.3

X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = µ_t1[1] + rnorm(1)
Y[1] = µ_t2[1] + rnorm(1)
Z[1] = µ_t3[1] + rnorm(1)

e_1 = rnorm(T, 0, 1)
e_2 = rnorm(T, 0, 0.3)
e_3 = rnorm(T, 0, 0.6)

for (t in 2:T) {
  X[t] = a1 * (X[t - 1] - µ_t1[t]) + b1 * (Y[t - 1] - µ_t1[t]) +
    c1 * (Z[t - 1] - µ_t1[t]) + µ_t1[t] + e_1[t]
  
  Y[t] = a2 * (X[t - 1] - µ_t2[t]) + b2 * (Y[t - 1] - µ_t2[t]) +
    c2 * (Z[t - 1] - µ_t2[t]) + µ_t2[t] + e_2[t]
  
  Z[t] = a3 * (X[t - 1] - µ_t3[t]) + b3 * (Y[t - 1] - µ_t3[t]) +
    c3 * (Z[t - 1] - µ_t3[t]) + µ_t3[t] + e_3[t]
}

VAR_data = ts(cbind(X, Y, Z))
ts.plot(VAR_data, col = 1:3, rank)

VARselect(VAR_data, lag.max = 10, type = "trend")$selection

j_test = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(j_test)

α = j_test@W[,c(1:2)] # alpha is the loading matrix
β = j_test1@V[1:3,c(1:2)] # beta is eigenvectors matrix, both columns since r = 3
Π = α %*% t(β)
Π + diag(3)

# Rank 1
set.seed(0)
T = 200

rank = 1

trend1 = cumsum(rnorm(T, 0.1, 1))
trend2 = cumsum(rnorm(T, 0.05, 1))
trend3 = cumsum(rnorm(T, 0.02, 1))

if (rank == 0) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = trend3
} else if (rank == 1) {
  µ_t1 = trend1
  µ_t2 = trend1 + rnorm(T, 0, 0.1)
  µ_t3 = trend2
} else if (rank == 2) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = 0.5 * (trend1 + trend2)
}

a1 = 0.5; b1 = 0.2; c1 = 0.1
a2 = 0.1; b2 = 0.4; c2 = 0.2
a3 = 0.2; b3 = 0.1; c3 = 0.3

X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = µ_t1[1] + rnorm(1)
Y[1] = µ_t2[1] + rnorm(1)
Z[1] = µ_t3[1] + rnorm(1)

e_1 = rnorm(T, 0, 1)
e_2 = rnorm(T, 0, 0.3)
e_3 = rnorm(T, 0, 0.6)

for (t in 2:T) {
  X[t] = a1 * (X[t - 1] - µ_t1[t]) + b1 * (Y[t - 1] - µ_t1[t]) +
    c1 * (Z[t - 1] - µ_t1[t]) + µ_t1[t] + e_1[t]
  
  Y[t] = a2 * (X[t - 1] - µ_t2[t]) + b2 * (Y[t - 1] - µ_t2[t]) +
    c2 * (Z[t - 1] - µ_t2[t]) + µ_t2[t] + e_2[t]
  
  Z[t] = a3 * (X[t - 1] - µ_t3[t]) + b3 * (Y[t - 1] - µ_t3[t]) +
    c3 * (Z[t - 1] - µ_t3[t]) + µ_t3[t] + e_3[t]
}

VAR_data = ts(cbind(X, Y, Z))
ts.plot(VAR_data, col = 1:3, rank)

VARselect(VAR_data, lag.max = 10, type = "trend")$selection

j_test = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(j_test)


# Rank 2
set.seed(0)
T = 200

rank = 2

trend1 = cumsum(rnorm(T, 0.1, 1))
trend2 = cumsum(rnorm(T, 0.05, 1))
trend3 = cumsum(rnorm(T, 0.02, 1))

if (rank == 0) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = trend3
} else if (rank == 1) {
  µ_t1 = trend1
  µ_t2 = trend1 + rnorm(T, 0, 0.1)
  µ_t3 = trend2
} else if (rank == 2) {
  µ_t1 = trend1
  µ_t2 = trend2
  µ_t3 = 0.5 * (trend1 + trend2)
}

a1 = 0.5; b1 = 0.2; c1 = 0.1
a2 = 0.1; b2 = 0.4; c2 = 0.2
a3 = 0.2; b3 = 0.1; c3 = 0.3

X = numeric(T)
Y = numeric(T)
Z = numeric(T)

X[1] = µ_t1[1] + rnorm(1)
Y[1] = µ_t2[1] + rnorm(1)
Z[1] = µ_t3[1] + rnorm(1)

e_1 = rnorm(T, 0, 1)
e_2 = rnorm(T, 0, 0.3)
e_3 = rnorm(T, 0, 0.6)

for (t in 2:T) {
  X[t] = a1 * (X[t - 1] - µ_t1[t]) + b1 * (Y[t - 1] - µ_t1[t]) +
    c1 * (Z[t - 1] - µ_t1[t]) + µ_t1[t] + e_1[t]
  
  Y[t] = a2 * (X[t - 1] - µ_t2[t]) + b2 * (Y[t - 1] - µ_t2[t]) +
    c2 * (Z[t - 1] - µ_t2[t]) + µ_t2[t] + e_2[t]
  
  Z[t] = a3 * (X[t - 1] - µ_t3[t]) + b3 * (Y[t - 1] - µ_t3[t]) +
    c3 * (Z[t - 1] - µ_t3[t]) + µ_t3[t] + e_3[t]
}

VAR_data = ts(cbind(X, Y, Z))
ts.plot(VAR_data, col = 1:3, rank)

VARselect(VAR_data, lag.max = 10, type = "const")$selection

j_test = ca.jo(VAR_data, type = "trace", ecdet = "const", K = 2)
summary(j_test)
