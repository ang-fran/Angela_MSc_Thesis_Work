rm(list = ls())
library(tseries)
library(vars)
library(urca)
library(forecast)

# ---- no mean ----
T = 2000
n = 3 # Number of variables
set.seed(0)

Var_vector = matrix(0, nrow = T, ncol = 3)
colnames(Var_vector) = c("X", "Y", "Z")

# Initial values
Var_vector[1, ] = rnorm(3)

# Rank 0
Π = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Doesn't look stationary for all time series

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff) # Only first diff
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,3]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r = 0

# rank 1
α = matrix(c(1, 2, 3), nrow = 3)
β = matrix(c(1, 0, -1), nrow = 3)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)
ts.plot(VAR_data[,1], col = 1)
ts.plot(VAR_data[,2], col = 2)
ts.plot(VAR_data[,3], col = 3)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Not looking so stationary eitherrr

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,2]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r <= 1
jotest@W[,1] %*% t(jotest@V[,1])


# Rank 2
α = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3, byrow = F)
β = matrix(c(1, 0, -1, 0, 1, -1), nrow = 3, byrow = F)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + (Π %*% Var_vector[t - 1, ]) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# This doesn't look stationary
# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
# Doesn't look stationary
d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2)) # Much better
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)
# Now we're talking

ts.plot(d_data[,2]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r <= 2
jotest@W[,-3] %*% t(jotest@V[,-3])
# Close enough to Π

# Rank 3: Has to be a VAR where eigenvalues of A are all < 1
# My A here is essentially Π since this is lag 1
α = matrix(c(-0.5, -0.5, -0.5, 0.55, 0, -0.5, 0.3, -0.5
             , 0.3), nrow = 3, byrow = F)
β = matrix(c(0.5, 0.5, 0.5, -0.55, 0, 0.5, -0.3, 0.5
             , -0.3), nrow = 3, byrow = F)
Π = α %*% t(β)
eigen(Π + diag(3))

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + (Π %*% Var_vector[t - 1, ]) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)
ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Looks pretty stationary - as expected
ts.plot(VAR_data[,1], col = 1)
ts.plot(VAR_data[,2], col = 2)
ts.plot(VAR_data[,3], col = 3)
# Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly rejecting r <= 2
jotest@W %*% t(jotest@V)
# Close enough to Π

# ---- Does saying there is a mean affect me carrying out the model with no mean? ----
T = 2000
n = 3 # Number of variables
set.seed(0)

Var_vector = matrix(0, nrow = T, ncol = 3)
colnames(Var_vector) = c("X", "Y", "Z")

# Initial values
Var_vector[1, ] = rnorm(3)

# Rank 0
Π = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Doesn't look stationary for all time series

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff) # Only first diff
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,3]) # Just checking

VARselect(VAR_data)$selection

# Case for constant mean
jotest = ca.jo(VAR_data, type = "trace", ecdet = "const", K = 2)
summary(jotest)
# Still clearly fail to reject r = 0, and results are as expected

# Case for increasing(trend) mean
jotest1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(jotest1)
# Still clearly fail to reject r = 0, and results are still expected


# rank 1
α = matrix(c(1, 2, 3), nrow = 3)
β = matrix(c(1, 0, -1), nrow = 3)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)
ts.plot(VAR_data[,1], col = 1)
ts.plot(VAR_data[,2], col = 2)
ts.plot(VAR_data[,3], col = 3)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Not looking so stationary eitherrr

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,2]) # Just checking

VARselect(VAR_data)$selection

# Case for constant mean
jotest = ca.jo(VAR_data, type = "trace", ecdet = "const", K = 2)
summary(jotest)
# Clearly fail to reject r <= 1, and alpha & beta are still pretty
# close to given inputs; Constant mean here is negligble

# Case for increasing(trend) mean
jotest1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(jotest1)
# Same thing here, yippee
jotest@W[,1] %*% t(jotest@V[,1])
jotest1@W[,1] %*% t(jotest1@V[,1])


# Rank 2
α = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3, byrow = F)
β = matrix(c(1, 0, -1, 0, 1, -1), nrow = 3, byrow = F)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + (Π %*% Var_vector[t - 1, ]) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# This doesn't look stationary
# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
# Doesn't look stationary
d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2)) # Much better
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)
# Now we're talking

ts.plot(d_data[,2]) # Just checking

VARselect(VAR_data)$selection

# Case for constant mean
jotest = ca.jo(VAR_data, type = "trace", ecdet = "const", K = 2)
summary(jotest)
# Still clearly fail to reject r = 2, and results are as expected

# Case for increasing(trend) mean
jotest1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(jotest1)
# Yupp

# Clearly fail to reject r <= 2
jotest@W[,-3] %*% t(jotest@V[,-3])
jotest1@W[,-3] %*% t(jotest1@V[,-3])
# Both close enough to Π

# Rank 3: Has to be a VAR where eigenvalues of A are all < 1
# My A here is essentially Π since this is lag 1
α = matrix(c(-0.5, -0.5, -0.5, 0.55, 0, -0.5, 0.3, -0.5
             , 0.3), nrow = 3, byrow = F)
β = matrix(c(0.5, 0.5, 0.5, -0.55, 0, 0.5, -0.3, 0.5
             , -0.3), nrow = 3, byrow = F)
Π = α %*% t(β)
eigen(Π + diag(3))

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = Var_vector[t - 1, ] + (Π %*% Var_vector[t - 1, ]) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)
ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Looks pretty stationary - as expected
ts.plot(VAR_data[,1], col = 1)
ts.plot(VAR_data[,2], col = 2)
ts.plot(VAR_data[,3], col = 3)
# Just checking

VARselect(VAR_data)$selection

# Case for constant mean
jotest = ca.jo(VAR_data, type = "trace", ecdet = "const", K = 2)
summary(jotest)
# Clearly rejecting r <= 2, and weights of constant is practically zero

# Case for increasing(trend) mean
jotest1 = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(jotest1)
# Yeah same
jotest@W %*% t(jotest@V)
jotest1@W %*% t(jotest1@V)
# Close enough to Π


# ---- With constant mean ----
T = 2000
n = 3 # Number of variables
set.seed(0)

Var_vector = matrix(0, nrow = T, ncol = 3)
colnames(Var_vector) = c("X", "Y", "Z")

# Initial values
Var_vector[1, ] = rnorm(3)

# Rank 0
Π = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu = as.matrix(rep(2, times = T*n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Definitely not stationary - to be expected from rank 0

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)
ts.plot(d_data[,3]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r = 0

# rank 1
α = matrix(c(1, 2, 3), nrow = 3)
β = matrix(c(1, 0, -1), nrow = 3)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu = as.matrix(rep(2, times = T*n), nrow = T, ncol = n)
# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Still looks nonstationary
# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,3]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r <= 1
jotest@W[,1] %*% t(jotest@V[,1])
# Similar results as no mean

# Rank 2
α = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3, byrow = F)
β = matrix(c(1, 0, -1, 0, 1, -1), nrow = 3, byrow = F)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu = as.matrix(rep(2, times = T*n), nrow = T, ncol = n)
# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Same here - non stationary :P

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6) # This is not doing for me

d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2)) # Much better
ts.plot(d_data, col = 4:6) # Better
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,1]) # Just checking

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r <= 2
jotest@W[,-3] %*% t(jotest@V[,-3])
# Close enough to Π

# Rank 3: Has to be a VAR where eigenvalues of A are all < 1
# My A here is essentially Π since this is lag 1
α = matrix(c(-0.5, -0.5, -0.5, 0.55, 0, -0.5, 0.3, -0.5
             , 0.3), nrow = 3, byrow = F)
β = matrix(c(0.5, 0.5, 0.5, -0.55, 0, 0.5, -0.3, 0.5
             , -0.3), nrow = 3, byrow = F)
Π = α %*% t(β)
eigen(Π + diag(3))

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu = as.matrix(rep(2, times = T*n), nrow = T, ncol = n)
# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Looks stationary

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly rejecting r <= 2
jotest@W %*% t(jotest@V)
# Close enough to Π


# ---- Increasing mean: Work in progress ----
T = 2000
n = 3 # Number of variables
set.seed(0)

Var_vector = matrix(0, nrow = T, ncol = 3)
colnames(Var_vector) = c("X", "Y", "Z")

# Initial values
Var_vector[1, ] = rnorm(3)

# Rank 0
Π = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
Π_plus = matrix(0, n, n+1)

mu0 = c(1, 2, 3)
mu1 = c(0.025, 0.05, 0.0015)
mu_t = matrix(NA, nrow = T, ncol = n)

ν = -Π %*% (matrix(mu0, ncol = 1)) + (matrix(mu1, ncol = 1))

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

# Simulation
for (t in 2:T) {
  Var_lag = Var_vector[t - 1,]
  x = c(Var_lag, t - 1)
  Var_vector[t, ] = ν + Var_vector[t - 1, ] + 
    (Π_plus %*% x) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Definitely not stationary - to be expected from rank 0 and the fact
# an increasing mean has been added

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6) # Sharp

VARselect(VAR_data, lag.max = 20)$selection # Saying lag of 1/2

jotest = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 2)
summary(jotest)
# Clearly fail to reject r = 0
jotest1 = ca.jo(VAR_data, type = 'trace', ecdet = 'none', K = 2)
summary(jotest1) # Still failing to reject

# rank 1
α = matrix(c(1, 2, 3), nrow = 3)
β = matrix(c(1, 0, -1), nrow = 3)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu0 = c(1, 2, 3)
mu1 = c(0.025, 0.05, 0.0015)
mu_t = matrix(NA, nrow = T, ncol = n)
mu_t[1,] = mu0 + mu1

η = -t(β)%*%mu1

Π_plus = α %*% cbind(t(β), η)

ν = -Π %*% (matrix(mu0, ncol = 1)) + (matrix(mu1, ncol = 1))

# Simulation
for (t in 2:T) {
  Var_lag = Var_vector[t - 1,]
  x = c(Var_lag, t - 1)
  Var_vector[t, ] = ν + Var_vector[t - 1, ] + 
    (Π_plus %*% x) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Not stationary

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6) # Not doing it easily

ts.plot(d_data[,2]) # Just checking
adf.test(d_data[,3])
adf.test(VAR_data[,1])

VARselect(VAR_data, lag.max = 20)$selection # Saying lag of 14,8,7

jotest = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 8)
summary(jotest)
# Clearly fail to reject r <= 1
jotest@W[,1] %*% t(jotest@V[,1])

# Me saying no trend/mean
jotest1 = ca.jo(VAR_data, type = 'trace', ecdet = 'none', K = 2) # Lag 8 shows failing to reject 
summary(jotest1) # Still failing to reject
jotest1@W[,1] %*% t(jotest1@V[,1]) # Def not the same as Π+

# Similar results as no mean

# Rank 2
α = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3, byrow = F)
β = matrix(c(1, 0, -1, 0, 1, -1), nrow = 3, byrow = F)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu0 = c(1, 2, 3)
mu1 = c(0.025, 0.05, 0.0015)
mu_t = matrix(NA, nrow = T, ncol = n)
mu_t[1,] = mu0 + mu1

η = -t(β)%*%mu1

Π_plus = α %*% cbind(t(β), η)

ν = -Π %*% (matrix(mu0, ncol = 1)) + (matrix(mu1, ncol = 1))

# Simulation
for (t in 2:T) {
  Var_lag = Var_vector[t - 1,]
  x = c(Var_lag, t - 1)
  Var_vector[t, ] = ν + Var_vector[t - 1, ] + 
    (Π_plus %*% x) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Not stationary obviously

# Checking my individual time series are stationary by themselves
# or by differencing
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6) # Not doing it easily

d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2))
ts.plot(d_data, col = 4:6)
# Cool
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,2]) # Just checking
adf.test(d_data[,3])
adf.test(VAR_data[,2])

VARselect(VAR_data)$selection # 10,8,6,10

jotest = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 9)
summary(jotest)
# Clearly fail to reject r <= 2
jotest@W[,-3] %*% t(jotest@V[,-3])
# Close enough to Π+: nice!

# Rank 3: Has to be a VAR where eigenvalues of A are all < 1
# My A here is essentially Π since this is lag 1
α = matrix(c(-0.5, -0.5, -0.5, 0.55, 0, -0.5, 0.3, -0.5
             , 0.3), nrow = 3, byrow = F)
β = matrix(c(0.5, 0.5, 0.5, -0.55, 0, 0.5, -0.3, 0.5
             , -0.3), nrow = 3, byrow = F)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu0 = c(1, 2, 3)
mu1 = c(0.025, 0.05, 0.0015)
mu_t = matrix(NA, nrow = T, ncol = n)
mu_t[1,] = mu0 + mu1

η = -t(β)%*%mu1

Π_plus = α %*% cbind(t(β), η)

ν = -Π %*% (matrix(mu0, ncol = 1)) + (matrix(mu1, ncol = 1))

# Simulation
for (t in 2:T) {
  Var_lag = Var_vector[t - 1,]
  x = c(Var_lag, t - 1)
  Var_vector[t, ] = ν + Var_vector[t - 1, ] + 
    (Π_plus %*% x) + e_t[t,]
}

# Convert to time series
VAR_data = ts(Var_vector)

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)
# Not stationary but expected due to trend
d_data = apply(VAR_data, 2, diff)
ts.plot(d_data, col = 4:6) 
# Nice
VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "trend", K = 6)
summary(jotest)
# Clearly rejecting r <= 2
jotest@W %*% t(jotest@V)
# Close enough to Π_plus
