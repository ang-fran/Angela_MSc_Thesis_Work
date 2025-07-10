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

ts.plot(d_data[,1]) # Just checking

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
e_t = matrix(rnorm(T * n), nrow = T, ncol = n)

mu1 = cumsum(abs(rnorm(T)))
mu2 = cumsum(abs(rnorm(T)))
mu3 = cumsum(abs(rnorm(T)))
mu = cbind(mu1, mu2, mu3)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
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
ts.plot(d_data, col = 4:6) # Not doing it easily

d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2)) # Much better
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,3]) # Just checking

VARselect(VAR_data)$selection # Saying lag of 5-6

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 5)
summary(jotest)
# Clearly fail to reject r = 0

# rank 1
α = matrix(c(1, 2, 3), nrow = 3)
β = matrix(c(1, 0, -1), nrow = 3)
Π = α %*% t(β)

e_t = matrix(rnorm(T * n), nrow = T, ncol = n)
mu1 = cumsum(abs(rnorm(T)))
mu2 = cumsum(abs(rnorm(T)))
mu3 = cumsum(abs(rnorm(T)))
mu = cbind(mu1, mu2, mu3)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
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

d_data = apply(VAR_data, 2, function(x) diff(x, differences = 2)) # Much better
ts.plot(d_data, col = 4:6)
legend("topleft", legend = colnames(Var_vector), col = 4:6, lty = 1)

ts.plot(d_data[,2]) # Just checking
adf.test(d_data[,3])
adf.test(VAR_data[,1])

VARselect(VAR_data)$selection # Saying lag of 5-6

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
# Had to scale this down
mu1 = 0.1 * cumsum(abs(rnorm(T)))
mu2 = 0.1 * cumsum(abs(rnorm(T)))
mu3 = 0.1 * cumsum(abs(rnorm(T)))
mu = cbind(mu1, mu2, mu3)

# Simulation
for (t in 2:T) {
  Var_vector[t, ] = mu[t, ] + Var_vector[t - 1, ] + Π %*% Var_vector[t - 1, ] + e_t[t,]
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

VARselect(VAR_data)$selection # 10,9,6,10

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly fail to reject r <= 2
jotest@W[,-3] %*% t(jotest@V[,-3])
# Close enough to Π

################### 
# Come back to this to edit
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
