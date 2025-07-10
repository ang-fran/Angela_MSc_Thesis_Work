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

ts.plot(VAR_data, col = 1:3, main = "Simulated VECM Series", ylab = "Value")
legend("topleft", legend = colnames(Var_vector), col = 1:3, lty = 1)

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

VARselect(VAR_data)$selection

jotest = ca.jo(VAR_data, type = "trace", ecdet = "none", K = 2)
summary(jotest)
# Clearly rejecting r <= 2
jotest@W %*% t(jotest@V)
# Close enough to Π


# ---- Increasing mean ----
