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