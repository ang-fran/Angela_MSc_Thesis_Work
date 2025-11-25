rm(list = ls())
#install.packages('mvtnorm')
library(mvtnorm)

set.seed(0)

T = 200 
M = 2 

P = matrix(c(0.75, 0.25, 0.65, 0.35), 2, 2, byrow = TRUE)
A = matrix(c(0.5,0.1,0.2,0.4),2,2)
n = 2
mu_list = list(c(0, 0), #s1
               c(3, -2) # s2
               )
Sigma = matrix(c(1,0.2, 0.2,1), 2, 2)
s = numeric(T)
y1 = matrix(0, T, n)
y2 = matrix(0, T, n)
s[1] = sample(1:M, 1)

for(t in 2:T){
  s[t] = sample(1:M, 1, prob = P[s[t-1],])
  } 
mu_t = matrix(0, T, 2)

for(t in 1:T){
  mu_t[t,] = mu_list[[s[t]]]
}

# MSI
y1[1,] = y2[1,] = mu_list[[s[1]]] + rmvnorm(1, rep(0,n), Sigma)

for(t in 2:T){
  y1[t,] = mu_list[[s[t]]] + as.numeric(A %*% y1[t-1,]) +
    rmvnorm(1, mean = rep(0,2), sigma = Sigma)
}
y1

# MSM
for(t in 2:T){
  y2[t,] = mu_list[[s[t]]] + as.numeric(A %*% (y2[t-1,] - mu_list[[s[t-1]]])) +
    rmvnorm(1, mean = rep(0,2), sigma = Sigma)
}
y2

# Estimating AR matrix
# MSI
Y = y1[2:T, ]
X = y1[1:(T-1), ]
Mu = mu_t[2:T, ]  # matrix of mu at each t
Yc = Y - Mu

A_hat = solve(t(X) %*% X) %*% t(X) %*% Yc
A_hat


# MSM
Y = y2[2:T, ]
X = y2[1:(T-1), ] - mu_t[1:(T-1), ]
Mu = mu_t[2:T, ]
Yc = Y - Mu

A_hat_MSM = solve(t(X) %*% X) %*% t(X) %*% Yc
A_hat_MSM

