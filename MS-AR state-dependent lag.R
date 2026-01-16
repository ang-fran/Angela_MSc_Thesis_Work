rm(list = ls())
set.seed(0)

# ---- Parameters ----
T = 200
M = 2
alpha_mu = c(5, -5)
alpha_sigma = c(0.5, 0.6)
alpha_phi1 = c(0.8, 0.5)
alpha_phi2 = c(0, 0.4)
P = matrix(c(0.75, 0.25, 0.7, 0.3), 2, 2, byrow = TRUE)

# ---- Simulate states and series ----
s = numeric(T)
y = numeric(T)
s[1] = sample(1:M, 1)
y[1] = rnorm(1, mean = alpha_mu[s[1]], sd = alpha_sigma[s[1]])
s[2] = sample(1:M, 1)
y[2] = rnorm(1, mean = alpha_mu[s[2]], sd = alpha_sigma[s[2]])

for(t in 3:T){
  s[t] = sample(1:M, 1, prob = P[s[t-1],])
  y[t] = alpha_mu[s[t]] + alpha_phi1[s[t]] * (y[t-1] - alpha_mu[s[t-1]]) +
    (s[t] - 1)*(alpha_phi2[s[t]] * (y[t-2] - alpha_mu[s[t-2]])) +
    rnorm(1, 0, alpha_sigma[s[t]])
}

# ---- Initial guesses ----
alpha_mu_curr    = c(4.5, -4.5)
alpha_sigma_curr = c(0.6, 0.8)
alpha_phi1_curr  = c(0.7, 0.3)
alpha_phi2_curr  = c(0.2, 0.1)
P_hat_curr       = matrix(c(0.8, 0.2, 0.6, 0.4), nrow = 2, byrow = T)

# ---- Conditional density function ----
compute_f_y_cond_st = function(y, alpha_mu, alpha_sigma,
                                alpha_phi1, alpha_phi2, M, T){
  # Inputs
  # alpha_mu: Mean
  # alpha_phi's: AR terms
  # alpha_sigma: Variance
  f_y_cond_st = array(0, dim = c(M, M, M, T))
  # [s_t, s_{t-1}, s_{t-2}, time]
  
  for(t in 1:T){
    for(k in 1:M){
      
      sigma2 = alpha_sigma[k]^2
      norm_c = 1 / sqrt(2*pi*sigma2)
      
      if(t == 1){
        resid = y[t] - alpha_mu[k]
        f_y_cond_st[k,,,t] = norm_c * exp(-(resid^2)/(2*sigma2))
      }
      
      if(t == 2){
        for(l in 1:M){
          mean_t = alpha_mu[k] +
            alpha_phi1[k] * (y[t-1] - alpha_mu[l])
          resid = y[t] - mean_t
          f_y_cond_st[k,l,,t] = norm_c * exp(-(resid^2)/(2*sigma2))
        }
      }
      
      if(t >= 3){
        for(l in 1:M){
          for(m in 1:M){
            mean_t = alpha_mu[k] +
              (alpha_phi1[k] * (y[t-1] - alpha_mu[l])) +
              ((k-1) * alpha_phi2[k] * (y[t-2] - alpha_mu[m]))
            resid = y[t] - mean_t
            f_y_cond_st[k,l,m,t] = norm_c * exp(-(resid^2)/(2*sigma2))
          }
        }
      }
    }
  }
  return(f_y_cond_st)
}

# CHECK
A = compute_f_y_cond_st(y, alpha_mu, alpha_sigma,
                    alpha_phi1, alpha_phi2, M, T)



# ---- E-step ----

msar_e_step = function(y, P_hat, alpha_mu, alpha_sigma, alpha_phi1, alpha_phi2, M, T){
  # Step 2: Filtering
  p_s1 = rep(1/M, M)
  
  # t = 1: Probability of s1 given y1
  # f_y_cond_st is [s_t, s_t-1, s_t-2, t]
  p_s1_given_y1 = f_y_cond_st[, 1, 1, 1] * p_s1 
  p_s1_given_y1 = p_s1_given_y1 / sum(p_s1_given_y1)
  
  p_st_stminus1_given_y = vector("list", T)
  
  # t = 2: Joint P(s2, s1 | Y2)
  # Density depends on (s2, s1)
  joint2 = matrix(0, M, M)
  for(i in 1:M){ # s1
    for(j in 1:M){ # s2
      # f is [k=s_t, i=s_t-1, m=s_t-2, t]
      joint2[i,k] = f_y_cond_st[k, i, 1, 2] * P_hat[i, k] * p_s1_given_y1[i]
    }
  }
  p_st_stminus1_given_y[[2]] = joint2 / sum(joint2)
  
  # t >= 3: Joint P(st, st-1 | Yt)
  for(t in 3:T){
    joint_t = matrix(0, M, M)
    for(l in 1:M){   # l = s_{t-1}
      for(k in 1:M){ # k = s_{t}
        
        # For MS-AR(2), y_t depends on s_{t-2}. 
        # We must sum over m = s_{t-2}
        sum_over_m = 0
        for(m in 1:M){
          # prob = P(st | st-1) * P(st-1, st-2 | Y_{t-1})
          prob_prev = p_st_stminus1_given_y[[t-1]][m, l] 
          
          # density = f(yt | st=k, st-1=l, st-2=m)
          dens = f_y_cond_st[k, l, m, t]
          
          sum_over_m = sum_over_m + (dens * P_hat[l, k] * prob_prev)
        }
        joint_t[l, k] = sum_over_m
      }
    }
    p_st_stminus1_given_y[[t]] = joint_t / sum(joint_t)
  }
  
  # Step 3: Smoothing
  # Initialize storage for Smoothed Joint State Probabilities P(s_t, s_{t-1} | Y_T)
  p_st_stminus1_given_YT <- vector("list", T)
  
  # Step 3b: At t = T, smoothed is identical to filtered
  p_st_stminus1_given_YT[[T]] <- p_st_stminus1_given_y[[T]]
  
  # Step 3d: Iterate backward for all t from T-1 down to 2
  for (t in (T-1):2) {
    
    # Matrix to hold the 4 smoothed joint valuations for this specific t
    smoothed_t_matrix <- matrix(0, M, M)
    
    for (l in 1:M) {     # Valuation for s_{t-1}
      for (k in 1:M) {   # Valuation for s_{t}
        
        # --- Step 3a: Initialization (tau = t + 1) ---
        # P(s_{t+1}, s_t, s_{t-1} | Y_{t+1})
        # We calculate this for all possible current states 'k_next' at t+1
        joint_tau <- numeric(M) 
        for (k_next in 1:M) {
          # Density at t+1 depends on s_{t+1}, s_t, s_{t-1}
          dens_next <- f_y_cond_st[k_next, k, l, t+1]
          
          # Num = f(y_{t+1}|s) * P(s_{t+1}|s_t) * P(s_t, s_{t-1}|Y_t)
          num <- dens_next * P_hat[k, k_next] * p_st_stminus1_given_y[[t]][l, k]
          
          # Normalize by marginal likelihood of y_{t+1}
          joint_tau[k_next] <- num / sum(f_y_cond_st[,, , t+1] * p_st_stminus1_given_y[[t]]) # Denom from 2b
        }
        
        # --- Step 3a: Sequential Calculation (tau = t + 2 ... T) ---
        # If t+1 is already T, we jump to the final summation
        if ((t + 1) < T) {
          # 'current_joint' tracks P(s_{tau}, s_{tau-1}, s_t, s_{t-1} | Y_{tau})
          current_joint <- matrix(0, M, M) # Rows: s_{tau-1}, Cols: s_{tau}
          
          # Initialize the first matrix for tau = t+2
          # Here s_{tau-2} is fixed as the 'k' we are currently looping over
          for(s_tau in 1:M) {
            for(s_prev in 1:M) {
              dens_tau <- f_y_cond_st[s_tau, s_prev, k, t+2]
              current_joint[s_prev, s_tau] <- dens_tau * P_hat[s_prev, s_tau] * joint_tau[s_prev]
            }
          }
          current_joint <- current_joint / sum(current_joint) # Normalize per step
          
          # Continue the chain to T
          if ((t + 2) < T) {
            for (tau in (t + 3):T) {
              next_joint <- matrix(0, M, M)
              for (s_tau in 1:M) {
                for (s_prev in 1:M) {
                  # Sum over s_{tau-2}
                  sum_m <- sum(f_y_cond_st[s_tau, s_prev, , tau] * P_hat[s_prev, s_tau] * current_joint[, s_prev])
                  next_joint[s_prev, s_tau] <- sum_m
                }
              }
              current_joint <- next_joint / sum(next_joint)
            }
          }
          # Step 3b: Marginalize out terminal states s_T, s_{T-1}
          smoothed_t_matrix[l, k] <- sum(current_joint)
        } else {
          # If tau=t+1=T, the smoothed joint is just the sum of our initialization
          smoothed_t_matrix[l, k] <- sum(joint_tau)
        }
      }
    }
    # Store results for time t
    p_st_stminus1_given_YT[[t]] <- smoothed_t_matrix / sum(smoothed_t_matrix)
  }
  
  # Step 4: marginalization
  # Initialize storage for P(s_t = k | Y_T)
  p_st_given_YT <- matrix(0, T, M)
  
  for (t in 2:T) {
    # Sum over s_{t-1} to get marginal for s_t
    # In our matrix, l (s_{t-1}) is rows, k (s_t) is columns
    p_st_given_YT[t, ] <- colSums(p_st_stminus1_given_YT[[t]])
  }
  
  return(list(
    smoothed_joint = p_st_stminus1_given_YT,
    smoothed_marginal = p_st_given_YT,
    filtered_joint = p_st_stminus1_given_y
  ))
}


msar_e_step = function(y, P_hat, alpha_mu, alpha_sigma, alpha_phi1, alpha_phi2, M, T){
  
  # --- Step 1: Conditional Densities (f) ---
  f_y_cond_st = array(0, dim = c(M, M, M, T))
  for(t in 1:T){
    for(k in 1:M){ # s_t
      sigma2 = alpha_sigma[k]^2
      norm_c = 1 / sqrt(2*pi*sigma2)
      
      if(t == 1){
        resid = y[t] - alpha_mu[k]
        f_y_cond_st[k,,,t] = norm_c * exp(-(resid^2)/(2*sigma2))
      } else if(t == 2){
        for(l in 1:M){ # s_{t-1}
          mean_t = alpha_mu[k] + alpha_phi1[k] * (y[t-1] - alpha_mu[l])
          resid = y[t] - mean_t
          f_y_cond_st[k,l,,t] = norm_c * exp(-(resid^2)/(2*sigma2))
        }
      } else { # t >= 3
        for(l in 1:M){ # s_{t-1}
          for(m in 1:M){ # s_{t-2}
            # State 2 (k=2) logic for AR(2)
            mean_t = alpha_mu[k] + alpha_phi1[k] * (y[t-1] - alpha_mu[l]) +
              ((k-1) * alpha_phi2[k] * (y[t-2] - alpha_mu[m]))
            resid = y[t] - mean_t
            f_y_cond_st[k,l,m,t] = norm_c * exp(-(resid^2)/(2*sigma2))
          }
        }
      }
    }
  }
  
  # --- Step 2: Filtering Pass (Forward) ---
  p_s1 = rep(1/M, M)
  p_s1_given_y1 = f_y_cond_st[, 1, 1, 1] * p_s1 
  p_s1_given_y1 = p_s1_given_y1 / sum(p_s1_given_y1)
  
  p_st_stminus1_given_y = vector("list", T)
  
  # t = 2: Joint P(s2, s1 | Y2)
  joint2 = matrix(0, M, M)
  for(i in 1:M){ # s1
    for(k in 1:M){ # s2
      joint2[i,k] = f_y_cond_st[k, i, 1, 2] * P_hat[i, k] * p_s1_given_y1[i]
    }
  }
  p_st_stminus1_given_y[[2]] = joint2 / sum(joint2)
  
  # t >= 3: Joint P(st, st-1 | Yt)
  for(t in 3:T){
    joint_t = matrix(0, M, M)
    for(l in 1:M){   # l = s_{t-1}
      for(k in 1:M){ # k = s_{t}
        sum_over_m = 0
        for(m in 1:M){ # m = s_{t-2}
          prob_prev = p_st_stminus1_given_y[[t-1]][m, l] 
          dens = f_y_cond_st[k, l, m, t]
          sum_over_m = sum_over_m + (dens * P_hat[l, k] * prob_prev)
        }
        joint_t[l, k] = sum_over_m
      }
    }
    p_st_stminus1_given_y[[t]] = joint_t / sum(joint_t)
  }
  
  # --- Step 3: Smoothing Pass (Forward Recursion) ---
  p_st_stminus1_given_YT <- vector("list", T)
  p_st_stminus1_given_YT[[T]] <- p_st_stminus1_given_y[[T]] #
  
  for (t in (T-1):2) {
    smoothed_t_matrix <- matrix(0, M, M)
    
    for (l in 1:M) {     # s_{t-1}
      for (k in 1:M) {   # s_{t}
        
        # FIX: Denom for tau = t + 1
        denom_t1 = 0
        for(s_t1 in 1:M) for(s_t in 1:M) for(s_t_minus1 in 1:M) {
          denom_t1 = denom_t1 + (f_y_cond_st[s_t1, s_t, s_t_minus1, t+1] * P_hat[s_t, s_t1] * p_st_stminus1_given_y[[t]][s_t_minus1, s_t])
        }
        
        # Initialization
        joint_tau <- numeric(M) 
        for (k_next in 1:M) {
          num <- f_y_cond_st[k_next, k, l, t+1] * P_hat[k, k_next] * p_st_stminus1_given_y[[t]][l, k]
          joint_tau[k_next] <- num / denom_t1
        }
        
        if ((t + 1) < T) {
          current_joint <- matrix(0, M, M)
          # FIX: Denom for tau = t + 2
          denom_t2 = 0
          for(s_t2 in 1:M) for(s_t1 in 1:M) for(s_t0 in 1:M) {
            denom_t2 = denom_t2 + (f_y_cond_st[s_t2, s_t1, s_t0, t+2] * P_hat[s_t1, s_t2] * p_st_stminus1_given_y[[t+1]][s_t0, s_t1])
          }
          
          for(s_tau in 1:M) {
            for(s_prev in 1:M) {
              current_joint[s_prev, s_tau] <- (f_y_cond_st[s_tau, s_prev, k, t+2] * P_hat[s_prev, s_tau] * joint_tau[s_prev]) / denom_t2
            }
          }
          
          if ((t + 2) < T) {
            for (tau in (t + 3):T) {
              next_joint <- matrix(0, M, M)

              denom_tau = 0
              for(s_c in 1:M) for(s_p in 1:M) for(s_pp in 1:M) {
                denom_tau = denom_tau + (f_y_cond_st[s_c, s_p, s_pp, tau] * P_hat[s_p, s_c] * p_st_stminus1_given_y[[tau-1]][s_pp, s_p])
              }
              for (s_tau in 1:M) {
                for (s_prev in 1:M) {
                  sum_m <- sum(f_y_cond_st[s_tau, s_prev, , tau] * P_hat[s_prev, s_tau] * current_joint[, s_prev])
                  next_joint[s_prev, s_tau] <- sum_m / denom_tau
                }
              }
              current_joint <- next_joint
            }
          }
          smoothed_t_matrix[l, k] <- sum(current_joint) #
        } else {
          smoothed_t_matrix[l, k] <- sum(joint_tau)
        }
      }
    }
    p_st_stminus1_given_YT[[t]] <- smoothed_t_matrix / sum(smoothed_t_matrix)
  }
  
  # --- Step 4: Marginalization ---
  p_st_given_YT <- matrix(0, T, M)
  p_st_given_YT[1, ] <- rowSums(p_st_stminus1_given_YT[[2]]) # Marginalize s1 from P(s2, s1)
  for (t in 2:T) {
    p_st_given_YT[t, ] <- colSums(p_st_stminus1_given_YT[[t]]) #
  }
  
  return(list(smoothed_marginal = p_st_given_YT, smoothed_joint = p_st_stminus1_given_YT))
}

B = msar_e_step(y, P_hat_curr, alpha_mu_curr, alpha_sigma_curr, alpha_phi1_curr, alpha_phi2_curr, M, T)

# ---- M-step ----
msar_m_step = function(y, p_st_stminus1_given_YT, alpha_mu, alpha_phi1, alpha_phi2, M, T) {
  
  # --- 1. Define Gamma (Smoothed Marginals) ---
  # gamma[t, k] corresponds to your gamma_t(k) logic
  gamma <- matrix(0, T, M)
  # t=1: Marginalize s1 from the joint at t=2
  gamma[1, ] <- rowSums(p_st_stminus1_given_YT[[2]])
  # t > 1: Marginalize s_t from the joint at t
  for (t in 2:T) {
    gamma[t, ] <- colSums(p_st_stminus1_given_YT[[t]])
  }
  
  # --- 2. Calculate Expected Lagged Means ---
  # These account for the mu_l and mu_m terms in your derivations
  exp_mu_lag = numeric(T)
  for(t in 2:T) {
    # Expected value of previous mean using marginal for s_{t-1}
    exp_mu_lag[t-1] = sum(alpha_mu * rowSums(p_st_stminus1_given_YT[[t]]))
  }
  
  # Initialize containers for updated parameters
  new_mu = numeric(M)
  new_phi1 = numeric(M)
  new_phi2 = numeric(M)
  new_sigma2 = numeric(M)
  
  for(k in 1:M) {
    
    # --- 3. Update Phi 1 & 2 (Iterative Conditional Logic) ---
    num_p1 = 0; den_p1 = 0
    num_p2 = 0; den_p2 = 0
    
    for(t in 3:T) {
      x1 = y[t-1] - exp_mu_lag[t-1]
      x2 = y[t-2] - exp_mu_lag[t-2]
      
      # Phi 1k: (y_t - mu_k - phi_2k * x2) / x1
      term_p1 = y[t] - alpha_mu[k] - ((k-1) * alpha_phi2[k] * x2)
      num_p1 = num_p1 + gamma[t, k] * x1 * term_p1
      den_p1 = den_p1 + gamma[t, k] * (x1^2)
      
      # Phi 2k (State 2 only): (y_t - mu_k - phi_1k * x1) / x2
      if(k == 2) {
        term_p2 = y[t] - alpha_mu[k] - (alpha_phi1[k] * x1)
        num_p2 = num_p2 + gamma[t, k] * x2 * term_p2
        den_p2 = den_p2 + gamma[t, k] * (x2^2)
      }
    }
    new_phi1[k] = num_p1 / den_p1
    new_phi2[k] = if(k == 2) num_p2 / den_p2 else 0
    
    # --- 4. Update Mu (Weighted Mean Adjustment) ---
    num_mu = 0; den_mu = 0
    for(t in 3:T) {
      ar_effect = new_phi1[k] * (y[t-1] - exp_mu_lag[t-1]) + 
        ((k-1) * new_phi2[k] * (y[t-2] - exp_mu_lag[t-2]))
      
      num_mu = num_mu + gamma[t, k] * (y[t] - ar_effect)
      den_mu = den_mu + gamma[t, k]
    }
    new_mu[k] = num_mu / den_mu
    
    # --- 5. Update Sigma^2 (Weighted Residual Variance) ---
    num_sig = 0
    for(t in 3:T) {
      eps = y[t] - new_mu[k] - 
        new_phi1[k] * (y[t-1] - exp_mu_lag[t-1]) - 
        ((k-1) * new_phi2[k] * (y[t-2] - exp_mu_lag[t-2]))
      num_sig = num_sig + gamma[t, k] * (eps^2)
    }
    new_sigma2[k] = num_sig / den_mu # den_mu is the sum of weights gamma_t(k)
  }
  
  # --- 6. Update Transition Matrix P ---
  new_P = matrix(0, M, M)
  for(i in 1:M) {
    # Denominator: P(s_{t-1} = i | Y_T)
    denom_p = 0
    for(t in 2:T) denom_p = denom_p + sum(p_st_stminus1_given_YT[[t]][i, ])
    
    for(j in 1:M) {
      # Numerator: P(s_t = j, s_{t-1} = i | Y_T)
      num_p = 0
      for(t in 2:T) num_p = num_p + p_st_stminus1_given_YT[[t]][i, j]
      new_P[i, j] = num_p / denom_p
    }
  }
  
  return(list(
    mu = new_mu,
    phi1 = new_phi1,
    phi2 = new_phi2,
    sigma = sqrt(new_sigma2),
    P = new_P
  ))
}

# 1. Run the M-step and save the output to 'm_out'
# Correct call: Pass the list of joint matrices, not the marginal matrix
m_out <- msar_m_step(y, B$smoothed_joint, alpha_mu, alpha_phi1, alpha_phi2, M, T)
# 2. Perform the three basic checks:

# A. Probability Check: Transition rows must sum to 1
print("--- Transition Matrix Row Sums (Check for 1.0) ---")
print(rowSums(m_out$P))

# B. Parameter Check: Sigma must be positive
print("--- Updated Sigmas ---")
print(m_out$sigma)

# C. Value Check: New means for each regime
print("--- Updated Mus ---")
print(m_out$mu)



# ---- EM Loop ----
max_iter = 3000
tol = 1e-6
eps = 1e-12
loglik_vec <- numeric(max_iter)
for (iter in 1:max_iter) {
  
  # 1. E-STEP: Get Smoothed Probabilities
  # This function must return 'smoothed_joint' and 'loglik'
  E_out <- msar_e_step(y, P_hat_curr, alpha_mu_curr, alpha_sigma_curr, alpha_phi1_curr, alpha_phi2_curr, M, T)
  
  # Store log-likelihood for convergence check
  loglik_vec[iter] <- E_out$loglik 
  
  # 2. M-STEP: Update Parameters
  # Pass the list of joint probabilities to the M-step function
  M_out <- msar_m_step(y, E_out$smoothed_joint, alpha_mu_curr, alpha_phi1_curr, alpha_phi2_curr, M, T)
  
  # 3. UPDATE PARAMETERS FOR NEXT ITERATION
  alpha_mu_curr    <- M_out$mu    #
  alpha_phi1_curr  <- M_out$phi1  #
  alpha_phi2_curr  <- M_out$phi2  #
  alpha_sigma_curr <- M_out$sigma #
  P_hat_curr       <- M_out$P     #
  
  # 4. CONVERGENCE CHECK
  if (iter > 1) {
    change <- loglik_vec[iter] - loglik_vec[iter-1]
    cat(sprintf("Iteration %d | Log-Likelihood: %.4f | Change: %.6f\n", iter, loglik_vec[iter], change))
    
    if (abs(change) < tol) {
      cat("EM Algorithm converged.\n")
      break
    }
  }
}

# Final Results
final_params <- list(mu=alpha_mu, phi1=alpha_phi1, phi2=alpha_phi2, sigma=alpha_sigma, P=P_hat)
























max_iter = 3000
tol = 1e-6
eps = 1e-12
loglik_history = numeric(max_iter)

for(iter in 1:max_iter){
  estep = msar_e_step(y, P_hat_curr, alpha_mu_curr, alpha_sigma_curr,
                      alpha_phi1_curr, alpha_phi2_curr, M, T)
  
  smoothed_marg  = estep$smoothed_marginals
  smoothed_joint = estep$smoothed
  f_y_cond_st    = estep$f_y_cond_st
  
  # Log-likelihood
  loglik = log(sum(estep$p_s1 * diag(f_y_cond_st[,,1])) + eps)
  for(t in 2:T){
    for(i in 1:M){
      for(j in 1:M){
        loglik = loglik + smoothed_joint[t,i,j] *
          (log(P_hat_curr[i,j] + eps) + log(f_y_cond_st[j,i,t] + eps))
      }
    }
  }
  loglik_history[iter] = loglik
  
  # M-step
  mstep = msar_m_step(y, smoothed_marg, smoothed_joint, alpha_mu_curr)
  alpha_mu_curr    = mstep$alpha_mu
  alpha_sigma_curr = mstep$alpha_sigma
  alpha_phi1_curr  = mstep$alpha_phi1
  alpha_phi2_curr  = mstep$alpha_phi2
  P_hat_curr       = mstep$P_hat
  
  # Convergence
  if(iter > 1 && abs(loglik_history[iter] - loglik_history[iter-1]) < tol){
    loglik_history = loglik_history[1:iter]
    cat("Converged at iteration", iter, "\n")
    break
  }
  
  if(iter %% 50 == 0){
    cat("Iter:", iter, "LogLik:", round(loglik,4),
        "Regime probs:", round(colMeans(smoothed_marg),3), "\n")
  }
}

# ---- Final results ----
alpha_mu_curr
alpha_sigma_curr
alpha_phi1_curr
alpha_phi2_curr
P_hat_curr
plot(loglik_history, type='b', main="Log-likelihood EM", xlab="Iteration", ylab="LogLik")
