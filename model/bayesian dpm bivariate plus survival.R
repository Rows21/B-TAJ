varying_coef_code <- nimbleCode({
  
  # Coefficients
  gamma_amy[1] ~ dnorm(50, sd = 10)
  gamma_tau[1] ~ dnorm(50, sd = 10)
  
  for(p in 2:P){

    gamma_amy[p] ~ dnorm(0, 1/10)
    gamma_tau[p] ~ dnorm(0, 1/10)
    
    gamma_d[p] <- gamma_tau[p] - gamma_amy[p]

  }
  
  
  for(p in 1:3){

    gamma_weib[p] ~ dnorm(0, 1/100)
    
  }
  
  
    rho ~ dunif(-1, 1)
    sigma_pen[1] ~ dgamma(0.001, 0.001)
    sigma_pen[2] ~ dgamma(0.001, 0.001)
    
    theta_amy[1] ~ dnorm(-9, sd = .5)
    theta_tau[1] ~ dnorm(-9, sd = .5)
    
    shape_weib ~ dhalfflat()
    
    for(k in 2:K){
      
      lambda_1[k] ~ dexp(1)
      lambda_2[k] ~ dexp(1)
      
      theta_amy[k] ~ dnorm(theta_amy[k - 1], sigma_pen[1]/lambda_1[k])
      theta_tau[k] ~ dnorm(theta_tau[k - 1], sigma_pen[2]/lambda_2[k])
    }
    

    sigma_y_amy ~ T(dnorm(0, 1), 0, )
    sigma_y_tau ~ T(dnorm(0, 1), 0, )
        
    sigma_a_amy ~ T(dnorm(5, sd = 2), 0, )
    sigma_a_tau ~ T(dnorm(5, sd = 2), 0, )
    

    idx_bound[1:3] <- c(1, n_pred, n_pred)
    idx <- middle_index
    a_age <- middle_index
  
    
  
    r_pred_amy[1:n_pred] <- exp(B_d[1:n_pred, 1:K] %*% theta_amy[1:K])
    f_pred_amy0[1:n_pred] <- cumulative_trapz(r_pred_amy[1:n_pred], step)
    f_pred_amy[1:n_pred] <- f_pred_amy0[1:n_pred] - f_pred_amy0[middle_index] + amy_thresh
    
    r_pred_tau[1:n_pred] <- exp(B_d[1:n_pred, 1:K] %*% theta_tau[1:K])
    f_pred_tau0[1:n_pred] <- cumulative_trapz(r_pred_tau[1:n_pred], step)
    f_pred_tau[1:n_pred] <- f_pred_tau0[1:n_pred] - f_pred_tau0[middle_index] + tau_thresh
    
    
    cov_mat[1, 1] <- sigma_a_amy ^ 2
    cov_mat[2, 2] <- sigma_a_tau ^ 2
    cov_mat[1, 2] <- rho * sigma_a_amy * sigma_a_tau
    cov_mat[2, 1] <- cov_mat[1, 2]
    
    # L[1, 1] <- sigma_a_amy
    # L[1, 2] <- 0
    # L[2, 1] <- rho * sigma_a_tau
    # L[2, 2] <- sigma_a_tau * sqrt(1 - rho^2)
    
  for(n in 1:n_subjects_amy){
    
    xg[1, n] <- inprod(Xsub_amy[n, 1:P], gamma_amy[1:P])
    xg[2, n] <- inprod(Xsub_amy[n, 1:P], gamma_tau[1:P])
    
    # z[1, n] ~ dnorm(0, 1)
    # z[2, n] ~ dnorm(0, 1)
    
    #alpha[1:2, n] <- xg[1:2, n] + L[1:2, 1:2] %*% z[1:2, n]
    
    alpha[1:2, n] ~ dmnorm(xg[1:2, n], cov = cov_mat[1:2, 1:2])  
    
    constraint_alpha_amy[n] ~ dconstraint(alpha[1, n] > (max_t_amy[n] - (idx_bound[2] - middle_index)  * step) & alpha[1, n] < (min_t_amy[n] + (middle_index - idx_bound[1]) * step))
    constraint_alpha_tau[n] ~ dconstraint(alpha[2, n] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha[2, n] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
    #constraint_alpha_tau[n] ~ dconstraint(alpha[2, tau_amy_link[n]] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha[2, tau_amy_link[n]] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
  }
    
  # for(n in 1:n_subjects_tau){
  #   alpha_tau[n] <- alpha[2, tau_amy_link[n]]
  #   constraint_alpha_tau[n] ~ dconstraint(alpha_tau[n] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha_tau[n] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
  #   #constraint_alpha_tau[n] ~ dconstraint(alpha[2, tau_amy_link[n]] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha[2, tau_amy_link[n]] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
  # 
  # }
    
  #   
  # for(n in 1:n_subjects_tau){
  #     
  #   #xg_tau[n] <- inprod(Xsub_tau[n, 1:P], gamma_tau[1:P]) + beta_amy_age * (alpha_amy[tau_amy_link[n]] - xg_amy[tau_amy_link[n]]) / 10 # need a common index, a link
  #   xg_tau[n] <- inprod(Xsub_tau[n, 1:P], gamma_tau[1:P])
  #   #alpha_tau[n] ~ T(dt(df = 5, mu = xg_tau[n], sigma = sigma_a_tau), -Inf, alpha_tau_pos[n])
  #   #alpha_tau[n] ~ dt(df = 5, mu = xg_tau[n], sigma = sigma_a_tau)
  #   
  #   cond_mu[n] <- xg_tau[n] + rho * sigma_a_tau * (alpha_amy[n] - xg_amy[n]) / sigma_a_amy
  #   alpha_tau[n] ~ dnorm(cond_mu[n] , sd = sqrt((1 - rho^2) * sigma_a_tau ^ 2))
  #   
  #   constraint_alpha_tau[n] ~ dconstraint(alpha_tau[n] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha_tau[n] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
  #   alpha_d[n] <- alpha_tau[n] - alpha_amy[tau_amy_link[n]]
  #   
  #   }
  
    
    for(n in 1:n_subjects_surv){

      censored[n] ~ dinterval(tt[n], c[n])
      #t[n] ~ dweib(shape_weib, lambda[n])
      tt[n] ~ dmyweib3(shape_weib, lambda[n], alpha[2, tau_surv_link[n]])
      #t[n] ~ dmyweib3(shape_weib, lambda[n], 0)
      #lambda[n] <- exp(inprod(Xsub_surv[n, 1:2], gamma_weib[1:2]))
      #lambda[n] <- exp(inprod(Xsub_surv[n, 1:P], gamma_weib[1:P]) +  alpha_tau[tau_surv_link[n]] * gamma_weib[(P + 1)])
      lambda[n] <- exp(gamma_weib[1] + Xsub_surv[n, 2] * gamma_weib[2] +  alpha[2, tau_surv_link[n]] * gamma_weib[3])
    }

  
    
  for(i in 1:N_amy){
    
    # age since onset
    dd_amy[i] <- (t_amy[i] - alpha[1, subject_amy[i]])  
    row_idx_amy[i] <- (dd_amy[i]/step + idx) 
    safe_idx_amy[i] <- round(min(max(idx_bound[1], row_idx_amy[i]), idx_bound[2]))

    # # No covariate
    f_val_amy[i] <- f_pred_amy[safe_idx_amy[i]]
    y_pred_amy[i] ~ dnorm(f_val_amy[i], sd = sigma_y_amy)
    
  }
    
  for(i in 1:N_tau){
      
    # age since onset
    dd_tau[i] <- (t_tau[i] - alpha[2, subject_tau[i]])  
    row_idx_tau[i] <- (dd_tau[i]/step + idx) 
    safe_idx_tau[i] <- round(min(max(idx_bound[1], row_idx_tau[i]), idx_bound[2]))
    
    # # No covariate
    f_val_tau[i] <- f_pred_tau[safe_idx_tau[i]]
    y_pred_tau[i] ~ dnorm(f_val_tau[i], sd = sigma_y_tau)
      
  }
    
    
})



