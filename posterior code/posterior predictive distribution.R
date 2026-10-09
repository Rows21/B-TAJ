### Extract posterior distributions from fitted model
#source('Real Data Analysis/real data mcmc_posterior.R')


# ### As a reminder, we need to pull in outcomes data and other constants for unobserved data
# 
# simudata$data$y_pred_amy
# simudata$constants$N_amy
# 
# # MCMC samples
# f_pred_amy <- (samples[,grep("f_pred_amy", colnames(samples))]) 
# f_pred_tau <- (samples[,grep("f_pred_tau", colnames(samples))]) 
# 
# gamma_amy <- samples[, grepl("^gamma_amy(\\[|$|_)", colnames(samples)), drop = FALSE]
# gamma_tau <- (samples[, grepl("^gamma_tau(\\[|$|_)", colnames(samples)), drop = FALSE]) 
# 
# sigma_a_amy <- (samples[, grepl("sigma_a_amy", colnames(samples)), drop = FALSE])
# sigma_a_tau <- (samples[, grepl("sigma_a_tau", colnames(samples)), drop = FALSE])
# 
# sigma_y_amy <- (samples[, grepl("sigma_y_amy", colnames(samples)), drop = FALSE])
# sigma_y_tau <- (samples[, grepl("sigma_y_tau", colnames(samples)), drop = FALSE])
# 
# rho <- (samples[, grepl("rho", colnames(samples)), drop = FALSE])
# 
# simudata$data <- c(
#   simudata$data,
#   list(
#     f_pred_amy = (f_pred_amy),
#     f_pred_tau = (f_pred_tau),
#     
#     gamma_amy = (gamma_amy),
#     gamma_tau = (gamma_tau),
#     
#     sigma_a_amy = c(sigma_a_amy),
#     sigma_a_tau = c(sigma_a_tau),
#     
#     sigma_y_amy = c(sigma_y_amy),
#     sigma_y_tau = c(sigma_y_tau),
#     rho = c(rho)
#   )
# )
# 
# simudata$constants$n_post <- simudata$constants$n_subjects_amy
# simudata$constants$n_mcmc <- length(rho)

posterior_predictions_code <- nimbleCode({
  
  idx_bound[1:3] <- c(1, n_pred, n_pred)
  idx <- middle_index
  a_age <- middle_index
  
  u ~ dunif(0, 1)
  
  gamma_amy_post[1:P] <- random_row(gamma_amy[1:n_mcmc, 1:P]) + u - u 
  gamma_tau_post[1:P] <- random_row(gamma_tau[1:n_mcmc, 1:P]) + u - u
  
  f_pred_amy_post[1:201] <- random_row(f_pred_amy[1:n_mcmc, 1:201]) + u - u
  f_pred_tau_post[1:201] <- random_row(f_pred_tau[1:n_mcmc, 1:201]) + u - u
    
  sigma_a_amy_post <- random_scalar(sigma_a_amy[1:n_mcmc]) + u - u
  sigma_a_tau_post <- random_scalar(sigma_a_tau[1:n_mcmc]) + u - u
  
  sigma_y_amy_post <- random_scalar(sigma_y_amy[1:n_mcmc]) + u - u
  sigma_y_tau_post <- random_scalar(sigma_y_tau[1:n_mcmc]) + u - u
  
  rho_post <- random_scalar(rho[1:n_mcmc]) + u - u
  
  
  cov_mat[1, 1] <- sigma_a_amy_post ^ 2
  cov_mat[2, 2] <- sigma_a_tau_post ^ 2
  cov_mat[1, 2] <- rho_post * sigma_a_amy_post * sigma_a_tau_post
  cov_mat[2, 1] <- cov_mat[1, 2]
  
  for(n in 1:n_post){
    
    xg[1, n] <- inprod(Xsub_amy[n, 1:P], gamma_amy_post[1:P])
    xg[2, n] <- inprod(Xsub_amy[n, 1:P], gamma_tau_post[1:P])
    
    alpha[1:2, n] ~ dmnorm(xg[1:2, n], cov = cov_mat[1:2, 1:2])  
    
    constraint_alpha_amy[n] ~ dconstraint(alpha[1, n] > (max_t_amy[n] - (idx_bound[2] - middle_index)  * step) & alpha[1, n] < (min_t_amy[n] + (middle_index - idx_bound[1]) * step))
    constraint_alpha_tau[n] ~ dconstraint(alpha[2, n] > (max_t_tau[n] - (idx_bound[2] - middle_index)  * step) & alpha[2, n] < (min_t_tau[n] + (middle_index - idx_bound[1]) * step))
  }
  
  
  for(i in 1:N_amy){
    
    # age since onset
    dd_amy[i] <- (t_amy[i] - alpha[1, subject_amy[i]])  
    row_idx_amy[i] <- (dd_amy[i]/step + idx) 
    safe_idx_amy[i] <- round(min(max(idx_bound[1], row_idx_amy[i]), idx_bound[2]))
    
    # # No covariate
    f_val_amy[i] <- f_pred_amy_post[safe_idx_amy[i]]
    y_pred_amy[i] ~ dnorm(f_val_amy[i], sd = sigma_y_amy_post)
    
  }
  
  for(i in 1:N_tau){
    
    # age since onset
    dd_tau[i] <- (t_tau[i] - alpha[2, tau_amy_link[subject_tau[i]]])  
    row_idx_tau[i] <- (dd_tau[i]/step + idx) 
    safe_idx_tau[i] <- round(min(max(idx_bound[1], row_idx_tau[i]), idx_bound[2]))
    
    # # No covariate
    f_val_tau[i] <- f_pred_tau_post[safe_idx_tau[i]]
    y_pred_tau[i] ~ dnorm(f_val_tau[i], sd = sigma_y_tau_post)
    
  }
  
  
})

random_row <- nimbleFunction(
  run = function(x = double(2)) {
    returnType(double(1))
    
    nrow_x <- dim(x)[1]
    ncol_x <- dim(x)[2]
    
    i <- floor(runif(1, 0, 1) * nrow_x) + 1
    
    out <- numeric(ncol_x, init = FALSE)
    for(j in 1:ncol_x) {
      out[j] <- x[i, j]
    }
    
    return(out)
  }
)

random_scalar <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(0))
    n <- length(x)
    j <- floor(runif(1) * n) + 1
    return(x[j])
  }
)

always_accept_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Extract the node name
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    # 1. Propose a new value (e.g., using a normal distribution)
    current_val <- model[[target]]
    proposal <- rnorm(1, mean = current_val, sd = 1) 
    
    # 2. Force the update into the model
    model[[target]] <<- proposal
    
    # 3. Synchronize the model and dependent nodes
    calculate(model, calcNodes)
    
    # 4. Copy the new state into mvSaved (the MCMC storage)
    # Since we ALWAYS accept, we always copy the current model state
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)
    
    
    




