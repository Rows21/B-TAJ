## ---- nimble simulation model 
# ---- NIMBLE Model Code ----
varying_coef_code_sim <- nimbleCode({
  
  sigma ~ dunif(0, 1)
  tau ~ T(dnorm(0, 1), 0, )  # Global horseshoe scale (half-Cauchy)
  
  beta_amy_age ~ dnorm(0, sd = 1)
  
  rho <- 0.3
  shape_weib <- 1.84 
  gamma_weib[1] <- 4.765
  gamma_weib[2] <- 0.062
  gamma_weib[3] <- -0.039

  gamma_amy[1] ~ dnorm(50, sd = 1)
  gamma_tau[1] ~ dnorm(55, sd = 1)
  
  for (p in 2:P) {
    gamma_amy[p] ~ dnorm(0, 1/10)
    gamma_tau[p] ~ dnorm(0, 1/10)
    
  }
    
    idx_bound[1:3] <- c(1, n_pred, n_pred)
    idx <- n_pred/2
    a_age <- n_pred/2
  
    
    r_pred_int[1:n_pred] <- exp(B_d[1:nd_pred, 1:K] %*% theta[1:K, 1])
    f_pred_int0[1:n_pred] <- cumulative_trapz(r_pred_int[1:n_pred], step)
    f_pred[1:n_pred] <- f_pred_int0[1:n_pred] - f_pred_int0[n_pred/2] + a_thresh
    
    for(n in 1:n_subjects_amy){
      
      xg[1, n] <- inprod(Xsub[n, 1:P], gamma_amy[1:P])
      xg[2, n] <- inprod(Xsub[n, 1:P], gamma_tau[1:P])
      #alpha_amy[n] ~ dt(df = 5, mu = xg_amy[n], sigma = sa)
      #alpha_amy[n] ~ dnorm(xg_amy[n], sd = sa)
      for(i in 1:n_pred){
        y_pred_amy[i, n] ~ dnorm(f_pred[i], sd = sy)
      }
    
      upper_range_a[n] <- (idx_bound[2] - range  - idx_bound[1] + 1) # tune this later to be subject-specific
      d_start_a[n] <- idx_bound[1] + rcat(1, prob = rep(1 /upper_range_a[n], upper_range_a[n])) - 1
    }
    
    for(n in 1:n_subjects_tau){
      
      #xg_tau[n] <- inprod(Xsub[tau_amy_link[n], 1:P], gamma_tau[1:P]) + beta_amy_age * (alpha_amy[tau_amy_link[n]] - xg_amy[tau_amy_link[n]]) / 10 # need a common index, a link
    
      #alpha_tau[n] ~ dt(df = 5, mu = xg_tau[n], sigma = sa)
      #cond_mu[n] <- xg_tau[n] + rho * sa * (alpha_amy[n] - xg_amy[n]) / sa
      #alpha_tau[n] ~ dnorm(cond_mu[n] , sd = sqrt((1 - rho^2) * sa ^ 2))

      for(i in 1:n_pred){
        y_pred_tau[i, n] ~ dnorm(f_pred[i], sd = sy)
      }
      upper_range_t[n] <- (idx_bound[2] - range  - idx_bound[1] + 1) # tune this later to be subject-specific
      d_start_t[n] <- idx_bound[1] + rcat(1, prob = rep(1 /upper_range_t[n], upper_range_t[n])) - 1
    }
  
    cov_mat[1, 1] <- sa ^ 2
    cov_mat[2, 2] <- sa ^ 2
    cov_mat[1, 2] <- rho * sa * sa
    cov_mat[2, 1] <- cov_mat[1, 2]
    
    for(n in 1:n_subjects_amy){
     
      alpha[1:2, n] ~ dmnorm(xg[1:2, n], cov = cov_mat[1:2, 1:2])  
       
    }
    
  for(i in 1:N_amy){
      
    d_amy[i] <- d_start_a[subject_amy[i]] + rcat(1, prob = rep(1/range, range)) - 1
    y_sparse_amy[i] <- y_pred_amy[d_amy[i], subject_amy[i]]
    t_amy[i] <- ((d_amy[i] - idx) * step) + alpha[1, subject_amy[i]]
    dd_amy[i] <- (t_amy[i] - alpha[1, subject_amy[i]])
  }
  
  for(i in 1:N_tau){ 
      
    d_tau[i] <- d_start_t[subject_tau_dense[i]] + rcat(1, prob = rep(1/range, range)) - 1
    y_sparse_tau[i] <- y_pred_tau[d_tau[i], subject_tau_dense[i]]
    t_tau[i] <- ((d_tau[i] - idx) * step) + alpha[2, subject_tau[i]]
    dd_tau[i] <- (t_tau[i] - alpha[2, subject_tau[i]])

    tt[i] ~ dmyweib3(shape_weib, lambda[i], alpha[2, subject_tau[i]])
    lambda[i] <- exp(gamma_weib[1] + Xsub[subject_tau[i], 2] * gamma_weib[2] +  alpha[2, subject_tau[i]] * gamma_weib[3])

    
  }
    
    
})

## ---- varying coef model compile
generate_vc <- function(varying_coef_code_sim,
                        data_1,
                        scenario = 1,
                        sb = NULL,
                        iter = NULL,
                        save_data = T,
                        root = NULL
){
  
  constants_simu <- data_1$constants
  K <- constants_simu$K
  P <- constants_simu$P
  n <- constants_simu$N_amy
  n_subjects <- constants_simu$n_subjects
  n_subjects_amy <- constants_simu$n_subjects_amy
  n_subjects_tau <- constants_simu$n_subjects_tau
  subject_id <- constants_simu$subject
  tau_amy_link <- constants_simu$tau_amy_link
  subject_tau <- constants_simu$subject_tau # subject_id[subject_id %in% tau_amy_link]
  
  # set up simulation model
  modelsim <- nimbleModel(varying_coef_code_sim, constants = constants_simu, check = FALSE, calculate = FALSE)
  
  # print(0)
  cModelsim <- compileNimble(modelsim)
  cModelsim$simulate()

  data_t <- data.frame(subject = subject_id, 
                       t_amy = cModelsim$t_amy,  
                       t_tau = NA)
  data_t$t_tau[subject_id %in% tau_amy_link] <- cModelsim$t_tau
  
  per_subject_amy <- data_t %>%
    group_by(subject) %>%
    summarise(t_min_amy = min(t_amy, na.rm = TRUE),
              t_max_amy = max(t_amy, na.rm = TRUE),
              mean_t_amy = mean(t_amy, na.rm = TRUE),
              .groups = "drop")
  
  per_subject_tau <- data_t %>% 
    drop_na(t_tau) %>%
    group_by(subject) %>%
    summarise(t_min_tau = min(t_tau, na.rm = TRUE),
              t_max_tau = max(t_tau, na.rm = TRUE),
              mean_t_tau = mean(t_tau, na.rm = TRUE),
              .groups = "drop")

  # data processing for simulation
  data <- list(y_pred_amy = cModelsim$y_sparse_amy, 
               y_pred_tau = cModelsim$y_sparse_tau)
  data$constraint_alpha_amy <- rep(1, n_subjects)
  data$constraint_alpha_tau <- rep(1, length(tau_amy_link))#rep(1, n_subjects)

  
  constants <- list(
    N_amy = n,          # total number of observations
    N_tau = length(subject_tau),
    
    n_subjects_amy = n_subjects_amy, 
    n_subjects_tau = n_subjects_tau, 
  
    subject_amy = subject_id,  
    subject_tau = subject_tau,  
    
    Xsub_amy = constants_simu$Xsub,      
    Xsub_tau = constants_simu$Xsub[tau_amy_link,],
    
    amy_thresh = constants_simu$a_thresh, # transform from CL to standardized SUVR scale (linear transformation)
    tau_thresh = constants_simu$a_thresh,
    tau_amy_link = tau_amy_link,
    
    K = K,          # number of spline basis functions
    #Kd = constants_simu$Kd,

    P = P,          # covariate dimension
    #B_pred = constants_simu$B_pred, # spline basis for predictions
    B_d = constants_simu$B_d,
    #B_dc = constants_simu$B_dc,
    #d_pred = constants_simu$d_pred,
    #nd_pred = constants_simu$nd_pred,
    n_pred = constants_simu$n_pred, # length of prediction interval
    #nv_pred = length(constants_simu$v_pred),
    #v_pred = constants_simu$v_pred,        # dense value grid
    
    max_t_amy = per_subject_amy$t_max_amy,
    min_t_amy = per_subject_amy$t_min_amy,
    t_amy = cModelsim$t_amy, # age mimic
    
    max_t_tau = per_subject_tau$t_max_tau,#[tau_amy_link],
    min_t_tau = per_subject_tau$t_min_tau,#[tau_amy_link],
    t_tau = cModelsim$t_tau, # age mimic
  
    step = constants_simu$step,
    #target_threshold = constants_simu$target_threshold,
    #a_thresh = constants_simu$a_thresh,
    middle_index = round(constants_simu$n_pred/2)
    #mid_point = constants_simu$target_threshold,
    #vec_mid_index = which.min(abs(constants_simu$v_pred - constants_simu$target_threshold)),
    #B_cov = constants_simu$B_cov,
    #Kc = constants_simu$Kc,
    #phi = constants_simu$phi,
    #cov_grid = constants_simu$cov_grid,
    #nc_pred = constants_simu$nc_pred
  )
  
  alpha_inits <- matrix(NA, nrow = 2, ncol = length(per_subject_amy$mean_t_amy))
  alpha_inits[1, ] <- per_subject_amy$mean_t_amy
  alpha_inits[2, ] <- NA
  alpha_inits[2, tau_amy_link] <-  per_subject_tau$mean_t_tau
  alpha_inits[2, -tau_amy_link] <-  per_subject_amy$mean_t_amy[-tau_amy_link]
    
  inits <- list(theta_amy = seq(-9, -3, length = K), 
                theta_tau = seq(-9, -3, length = K),
                alpha = alpha_inits) #
  return(list(data = data, 
              constants = constants, 
              inits = inits, 
              cModelsim = cModelsim))
  
}
