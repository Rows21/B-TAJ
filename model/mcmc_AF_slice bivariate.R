#library(nimbleHMC)

mcchains <- function(
    varying_coef_code,
    simudata,
    monitors = c("theta", "beta_pred", "r_pred", "f_pred", "y_pred", "row_idx", "idx","safe_idx", "idx_bound",
                 "dd", "gamma_r", "alpha_i", "xg", "f_val", "sigma_a", "sigma_y"),
    niter = 200,
    nburnin = 100, 
    thin = 5,
    nchains = 1,
    fullres = TRUE){
  
  s1 <- Sys.time()
  model <- nimbleModel(varying_coef_code, 
                       constants = simudata$constants, 
                       data = simudata$data, 
                       check = FALSE, calculate = FALSE, 
                       inits = simudata$inits)
  #model$initializeInfo()
  ## ---- mcmc ----
  cModel <- compileNimble(model)
  conf <- configureMCMC(model)
  conf$addMonitors(monitors)
  # conf$removeSampler("theta_amy") 
  # conf$removeSampler("theta_tau") 
  # conf$removeSampler("gamma_amy") 
  # conf$removeSampler("gamma_tau") 
  # conf$removeSampler("beta_amy_age") 
  # 
  # conf$addSampler(target = c("gamma_amy[1:3]", "gamma_tau[1:3]", "beta_amy_age"), type = 'AF_slice')
  # conf$addSampler(target = c("theta_amy[1:12]"), type = 'AF_slice')
  # conf$addSampler(target = c("theta_tau[1:12]"), type = 'AF_slice')
  conf$removeSampler("theta_amy")
  conf$removeSampler("theta_tau")
  conf$removeSampler("gamma_amy")
  conf$removeSampler("gamma_tau")
  #conf$removeSampler("beta_amy_age")
  #conf$removeSampler("shape_weib")
  conf$removeSampler("gamma_weib")
  conf$removeSampler("alpha")
  conf$removeSampler("rho")
  conf$removeSampler("sigma_a_tau")
  conf$removeSampler("sigma_a_amy")
  
  conf$addSampler(target = c("gamma_amy[1:3]", "gamma_tau[1:3]", 'rho'), type = 'AF_slice')
  conf$addSampler(target = c("sigma_a_amy", "sigma_a_tau"), type = 'AF_slice')
  #conf$addSampler(target = c("gamma_amy[1:3]"), type = 'AF_slice')
  #conf$addSampler(target = c("gamma_tau[1:3]"), type = 'AF_slice')
  conf$addSampler(target = c("theta_amy[1:10]"), type = 'AF_slice')
  conf$addSampler(target = c("theta_tau[1:10]"), type = 'AF_slice')
  #conf$addSampler(target = c("gamma_weib[1:3]", "shape_weib"), type = 'AF_slice')
  conf$addSampler(target = c("gamma_weib[1:3]"), type = 'AF_slice')

for(n in 1:1000) {
  conf$addSampler(target = paste0('alpha[1:2, ', n, ']'),
                  type = 'AF_slice')
}


  mcmc <- buildMCMC(conf)
  cMCMC <- compileNimble(mcmc, project = model)
  samples <- runMCMC(cMCMC, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)
  s2 <- Sys.time()
  #samples <- runMCMC(cMCMC, niter = 500, nburnin = 0, thin = 5, nchains = 4)
  ## ---- storage ----
  if(fullres == T){
    return(list(samples = samples, 
                time = as.numeric(s2-s1),
                mcmc = mcmc,
                conf = conf,
                cMCMC = cMCMC
                ))
  }else{
    conf$addMonitors(c("theta", "beta_pred","r_pred", "f_pred", "y_pred","sigma_y", "f_val", "bp", "lambda"))
    
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
    # return()
  }
}

