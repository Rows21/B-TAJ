
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
  
  ## ---- mcmc ----
  cModel <- compileNimble(model)
  conf <- configureMCMC(model)
  conf$addMonitors(monitors)
  conf$removeSampler("u")
  conf$addSampler(target = "u", type = "always_accept_sampler")
  conf$removeSampler("alpha")
  
  for(n in 1:simudata$constants$n_post) {
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