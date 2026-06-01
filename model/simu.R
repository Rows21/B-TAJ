# path <- c("/Users/aaron.scheffler/Downloads/Bayesian DPM Model Simulation new")
# path <- c("M:/Bayesian DPM Model/Bayesian DPM Model/Bayesian DPM Model Simulation new")
# setwd(path)
# library(matrixStats)
# source('data_gen_new.R')
# source('dweibull3.R')
# source('data_compile bivariate.R')
# source('bayesian dpm bivariate plus survival.R')
# source('mcmc_AF_slice bivariate.R')
# source('cumulative_trapz.R')
# source('simu.R')

simu <- function(varying_coef_code_sim,
                       varying_coef_code,
                       outer_iter = 1,
                       # data generate control
                       scenario = 1,
                       n_subjects = 300,
                       sa = 2.5,
                       sy = 0.07,
                       # mcmc control
                       niter = 1000,
                       nburnin = 800,
                       nchains = 1,
                       # saving
                       save_res = TRUE,
                       fullres = TRUE,
                       seeds = NULL){
  
  ## control parameters
  # monitors = c("theta_amy", "alpha_amy", "f_val_amy", "f_pred_amy", "r_pred_amy", "row_idx_amy", "gamma_amy", "alpha_amy", "xg_amy", "sigma_a_amy", "sigma_y_amy",
  #              "theta_tau", "alpha_tau", "f_val_tau", "f_pred_tau","r_pred_tau", "row_idx_tau", "gamma_tau", "alpha_tau", "xg_tau", "beta_amy_age",
  #              "gamma_d", "alpha_d", "gamma_weib", "shape_weib", "lambda", "tt")
  monitors = c("alpha", "xg",
               "f_pred_amy", "f_pred_tau", "r_pred_amy", "r_pred_tau",
               "gamma_amy", "gamma_tau", "gamma_weib", "shape_weib",
               "sigma_a_amy", "sigma_a_tau", 
               "sigma_y_amy", "sigma_y_tau", "rho")
  # hard set P
  # if(scenario == 1){P = 4}else{P = 3}
  P <- 3
  seed_list <- rep(NA, outer_iter)
  for(i in 1:outer_iter){
    
    #seed <- 59724504

    seed <- sample.int(.Machine$integer.max,1)
    
    seed_list[i] <- seed
    
    set.seed(seed)
    data_1 <- generate_d1(n_subjects = n_subjects,
                          scenario = 1,
                          K = 8, 
                          Kd = 8,
                          sa = sa,
                          step = 0.25,
                          sy = sy,
                          start_index = 1,
                          end_index = 250,
                          target_threshold = 1, 
                          mid_point = 1, 
                          a_thresh = 1,
                          C = NULL)
    
    simudata <- generate_vc(varying_coef_code_sim, data_1, 
                            scenario = scenario, #iter = i,
                            save_data = save_data, root = root)
    
    # construct survival 
    subject_id <- data_1$constants$subject_tau
    subject_tau <- simudata$constants$subject_tau
    cModelsim <- simudata$cModelsim
    
    data_surv <- data.frame(subject_id, cModelsim$y_sparse_tau, cModelsim$tt) |> 
      filter(cModelsim.y_sparse_tau >= 1, subject_id %in% subject_tau) |> #filter(cModelsim.t > 5) |> 
      group_by(subject_id) |>
      slice_min(cModelsim.tt, n=1, with_ties = F) |>
      ungroup()

    data_surv$censored <- rbinom(nrow(data_surv),1, prob=0.3)
    data_surv$thresh <- 0
    data_surv$thresh[which(data_surv$cModelsim.y_sparse_tau >= 1)] <- 1
    data_surv$t <- data_surv$cModelsim.tt
    data_surv$c <- data_surv$cModelsim.tt
    data_surv$t[data_surv$censored == 1] <- NA
    data_surv$c[data_surv$censored == 0] <- Inf
    
    
    n_subjects <- data_1$constants$n_subjects
    alpha_tau_pos <- rep(Inf, n_subjects)
    alpha_tau_pos[which(1:n_subjects %in% data_surv$subject_id)] <- data_surv$cModelsim.tt
    simudata$constants$alpha_tau_pos = alpha_tau_pos[simudata$constants$tau_amy_link]
    
    simudata$constants$n_subjects_surv <- length(data_surv$subject_id)
    simudata$constants$Xsub_surv <- simudata$constants$Xsub_amy[data_surv$subject_id,]
    simudata$constants$tau_surv_link <- data_surv$subject_id
    simudata$constants$c <- data_surv$c
    
    tinit <- rep(NA, length(data_surv$subject_id))
    tinit[data_surv$censored == 1] <- data_surv$cModelsim.tt[data_surv$censored == 1] + 1
    
    simudata$inits$tt <- tinit
    simudata$inits$shape_weib <- 1
    simudata$inits$gamma_weib <- c(log(median(data_surv$cModelsim.tt, na.rm = TRUE)), rep(0, 2))
    
    simudata$inits$gamma_amy <- c(50, 0, 0)
    simudata$inits$gamma_tau <- c(55, 0, 0)
    
    # Need to truncate to respect the constraints of the integral
    simudata$inits$alpha[2, data_surv$subject_id] <- pmin(
      data_surv$cModelsim.tt - 1, 
      simudata$inits$alpha[2, data_surv$subject_id]
    )

    # Initialize tau to be same dimension as amyloid 
    
    simudata$data$censored <- data_surv$censored
    simudata$data$tt <- data_surv$t

    max_t_tau <- rep(NA, simudata$constants$n_subjects_amy)
    min_t_tau <- rep(NA, simudata$constants$n_subjects_amy)
    
    max_t_tau[simudata$constants$tau_amy_link] <- simudata$constants$max_t_tau
    min_t_tau[simudata$constants$tau_amy_link] <- simudata$constants$min_t_tau
    
    max_t_tau[-simudata$constants$tau_amy_link] <- 0
    min_t_tau[-simudata$constants$tau_amy_link] <- 120
    
    simudata$constants$max_t_tau <- max_t_tau
    simudata$constants$min_t_tau <- min_t_tau
    
    simudata$data$constraint_alpha_tau <- simudata$data$constraint_alpha_amy
    
    res <- mcchains(varying_coef_code,
                    simudata,
                    monitors = monitors,
                    niter = niter,
                    nburnin = nburnin, 
                    thin = 1, 
                    nchains = 1,
                    fullres = TRUE) # fullres fixed
    samples <- res$samples
    
    #save(simudata, file = paste0("M:/Bayesian DPM Model new/Bayesian DPM Model/Simulation/Data/simu_data",i,".Rdata"))
    #save(samples, file = paste0("M:/Bayesian DPM Model new/Bayesian DPM Model/Simulation/Sample/samples",i,".Rdata"))
    #return(list(samples = samples, simudata = simudata))
    
    # if(nchains > 1){
    #   return(samples)
    # }
    
    chains <- list()
    # save iteration output
    for(nm in monitors){
      theta_raw <- samples[,grep(nm, colnames(samples))]
        chains[[nm]] <- theta_raw
    }

    if(save_res == T){
      root_res <- file.path('results',paste0('s',scenario,'_',n_subjects,sa*10,sy*100),paste0('iter',i))
      # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
      if (!dir.exists("results")) dir.create("results")
      if (!dir.exists(file.path('results',paste0('s',scenario,'_',n_subjects,sa*10,sy*100)))) dir.create(file.path('results',paste0('s',scenario,'_',n_subjects,sa*10,sy*100)))
      # file.path('results',paste0('s',scenario,'_',n_subjects,sa,sy*100))
      if (!dir.exists(root_res)) dir.create(root_res)
      for(nm in monitors){
        write.csv(chains[[nm]], file.path(root_res, paste0(nm, ".csv")), row.names = F)
      }
    }

    if(i == 1){

        f_pred <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$nd_pred)
        f_pred[1,] <- cModelsim$f_pred
        
        alpha_amy <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$n_subjects)
        alpha_amy[1,] <- cModelsim$alpha[1,]
        alpha_tau <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$n_subjects/2)
        alpha_tau[1,] <- cModelsim$alpha[2,simudata$constants$tau_amy_link]
        
        gamma_amy <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$P)
        gamma_amy[1,] <- cModelsim$gamma_amy
        gamma_tau <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$P)
        gamma_tau[1,] <- cModelsim$gamma_tau
        
        #y_sparse_amy[1,] <- simudata$cModelsim$y_sparse_amy
        
        
        tau_amy_link <- matrix(NA, nrow = outer_iter, ncol = data_1$constants$n_subjects/2)
        tau_amy_link[1,] <- simudata$constants$tau_amy_link
        t <- rep(NA, outer_iter)
        t[1] <- res$time
        
        # y_pred_amy <- matrix(NA, nrow = outer_iter, ncol = 6 * data_1$constants$n_subjects)
        # y_pred_amy[1,1:length(cModelsim$y_sparse_amy)] <- cModelsim$y_sparse_amy
        # y_pred_amy_est <- matrix(NA, nrow = outer_iter, ncol = 6 * data_1$constants$n_subjects)
        # y_pred_amy_est[1,1:length(cModelsim$y_sparse_amy)] <- chains$f_val_amy
        # 
        # y_pred_tau <- matrix(NA, nrow = outer_iter, ncol = 6 * data_1$constants$n_subjects)
        # y_pred_tau[1,1:length(cModelsim$y_sparse_tau)] <- cModelsim$y_sparse_tau
        # y_pred_tau_est <- matrix(NA, nrow = outer_iter, ncol = 6 * data_1$constants$n_subjects)
        # y_pred_tau_est[1,1:length(cModelsim$y_sparse_tau)] <- chains$f_val_tau

    }else{

      f_pred[i,] <- cModelsim$f_pred

      alpha_amy[i,] <- cModelsim$alpha[1,]
      alpha_tau[i,] <- cModelsim$alpha[2,simudata$constants$tau_amy_link]

      gamma_amy[i,] <- cModelsim$gamma_amy
      gamma_tau[i,] <- cModelsim$gamma_tau
      
      tau_amy_link[i,] <- simudata$constants$tau_amy_link

      t[i] <- res$time
    }
  }
  results <- list(
    # samples = samples,
    # simudata = simudata,
    # cModelsim = cModelsim,
    seed_list = seed_list,
    tau_amy_link = tau_amy_link,
    
    f_pred = f_pred,
    
    alpha_amy = alpha_amy,
    alpha_tau = alpha_tau,
    
    gamma_amy = gamma_amy,
    gamma_tau = gamma_tau,
    
    

    t = t
  )
  
  if(save_res == T){
    root_gt <- file.path("results",paste0('s',scenario,'_',n_subjects,sa*10,sy*100),"gt")
    if (!dir.exists(root_gt)){
      dir.create(root_gt)
    }
    lapply(names(results[1:length(results)]), function(nm) {
      write.csv(as.data.frame(results[[nm]]),
                file = file.path(root_gt, paste0(nm, ".csv")),
                row.names = FALSE)
    })
  }
  if(save_res == 'T'){return(print("finish!"))}else{
    return(list(
      simudata = simudata,
      res = res
           ))
  }
  
}