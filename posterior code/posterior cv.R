# Define initial values for parameters
library(nimble)
library(matrixStats)

setwd('M:/Bayesian DPM Model/Bayesian DPM Model/')
setwd("E:/Github/B-TAJ/posterior code/")
source("real data bivariate plus survival data cv.R")
source("Real Data Analysis/Data_preprocess_amy.R")
source("Real Data Analysis/Data_preprocess_tau.R")
source("Real Data Analysis/real data bivariate plus survival data.R")
source('real data mcmc_posterior.R')
source('posterior predictive distribution.R')
set.seed(123)


id <- 1:1900
fold <- sample(rep(1:5, length.out = length(id)))
fold_list <- split(id, fold)
## Tutorial
## pred
test_id  <- fold_list[[5]]
train_id <- sort(as.numeric(unlist(fold_list[-5])))
load(file = "M:/Bayesian DPM Model/Bayesian DPM Model/Real Data Analysis/real_data_cv5.Rdata")
load(file = "M:/Bayesian DPM Model/Bayesian DPM Model/Real Data Analysis/model_out5.Rdata")
pred_listf <- realdata_pred(X_amy, X_tau, test_id, simudata$constants$tau_amy_link, dta_surv_df, res$samples, simudata_cv, biomarker = 'amy', type = 'forward')
posterior_res_forward <- mcchains(posterior_predictions_code,
                                  pred_listf$simudata_cv,
                                  monitors = c("alpha", "gamma_amy_post", "rho_post"),
                                  niter = 5000,
                                  nburnin = 3000,
                                  thin = 1,
                                  nchains = 1,
                                  fullres = TRUE) # fullres fixed

pred_listb <- realdata_pred(X_amy, X_tau, test_id, simudata$constants$tau_amy_link, dta_surv_df, res$samples, simudata_cv, biomarker = 'amy', type = 'backward')
posterior_res_backward <- mcchains(posterior_predictions_code,
                                   pred_listb$simudata_cv,
                                   monitors = c("alpha", "gamma_amy_post", "rho_post"),
                                   niter = 5000,
                                   nburnin = 3000,
                                   thin = 1,
                                   nchains = 1,
                                   fullres = TRUE) # fullres fixed


post_alpha_amyf <- posterior_res_forward$samples[,grep("alpha", colnames(posterior_res_forward$samples))][, seq(1, pred_listf$simudata_pred$constants$n_subjects_amy*2, by = 2)]
post_alpha_amyb <- posterior_res_backward$samples[,grep("alpha", colnames(posterior_res_backward$samples))][, seq(1, pred_listb$simudata_pred$constants$n_subjects_amy*2, by = 2)]
plot(colMedians(post_alpha_amyf), colMedians(post_alpha_amyb)) 

pred_listf <- realdata_pred(X_amy, X_tau, test_id, simudata$constants$tau_amy_link, dta_surv_df, res$samples, simudata_cv, biomarker = 'tau', type = 'forward')
posterior_res_forward <- mcchains(posterior_predictions_code,
                                  pred_listf$simudata_cv,
                                  monitors = c("alpha", "gamma_amy_post", "rho_post"),
                                  niter = 5000,
                                  nburnin = 3000,
                                  thin = 1,
                                  nchains = 1,
                                  fullres = TRUE) # fullres fixed

pred_listb <- realdata_pred(X_amy, X_tau, test_id, simudata$constants$tau_amy_link, dta_surv_df, res$samples, simudata_cv, biomarker = 'tau', type = 'backward')
posterior_res_backward <- mcchains(posterior_predictions_code,
                                   pred_listb$simudata_cv,
                                   monitors = c("alpha", "gamma_amy_post", "rho_post"),
                                   niter = 5000,
                                   nburnin = 3000,
                                   thin = 1,
                                   nchains = 1,
                                   fullres = TRUE) # fullres fixed


post_alpha_tauf <- posterior_res_forward$samples[,grep("alpha", colnames(posterior_res_forward$samples))][, seq(2, pred_listf$simudata_pred$constants$n_subjects_amy*2, by = 2)][,pred_listf$simudata_pred$constants$tau_amy_link]
post_alpha_taub <- posterior_res_backward$samples[,grep("alpha", colnames(posterior_res_backward$samples))][, seq(2, pred_listb$simudata_pred$constants$n_subjects_amy*2, by = 2)][,pred_listb$simudata_pred$constants$tau_amy_link]
plot(colMedians(post_alpha_taub), colMedians(post_alpha_tauf)) 


## 5 fold cv 
for(cv in 1:5){
  test_id  <- fold_list[[cv]]
  train_id <- sort(as.numeric(unlist(fold_list[-cv])))
  
  simudata_cv <- realdata_cv(X_amy, X_tau, train_id, simudata$constants$tau_amy_link, dta_surv_df)
  save(simudata_cv, file = paste0("E:/Github/data/real_data_cv",cv,".Rdata"))
  # Constraints
  max_t_tau <- rep(NA, simudata_cv$constants$n_subjects_amy)
  min_t_tau <- rep(NA, simudata_cv$constants$n_subjects_amy)
  
  max_t_tau[simudata_cv$constants$tau_amy_link] <- simudata_cv$constants$max_t_tau
  min_t_tau[simudata_cv$constants$tau_amy_link] <- simudata_cv$constants$min_t_tau
  
  max_t_tau[-simudata_cv$constants$tau_amy_link] <- 0
  min_t_tau[-simudata_cv$constants$tau_amy_link] <- 140
  
  simudata_cv$constants$max_t_tau <- max_t_tau
  simudata_cv$constants$min_t_tau <- min_t_tau
  
  simudata_cv$data$constraint_alpha_tau <- simudata_cv$data$constraint_alpha_amy
  
  # Alpha inits
  simudata_cv$inits$alpha <- matrix(NA, nrow = 2, ncol = simudata_cv$constants$n_subjects_amy)
  simudata_cv$inits$alpha[1, ] <- simudata_cv$inits$alpha_amy
  simudata_cv$inits$alpha[2, simudata_cv$constants$tau_amy_link] <- simudata_cv$inits$alpha_tau
  simudata_cv$inits$alpha[2, -simudata_cv$constants$tau_amy_link] <- simudata_cv$inits$alpha_amy[-simudata_cv$constants$tau_amy_link]
  
  # Survival inits
  # Needed?
  
  # Run model
  res <- mcchains(varying_coef_code,
                  simudata_cv,
                  monitors = c("theta_amy", "f_val_amy", "f_pred_amy", "r_pred_amy", "row_idx_amy", "gamma_amy", "alpha",
                               "theta_tau", "f_val_tau", "f_pred_tau","r_pred_tau", "row_idx_tau", "gamma_tau", "xg", 
                               "gamma_d", "gamma_weib", "shape_weib", "sigma_y_amy", "sigma_a_amy", "sigma_y_tau", "sigma_a_tau", "rho"),
                  niter = 5000,
                  nburnin = 3000,
                  thin = 1,
                  nchains = 1,
                  fullres = TRUE) # fullres fixed
  save(res, file = paste0("E:/Github/data/model_out",cv,".Rdata"))
}


simudata_cv <- realdata_cv(X_amy, X_tau, train_id, simudata$constants$tau_amy_link, dta_surv_df)




post_gamma_amy <- posterior_res$samples[,grep("gamma_amy_post", colnames(posterior_res$samples))]


plot_overlap_density <- function(alpha_amy, post_alpha_amy, j,
                                 fill = FALSE) {
  
  x1 <- alpha_amy[, j]
  x2 <- post_alpha_amy[, j]
  
  d1 <- density(x1, na.rm = TRUE)
  d2 <- density(x2, na.rm = TRUE)
  
  # Common x and y limits so both densities are fully visible
  xlim <- range(c(d1$x, d2$x))
  ylim <- range(c(d1$y, d2$y))
  
  plot(d1,
       xlim = xlim,
       ylim = ylim,
       lwd = 2,
       col = "blue",
       main = paste("Density comparison (column", j, ")"),
       xlab = "Value")
  
  if(fill) {
    polygon(d1, col = rgb(0, 0, 1, 0.3), border = NA)
    polygon(d2, col = rgb(1, 0, 0, 0.3), border = NA)
  }
  
  lines(d1, lwd = 2, col = "blue")
  lines(d2, lwd = 2, col = "red", lty = 2)
  
  legend("topright",
         legend = c("alpha_amy", "post_alpha_amy"),
         col = c("blue", "red"),
         lwd = 2,
         lty = c(1, 2),
         bty = "n")
}


plot(alpha_amy, colMedians(post_alpha_amy))


plot_overlap_density(gamma_amy, post_gamma_amy, j = 1)

# with filled densities
plot_overlap_density(alpha_amy, post_alpha_amy, j = 115, fill = TRUE)
