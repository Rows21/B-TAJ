realdata_cv <- function(
    X_amy,
    X_tau,
    train_id,
    tau_amy_link,
    dta_surv_df,
    step = 1 / 4,
    start_index = 1,
    end_index = 201,
    middle_index = 81,
    lower_limit = -20,
    upper_limit = 30
    )
{
  # ---- Parameters ----
  #load(file = "M:/Bayesian DPM Model/Bayesian DPM Model/Real Data Analysis/real_data.Rdata")
  # step = 1 / 4
  # start_index = 1
  # end_index = 201
  # middle_index = 81
  # lower_limit = -20
  # upper_limit = 30
  n_pred <- length(start_index:end_index)
  

  
  # Number of subjects and indexing
  X_amy_cv <- X_amy[train_id,]
  rownames(X_amy_cv) <- NULL
  X_tau_cv <- X_tau[which(tau_amy_link %in% train_id),]
  rownames(X_tau_cv) <- NULL
  dta_amy_cv <- dta_amy[which(dta_amy$RID_num %in% train_id),]
  dta_amy_cv$RID_num <- as.integer(factor(dta_amy_cv$RID))
  rownames(dta_amy_cv) <- NULL
  dta_tau_cv <- dta_tau[which(dta_tau$RID %in% dta_amy_cv$RID),]
  dta_tau_cv$RID_num <- as.integer(factor(dta_tau_cv$RID))
  rownames(dta_tau_cv) <- NULL
  tau_amy_link <- match(unique(dta_tau_cv$RID), unique(dta_amy_cv$RID))
  
  n_subjects_amy <- dim(X_amy_cv)[1]
  n_subjects_tau <- dim(X_tau_cv)[1]
    
  # Number of observations
  N_amy <- dim(dta_amy_cv)[1]
  N_tau <- dim(dta_tau_cv)[1]
  
  dta_surv_df_cv <- dta_surv_df[which(dta_surv_df$RID %in% intersect(unique(dta_tau_cv$RID),dta_surv_df$RID)),]
  rownames(dta_surv_df_cv) <- NULL
  n_subjects_surv <- nrow(dta_surv_df_cv)
  tau_surv_link <- which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)
  
  alpha_tau_pos <- rep(Inf, length(unique(dta_tau_cv$RID)))
  alpha_tau_pos[which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)] <- dta_surv_df_cv$time
  
  # Covariate dimension
  P <- dim(X_amy)[2]
  
  # Min and max age for tau and amyloid observations
  age_amy <- dta_amy_cv[,c('RID',"age_at_scan")] %>%
    group_by(RID) %>%
    summarise(t_min = min(age_at_scan, na.rm = TRUE),
              t_max = max(age_at_scan, na.rm = TRUE), 
              mean_t = mean(age_at_scan, na.rm = TRUE),
              .groups = "drop")
  
  age_tau <- dta_tau_cv[,c('RID',"age_at_scan")] %>%
    group_by(RID) %>%
    summarise(t_min = min(age_at_scan, na.rm = TRUE),
              t_max = max(age_at_scan, na.rm = TRUE), 
              mean_t = mean(age_at_scan, na.rm = TRUE),
              .groups = "drop")
  
  K <- 8
  
  # Data
  data <- list(y_pred_amy = dta_amy_cv$SUVR, y_pred_tau = dta_tau_cv$mesial_temporal_tau_pet_metaroi,
               t = dta_surv_df_cv$t, censored = dta_surv_df_cv$censored)
  data$constraint_alpha_amy <- rep(1, n_subjects_amy)
  data$constraint_alpha_tau <- rep(1, n_subjects_tau)
  
  # Spline basis
  sbBd <- function(K = 10, lower_limit = 0.5, upper_limit = 2, step = step){
    
    knots = seq(lower_limit, upper_limit, length = K)
    knts <- seq(lower_limit, upper_limit, length = K)[-c(1, K)]
    rtt <- seq(lower_limit, upper_limit, length = K)[c(1,K)]
    K <- K + 2
    # if(!is.null(v)){
    #   B <- bSpline(v, intercept = TRUE, knots = knts, Boundary.knots = rtt)
    # }else{
    #   B <- NULL
    # }
    
    # Generate spline basis for prediction
    v_pred <- seq(lower_limit, upper_limit, by = step)
    B_pred <- bSpline(v_pred, intercept = TRUE, knots = knts, Boundary.knots = rtt)
    
    return(list(B_pred = B_pred, v_pred = v_pred, knots = knots))
  }
  res_d <- sbBd(K = K, lower_limit = lower_limit, upper_limit = upper_limit, step = step)
  
  K <- K + 2
  
  # ---- Constants, Data, Inits ----
  constants <- list(
    N_amy = N_amy,
    N_tau = N_tau,
    
    n_subjects_amy = n_subjects_amy,
    n_subjects_tau = n_subjects_tau,
    n_subjects_surv = n_subjects_surv,
    
    K = K,          # number of spline basis functions
    P = P,
    B_d = res_d$B_pred, 
    
    n_pred = n_pred, # length of prediction interval
    
    subject_amy = dta_amy_cv$RID_num, 
    subject_tau = dta_tau_cv$RID_num, 
    
    Xsub_amy = X_amy_cv,
    Xsub_tau = X_tau_cv,
    Xsub_surv = X_tau_cv[tau_surv_link, ],
    
    max_t_amy = age_amy$t_max,
    min_t_amy =  age_amy$t_min,
    
    max_t_tau = age_tau$t_max,
    min_t_tau = age_tau$t_min,
    
    alpha_tau_pos = alpha_tau_pos,
  
    tau_amy_link = tau_amy_link,
    tau_surv_link = tau_surv_link,
    
    t_amy = dta_amy_cv$age_at_scan, 
    t_tau = dta_tau_cv$age_at_scan, # age mimic
    
    step = step,
    
    amy_thresh = (24 + 102.35)/169.27, # transform from CL to standardized SUVR scale (linear transformation)
    tau_thresh = 1.41,
      
    middle_index = middle_index,
    c = dta_surv_df_cv$c)
  
  # Inits
  
  theta <- seq(-9, -3, length = K)
  tinit <- rep(NA, n_subjects_surv)
  tinit[dta_surv_df_cv$censored == 1] <- dta_surv_df_cv$time[dta_surv_df_cv$censored == 1] + 10
  inits <- list(theta_amy = theta, theta_tau = theta,
                alpha_amy = age_amy$mean_t, alpha_tau = age_tau$mean_t,
                t = tinit, shape_weib = 1, gamma_weib = c(log(median(dta_surv_df_cv$time, na.rm = TRUE)), rep(0, P)))
  
  inits$alpha_tau[which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)] <- dta_surv_df_cv$time - 5
  
  
  # Format data

  simudata_cv <- list(constants = constants, data = data, inits = inits)
  return(simudata_cv)
  #save(simudata_cv, file = "~/Downloads/Bayesian DPM Code/Real Data Analysis/real_data_cv.Rdata")
}


realdata_pred <- function(
    X_amy,
    X_tau,
    test_id,
    tau_amy_link,
    dta_surv_df,
    samples,
    simudata_cv,
    biomarker,
    type,
    step = 1 / 4,
    start_index = 1,
    end_index = 201,
    middle_index = 81,
    lower_limit = -20,
    upper_limit = 30,
    threshold = 24
)
{
  
  # biomarker = 'tau'
  # type = 'forward'
  # step = 1 / 4
  # start_index = 1
  # end_index = 201
  # middle_index = 81
  # lower_limit = -20
  # upper_limit = 30
  # tau_amy_link = simudata$constants$tau_amy_link
# filter the idx with backward and forward
  n_pred <- length(start_index:end_index)
if(biomarker == 'amy'){
df_filtered_amy <- dta_amy |>
  dplyr::group_by(RID) |>
  dplyr::filter(
    any(CENTILOIDS < threshold, na.rm = TRUE) &
      any(CENTILOIDS > threshold, na.rm = TRUE)
  ) |>
  dplyr::ungroup()

  pred_id <- intersect(test_id, df_filtered_amy$RID_num)
  X_amy_cv <- X_amy[pred_id,]
  rownames(X_amy_cv) <- NULL
  X_tau_cv <- X_tau[which(tau_amy_link %in% pred_id),]
  rownames(X_tau_cv) <- NULL
  dta_amy_cv <- dta_amy[which(dta_amy$RID_num %in% pred_id),]
  
  if(type == 'forward'){
    dta_amy_cv <-  dta_amy_cv |> dplyr::filter(CENTILOIDS < threshold)
  }else if(type == 'backward'){
    dta_amy_cv <-  dta_amy_cv |> dplyr::filter(CENTILOIDS > threshold)
    }
  dta_amy_cv$RID_num <- as.integer(factor(dta_amy_cv$RID))
  rownames(dta_amy_cv) <- NULL
  dta_tau_cv <- dta_tau[which(dta_tau$RID %in% dta_amy_cv$RID),]
  dta_tau_cv$RID_num <- as.integer(factor(dta_tau_cv$RID))
  
}else if(biomarker == 'tau'){
  df_filtered_tau <- dta_tau |>
    dplyr::group_by(RID) |>
    dplyr::filter(
      any(SUVR < threshold, na.rm = TRUE) &
        any(SUVR > threshold, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  df_filtered_amy <- dta_amy[which(dta_amy$RID %in% unique(df_filtered_tau$RID)),] 
  pred_id <- intersect(test_id, df_filtered_amy$RID_num)
  X_amy_cv <- X_amy[pred_id,]
  rownames(X_amy_cv) <- NULL
  X_tau_cv <- X_tau[which(tau_amy_link %in% pred_id),]
  rownames(X_tau_cv) <- NULL
  dta_amy_cv <- dta_amy[which(dta_amy$RID_num %in% pred_id),]
  dta_amy_cv$RID_num <- as.integer(factor(dta_amy_cv$RID))
  rownames(dta_amy_cv) <- NULL
  dta_tau_cv <- dta_tau[which(dta_tau$RID %in% dta_amy_cv$RID),]
  if(type == 'forward'){
    dta_tau_cv <-  dta_tau_cv |> dplyr::filter(SUVR < threshold)
  }else if(type == 'backward'){
    dta_tau_cv <-  dta_tau_cv |> dplyr::filter(SUVR > threshold)
  }
  dta_tau_cv$RID_num <- as.integer(factor(dta_tau_cv$RID))
}

rownames(dta_tau_cv) <- NULL
tau_amy_link <- match(unique(dta_tau_cv$RID), unique(dta_amy_cv$RID))

n_subjects_amy <- dim(X_amy_cv)[1]
n_subjects_tau <- dim(X_tau_cv)[1]

# Number of observations
N_amy <- dim(dta_amy_cv)[1]
N_tau <- dim(dta_tau_cv)[1]

dta_surv_df_cv <- dta_surv_df[which(dta_surv_df$RID %in% intersect(unique(dta_tau_cv$RID),dta_surv_df$RID)),]
rownames(dta_surv_df_cv) <- NULL
n_subjects_surv <- nrow(dta_surv_df_cv)
tau_surv_link <- which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)

alpha_tau_pos <- rep(Inf, length(unique(dta_tau_cv$RID)))
alpha_tau_pos[which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)] <- dta_surv_df_cv$time

# Covariate dimension
P <- dim(X_amy)[2]

# Min and max age for tau and amyloid observations
age_amy <- dta_amy_cv[,c('RID',"age_at_scan")] %>%
  group_by(RID) %>%
  summarise(t_min = min(age_at_scan, na.rm = TRUE),
            t_max = max(age_at_scan, na.rm = TRUE), 
            mean_t = mean(age_at_scan, na.rm = TRUE),
            .groups = "drop")

age_tau <- dta_tau_cv[,c('RID',"age_at_scan")] %>%
  group_by(RID) %>%
  summarise(t_min = min(age_at_scan, na.rm = TRUE),
            t_max = max(age_at_scan, na.rm = TRUE), 
            mean_t = mean(age_at_scan, na.rm = TRUE),
            .groups = "drop")

K <- 8

# Data
data <- list(y_pred_amy = dta_amy_cv$SUVR, y_pred_tau = dta_tau_cv$mesial_temporal_tau_pet_metaroi,
             t = dta_surv_df_cv$t, censored = dta_surv_df_cv$censored)
data$constraint_alpha_amy <- rep(1, n_subjects_amy)
data$constraint_alpha_tau <- rep(1, n_subjects_tau)

# Spline basis
sbBd <- function(K = 10, lower_limit = 0.5, upper_limit = 2, step = step){
  
  knots = seq(lower_limit, upper_limit, length = K)
  knts <- seq(lower_limit, upper_limit, length = K)[-c(1, K)]
  rtt <- seq(lower_limit, upper_limit, length = K)[c(1,K)]
  K <- K + 2
  # if(!is.null(v)){
  #   B <- bSpline(v, intercept = TRUE, knots = knts, Boundary.knots = rtt)
  # }else{
  #   B <- NULL
  # }
  
  # Generate spline basis for prediction
  v_pred <- seq(lower_limit, upper_limit, by = step)
  B_pred <- bSpline(v_pred, intercept = TRUE, knots = knts, Boundary.knots = rtt)
  
  return(list(B_pred = B_pred, v_pred = v_pred, knots = knots))
}
res_d <- sbBd(K = K, lower_limit = lower_limit, upper_limit = upper_limit, step = step)

K <- K + 2

# ---- Constants, Data, Inits ----
constants <- list(
  N_amy = N_amy,
  N_tau = N_tau,
  
  n_subjects_amy = n_subjects_amy,
  n_subjects_tau = n_subjects_tau,
  n_subjects_surv = n_subjects_surv,
  
  K = K,          # number of spline basis functions
  P = P,
  B_d = res_d$B_pred, 
  
  n_pred = n_pred, # length of prediction interval
  
  subject_amy = dta_amy_cv$RID_num, 
  subject_tau = dta_tau_cv$RID_num, 
  
  Xsub_amy = X_amy_cv,
  Xsub_tau = X_tau_cv,
  Xsub_surv = X_tau_cv[tau_surv_link, ],
  
  max_t_amy = age_amy$t_max,
  min_t_amy =  age_amy$t_min,
  
  max_t_tau = age_tau$t_max,
  min_t_tau = age_tau$t_min,
  
  alpha_tau_pos = alpha_tau_pos,
  
  tau_amy_link = tau_amy_link,
  tau_surv_link = tau_surv_link,
  
  t_amy = dta_amy_cv$age_at_scan, 
  t_tau = dta_tau_cv$age_at_scan, # age mimic
  
  step = step,
  
  amy_thresh = (24 + 102.35)/169.27, # transform from CL to standardized SUVR scale (linear transformation)
  tau_thresh = 1.41,
  
  middle_index = middle_index,
  c = dta_surv_df_cv$c)

# Inits

theta <- seq(-9, -3, length = K)
tinit <- rep(NA, n_subjects_surv)
tinit[dta_surv_df_cv$censored == 1] <- dta_surv_df_cv$time[dta_surv_df_cv$censored == 1] + 10

inits <- list(theta_amy = theta, theta_tau = theta,
              alpha_amy = age_amy$mean_t, alpha_tau = age_tau$mean_t,
              t = tinit, shape_weib = 1, gamma_weib = c(log(median(dta_surv_df_cv$time, na.rm = TRUE)), rep(0, P)))

inits$alpha_tau[which(unique(dta_tau_cv$RID) %in% dta_surv_df_cv$RID)] <- dta_surv_df_cv$time - 5


# Format data

simudata_pred <- list(constants = constants, data = data, inits = inits)

#samples <- res$samples
f_pred_amy <- (samples[,grep("f_pred_amy", colnames(samples))]) 
f_pred_tau <- (samples[,grep("f_pred_tau", colnames(samples))]) 

gamma_amy <- samples[, grepl("^gamma_amy(\\[|$|_)", colnames(samples)), drop = FALSE]
gamma_tau <- (samples[, grepl("^gamma_tau(\\[|$|_)", colnames(samples)), drop = FALSE]) 

sigma_a_amy <- (samples[, grepl("sigma_a_amy", colnames(samples)), drop = FALSE])
sigma_a_tau <- (samples[, grepl("sigma_a_tau", colnames(samples)), drop = FALSE])

sigma_y_amy <- (samples[, grepl("sigma_y_amy", colnames(samples)), drop = FALSE])
sigma_y_tau <- (samples[, grepl("sigma_y_tau", colnames(samples)), drop = FALSE])

rho <- (samples[, grepl("rho", colnames(samples)), drop = FALSE])

simudata_cv$data <- c(
  simudata_cv$data,
  list(
    f_pred_amy = (f_pred_amy),
    f_pred_tau = (f_pred_tau),
    
    gamma_amy = (gamma_amy),
    gamma_tau = (gamma_tau),
    
    sigma_a_amy = c(sigma_a_amy),
    sigma_a_tau = c(sigma_a_tau),
    
    sigma_y_amy = c(sigma_y_amy),
    sigma_y_tau = c(sigma_y_tau),
    rho = c(rho)
  )
)

# Constants
simudata_cv$constants$n_mcmc <- length(rho)
simudata_cv$constants$n_post <- simudata_pred$constants$n_subjects_amy
simudata_cv$constants$N_amy  <- simudata_pred$constants$N_amy
simudata_cv$constants$N_tau  <- simudata_pred$constants$N_tau

# Data
simudata_cv$data$y_pred_amy <- simudata_pred$data$y_pred_amy 
simudata_cv$data$y_pred_tau <- simudata_pred$data$y_pred_tau

simudata_cv$data$t_amy <- simudata_pred$data$t_amy
simudata_cv$data$t_tau <- simudata_pred$data$t_tau

simudata_cv$constants$Xsub_amy  <- simudata_pred$constants$Xsub_amy 

# Constraint
simudata_cv$constants$subject_amy  <- simudata_pred$constants$subject_amy 
simudata_cv$constants$subject_tau  <- simudata_pred$constants$subject_tau
simudata_cv$constants$tau_amy_link <- simudata_pred$constants$tau_amy_link
simudata_cv$data$constraint_alpha_amy <- rep(1, simudata_pred$constants$n_subjects_amy)
simudata_cv$data$constraint_alpha_tau <- rep(1, simudata_pred$constants$n_subjects_amy)

# Subject-specific boundaries for the dconstraint loops
simudata_cv$constants$max_t_amy <- simudata_pred$constants$max_t_amy
simudata_cv$constants$min_t_amy <- simudata_pred$constants$min_t_amy

max_t_tau <- rep(NA, simudata_pred$constants$n_subjects_amy)
min_t_tau <- rep(NA, simudata_pred$constants$n_subjects_amy)

max_t_tau[simudata_pred$constants$tau_amy_link] <- simudata_pred$constants$max_t_tau
min_t_tau[simudata_pred$constants$tau_amy_link] <- simudata_pred$constants$min_t_tau

max_t_tau[-simudata_pred$constants$tau_amy_link] <- 0
min_t_tau[-simudata_pred$constants$tau_amy_link] <- 140

simudata_cv$constants$max_t_tau <- max_t_tau
simudata_cv$constants$min_t_tau <- min_t_tau

# inits
simudata_cv$inits$alpha <- matrix(NA, nrow = 2, ncol = simudata_pred$constants$n_subjects_amy)
simudata_cv$inits$alpha[1, ] <- simudata_pred$inits$alpha_amy
simudata_cv$inits$alpha[2, simudata_pred$constants$tau_amy_link] <- simudata_pred$inits$alpha_tau
simudata_cv$inits$alpha[2, -simudata_pred$constants$tau_amy_link] <- simudata_pred$inits$alpha_amy[-simudata_pred$constants$tau_amy_link]

return(list(simudata_cv = simudata_cv,simudata_pred = simudata_pred, pred_id = pred_id))
}

