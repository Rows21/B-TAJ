library(nimble)
library(splines2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mgcv)
library(MASS)

## ---- Data 1 & 2 Generation ----
# ---- function wrap ----
# ---- spline basis ----
sbB <- function(K = 10, lower_limit = 0.5, upper_limit = 2, v = NULL){
  
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
  v_pred <- seq(lower_limit, upper_limit, length.out = 50)
  B_pred <- bSpline(v_pred, intercept = TRUE, knots = knts, Boundary.knots = rtt)
  
  return(list(B_pred = B_pred, v_pred = v_pred, knots = knots))
}

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


#v_pred <- seq(lower_limit, upper_limit, length.out = 50)
sbB_mgcv <- function(K = 10, lower_limit = 0.5, upper_limit = 2, v = NULL) {
  
  n <- 50
  if (is.null(v)) v <- seq(lower_limit, upper_limit, length.out = n)
  
  smooth_obj <- mgcv::smoothCon(
    s(v, bs = "cr", k = K), 
    data = data.frame(v = v),
    absorb.cons = TRUE,
    diagonal.penalty = TRUE
  ) 
  
  B_pred <- matrix(1, nrow = n, ncol = K)
  B_pred[,2:K] <- smooth_obj[[1]]$X
  v_pred <- v
  
  S <- smooth_obj$S[[1]]
  
  return(list(B_pred = B_pred, v_pred = v_pred, S = S, smooth = smooth_obj))
}

# ---- generate ----
generate_d1 <- function(n_subjects = 300, 
                        P = 3, 
                        K = 10,
                        Kd = 5,
                        scenario = 1,
                        start_index = 1, 
                        end_index = 100, 
                        range = 30, 
                        step = step,
                        lower_limit = 0.5, 
                        upper_limit = 2,
                        sg_r = 10, # sigma_gamma_r
                        #sb = 0.2, # sigma_b[n]
                        #sr = 0.1, # sigma_r
                        sa = 4,   # sigma_alpha_i
                        sy = 0.1, # sigma_y_pred,
                        target_threshold = target_threshold,
                        mid_point = mid_point,
                        a_thresh = a_thresh,
                        C = "cr"
                          ){

  # Assign 2–6 observations per subject
  obs_per_subject <- sample(2:6, size = n_subjects, prob = c(0.44, .27, .13, .1, .06), replace = TRUE)
  # Total number of observations
  n <- sum(obs_per_subject)
  n_pred <- length(start_index:end_index)
  
  # Subject ID for each row
  subject_id <- rep(1:n_subjects, times = obs_per_subject)

  # Generate covariates
  if(scenario == 1){
    x1 <- sample(c(0, 1), n_subjects, replace = TRUE)
    x2 <- sample(c(0, 1), n_subjects, replace = TRUE)
    #x3 <- subject <- rnorm(n_subjects, 0, 1)
    #x3 <- x3 - mean(x3)
    int <- rep(1, n_subjects)
    #X_subject <- cbind(int,x1, x2, x3)
    X_subject <- cbind(int,x1, x2)
    #X_subject_S <- cbind(int,x1, x2, scale(x3))
  }else{
    X_subject <- (matrix(runif(n_subjects * P, -0.5, 5), nrow = n_subjects, ncol = P)) 
    #X_subject_S <- scale(X_subject)
  }
  

  
  X <- X_subject[subject_id, ]  # shape: n x r

  
  # Time: random per observation
  # Create a "cluster center" for each subject
  time_center <- runif(n_subjects, 0.1, 0.9)  # avoid edges
  
  # Generate time points clustered around each subject's center
  v <- unlist(lapply(1:n_subjects, function(s) {
    rnorm(obs_per_subject[s], mean = time_center[s], sd = 0.05)
  }))
  
  # Truncate to [0, 1] to keep within bounds
  v <- pmin(pmax(v, 0), 1)
  
  res_sb <- sbB(K = K, lower_limit = lower_limit, upper_limit = upper_limit)
  
  res_d <- sbBd(K = K, lower_limit = (1 - n_pred/2) * step, upper_limit = (n_pred - n_pred/2) * step, step = step)
  res_dc <- sbBd(K = Kd, lower_limit = (1 - n_pred/2) * step, upper_limit = (n_pred - n_pred/2) * step, step = step)

  K <- K+2
  Kd <- Kd + 2
  # generate theta true
  theta_true <- matrix(NA, nrow = K, ncol = P)
  
  for (p in 1:P) {
    # Create increasing-then-flat pattern
    increase <- c(-9, seq(-9, -6/p, length = 2))
    flat2 <- rep(-3.8/p, K - (length(c(increase)) + 2))
    theta_true[,p] <- c(increase, flat2, -4, -4)
  }

  theta_true[, 3] <- 0.0
  theta_true[, 2] <- (theta_true[, 2]+9)* 1e-1 + 0.05
  #theta_true[, 4] <- theta_true[, 2]/2
  theta_true[1,2:3] <- rep(0,2)
  
  theta_true[,2:P] <- 0
  
  # for (p in 1:P) {
  #   # Create increasing-then-flat pattern
  #   increase <- c(0, seq(0, .06/p, length = 4))
  #   flat2 <- rep(.06/p, 12 - (length(c(increase)) + 2))
  #   theta_true[,p ] <- c(increase, flat2, 0.01, .005)
  # }
  # 
  # theta_true[, 3] <- 0
  # theta_true[, 2] <- theta_true[, 2] * 1.25
  # theta_true[, 4] <- -theta_true[, 2]/3
  
  
  #tau_amy_link <- sort(sample(n_subjects, n_subjects))
  tau_amy_link <- sort(sample(n_subjects, n_subjects/2))
  subject_tau <- subject_id[subject_id %in% tau_amy_link]
  
  constants <- list(
    N = n,          # total number of observations
    N_amy = n,
    N_tau = length(subject_tau),
    
    n_subjects = n_subjects, # number of unique subjects
    n_subjects_amy = n_subjects,
    #n_subjects_tau = n_subjects,
    n_subjects_tau = n_subjects/2,
    
    K = K,          # number of spline basis functions
    P = P,          # covariate dimension
    #B_pred = res_sb$B_pred, # spline basis for predictions
    B_d = res_d$B_pred,
    #B_dc = res_dc$B_pred,
    #Kd = Kd,
    n_pred = n_pred, # length of prediction interval
    #v_pred = res_sb$v_pred,        # dense value grid
    #d_pred = res_d$v_pred,
    nd_pred = length(res_d$v_pred),
    #nv_pred = length(res_sb$v_pred),
    subject = subject_id,  # integer vector of subject IDs (length N)
    subject_amy = subject_id,
    subject_tau = subject_tau,
    subject_tau_dense = dense_rank(subject_tau),
    tau_amy_link = tau_amy_link,
    
    Xsub = X_subject,        # matrix of subject-level covariates (nSubjects × P)
    theta = theta_true,
    
    
    ## ---- parameters ----
    #sg_r = sg_r,
    sa = sa,
    sy = sy,
    step = step,
    range = range,
    #mid_point = mid_point,
    #knots = res_sb$knots,
    #target_threshold = target_threshold,
    a_thresh = a_thresh
  )
  
  theta <- matrix(rnorm(K * P, 0, .1), nrow = K, ncol = P)
  if(scenario == 1){theta[,1] <- c(seq(-9,-3,by = 1),rep(-3,3))}
  
  
  if(!is.null(C)){
    K_new <- 10
    constants$K <- K_new
    v_pred <- res_sb$v_pred
    data <- data.frame(v_pred = v_pred)
    smcon <- smoothCon(s(v_pred, bs=C, k=K_new), data = data, 
                       absorb.cons = FALSE, diagonal.penalty = TRUE)
    C <- smcon[[1]]$X
    
    randeff = mgcv::smooth2random(smcon[[1]], vnames = names(data), type = 2)
    X_mat_r = randeff$rand$Xr
    X_mat_f = randeff$Xf
    K0 <- dim(X_mat_f)[2]
    K1 <- dim(X_mat_r)[2]
    constants$X_mat_f <- C[,(K_new - K0 + 1):K_new]
    constants$X_mat_r <- C[,1:K1]
    constants[['K0']] <- K0
    constants[['K1']] <- K1
    
    f_old <- res_sb$B_pred %*% theta_true
    phi <- solve(t(C) %*% C, t(C) %*% f_old)
    constants$B_pred <- C
    theta_true <- phi
    constants$theta <- phi
    theta <- ginv(C) %*% res_sb$B_pred %*% theta
    # theta[,1] <- c(seq(-9,-0.5,length.out = 3), rep(-0.5,K_new - 3))
  }
  
  # if(scenario == 1){
  # }else{
  # }
  
  #theta_true <- log(theta_true)
  # constants
  
  return(list(theta_true = theta_true,
              constants = constants,
              theta = theta))
}
# --------


