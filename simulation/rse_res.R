library(matrixStats)
library(nimble)
library(ggthemr)

cover_ci <- function(samples, truth, lower_q = 0.025, upper_q = 0.975) {
  
  lower <- apply(samples, 2, quantile, probs = lower_q, na.rm = TRUE)
  upper <- apply(samples, 2, quantile, probs = upper_q, na.rm = TRUE)
  
  inside <- truth >= lower & truth <= upper
  
  list(
    count = sum(inside),
    coverage = mean(inside),
    lower = lower,
    upper = upper
  )
}

cover_ci196 <- function(samples, truth, lower_q = 0.025, upper_q = 0.975, scale = 1.96) {
  
  lower <- apply(samples, 2, quantile, probs = lower_q, na.rm = TRUE)
  upper <- apply(samples, 2, quantile, probs = upper_q, na.rm = TRUE)
  
  mid <- apply(samples, 2, quantile, probs = 0.5, na.rm = TRUE)
  half_width <- (upper - lower) / 2
  
  lower_ext <- mid - scale * half_width
  upper_ext <- mid + scale * half_width
  
  inside <- truth >= lower_ext & truth <= upper_ext
  
  list(
    count = sum(inside),
    coverage = mean(inside),
    lower = lower_ext,
    upper = upper_ext,
    original_lower = lower,
    original_upper = upper
  )
}



df_plot_a <- list()
df_plot_t <- list()
# load ground truth path
root <- "M:/Bayesian DPM Model/Bayesian DPM Model/Bayesian DPM Model Simulation new/arch/results"
n_subjects <- 1000
# sa <- 8
# sy <- 0.05
k <- 1
for(sa in c(4,8)){
  for(sy in c(0.05,0.1)){
    sc <- paste0("s1_",n_subjects,sa*10,sy*100)
    path_gt <- file.path(root,sc,"gt")
    
    ## ====
    f_pred_gt <- read.csv(file.path(path_gt, "f_pred.csv"))
    
    outer_iter <- 100
    
    f_pred_amy_est_mat <- matrix(NA, nrow = outer_iter, ncol = 250)
    f_pred_tau_est_mat <- matrix(NA, nrow = outer_iter, ncol = 250)
    
    coverage <- list()
    coverage[['f_pred_amy']] <- rep(NA,100)
    coverage[['f_pred_tau']] <- rep(NA,100)
    
    
    # for the i th iter
    for(i in 1:outer_iter){
      
      path_iter <- file.path(root,sc,paste0("iter",i))
      
      f_pred <- as.numeric(f_pred_gt[i,])
      f_pred_amy_samples <- as.matrix(read.csv(file.path(path_iter, "f_pred_amy.csv")))
      f_pred_amy_est <- colMedians(f_pred_amy_samples[,])
      f_pred_amy_est_mat[i,] <- f_pred_amy_est
      cover_res <- cover_ci196(f_pred_amy_samples, f_pred)
      coverage[['f_pred_amy']][i] <- cover_res[[2]]
      
      f_pred_tau_samples <- as.matrix(read.csv(file.path(path_iter, "f_pred_tau.csv")))
      f_pred_tau_est <- colMedians(f_pred_tau_samples[,])
      f_pred_tau_est_mat[i,] <- f_pred_tau_est
      cover_res <- cover_ci196(f_pred_tau_samples, f_pred)
      coverage[['f_pred_tau']][i] <- cover_res[[2]]
    
    }
    
    
    rse_amy <- rep(NA,100)
    for(i in 1:100){
      rse_amy[i] <- sum((f_pred_amy_est_mat[i,] - f_pred_gt[i,])^2) / sum(as.numeric(f_pred_gt[i,])^2)
    }
    
    rse_tau <- rep(NA,100)
    for(i in 1:100){
      rse_tau[i] <- sum((f_pred_tau_est_mat[i,] - f_pred_gt[i,])^2) / sum(as.numeric(f_pred_gt[i,])^2)
    }
    mean(rse_amy)
    mean(rse_tau)
    
    mean(coverage$f_pred_amy)
    mean(coverage$f_pred_tau)
    x <- 1:250
    q10 <- apply(f_pred_amy_est_mat, 2, quantile, 0.10)
    med <- apply(f_pred_amy_est_mat, 2, median)
    q90 <- apply(f_pred_amy_est_mat, 2, quantile, 0.90)
    gt <- as.numeric(f_pred_gt[1,])
    df_plot_amy <- data.frame(
      x = x,
      gt = gt,
      q10 = q10,
      q90 = q90,
      med = med,
      ss = paste0("(",sa,", ",sy,")")
    )
    
    
    x <- 1:250
    q10 <- apply(f_pred_tau_est_mat, 2, quantile, 0.10)
    med <- apply(f_pred_tau_est_mat, 2, median)
    q90 <- apply(f_pred_tau_est_mat, 2, quantile, 0.90)
    gt <- as.numeric(f_pred_gt[1,])
    
    # plot(x, gt, type = "l", lwd = 2, col = "black", ylim = range(c(gt, q10, q90)))
    # lines(x, med, lwd = 2, col = "blue")
    # lines(x, q10, lty = 2, col = "blue")
    # lines(x, q90, lty = 3, col = "blue")
    
    df_plot_tau <- data.frame(
      x = x,
      gt = gt,
      q10 = q10,
      q90 = q90,
      med = med,
      ss = paste0("(",sa,", ",sy,")")
    )
    df_plot_a[[k]] <- df_plot_amy
    df_plot_t[[k]] <- df_plot_tau
    k <- k+1
  }}

df_plot1 <- rbind(df_plot_a[[1]],df_plot_a[[2]],
                 df_plot_a[[3]],df_plot_a[[4]])
df_plot1['type'] <- '1'
df_plot2 <- rbind(df_plot_t[[1]],df_plot_t[[2]],
                  df_plot_t[[3]],df_plot_t[[4]])
df_plot2['type'] <- '2'
df_plot <- rbind(df_plot1, df_plot2)
library(ggplot2)
ggthemr('fresh')
ggplot(df_plot, aes(x = x)) +
  
  # # 不确定性区间（10–90）
  # geom_ribbon(aes(ymin = q10, ymax = q90),
  #             fill = "steelblue", alpha = 0.2) +
  # ground truth
  geom_line(aes(y = gt), color = "black",
            linewidth = 1) +
  
  # quantile lines
  geom_line(aes(y = q10, linetype = "10th", color = type), linewidth = 0.5) +
  geom_line(aes(y = med, linetype = "median", color = type), linewidth = 0.5) +
  geom_line(aes(y = q90, linetype = "90th", color = type), linewidth = 0.5) +
  
  facet_grid(type~ss)+
  scale_linetype_manual(
    values = c(
      "10th" = "dashed",
      "median" = "solid",
      "90th" = "dotdash"
    )
  ) +
  scale_color_manual(values = c("#E41A1C", 
                                "#E41A1C"
  ))+
  
  labs(
    x = "Scenario",
    y = "Clinical marker",
    linetype = "",
    color = ""
  ) +
  # theme_minimal(base_size = 10)+
  theme(
    legend.position = "bottom",

    axis.title = element_text(size = 30),
    axis.text = element_text(size = 10),

    strip.text = element_text(size = 20),

    legend.text = element_text(size = 20),
    legend.title = element_text(size = 10)
  )
ggsave(paste0('M:/Bayesian DPM Model/Bayesian DPM Model/Bayesian DPM Model Simulation new/','curve.png'),
       width = 19,
       height = 10)
