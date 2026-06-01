library(matrixStats)
library(nimble)

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

# load ground truth path
root <- "M:/Bayesian DPM Model/Bayesian DPM Model/Bayesian DPM Model Simulation new/arch/results"
n_subjects <- 1000
sa <- 8
sy <- 0.05
sc <- paste0("s1_",n_subjects,sa*10,sy*100)
path_gt <- file.path(root,sc,"gt")

## ====
alpha_amy_gt <- read.csv(file.path(path_gt, "alpha_amy.csv"))
alpha_tau_gt <- read.csv(file.path(path_gt, "alpha_tau.csv"))

f_pred_gt <- read.csv(file.path(path_gt, "f_pred.csv"))

gamma_amy_gt <- read.csv(file.path(path_gt, "gamma_amy.csv"))
gamma_tau_gt <- read.csv(file.path(path_gt, "gamma_tau.csv"))

tau_amy_link <- read.csv(file.path(path_gt, "tau_amy_link.csv"))

outer_iter <- 100
alpha_amy_est_mat <- matrix(NA, nrow = outer_iter, ncol = n_subjects)
alpha_tau_est_mat <- matrix(NA, nrow = outer_iter, ncol = n_subjects/2)

f_pred_amy_est_mat <- matrix(NA, nrow = outer_iter, ncol = 250)
f_pred_tau_est_mat <- matrix(NA, nrow = outer_iter, ncol = 250)

gamma_amy_est_mat <- matrix(NA, nrow = outer_iter, ncol = 3)
gamma_tau_est_mat <- matrix(NA, nrow = outer_iter, ncol = 3)

shape_weib_est_vec <- rep(NA, outer_iter)
gamma_weib_est_mat <- matrix(NA, nrow = outer_iter, ncol = 3)

sa_amy_vec <- rep(NA, outer_iter)
sa_tau_vec <- rep(NA, outer_iter)
sy_amy_vec <- rep(NA, outer_iter)
sy_tau_vec <- rep(NA, outer_iter)

# coverage
coverage <- list()
coverage[['alpha_amy']] <- rep(NA,100)
coverage[['alpha_tau']] <- rep(NA,100)

coverage[['f_pred_amy']] <- rep(NA,100)
coverage[['f_pred_tau']] <- rep(NA,100)

coverage[['gamma_amy']] <- rep(NA,100)
coverage[['gamma_tau']] <- rep(NA,100)

coverage[['gamma_weib']] <- rep(NA,100)
coverage[['shape_weib']] <- rep(NA,100)

# for the i th iter
for(i in 1:outer_iter){

path_iter <- file.path(root,sc,paste0("iter",i))

alpha_amy <- as.numeric(alpha_amy_gt[i,])
alpha_samples <- as.matrix(read.csv(file.path(path_iter, "alpha.csv")))
alpha_amy_samples <- alpha_samples[,seq(1, 2000, by = 2)]
alpha_amy_est <- colMedians(alpha_amy_samples)
alpha_amy_est_mat[i,] <- as.numeric(alpha_amy_est)

cover_res <- cover_ci(alpha_amy_samples, alpha_amy)
coverage[['alpha_amy']][i] <- cover_res[[2]]
# plot(alpha_amy, alpha_amy_est)
# abline(0, 1)
# cor(alpha_amy, alpha_amy_est)
# mean(((alpha_amy_est - alpha_amy)^2))

alpha_tau <- as.numeric(alpha_tau_gt[i,])
alpha_tau_samples <- alpha_samples[,seq(2, 2000, by = 2)]
alpha_tau_samples <- alpha_tau_samples[,as.numeric(tau_amy_link[i,])]
alpha_tau_est <- colMedians(alpha_tau_samples)
alpha_tau_est_mat[i,] <- alpha_tau_est

cover_res <- cover_ci(alpha_tau_samples, alpha_tau)
coverage[['alpha_tau']][i] <- cover_res[[2]]

f_pred <- as.numeric(f_pred_gt[i,])
f_pred_amy_samples <- as.matrix(read.csv(file.path(path_iter, "f_pred_amy.csv")))
f_pred_amy_est <- colMedians(f_pred_amy_samples[,])
f_pred_amy_est_mat[i,] <- f_pred_amy_est

cover_res <- cover_ci(f_pred_amy_samples, f_pred)
coverage[['f_pred_amy']][i] <- cover_res[[2]]
# plot(f_pred, type ='l')
# lines(f_pred_amy_est, lty = 2)

f_pred_tau_samples <- as.matrix(read.csv(file.path(path_iter, "f_pred_tau.csv")))
f_pred_tau_est <- colMedians(f_pred_tau_samples[,])
f_pred_tau_est_mat[i,] <- f_pred_tau_est
cover_res <- cover_ci(f_pred_tau_samples, f_pred)
coverage[['f_pred_tau']][i] <- cover_res[[2]]
# plot(f_pred, type ='l')
# lines(f_pred_tau_est, lty = 2)

gamma_amy <- as.numeric(gamma_amy_gt[i,])
gamma_amy_samples <- as.matrix(read.csv(file.path(path_iter, "gamma_amy.csv")))
gamma_amy_est <- colMedians(gamma_amy_samples[,])
gamma_amy_est_mat[i,] <- gamma_amy_est

cover_res <- cover_ci(gamma_amy_samples, gamma_amy)
coverage[['gamma_amy']][i] <- cover_res[[2]]

# plot(gamma_amy, gamma_amy_est)
# abline(0, 1)
# cor(gamma_amy, gamma_amy_est)
# mean(((gamma_amy_est - gamma_amy)^2))

gamma_tau <- as.numeric(gamma_tau_gt[i,])
gamma_tau_samples <- as.matrix(read.csv(file.path(path_iter, "gamma_tau.csv")))
gamma_tau_est <- colMedians(gamma_tau_samples[,])
gamma_tau_est_mat[i,] <- gamma_tau_est

cover_res <- cover_ci(gamma_tau_samples, gamma_tau)
coverage[['gamma_tau']][i] <- cover_res[[2]]
# plot(gamma_tau, gamma_tau_est)
# abline(0, 1)
# cor(gamma_tau, gamma_tau_est)
# mean(((gamma_tau_est - gamma_tau)^2))

gamma_weib <- c(4.765,0.062,-0.039)
gamma_weib_samples <- as.matrix(read.csv(file.path(path_iter, "gamma_weib.csv")))
gamma_weib_est <- colMedians(gamma_weib_samples[,])
gamma_weib_est_mat[i,] <- gamma_weib_est

cover_res <- cover_ci(gamma_weib_samples[,], gamma_weib)
coverage[['gamma_weib']][i] <- cover_res[[2]]

shape_weib <- 1.84
shape_weib_samples <- as.vector(read.csv(file.path(path_iter, "shape_weib.csv")))
shape_weib_est <- median(shape_weib_samples$x[])
shape_weib_est_vec[i] <- shape_weib_est

q1 <- quantile(shape_weib_samples$x[],prob = 0.025)
q2 <- quantile(shape_weib_samples$x[],prob = 0.975)
inside <- shape_weib >= q1 & shape_weib <= q2
coverage[['shape_weib']][i] <- mean(inside)

sa_amy <- as.vector(read.csv(file.path(path_iter, "sigma_a_amy.csv")))
sa_amy_vec[i] <- mean(sa_amy$x[])
#plot(sa_amy$x, type  = 'l')
sa_tau <- as.vector(read.csv(file.path(path_iter, "sigma_a_tau.csv")))
sa_tau_vec[i] <- mean(sa_tau$x[])
#plot(sa_tau$x, type  = 'l')

sy_amy <- as.vector(read.csv(file.path(path_iter, "sigma_y_amy.csv")))
sy_amy_vec[i] <- mean(sy_amy$x[])
#plot(sy_amy$x, type  = 'l')
sy_tau <- as.vector(read.csv(file.path(path_iter, "sigma_y_tau.csv")))
sy_tau_vec[i] <- mean(sy_tau$x[])
#plot(sy_tau$x, type  = 'l')
}

aa <- (alpha_amy_est_mat - alpha_amy_gt)
gamma_weib_gt <- t(matrix(c(4.765,0.062,-0.039), ncol = outer_iter, nrow = 3))

mean(as.matrix((alpha_amy_est_mat - alpha_amy_gt)^2))
mean(as.matrix((alpha_tau_est_mat - alpha_tau_gt[1:(dim(alpha_tau_gt)[2]/2)])^2))
mean(as.matrix((gamma_amy_est_mat - gamma_amy_gt)^2))
mean(as.matrix((gamma_tau_est_mat - gamma_tau_gt)^2))
mean(as.matrix((gamma_weib_est_mat - gamma_weib_gt)^2))
mean((shape_weib_est_vec - shape_weib)^2)
mean((sa_amy_vec - sa)^2)
mean((sa_tau_vec - sa)^2)
mean((sy_amy_vec - sy)^2)
mean((sy_tau_vec - sy)^2)

mean(as.matrix((alpha_amy_est_mat - alpha_amy_gt)^2/alpha_amy_gt))
mean(as.matrix((alpha_tau_est_mat - alpha_tau_gt[1:(dim(alpha_tau_gt)[2]/2)])^2/alpha_tau_gt[1:(dim(alpha_tau_gt)[2]/2)]))
mean(as.matrix((gamma_amy_est_mat - gamma_amy_gt)^2/abs(gamma_amy_gt)))
mean(as.matrix((gamma_tau_est_mat - gamma_tau_gt)^2/abs(gamma_tau_gt)))
mean(as.matrix((gamma_weib_est_mat - gamma_weib_gt)^2/abs(gamma_weib_gt)))
mean((shape_weib_est_vec - shape_weib)^2/abs(shape_weib))
mean((sa_amy_vec - sa)^2/sa)
mean((sa_tau_vec - sa)^2/sa)
mean((sy_amy_vec - sy)^2/sy)
mean((sy_tau_vec - sy)^2/sy)

sd(as.matrix((alpha_amy_est_mat - alpha_amy_gt)^2/alpha_amy_gt))/sqrt(n_subjects)
sd(as.matrix((alpha_tau_est_mat - alpha_tau_gt[1:(dim(alpha_tau_gt)[2]/2)])^2/alpha_tau_gt[1:(dim(alpha_tau_gt)[2]/2)]))/sqrt(n_subjects/2)
sd(as.matrix((gamma_amy_est_mat - gamma_amy_gt)^2/abs(gamma_amy_gt)))/sqrt(3)
sd(as.matrix((gamma_tau_est_mat - gamma_tau_gt)^2/abs(gamma_tau_gt)))/sqrt(3)
sd(as.matrix((gamma_weib_est_mat - gamma_weib_gt)^2))/sqrt(3)



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

mean(coverage$alpha_amy)
mean(coverage$alpha_tau)
mean(coverage$gamma_amy)
mean(coverage$gamma_tau)
mean(coverage$f_pred_amy)
mean(coverage$f_pred_tau)
mean(coverage$shape_weib)
mean(coverage$gamma_weib)

x <- 1:250
q10 <- apply(f_pred_amy_est_mat, 2, quantile, 0.10)
med <- apply(f_pred_amy_est_mat, 2, median)
q90 <- apply(f_pred_amy_est_mat, 2, quantile, 0.90)
gt <- as.numeric(f_pred_gt[1,])

plot(x, gt, type = "l", lwd = 2, col = "black", ylim = range(c(gt, q10, q90)))
lines(x, med, lwd = 2, col = "blue")
lines(x, q10, lty = 2, col = "blue")
lines(x, q90, lty = 3, col = "blue")

x <- 1:250
q10 <- apply(f_pred_tau_est_mat, 2, quantile, 0.10)
med <- apply(f_pred_tau_est_mat, 2, median)
q90 <- apply(f_pred_tau_est_mat, 2, quantile, 0.90)
gt <- as.numeric(f_pred_gt[1,])

plot(x, gt, type = "l", lwd = 2, col = "black", ylim = range(c(gt, q10, q90)))
lines(x, med, lwd = 2, col = "blue")
lines(x, q10, lty = 2, col = "blue")
lines(x, q90, lty = 3, col = "blue")

df_plot <- data.frame(
  x = x,
  gt = gt,
  q10 = q10,
  q90 = q90,
  med = med
)
library(ggplot2)

ggplot(df_plot, aes(x = x)) +
  
  # 不确定性区间（10–90）
  geom_ribbon(aes(ymin = q10, ymax = q90),
              fill = "steelblue", alpha = 0.2) +
  
  # quantile lines
  geom_line(aes(y = q10, linetype = "10th"), linewidth = 0.4) +
  geom_line(aes(y = med, linetype = "median"), linewidth = 0.5) +
  geom_line(aes(y = q90, linetype = "90th"), linewidth = 0.4) +
  
  # ground truth
  geom_line(aes(y = gt, color = "Ground truth"),
            linewidth = 0.5) +
  
  scale_linetype_manual(
    values = c(
      "10th" = "dashed",
      "median" = "solid",
      "90th" = "dotdash"
    )
  ) +
  
  scale_color_manual(
    values = c("Ground truth" = "blue")
  ) +
  
  labs(
    x = "Index",
    y = "Value",
    linetype = "",
    color = ""
  ) +
  theme(legend.position = "bottom")+
  theme_minimal(base_size = 10)
