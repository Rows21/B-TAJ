# Minimal simulation -> single-subject posterior prediction workflow.
# Run this file from the repository root with R 4.6.1.

library(ggplot2)
library(matrixStats)

source("model/data_gen_new.R")
source("model/dweibull3.R")
source("model/data_compile bivariate.R")
source("model/bayesian dpm bivariate plus survival.R")
source("model/mcmc_AF_slice bivariate.R")
source("model/cumulative_trapz.R")
source("model/simu.R")

set.seed(20261008)
sim <- simu(
  varying_coef_code_sim = varying_coef_code_sim,
  varying_coef_code = varying_coef_code,
  n_subjects = 100,
  sa = 8,
  sy = 0.10,
  outer_iter = 1,
  niter = 2000,
  nburnin = 1000,
  nchains = 1,
  save_res = FALSE
)

simudata <- sim$simudata
samples <- sim$res$samples
cs <- simudata$constants
cm <- simudata$cModelsim

# Use one subject with both amyloid and tau observations.
subject_id <- intersect(unique(cs$subject_amy), cs$tau_amy_link)[1]
amy_idx <- which(cs$subject_amy == subject_id)
tau_idx <- which(cs$subject_tau == subject_id)
tau_number <- match(subject_id, cs$tau_amy_link)

if (!length(amy_idx) || !length(tau_idx)) {
  stop("The simulation did not produce a subject with both biomarkers.")
}

take <- function(pattern, draw_ids = seq_len(nrow(samples))) {
  z <- grep(pattern, colnames(samples))
  if (!length(z)) stop("Missing posterior monitor: ", pattern, call. = FALSE)
  samples[draw_ids, z, drop = FALSE]
}

draw_ids <- sample(seq_len(nrow(samples)), min(1000L, nrow(samples)))
f_amy <- take("^f_pred_amy\\[", draw_ids)
f_tau <- take("^f_pred_tau\\[", draw_ids)
gamma_amy <- take("^gamma_amy\\[", draw_ids)
gamma_tau <- take("^gamma_tau\\[", draw_ids)
alpha_names <- grep("^alpha\\[", colnames(samples))
if (length(alpha_names) == 0L) {
  stop("The simulation posterior does not contain alpha draws.", call. = FALSE)
}

get_scalar <- function(pattern, default) {
  z <- grep(pattern, colnames(samples))
  if (!length(z)) return(rep(default, length(draw_ids)))
  as.numeric(samples[draw_ids, z[1]])
}

sigma_a <- pmax(get_scalar("^sigma_a_amy$", 8), 1e-6)
sigma_y <- pmax(get_scalar("^sigma_y_amy$", 0.1), 1e-6)
rho <- pmax(pmin(get_scalar("^rho$", 0.3), 0.999), -0.999)
gamma_amy <- as.matrix(gamma_amy)
gamma_tau <- as.matrix(gamma_tau)
x <- as.numeric(cs$Xsub_amy[subject_id, ])
mu_amy <- drop(gamma_amy %*% x)
mu_tau <- drop(gamma_tau %*% x)

# This is the single-subject analogue of realdata_pred: retain the fitted
# posterior draw, plug in the subject's observations, and resample alpha.
grid_index <- function(age, alpha, n_grid, middle, step) {
  pmax(1L, pmin(n_grid, round((age - alpha) / step + middle)))
}

n_particles <- 1000L
alpha_post <- matrix(NA_real_, nrow = length(draw_ids), ncol = 2,
                     dimnames = list(NULL, c("alpha_amy", "alpha_tau")))
for (j in seq_along(draw_ids)) {
  z <- matrix(rnorm(2L * n_particles), ncol = 2L)
  z[, 2] <- rho[j] * z[, 1] + sqrt(1 - rho[j]^2) * z[, 2]
  particles <- cbind(mu_amy[j] + sigma_a[j] * z[, 1],
                     mu_tau[j] + sigma_a[j] * z[, 2])
  ia <- grid_index(cm$t_amy[amy_idx], particles[, 1], ncol(f_amy),
                   cs$middle_index, cs$step)
  it <- grid_index(cm$t_tau[tau_idx], particles[, 2], ncol(f_tau),
                   cs$middle_index, cs$step)
  log_w <- vapply(seq_len(n_particles), function(k) {
    sum(dnorm(cm$y_sparse_amy[amy_idx], f_amy[j, ia[k]], sigma_y[j],
              log = TRUE)) +
      sum(dnorm(cm$y_sparse_tau[tau_idx], f_tau[j, it[k]], sigma_y[j],
                log = TRUE))
  }, numeric(1))
  w <- exp(log_w - max(log_w))
  alpha_post[j, ] <- particles[sample.int(n_particles, 1L, prob = w / sum(w)), ]
}

# Posterior trajectory bands after plugging in the estimated alpha.
age_grid <- seq(min(c(cm$t_amy[amy_idx], cm$t_tau[tau_idx])) - 5,
                max(c(cm$t_amy[amy_idx], cm$t_tau[tau_idx])) + 5,
                length.out = 100)
trajectory <- function(f, alpha, y_name) {
  values <- vapply(seq_along(age_grid), function(k) {
    idx <- grid_index(age_grid[k], alpha, ncol(f), cs$middle_index, cs$step)
    f[cbind(seq_len(nrow(f)), idx)]
  }, numeric(nrow(f)))
  data.frame(age = age_grid,
             median = colMedians(values),
             lower = colQuantiles(values, probs = 0.025),
             upper = colQuantiles(values, probs = 0.975),
             biomarker = y_name)
}

amyloid_plot_data <- trajectory(f_amy, alpha_post[, 1], "Amyloid SUVR")
tau_plot_data <- trajectory(f_tau, alpha_post[, 2], "Tau PET")

observed <- rbind(
  data.frame(age = cm$t_amy[amy_idx], value = cm$y_sparse_amy[amy_idx],
             biomarker = "Amyloid SUVR"),
  data.frame(age = cm$t_tau[tau_idx], value = cm$y_sparse_tau[tau_idx],
             biomarker = "Tau PET")
)

true_trajectory <- rbind(
  data.frame(age = age_grid,
             value = cm$f_pred[grid_index(age_grid, cm$alpha[1, subject_id],
                                           length(cm$f_pred), cs$middle_index,
                                           cs$step)],
             biomarker = "Amyloid SUVR"),
  data.frame(age = age_grid,
             value = cm$f_pred[grid_index(age_grid, cm$alpha[2, subject_id],
                                           length(cm$f_pred), cs$middle_index,
                                           cs$step)],
             biomarker = "Tau PET")
)

plot_data <- rbind(
  transform(amyloid_plot_data, value = median),
  transform(tau_plot_data, value = median)
)

trajectory_plot <- ggplot() +
  geom_ribbon(data = rbind(amyloid_plot_data, tau_plot_data),
              aes(age, ymin = lower, ymax = upper, fill = biomarker),
              alpha = 0.20) +
  geom_line(data = true_trajectory,
            aes(age, value, color = biomarker), linetype = "dashed",
            linewidth = 0.8) +
  geom_line(data = plot_data,
            aes(age, value, color = biomarker), linewidth = 1) +
  geom_point(data = observed,
             aes(age, value, color = biomarker), size = 2) +
  facet_wrap(~ biomarker, scales = "free_y") +
  labs(title = paste("Simulation subject", subject_id),
       subtitle = "Points: observed data; dashed: simulation truth; band: posterior 95% interval",
       x = "Age (years)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

print(trajectory_plot)
saveRDS(list(subject_id = subject_id, observed = observed,
             alpha_truth = c(cm$alpha[1, subject_id], cm$alpha[2, subject_id]),
             alpha_post = alpha_post, trajectory_plot = trajectory_plot),
        file.path("new_subject_posterior_prediction",
                  "simulation_single_subject_prediction.rds"))
