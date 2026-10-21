##########################################################
# LOAD FUNCTIONS
##########################################################

path <- c("M:/Github/B-TAJ/model")
setwd(path)
library(matrixStats)
source('data_gen_new.R')
source('dweibull3.R')
source('data_compile bivariate.R')
source('bayesian dpm bivariate plus survival.R')
source('mcmc_AF_slice bivariate.R')
source('cumulative_trapz.R')
source('simu.R')

# path to generate data
path <- c("M:/Github/B-TAJ/")

##########################################################
# GENERATE DATA
##########################################################
for(n_subjects in c(1000)){
  for(sa in c(8)){
    for(sy in c(0.05)){
      res <- simu(varying_coef_code_sim, 
                  varying_coef_code,
                  n_subjects = n_subjects,
                  sa = sa,
                  sy = sy,
                  outer_iter = 1,
                  niter = 2000, nburnin = 1000, #thin = 1, 
                  nchains = 1,
                  save_res = T)
    }
  }
}

