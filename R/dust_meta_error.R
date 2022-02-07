# MSIRS2I2 dust METAPOP 2021

###################
## compile model ##
###################

library(odin.dust)
msirs_meta_dust <- odin.dust::odin_dust("models/dustmeta.R")

source("R/stochastic_init.R") # function for starting model at zoographic equilibrium

#############################
## user defined variables: ##
#############################

n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of subpops

# correction for balancing ex vs int foi

correction_ex <- matrix(c(2, 3, 3, 3, 2,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          2, 3, 3, 3, 2), nrow = 5, ncol = 5)

###################
## run the model ##
###################

###########################
###########################
#### test connectivity ####
###########################
###########################

n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta = par_grid_metapop_test$seasonality[1], N_0 = 40,
                      mu, N_age, n_r = n_r, n_c = n_c)

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, 
                                               beta = 0.4,
                                               gamma = gamma, 
                                               sigma = 1/90,
                                               sigma_m = sigma_m, 
                                               Ab_susc = 0.75,
                                               mAb_susc = mAb_susc, 
                                               reduced_shed = 1/90, 
                                               mu = mu, 
                                               N_0 = 40, 
                                               importation_rate = importation_rate, 
                                               imp_t = imp_t, 
                                               delta =  1, 
                                               connectivity = 0.001,
                                               correction_ex = correction_ex,
                                               S_ini_p = S_ini_p,
                                               foi_bg_usr = foi_bg_usr), 
                                   step = 0, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)