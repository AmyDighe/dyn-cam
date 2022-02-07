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


# number of age classes
N_age <- 49

# mean birth rate (per camel per day) 
alpha <- 0.000565 #0.0006 # 90% female * 50% reproductive age * 45.2% fecundity

# the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 
gamma <- 1/duration_infection

# the average rate of waning of mAbs
sigma_m <- 4.42/360 # from the catalytic work

# proportion of susceptibility experienced by calves with mAbs
mAb_susc <- 0 # default = 0

# the age dependent removal rate - balance birthrate
mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

mu <- vector(length = N_age)
mu[1:12] <- mu_1st_yr
mu[13:24] <- mu_2nd_yr
mu[25:36] <- mu_3rd_yr
mu[37:48] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

# importation rate for introducing infectious individuals
importation_rate <- 0

# background force of infection for first ten years
foi_bg_usr <- 0.000015
imp_t <- 1


# index for summing births for each age class
ind1 <- rep(0,12)
ind2 <- rep(0,12)

for(y in 2:13){
  ind1[y-1] <- 360 - ((y - 1) * 30) + 1 
  ind2[y-1] <- 360 - ((y - 2) * 30)
}

# repeating for 4 years to cover all age classes
ind1 <- rep(ind1, 4)
ind2 <- rep(ind2, 4)

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
S_ini_p <- stoch_init(alpha, delta = 1, N_0 = 40,
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
