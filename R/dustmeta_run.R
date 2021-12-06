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

N_age <- 49 # number of age classes
N_patch <- 25 # number of sub-populations

# input a value for birth rate (per camel per day) 
alpha <- 0.000565

# input a value for the baseline effective contact rate
beta <- 0.6

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14
gamma <- 1/duration_infection

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 60 # default = 
sigma <- 1/duration_immunity # default = 

# input a value for the average rate of waning of mAbs
sigma_m <- 4.42/365 # from the catalytic work

# input a value between 0 and 1 for susceptibility of individuals with mAbs/Abs
## 0 would mean mAbs/Abs afford complete protection
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from shedding in Alharbi 

# input values for the age dependent removal rate - balance birthrate
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

# input an initial population size
N_0 <- 1000000

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day
imp_t <- 1510000  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# set level of connectivity between patches
connectivity <- 0.05

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta, N_0, mu, N_age, N_patch)


###################
## run the model ##
###################

n_particles <- 2L

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, N_patch = N_patch, 
                                               alpha = alpha, beta = beta, gamma = gamma, 
                                               sigma = sigma, sigma_m = sigma_m, 
                                               Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                               reduced_shed = reduced_shed, mu = mu, N_0 = N_0, 
                                               importation_rate = importation_rate, imp_t = imp_t, 
                                               delta = delta, connectivity = connectivity, 
                                               S_ini_p = S_ini_p), 
                         step = 1, 
                         n_particles = n_particles, 
                         n_threads = 1L, 
                         seed = 1L) 
n_steps <- 360

msirs_model$run(10)

