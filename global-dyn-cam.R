# dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)
library(tictoc) # for benchmarking

##################################
# parameters that are not varied #
##################################

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

#############################
# stochastic initialisation #
#############################

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


###############
# R0 analysis #
###############

foi_df <- readRDS(file = "data/foi_df.rds")
foi_range <- readRDS(file = "data/foi_range.rds")

par_grid_R0 <- expand.grid(Ab_susc = c(0.25, 0.75, 1),
                           sigma = c(1/90, 1/30, 1),
                           red_shed = c(1/92, 1/2))

scenario <- vector(length = 18)
for(i in 1:18){
  
  if(par_grid_R0$sigma[i] < 1){
    scenario[i] <- paste("MSIRS2", par_grid_R0$Ab_susc[i], 
                         1/par_grid_R0$sigma[i], 1/par_grid_R0$red_shed[i],
                         sep = "_")
  } else {
    scenario[i] <- paste("MSIS2", par_grid_R0$Ab_susc[i], 
                         1/par_grid_R0$sigma[i], 1/par_grid_R0$red_shed[i],
                         sep = "_")
  }
}

beta_list <- vector(mode = "list", length = 18)
for(i in 1:18){
  beta_list[[i]] <- seq(0.08, 3.6, by = 0.08)
}
names(beta_list) <- scenario


REGION_1 <- c("South_Asia",rep("Africa", 3), rep("Middle_East", 4), 
              rep("Africa", 3), rep("Middle_East", 5), rep("South_Asia", 2), 
              rep("Africa", 3), "Middle_East", "Africa")


############################################
# setting up grid for persistence analyses #
############################################

waning <- 1/c(30, 30*3)
seasonality <- c(0, 0.5, 1)
beta <- c(0.25, 0.5, 1)
pop_og <- c(c(1,5) %o% 10^(3:7))
pop <- c(c(1,5) %o% 10^(3:6))
import_time <- c(151, 269)

# standard coarse single patch
par_grid <- expand.grid(waning = waning, 
                        seasonality = seasonality, 
                        beta = beta, 
                        pop = pop_og,
                        import_time = import_time)

# for finding point of bifurcation - fine beta 
par_grid_fine_beta <- expand.grid(waning = waning[2], 
                        seasonality = seasonality[2:3], 
                        beta = c(0.3, 0.35, 0.4, 0.45, 1.5), 
                        pop = pop,
                        import_time = import_time[1])

# for finding the critical community size - fine pop
par_grid_all_pop <- expand.grid(waning = waning, 
                               seasonality = seasonality, 
                               beta = c(0.25, 0.5, 1.0, 1.5), 
                               pop = c(c(1,2.5,5,7.5) %o% 10^(3:4), c(1,5) %o% 10^(5:6)),
                               import_time = import_time)

par_grid_fine_pop <- expand.grid(waning = waning[2], 
                                 seasonality = seasonality, 
                                 beta = c(0.25, 0.5, 1.0, 1.5), 
                                 pop = c(c(2.5,7.5) %o% 10^(3:5)),
                                 import_time = import_time)

# grid for meta-population model
par_grid_metapop <- expand.grid(waning = waning[2],
                        seasonality = seasonality,
                        beta = beta,
                        pop = pop_og/25,
                        import_time = import_time[1])
