# dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)
library(tictoc) # for benchmarking
library(rlist)

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

# background force of infection for first ten years
foi_bg_usr <- 0.000015
imp_t <- 1

e_t <- 360*15
m_t <- 360*25
p_t <- 360*35


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
REGION_1 <- c("South_Asia",rep("Africa", 3), rep("Middle_East", 4), 
              rep("Africa", 3), rep("Middle_East", 5), rep("South_Asia", 2), 
              rep("Africa", 3), "Middle_East", "Africa")

foi_df <- readRDS(file = "data/foi_df.rds")
foi_df$region <- REGION_1
foi_range <- readRDS(file = "data/foi_range.rds")

par_grid_R0 <- data.frame(Ab_susc = rep(c(0.25, 0.75, 1), times = 6),
                          sigma = rep(c(1/30, 1/90, 1/90, 1, 1, 1/30), times = 3),
                          red_shed = rep(c(1/90, 1/2, 1/4), each = 6))

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


beta_list <- list(MSIRS2_0.25_30_90 = seq(0.08, 3.6, by = 0.08),
                  MSIRS2_0.75_90_90 = seq(0.08, 3.6, by = 0.08),
                  MSIRS2_1_90_90 = seq(0.08, 3.6, by = 0.08),
                  MSIS2_0.25_1_90 = seq(0.08, 3.36, by = 0.08),
                  MSIS2_0.75_1_90 = seq(0.08, 3.2, by = 0.08),
                  MSIRS2_1_30_90 = seq(0.08, 3.2, by = 0.08),
                  MSIRS2_0.25_30_2 = seq(0.08, 0.72, by = 0.02),
                  MSIRS2_0.75_90_2 = seq(0.08, 0.64, by = 0.02),
                  MSIRS2_1_90_2 = seq(0.08, 0.64, by = 0.02),
                  MSIS2_0.25_1_2 = seq(0.08, 0.64, by = 0.02),
                  MSIS2_0.75_1_2 = seq(0.08, 0.32, by = 0.01),
                  MSIRS2_1_30_2 = seq(0.08, 0.4, by = 0.02),
                  MSIRS2_0.25_30_4 = seq(0.08, 1.2, by = 0.05),
                  MSIRS2_0.75_90_4 = seq(0.08, 1, by = 0.05),
                  MSIRS2_1_90_4 = seq(0.08, 1, by = 0.05),
                  MSIS2_0.25_1_4 = seq(0.08, 1, by = 0.05),
                  MSIS2_0.75_1_4 = seq(0.08, 0.48, by = 0.02),
                  MSIRS2_1_30_4 = seq(0.08, 0.6, by = 0.02))

############################################
# setting up grid for persistence analyses #
############################################

waning <- 1/c(90, 1)
susc <- c(0.75, 1)
shedding <- c("shedding_0.01" = 0.01, "shedding_0.25" = 0.25, "shedding_0.50" = 0.5)
beta <- list("shedding_0.01" = c(3.5, 7, 14, 21)/14,
             "shedding_0.25" = c(2.5, 4, 5.5, 7.5)/14,
             "shedding_0.50" = c(1.75, 2.5, 3.5, 4.5)/14)
seasonality <- c(0.5, 1)
pop <- c(c(1,2.5, 5, 7.5) %o% 10^(2:4), 1 %o% 10^(5:7))

# standard coarse single patch
par_grid_core_0.01 <- expand.grid(shedding = shedding["shedding_0.01"],
                             beta = beta[["shedding_0.01"]],
                             waning = waning[1], 
                             susc = susc[1],
                             seasonality = seasonality, 
                             pop = pop)
par_grid_core_0.25 <- expand.grid(shedding = shedding["shedding_0.25"],
                               beta = beta[["shedding_0.25"]],
                               waning = waning[1], 
                               susc = susc[1],
                               seasonality = seasonality, 
                               pop = pop)
par_grid_core_0.50 <- expand.grid(shedding = shedding["shedding_0.50"],
                               beta = beta[["shedding_0.50"]],
                               waning = waning[1], 
                               susc = susc[1],
                               seasonality = seasonality, 
                               pop = pop)
par_grid_core <- rbind(par_grid_core_0.01,
                       par_grid_core_0.25,
                       par_grid_core_0.50)

# for finding point of bifurcation - fine beta 

beta_extra <- list("shedding_0.01" = c(5)/14,
             "shedding_0.25" = c()/14,
             "shedding_0.50" = c()/14)

par_grid_fine_beta_0.01 <- expand.grid(shedding = shedding["shedding_0.01"],
                                       beta = beta_extra[["shedding_0.01"]],
                                       waning = waning[1], 
                                       susc = susc[1],
                                       seasonality = seasonality, 
                                       pop = pop)
par_grid_fine_beta_0.25 <- expand.grid(shedding = shedding["shedding_0.25"],
                                       beta = beta_extra[["shedding_0.25"]],
                                       waning = waning[1], 
                                       susc = susc[1],
                                       seasonality = seasonality, 
                                       pop = pop)
par_grid_fine_beta_0.50 <- expand.grid(shedding = shedding["shedding_0.50"],
                                       beta = beta_extra[["shedding_0.50"]],
                                       waning = waning[1], 
                                       susc = susc[1],
                                       seasonality = seasonality, 
                                       pop = pop)
par_grid_fine_beta <- rbind(par_grid_fine_beta_0.01,
                            par_grid_fine_beta_0.25,
                            par_grid_fine_beta_0.50)


# sensitivity analyses
par_grid_sens <- expand.grid(waning = waning[2])

# grid for meta-population model
par_grid_metapop <- par_grid_core %>% filter(seasonality == 1)
