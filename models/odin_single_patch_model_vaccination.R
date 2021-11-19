#############################
## Single Patch Model 2021 ##
#############################

# This simple stochastic model simulates MERS-CoV transmission in a single, homogenously mixed, 
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is approximated as 360 days 
# It is written using Odin
# It can then be run from a user-edited R script "odin_single_patch_run.R"

N_age <- 49 #number of age classes

##############################################################################################################
# RATES OF TRANSITION BETWEEN DISEASE STATES
# all rates are daily
# user-defined, defaults shown here in the parentheses
##############################################################################################################

################
## birth rate ##
################
alpha <- user(0.0006)

#############################
## infection rates S --> I ##
#############################
beta <- user(0.3) #base rate

# adjustment of beta values
mAb_susc <- user(0) # proportion of susceptibility experienced if maternal antibodies (mAbs) present
v_mAb_susc <- user(0) # proportion of susceptibility experienced if vaccinated AND maternal antibodies (mAbs) present
Ab_susc <- user(1) # proportion of susceptibility experienced if previously infected
v_susc <- user(1) # proportion of susceptibility experienced if vaccinated (and naive to natural infection)
v_Ab_susc <- user(1) # proportion of susceptibility experienced if vaccinated AND previously infected
reduced_shed <- user(1) # proportion of infectiousness seen in reinfections. default = no difference
v_reduced_shed <- user(1) # proportion of infectiousness seen in vaccinated first infections
v_shed <- user(1) # proportion of infectiousness seen in vaccinated reinfections

# frequency dependent rate of infection
update(rate_infection) <- (beta * sum(I[1:N_age]) + beta * reduced_shed * sum(I2[1:N_age]) + beta * v_shed * sum(vI[1:N_age]) + beta * v_reduced_shed * sum(vI2[1:N_age]))/sum(N[1:N_age])
rate_infection_vaccinated <- rate_infection * v_susc
rate_infection_mAb <- rate_infection * mAb_susc
rate_infection_mAb_vaccinated <- rate_infection * v_mAb_susc
rate_reinfection <- rate_infection * Ab_susc
rate_reinfection_vaccinated <- rate_infection * v_Ab_susc

#####################
## mortality rates ##
#####################
# user-defined age-dependent mortality rate
mu_1st_yr <- user(0.0005) # death rate for 1st year of life
mu_2nd_yr <- user(0.0005) # death rate for 2nd year of life
mu_3rd_yr <- user(0.00025) # death rate for 3rd year of life
mu_4th_yr <- user(0.00025) # death rate for 4th year of life
mu_adult_over_4 <- user(0.00025) # death rate in adulthood (>4 years)
# expand these across the age strata
mu[1:12] <- mu_1st_yr
mu[13:24] <- mu_2nd_yr
mu[25:36] <- mu_3rd_yr
mu[37:48] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

###########################
## recovery rate I --> R ## where R is non-infectious and completely immune to further infection
###########################
gamma <- user(1/14) # (/day) 
v_gamma <- user(1/14)

####################################################
## rate at which complete immunity wanes R --> S2 ## 
####################################################
sigma <- user(0.0005) # (/day) to be taken from catalytic work eventually
v_sigma <- user(0.0005)

###############################################################
## rate at which maternally-acquired immunity wanes Sm --> S ##
###############################################################
sigma_m <- user(0.006) # (/day) to be taken from catalytic work eventually
v_sigma_m <- user(0.006)

#########################
## rate of vaccination ##
#########################
vax_1m_3m <- user(0)
vax_4m_8m <- user(0)
vax_9m_12m <- user(0)
vax_1y_2y <- user(0)
vax_2y_3y <- user(0)
vax_3y_4y <- user(0)
vax_adult <- user(0)

vax[1:3] <-vax_1m_3m
vax[4:8] <- vax_4m_8m
vax[9:12] <- vax_9m_12m
vax[13:24] <- vax_1y_2y
vax[25:36] <- vax_2y_3y
vax[37:48] <- vax_3y_4y
vax[N_age] <- vax_adult

##################################################
## rate at which vaccine induced immunity wanes ##
##################################################
rho <- user(0.001) #(/day) to be varied as currently unknown

##############################################################################################################
# CONVERTING THESE RATES --> PROBABILITIES
# the above boil down to rates which are converted to probabilities below
##############################################################################################################
  ## for unvaccinated individuals
p_infection <- 1 - exp(-rate_infection)
p_infection_mAb <- 1 - exp(-rate_infection_mAb)
p_reinfection <- 1 - exp(-rate_reinfection)
p_mu[1:N_age] <- 1 - exp(-mu[i]) # prob death
p_gamma <- 1 - exp(-gamma) # prob recovery
p_sigma <- 1 - exp(-sigma) # prob waned
p_sigma_m <- 1 - exp(-sigma_m) #prob mAbs waned
p_vax[1:N_age] <- 1 - exp(-vax[i]) # prob become vaccinated
  ## for vaccinated individuals
p_v_infection <- 1 - exp(-rate_infection_vaccinated)
p_v_infection_mAb <- 1 - exp(-rate_infection_mAb_vaccinated)
p_v_reinfection <- 1 - exp(-rate_reinfection_vaccinated)
p_v_gamma <- 1 - exp(-v_gamma) # prob recovery if vaccinated
p_v_sigma <- 1 - exp(-v_sigma) # prob susceptible again if vaccinated
p_v_sigma_m <- 1 - exp(-v_sigma_m) #prob mAbs waned if vaccinated
p_rho <- 1 - exp(-rho) # prob vaccine induced immunity wanes

##############################################################################################################
# OUTFLOWS
# 1. infection, recovery, waning immunity and death
# 2. ageing
##############################################################################################################

#compartments are:
# Sm (and vSm for vaccinated individuals) = protected by mAbs
# S (vS) = susceptible
# I (vI) = infected and infectious
# R (vR) = recovered and completely immune
# S2 (vS2) = immunity has waned to some degree
# I2 (vI2) = infected and infectious for the 2nd+ time

# probability of leaving each compartment for any reason 
# (other than ageing which is dealt with later)
p_Sm[1:N_age] <- 1 - exp(- (sigma_m + rate_infection_mAb + mu[i] + vax[i])) 
p_S[1:N_age] <- 1 - exp(- (rate_infection + mu[i] + vax[i])) 
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i] + vax[i])) 
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i] + vax[i])) 
p_S2[1:N_age] <- 1 - exp(- (rate_reinfection + mu[i] + vax[i])) 
p_I2[1:N_age] <- 1 - exp(- (gamma + mu[i] + vax[i])) 
  ## for vaccinated individuals
p_vSm[1:N_age] <- 1 - exp(- (v_sigma_m + rate_infection_mAb + mu[i] + rho)) 
p_vS[1:N_age] <- 1 - exp(- (rate_infection + mu[i] + rho)) 
p_vI[1:N_age] <- 1 - exp(- (v_gamma + mu[i] + rho)) 
p_vR[1:N_age] <- 1 - exp(- (v_sigma + mu[i] + rho)) 
p_vS2[1:N_age] <- 1 - exp(- (rate_reinfection + mu[i] + rho)) 
p_vI2[1:N_age] <- 1 - exp(- (v_gamma + mu[i] + rho))

# outflows due to infection, recovery or death or vaccination
outflow_Sm[1:N_age] <- rbinom(Sm[i], prob = p_Sm[i])
outflow_S[1:N_age] <- rbinom(S[i], prob = p_S[i])
outflow_I[1:N_age] <- rbinom(I[i], prob = p_I[i])
outflow_R[1:N_age] <- rbinom(R[i], prob = p_R[i])
outflow_S2[1:N_age] <- rbinom(S2[i], prob = p_S2[i])
outflow_I2[1:N_age] <- rbinom(I2[i], prob = p_I2[i])
  ## for vaccinated individuals
outflow_vSm[1:N_age] <- rbinom(vSm[i], prob = p_vSm[i])
outflow_vS[1:N_age] <- rbinom(vS[i], prob = p_vS[i])
outflow_vI[1:N_age] <- rbinom(vI[i], prob = p_vI[i])
outflow_vR[1:N_age] <- rbinom(vR[i], prob = p_vR[i])
outflow_vS2[1:N_age] <- rbinom(vS2[i], prob = p_vS2[i])
outflow_vI2[1:N_age] <- rbinom(vI2[i], prob = p_vI2[i])


###############################################################################################################
# INFLOWS
# 1. where do these outflows go? 
# 2. births
# 3. importations
###############################################################################################################

################################################################
## new infections, recoveries, waning events and vaccinations ##
################################################################

#normalising the probabilities 
norm_p_sigma_m[1:N_age] <- p_sigma_m/(p_sigma_m + p_infection_mAb + p_mu[i] + p_vax[i])
norm_p_infection_mAb[1:N_age] <- p_infection_mAb/(p_sigma_m + p_infection_mAb + p_mu[i] + p_vax[i])
norm_p_infection[1:N_age] <- p_infection/(p_infection + p_mu[i] + p_vax[i])
norm_p_reinfection[1:N_age] <- p_reinfection/(p_reinfection + p_mu[i] + p_vax[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i] + p_vax[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i] + p_vax[i])

norm_p_vax_M[1:N_age] <- p_vax[i]/(p_sigma_m + p_infection_mAb + p_mu[i] + p_vax[i])
norm_p_vax_S[1:N_age] <- p_vax[i]/(p_infection + p_mu[i] + p_vax[i])
norm_p_vax_I[1:N_age] <- p_vax[i]/(p_gamma + p_mu[i] + p_vax[i])
norm_p_vax_R[1:N_age] <- p_vax[i]/(p_sigma + p_mu[i] + p_vax[i])
norm_p_vax_S2[1:N_age] <- p_vax[i]/(p_reinfection + p_mu[i] + p_vax[i])
norm_p_vax_I2[1:N_age] <- p_vax[i]/(p_gamma + p_mu[i] + p_vax[i])

norm_p_v_sigma_m[1:N_age] <- p_v_sigma_m/(p_v_sigma_m + p_infection_mAb + p_mu[i] + p_rho)
norm_p_v_infection_mAb[1:N_age] <- p_v_infection_mAb/(p_v_sigma_m + p_v_infection_mAb + p_mu[i] + p_rho)
norm_p_v_infection[1:N_age] <- p_v_infection/(p_v_infection + p_mu[i] + p_rho)
norm_p_v_reinfection[1:N_age] <- p_v_reinfection/(p_v_reinfection + p_mu[i] + p_rho)
norm_p_v_gamma[1:N_age] <- p_v_gamma/(p_v_gamma + p_mu[i] + p_rho)
norm_p_v_sigma[1:N_age] <- p_v_sigma/(p_v_sigma + p_mu[i] + p_rho)

norm_p_rho_vM[1:N_age] <- p_rho/(p_v_sigma_m + p_v_infection_mAb + p_mu[i] + p_rho)
norm_p_rho_vS[1:N_age] <- p_rho/(p_v_infection + p_mu[i] + p_rho)
norm_p_rho_vI[1:N_age] <- p_rho/(p_v_gamma + p_mu[i] + p_rho)
norm_p_rho_vR[1:N_age] <- p_rho/(p_v_sigma + p_mu[i] + p_rho)
norm_p_rho_vS2[1:N_age] <- p_rho/(p_v_reinfection + p_mu[i] + p_rho)
norm_p_rho_vI2[1:N_age] <- p_rho/(p_v_gamma + p_mu[i] + p_rho)

# number of new infections, vaccinations, recoveries and newly susceptible
new_waned_mAb[1:N_age] <- rbinom(outflow_Sm[i], prob = norm_p_sigma_m[i])
new_infections_mAb[1:N_age] <- rbinom(outflow_Sm[i], prob = norm_p_infection_mAb[i])
#update(new_infections[1:N_age]) <- rbinom(outflow_S[i], prob = norm_p_infection[i])
new_infections[1:N_age] <- rbinom(outflow_S[i], prob = norm_p_infection[i])
new_recoveries[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_gamma[i])
new_waned[1:N_age] <- rbinom(outflow_R[i], prob = norm_p_sigma[i])
#update(new_reinfections[1:N_age]) <- rbinom(outflow_S2[i], prob = norm_p_reinfection[i])
new_reinfections[1:N_age] <- rbinom(outflow_S2[i], prob = norm_p_reinfection[i])
new_recoveries_2[1:N_age] <- rbinom(outflow_I2[i], prob = norm_p_gamma[i])

new_vax_M[1:N_age] <- rbinom(outflow_Sm[i], prob = norm_p_vax_M[i])
new_vax_S[1:N_age] <- rbinom(outflow_S[i], prob = norm_p_vax_S[i])
new_vax_I[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_vax_I[i])
new_vax_R[1:N_age] <- rbinom(outflow_R[i], prob = norm_p_vax_R[i])
new_vax_S2[1:N_age] <- rbinom(outflow_S2[i], prob = norm_p_vax_S2[i])
new_vax_I2[1:N_age] <- rbinom(outflow_I2[i], prob = norm_p_vax_I2[i])

v_new_waned_mAb[1:N_age] <- rbinom(outflow_vSm[i], prob = norm_p_v_sigma_m[i])
v_new_infections_mAb[1:N_age] <- rbinom(outflow_vSm[i], prob = norm_p_v_infection_mAb[i])
#update(v_new_infections[1:N_age]) <- rbinom(outflow_vS[i], prob = norm_p_v_infection[i])
v_new_infections[1:N_age] <- rbinom(outflow_vS[i], prob = norm_p_v_infection[i])
v_new_recoveries[1:N_age] <- rbinom(outflow_vI[i], prob = norm_p_v_gamma[i])
v_new_waned[1:N_age] <- rbinom(outflow_vR[i], prob = norm_p_v_sigma[i])
#update(v_new_reinfections[1:N_age]) <- rbinom(outflow_vS2[i], prob = norm_p_v_reinfection[i])
v_new_reinfections[1:N_age] <- rbinom(outflow_vS2[i], prob = norm_p_v_reinfection[i])
v_new_recoveries_2[1:N_age] <- rbinom(outflow_vI2[i], prob = norm_p_v_gamma[i])

new_waned_vax_M[1:N_age] <- rbinom(outflow_vSm[i], prob = norm_p_rho_vM[i])
new_waned_vax_S[1:N_age] <- rbinom(outflow_vS[i], prob = norm_p_rho_vS[i])
new_waned_vax_I[1:N_age] <- rbinom(outflow_vI[i], prob = norm_p_rho_vI[i])
new_waned_vax_R[1:N_age] <- rbinom(outflow_vR[i], prob = norm_p_rho_vR[i])
new_waned_vax_S2[1:N_age] <- rbinom(outflow_vS2[i], prob = norm_p_rho_vS2[i])
new_waned_vax_I2[1:N_age] <- rbinom(outflow_vI2[i], prob = norm_p_rho_vI2[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)

birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial population size
new_births <- rpois(birth_rate) #per day

births_protected <- rbinom(new_births, prob = seropoz_A4) # protected by mAbs
births_not_protected <- new_births - births_protected #  NOT protected by mAbs

#########################
## importation process ##
#########################

importation_rate <- user(0.01) # should be proportional to population size?
imported_cases <- rpois(importation_rate) #per day
imp_t[] <- user() # a user defined time at which cases are imported

###############################################################################################################
# EQUATIONS for movement of individuals between disease states and age classes
# time-step = 1 day
###############################################################################################################

## STEP 1 - disease state changes
## all importations (whether using rate or pulse) occur into age class 25 (~2 years old)

# no need for tt %% here as no inflows
new_Sm[1] <- Sm[1] - outflow_Sm[1] + births_protected + new_waned_vax_M[1]
new_Sm[2:N_age] <- Sm[i] - outflow_Sm[i] + new_waned_vax_M[i]

new_S[1] <- S[1] - outflow_S[1] + new_waned_mAb[1] + births_not_protected + new_waned_vax_S[1] 
new_S[2:N_age] <- S[i] - outflow_S[i] + new_waned_mAb[i] + new_waned_vax_S[i]

new_I[1:24] <-  I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i] + new_waned_vax_I[i]
new_I[25] <- I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i] + imported_cases + new_waned_vax_I[i]
new_I[26:N_age] <-  I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i] + new_waned_vax_I[i]

new_R[1:N_age] <- R[i] - outflow_R[i] + new_recoveries[i] + new_recoveries_2[i] + new_waned_vax_R[i]
new_S2[1:N_age] <- S2[i] - outflow_S2[i] + new_waned[i] + new_waned_vax_S2[i]
new_I2[1:N_age] <- I2[i] - outflow_I2[i] + new_reinfections[i] + new_waned_vax_I2[i]

    ## for vaccinated individuals (no importation or births so only one line needed each)

new_vSm[1:N_age] <- vSm[i] - outflow_vSm[i] + new_vax_M[i]
new_vS[1:N_age] <- vS[i] - outflow_vS[i] + v_new_waned_mAb[i] + new_vax_S[i]
new_vI[1:N_age] <-  vI[i] - outflow_vI[i] + v_new_infections[i] + v_new_infections_mAb[i] + new_vax_I[i]
new_vR[1:N_age] <- vR[i] - outflow_vR[i] + v_new_recoveries[i] + v_new_recoveries_2[i]  + new_vax_R[i]
new_vS2[1:N_age] <- vS2[i] - outflow_vS2[i] + v_new_waned[i] + new_vax_S2[i]
new_vI2[1:N_age] <- vI2[i] - outflow_vI2[i] + v_new_reinfections[i]  + new_vax_I2[i]

## STEP 2 update with ageing

update(Sm[1]) <- if(tt %% 30 == 0) 0 else new_Sm[1]
update(Sm[2:48]) <- if(tt %% 30 == 0) new_Sm[i - 1] else new_Sm[i]
update(Sm[N_age]) <- if(tt %% 30 == 0) new_Sm[i - 1] + new_Sm[i] else new_Sm[i]

update(S[1]) <- if(tt %% 30 == 0) 0 else new_S[1]
update(S[2:48]) <- if(tt %% 30 == 0) new_S[i - 1] else new_S[i]
update(S[N_age]) <- if(tt %% 30 == 0) new_S[i - 1] + new_S[i] else new_S[i]

update(I[1]) <-  if(tt %% 30 == 0) 0 else new_I[1]
update(I[2:24]) <- if(tt %% 30 == 0) new_I[i - 1] else new_I[i]
update(I[25]) <- if(tt %% 30 == 0) new_I[i - 1] else if(tt == imp_t[1] || tt == imp_t[2] || tt == imp_t[3] || tt == imp_t[4] || tt == imp_t[5]) 5 + new_I[i] else new_I[i]
update(I[26:48]) <- if(tt %% 30 == 0) new_I[i - 1] else new_I[i]
update(I[N_age]) <- if(tt %% 30 == 0) new_I[i - 1] + new_I[i] else new_I[i]

update(R[1]) <- if(tt %% 30 == 0) 0 else new_R[1]
update(R[2:48]) <- if(tt %% 30 == 0) new_R[i - 1] else new_R[i]
update(R[N_age]) <- if(tt %% 30 == 0) new_R[i - 1] + new_R[i] else new_R[i]

update(S2[1]) <- if(tt %% 30 == 0) 0 else new_S2[1]
update(S2[2:48]) <- if(tt %% 30 == 0) new_S2[i - 1] else new_S2[i]
update(S2[N_age]) <- if(tt %% 30 == 0) new_S2[i - 1] + new_S2[i] else new_S2[i]

update(I2[1]) <- if(tt %% 30 == 0) 0 else new_I2[1]
update(I2[2:48]) <- if(tt %% 30 == 0) new_I2[i - 1] else new_I2[i]
update(I2[N_age]) <- if(tt %% 30 == 0) new_I2[i - 1] + new_I2[i] else new_I2[i]

    ## and for vaccinated individuals
update(vSm[1]) <- if(tt %% 30 == 0) 0 else new_vSm[1]
update(vSm[2:48]) <- if(tt %% 30 == 0) new_vSm[i - 1] else new_vSm[i]
update(vSm[N_age]) <- if(tt %% 30 == 0) new_vSm[i - 1] + new_vSm[i] else new_vSm[i]

update(vS[1]) <- if(tt %% 30 == 0) 0 else new_vS[1]
update(vS[2:48]) <- if(tt %% 30 == 0) new_vS[i - 1] else new_vS[i]
update(vS[N_age]) <- if(tt %% 30 == 0) new_vS[i - 1] + new_vS[i] else new_vS[i]

update(vI[1]) <-  if(tt %% 30 == 0) 0 else new_vI[1]
update(vI[2:48]) <- if(tt %% 30 == 0) new_vI[i - 1] else new_vI[i]
update(vI[N_age]) <- if(tt %% 30 == 0) new_vI[i - 1] + new_vI[i] else new_vI[i]

update(vR[1]) <- if(tt %% 30 == 0) 0 else new_vR[1]
update(vR[2:48]) <- if(tt %% 30 == 0) new_vR[i - 1] else new_vR[i]
update(vR[N_age]) <- if(tt %% 30 == 0) new_vR[i - 1] + new_vR[i] else new_vR[i]

update(vS2[1]) <- if(tt %% 30 == 0) 0 else new_vS2[1]
update(vS2[2:48]) <- if(tt %% 30 == 0) new_vS2[i - 1] else new_vS2[i]
update(vS2[N_age]) <- if(tt %% 30 == 0) new_vS2[i - 1] + new_vS2[i] else new_vS2[i]

update(vI2[1]) <- if(tt %% 30 == 0) 0 else new_vI2[1]
update(vI2[2:48]) <- if(tt %% 30 == 0) new_vI2[i - 1] else new_vI2[i]
update(vI2[N_age]) <- if(tt %% 30 == 0) new_vI2[i - 1] + new_vI2[i] else new_vI2[i]

#update(seroprevalence[1:N_age]) <- (I[i] + R[i] + S2[i] + I2[i]) / (Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i])

update(seropoz_A4) <- (sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]) + 
                         sum(vS2[N_age]) + sum(vI[N_age]) + sum(vI2[N_age]) + sum(vR[N_age]))/ 
  (sum(Sm[N_age]) + sum(S[N_age]) + sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]) + 
    sum(vSm[N_age]) + sum(vS[N_age]) + sum(vS2[N_age]) + sum(vI[N_age]) + sum(vI2[N_age]) + sum(vR[N_age]))

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

N[1:N_age] <- Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i] + vSm[i] + vS[i] + vI[i] + vR[i] + vS2[i] + vI2[i]

##################################################################################################################################
# initial conditions
# equilibrium solution approximated to speed up balancing of demography
# no equilibrium solution for infection at this point
##################################################################################################################################

## initial states
#initial(new_infections[1:N_age]) <- 0
#initial(new_reinfections[1:N_age]) <- 0
#initial(v_new_infections[1:N_age]) <- 0
#initial(v_new_reinfections[1:N_age]) <- 0
initial(rate_infection) <- 0
initial(Sm[1:N_age]) <- 0
initial(S[1:N_age]) <- S_ini_p[i]
initial(I[1:N_age]) <- 0
initial(R[1:N_age]) <- 0
initial(S2[1:N_age]) <- 0
initial(I2[1:N_age]) <- 0
initial(vSm[1:N_age]) <- 0
initial(vS[1:N_age]) <- 0
initial(vI[1:N_age]) <- 0
initial(vR[1:N_age]) <- 0
initial(vS2[1:N_age]) <- 0
initial(vI2[1:N_age]) <- 0

initial(tt) <- 1
#initial(N[1:N_age]) <- S_ini_p[i]
#initial(seroprevalence[1:N_age]) <- 0
initial(seropoz_A4) <- 0

## initial population size for use in birthrate

N_0 <- user(1000) # user-defined

## setting initial conditions using the equilibrium solution for age distribution

births_detr[1:360] <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * i / 360)))) # change to fixed N_0 = huge to avoid NaNs
births_det[1:360] <- rpois(births_detr[i]) 


## if we start the model with the equilibrium amount in each of the first 48 month-wide compartments,
## birthrate will be set to balance summed death rate of the equilibrium age distribution

a_max <- 32 ## estimated max number of years spent in the last age class before death (only 1% of the population above after 32 yrs)
## (calculated in "calculating_a_max.R using exp decay mod)
ind1[] <- user() # indexes of birth influx for demographic equilibrium solution
ind2[] <- user()

# setting initial number of camels at demographic equilibrium
S_ini[1] <- 0
S_ini[2:48] <- sum(births_det[ind1[i]:ind2[i]]) * exp(- (30 * sum(mu[1:(i - 1)])))

# special treatment for the final open-ended compartment 49
yearly_influx <- sum(births_det[1:360]) * exp(-(30 * (sum(mu[1:48])))) #annual inflow into the final age compartment
yr[1:a_max] <- i  # number of years of influx before max age expectancy reached
# camels remaining from each cohort to enter in the last 32 years influx that remain
cohort_remaining[1:a_max] <- yearly_influx * exp(- (360 * yr[i] * mu[N_age])) 
S_ini[N_age] <- sum(cohort_remaining[1:a_max])

# getting proportion of animals in each age-class, multiplying by N_0 and rounding to whole animals
S_ini_p[1:N_age] <- round((S_ini[i] / sum(S_ini[1:N_age])) * N_0)


################################################################################################################################

# OUTPUTS

################################################################################################################################

#########################################
## number of individuals in each state ##
#########################################


output(M) <- sum(Sm[1:N_age]) # individuals protected by maternal immunity
output(S_1) <- sum(S[1:N_age]) # susceptible individuals never infected
output(I_1) <- sum(I[1:N_age]) # individuals infectious for the 1st time
output(R_1) <- sum(R[1:N_age]) # individuals recovered from a 1st infection
output(S_2) <- sum(S2[1:N_age]) # susceptible individuals whose immunity has waned
output(I_2) <- sum(I2[1:N_age]) # individuals infectious for the 2nd+ time

output(vM) <- sum(vSm[1:N_age]) # individuals protected by maternal immunity
output(vS_1) <- sum(vS[1:N_age]) # susceptible individuals never infected
output(vI_1) <- sum(vI[1:N_age]) # individuals infectious for the 1st time
output(vR_1) <- sum(vR[1:N_age]) # individuals recovered from a 1st infection
output(vS_2) <- sum(vS2[1:N_age]) # susceptible individuals whose immunity has waned
output(vI_2) <- sum(vI2[1:N_age]) # individuals infectious for the 2nd+ time

output(N_pop[1:N_age]) <- N[i]

output(Stot) <- sum(S[1:N_age]) + sum(S2[1:N_age]) + sum(Sm[1:N_age]) # total number of susceptible individuals
output(Itot) <- sum(I[1:N_age]) + sum(I2[1:N_age]) # total number of infectious individuals
output(Vtot) <- sum(vSm[1:N_age]) + sum(vS[1:N_age]) + sum(vI[1:N_age]) + sum(vR[1:N_age]) + sum(vS2[1:N_age]) + sum(vI2[1:N_age])
output(Ntot) <- sum(N[1:N_age]) # total number of individuals

####################
## seroprevalence ##
####################
output(seropoz_A) <- 100 * (sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]))/(sum(Sm[25:N_age]) + sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]))
output(seropoz_J) <- 100 * (sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24])) / (sum(Sm[1:24]) + sum(S[1:24]) + sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]))
output(seropoz_tot) <- 100 * (sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]))/(sum(Sm[1:N_age]) + sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]))

###############
## incidence ##
###############
output(incidence_new_inf) <- sum(new_infections[1:N_age])
output(incidence_indig_inf) <- sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age])
output(total_incidence) <- (sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age]))
output(vax_incidence) <- sum(new_vax_M[1:N_age]) + sum(new_vax_S[1:N_age]) + sum(new_vax_I[1:N_age]) + sum(new_vax_R[1:N_age]) + sum(new_vax_S2[1:N_age]) + sum(new_vax_I2[1:N_age])

###########
## other ##
###########
output(inf) <- rate_infection
output(birthrate) <- birth_rate
output(births) <- new_births
output(importations) <- imported_cases
output(yy) <- yr[12]
#output(foi) <- beta * sum(I[1:N_age]) / sum(N[1:N_age]) # foi applying to first infections only


################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(vax) <- N_age
dim(p_mu) <- N_age
dim(p_vax) <- N_age
dim(p_Sm) <- N_age
dim(p_S) <- N_age
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- N_age
dim(p_I2) <- N_age
  ##vax
dim(p_vSm) <- N_age
dim(p_vS) <- N_age
dim(p_vI) <- N_age
dim(p_vR) <- N_age
dim(p_vS2) <- N_age
dim(p_vI2) <- N_age

dim(Sm) <- N_age
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S2) <- N_age
dim(I2) <- N_age
  ## vax
dim(vSm) <- N_age
dim(vS) <- N_age
dim(vI) <- N_age
dim(vR) <- N_age
dim(vS2) <- N_age
dim(vI2) <- N_age

dim(S_ini) <- N_age
dim(S_ini_p) <- N_age
dim(yr) <- a_max
dim(cohort_remaining) <- a_max
dim(N_pop) <- N_age

dim(norm_p_infection) <- N_age
dim(norm_p_infection_mAb) <- N_age
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- N_age
dim(norm_p_reinfection) <- N_age
  ##vax
dim(norm_p_v_infection) <- N_age
dim(norm_p_v_infection_mAb) <- N_age
dim(norm_p_v_gamma) <- N_age
dim(norm_p_v_sigma) <- N_age
dim(norm_p_v_sigma_m) <- N_age
dim(norm_p_v_reinfection) <- N_age
dim(norm_p_vax_M) <- N_age
dim(norm_p_vax_S) <- N_age
dim(norm_p_vax_I) <- N_age
dim(norm_p_vax_R) <- N_age
dim(norm_p_vax_S2) <- N_age
dim(norm_p_vax_I2) <- N_age
dim(norm_p_rho_vM) <- N_age
dim(norm_p_rho_vS) <- N_age
dim(norm_p_rho_vI) <- N_age
dim(norm_p_rho_vR) <- N_age
dim(norm_p_rho_vS2) <- N_age
dim(norm_p_rho_vI2) <- N_age

dim(new_waned_mAb) <- N_age
dim(new_infections) <- N_age
dim(new_infections_mAb) <- N_age
dim(new_recoveries) <- N_age
dim(new_waned) <- N_age
dim(new_reinfections) <- N_age
dim(new_recoveries_2) <- N_age
  ##vax
dim(v_new_waned_mAb) <- N_age
dim(v_new_infections) <- N_age
dim(v_new_infections_mAb) <- N_age
dim(v_new_recoveries) <- N_age
dim(v_new_waned) <- N_age
dim(v_new_reinfections) <- N_age
dim(v_new_recoveries_2) <- N_age
dim(new_vax_M) <- N_age
dim(new_vax_S) <- N_age
dim(new_vax_I) <- N_age
dim(new_vax_R) <- N_age
dim(new_vax_S2) <- N_age
dim(new_vax_I2) <- N_age
dim(new_waned_vax_M) <- N_age
dim(new_waned_vax_S) <- N_age
dim(new_waned_vax_I) <- N_age
dim(new_waned_vax_R) <- N_age
dim(new_waned_vax_S2) <- N_age
dim(new_waned_vax_I2) <- N_age

dim(outflow_Sm) <- N_age
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age
dim(outflow_S2) <- N_age
dim(outflow_I2) <- N_age
  ##vax
dim(outflow_vSm) <- N_age
dim(outflow_vS) <- N_age
dim(outflow_vI) <- N_age
dim(outflow_vR) <- N_age
dim(outflow_vS2) <- N_age
dim(outflow_vI2) <- N_age

dim(new_Sm) <- N_age
dim(new_S) <- N_age
dim(new_I) <- N_age
dim(new_R) <- N_age
dim(new_S2) <- N_age
dim(new_I2) <- N_age
  ##vax
dim(new_vSm) <- N_age
dim(new_vS) <- N_age
dim(new_vI) <- N_age
dim(new_vR) <- N_age
dim(new_vS2) <- N_age
dim(new_vI2) <- N_age

dim(births_det) <- 360
dim(births_detr) <- 360
#dim(seroprevalence) <- N_age
dim(N) <- N_age
dim(ind1) <- 48
dim(ind2) <- 48
dim(imp_t) <- 5
