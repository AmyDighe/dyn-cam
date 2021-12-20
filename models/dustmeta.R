####################################################
## Meta-population Model 2021 - odin.dust version ##
####################################################

# This simple stochastic model simulates MERS-CoV transmission in a structured population,
# where the total population is split intoa grid of 25 equally sized sub populations.
# sub populations can contribute to the foi in their neighbouring patches
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is approximated as 360 days 
# It is written using Odin
# It can then be run from a user-edited R script "odin_meta_population_run.R"

N_age <- user(49) #number of age classes
N_patch <- user(25) #number of subpopulations

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

# adjusted beta values
mAb_susc <- user(0) # proportion of susceptibility experienced if maternal antibodies (mAbs) present
Ab_susc <- user(1) # proportion of susceptibility experienced if previously infected
reduced_shed <- user(1) # proportion of shedding/infectiousness seen in reinfections. default = no difference

###############################
## meta-population structure ##
###############################

connectivity <- user(0.05) ##connection strength between patches (0<x<1)

# grid connections 5 x 5
# contribution of infections in neighbouring patches
## external_I1 from first-time infections in other patches
external_I1[1] <- (sum(I[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I[1:N_age,6])/sum(N[1:N_age,6]))
external_I1[2] <- (sum(I[1:N_age,1])/sum(N[1:N_age,1])) + (sum(I[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I[1:N_age,7])/sum(N[1:N_age,7]))
external_I1[3] <- (sum(I[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I[1:N_age,8])/sum(N[1:N_age,8]))
external_I1[4] <- (sum(I[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I[1:N_age,5])/sum(N[1:N_age,5])) + (sum(I[1:N_age,9])/sum(N[1:N_age,9]))
external_I1[5] <- (sum(I[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I[1:N_age,10])/sum(N[1:N_age,10]))
external_I1[6] <- (sum(I[1:N_age,1])/sum(N[1:N_age,1])) + (sum(I[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I[1:N_age,11])/sum(N[1:N_age,11]))
external_I1[7] <- (sum(I[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I[1:N_age,6])/sum(N[1:N_age,6])) + (sum(I[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I[1:N_age, 12])/sum(N[1:N_age,12]))
external_I1[8] <- (sum(I[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I[1:N_age, 13])/sum(N[1:N_age,13]))
external_I1[9] <- (sum(I[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I[1:N_age,10])/sum(N[1:N_age,10])) + (sum(I[1:N_age, 14])/sum(N[1:N_age,14]))
external_I1[10] <- (sum(I[1:N_age,5])/sum(N[1:N_age,5])) + (sum(I[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I[1:N_age,15])/sum(N[1:N_age,15]))
external_I1[11] <- (sum(I[1:N_age,6])/sum(N[1:N_age,6])) + (sum(I[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I[1:N_age,16])/sum(N[1:N_age,16]))
external_I1[12] <- (sum(I[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I[1:N_age,11])/sum(N[1:N_age,11])) + (sum(I[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I[1:N_age, 17])/sum(N[1:N_age,17]))
external_I1[13] <- (sum(I[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I[1:N_age, 18])/sum(N[1:N_age,18]))
external_I1[14] <- (sum(I[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I[1:N_age,15])/sum(N[1:N_age,15])) + (sum(I[1:N_age, 19])/sum(N[1:N_age,19]))
external_I1[15] <- (sum(I[1:N_age,10])/sum(N[1:N_age,10])) + (sum(I[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I[1:N_age,20])/sum(N[1:N_age,20]))
external_I1[16] <- (sum(I[1:N_age,11])/sum(N[1:N_age,11])) + (sum(I[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I[1:N_age,21])/sum(N[1:N_age,21]))
external_I1[17] <- (sum(I[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I[1:N_age,16])/sum(N[1:N_age,16])) + (sum(I[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I[1:N_age, 22])/sum(N[1:N_age,22]))
external_I1[18] <- (sum(I[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I[1:N_age, 23])/sum(N[1:N_age,23]))
external_I1[19] <- (sum(I[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I[1:N_age,20])/sum(N[1:N_age,20])) + (sum(I[1:N_age, 24])/sum(N[1:N_age,24]))
external_I1[20] <- (sum(I[1:N_age,15])/sum(N[1:N_age,15])) + (sum(I[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I[1:N_age,25])/sum(N[1:N_age,25]))
external_I1[21] <- (sum(I[1:N_age,16])/sum(N[1:N_age,16])) + (sum(I[1:N_age,22])/sum(N[1:N_age,22]))
external_I1[22] <- (sum(I[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I[1:N_age,21])/sum(N[1:N_age,21])) + (sum(I[1:N_age,23])/sum(N[1:N_age,23]))
external_I1[23] <- (sum(I[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I[1:N_age,22])/sum(N[1:N_age,22])) + (sum(I[1:N_age,24])/sum(N[1:N_age,24]))
external_I1[24] <- (sum(I[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I[1:N_age,23])/sum(N[1:N_age,23])) + (sum(I[1:N_age,25])/sum(N[1:N_age,25]))
external_I1[25] <- (sum(I[1:N_age,20])/sum(N[1:N_age,20])) + (sum(I[1:N_age,24])/sum(N[1:N_age,24]))

## FOI from second-plus-time infections in other patches
external_I2[1] <- (sum(I2[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I2[1:N_age,6])/sum(N[1:N_age,6]))
external_I2[2] <- (sum(I2[1:N_age,1])/sum(N[1:N_age,1])) + (sum(I2[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I2[1:N_age,7])/sum(N[1:N_age,7]))
external_I2[3] <- (sum(I2[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I2[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I2[1:N_age,8])/sum(N[1:N_age,8]))
external_I2[4] <- (sum(I2[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I2[1:N_age,5])/sum(N[1:N_age,5])) + (sum(I2[1:N_age,9])/sum(N[1:N_age,9]))
external_I2[5] <- (sum(I2[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I2[1:N_age,10])/sum(N[1:N_age,10]))
external_I2[6] <- (sum(I2[1:N_age,1])/sum(N[1:N_age,1])) + (sum(I2[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I2[1:N_age,11])/sum(N[1:N_age,11]))
external_I2[7] <- (sum(I2[1:N_age,2])/sum(N[1:N_age,2])) + (sum(I2[1:N_age,6])/sum(N[1:N_age,6])) + (sum(I2[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I2[1:N_age, 12])/sum(N[1:N_age,12]))
external_I2[8] <- (sum(I2[1:N_age,3])/sum(N[1:N_age,3])) + (sum(I2[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I2[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I2[1:N_age, 13])/sum(N[1:N_age,13]))
external_I2[9] <- (sum(I2[1:N_age,4])/sum(N[1:N_age,4])) + (sum(I2[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I2[1:N_age,10])/sum(N[1:N_age,10])) + (sum(I2[1:N_age, 14])/sum(N[1:N_age,14]))
external_I2[10] <- (sum(I2[1:N_age,5])/sum(N[1:N_age,5])) + (sum(I2[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I2[1:N_age,15])/sum(N[1:N_age,15]))
external_I2[11] <- (sum(I2[1:N_age,6])/sum(N[1:N_age,6])) + (sum(I2[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I2[1:N_age,16])/sum(N[1:N_age,16]))
external_I2[12] <- (sum(I2[1:N_age,7])/sum(N[1:N_age,7])) + (sum(I2[1:N_age,11])/sum(N[1:N_age,11])) + (sum(I2[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I2[1:N_age, 17])/sum(N[1:N_age,17]))
external_I2[13] <- (sum(I2[1:N_age,8])/sum(N[1:N_age,8])) + (sum(I2[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I2[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I2[1:N_age, 18])/sum(N[1:N_age,18]))
external_I2[14] <- (sum(I2[1:N_age,9])/sum(N[1:N_age,9])) + (sum(I2[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I2[1:N_age,15])/sum(N[1:N_age,15])) + (sum(I2[1:N_age, 19])/sum(N[1:N_age,19]))
external_I2[15] <- (sum(I2[1:N_age,10])/sum(N[1:N_age,10])) + (sum(I2[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I2[1:N_age,20])/sum(N[1:N_age,20]))
external_I2[16] <- (sum(I2[1:N_age,11])/sum(N[1:N_age,11])) + (sum(I2[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I2[1:N_age,21])/sum(N[1:N_age,21]))
external_I2[17] <- (sum(I2[1:N_age,12])/sum(N[1:N_age,12])) + (sum(I2[1:N_age,16])/sum(N[1:N_age,16])) + (sum(I2[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I2[1:N_age, 22])/sum(N[1:N_age,22]))
external_I2[18] <- (sum(I2[1:N_age,13])/sum(N[1:N_age,13])) + (sum(I2[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I2[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I2[1:N_age, 23])/sum(N[1:N_age,23]))
external_I2[19] <- (sum(I2[1:N_age,14])/sum(N[1:N_age,14])) + (sum(I2[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I2[1:N_age,20])/sum(N[1:N_age,20])) + (sum(I2[1:N_age, 24])/sum(N[1:N_age,24]))
external_I2[20] <- (sum(I2[1:N_age,15])/sum(N[1:N_age,15])) + (sum(I2[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I2[1:N_age,25])/sum(N[1:N_age,25]))
external_I2[21] <- (sum(I2[1:N_age,16])/sum(N[1:N_age,16])) + (sum(I2[1:N_age,22])/sum(N[1:N_age,22]))
external_I2[22] <- (sum(I2[1:N_age,17])/sum(N[1:N_age,17])) + (sum(I2[1:N_age,21])/sum(N[1:N_age,21])) + (sum(I2[1:N_age,23])/sum(N[1:N_age,23]))
external_I2[23] <- (sum(I2[1:N_age,18])/sum(N[1:N_age,18])) + (sum(I2[1:N_age,22])/sum(N[1:N_age,22])) + (sum(I2[1:N_age,24])/sum(N[1:N_age,24]))
external_I2[24] <- (sum(I2[1:N_age,19])/sum(N[1:N_age,19])) + (sum(I2[1:N_age,23])/sum(N[1:N_age,23])) + (sum(I2[1:N_age,25])/sum(N[1:N_age,25]))
external_I2[25] <- (sum(I2[1:N_age,20])/sum(N[1:N_age,20])) + (sum(I2[1:N_age,24])/sum(N[1:N_age,24]))

# frequency dependent rate of infection
# when reduced_shed = 1 I and I2 are essentially the same compartment
rate_infection[1:N_patch] <- beta * (sum(I[1:N_age, i]) / sum(N[1:N_age, i]))  + 
  beta * reduced_shed * (sum(I2[1:N_age, i]) / sum(N[1:N_age, i])) +
  beta * connectivity * external_I1[i] + 
  beta * reduced_shed * connectivity * external_I2[i]

rate_infection_mAb[1:N_patch] <- mAb_susc * rate_infection[i]
rate_reinfection[1:N_patch] <- Ab_susc * rate_infection[i]

#####################
## mortality rates ##
#####################
mu[] <- user() # user-defined age-dependent mortality rate

###########################
## recovery rate I --> R ## where R is non-infectious and completely immune to further infection
###########################
gamma <- user(0.05) # (/day) 

####################################################
## rate at which complete immunity wanes R --> S2 ## 
####################################################
sigma <- user(0.0005) # (/day) to be taken from catalytic work eventually

###############################################################
## rate at which maternally-acquired immunity wanes Sm --> S ##
###############################################################
sigma_m <- user(0.006) # (/day) to be taken from catalytic work eventually

##############################################################################################################
# CONVERTING THESE RATES --> PROBABILITIES
# the above boil down to 7 rates which are converted to probabilities below
##############################################################################################################
p_infection[] <- 1 - exp(-rate_infection[i])
p_infection_mAb[] <- 1 - exp(-rate_infection_mAb[i])
p_reinfection[] <- 1 - exp(-rate_reinfection[i])
p_mu[1:N_age] <- 1 - exp(-mu[i]) # prob death
p_gamma <- 1 - exp(-gamma) # prob recovery
p_sigma <- 1 - exp(-sigma) # prob waned
p_sigma_m <- 1 - exp(-sigma_m) #prob mAbs waned

##############################################################################################################
# OUTFLOWS
# 1. infection, recovery, waning immunity and death
# 2. ageing
##############################################################################################################

#compartments are:
# Sm = susceptible but protected by mAbs
# S = susceptible
# I = infected and infectious
# R = recovered and completely immune
# S2 = immunity has waned to some degree
# I2 = infected and infectious for the 2nd+ time

# probability of leaving each compartment for any reason 
# (other than ageing which is dealt with later)
p_Sm[1:N_age, ] <- 1 - exp(- (sigma_m + rate_infection_mAb[j] + mu[i])) 
p_S[1:N_age, ] <- 1 - exp(- (rate_infection[j] + mu[i])) 
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i])) 
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i])) 
p_S2[1:N_age, ] <- 1 - exp(- (rate_reinfection[j] + mu[i])) 
p_I2[1:N_age] <- 1 - exp(- (gamma + mu[i])) 

# outflows due to infection, recovery or death
outflow_Sm[1:N_age, ] <- rbinom(Sm[i,j], prob = p_Sm[i,j])
outflow_S[1:N_age, ] <- rbinom(S[i,j], prob = p_S[i,j])
outflow_I[1:N_age, ] <- rbinom(I[i,j], prob = p_I[i])
outflow_R[1:N_age, ] <- rbinom(R[i,j], prob = p_R[i])
outflow_S2[1:N_age, ] <- rbinom(S2[i,j], prob = p_S2[i,j])
outflow_I2[1:N_age, ] <- rbinom(I2[i,j], prob = p_I2[i])


###############################################################################################################
# INFLOWS
# 1. where do these outflows go? 
# 2. births
# 3. importations
###############################################################################################################

##################################################
## new infections, recoveries and waning events ##
##################################################

#normalising the probabilities 
norm_p_sigma_m[1:N_age, ] <- p_sigma_m/(p_sigma_m + p_infection_mAb[j] + p_mu[i])
norm_p_infection_mAb[1:N_age, ] <- p_infection_mAb[j]/(p_sigma_m + p_infection_mAb[j] + p_mu[i])
norm_p_infection[1:N_age, ] <- p_infection[j]/(p_infection[j] + p_mu[i])
norm_p_reinfection[1:N_age, ] <- p_reinfection[j]/(p_reinfection[j] + p_mu[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i])

# number of new infections, recoveries and newly susceptible
new_waned_mAb[1:N_age, ] <- rbinom(outflow_Sm[i,j], prob = norm_p_sigma_m[i,j])
new_infections_mAb[1:N_age, ] <- rbinom(outflow_Sm[i,j], prob = norm_p_infection_mAb[i,j])
new_infections[1:N_age, ] <- rbinom(outflow_S[i,j], prob = norm_p_infection[i,j])
new_recoveries[1:N_age, ] <- rbinom(outflow_I[i,j], prob = norm_p_gamma[i])
new_waned[1:N_age, ] <- rbinom(outflow_R[i,j], prob = norm_p_sigma[i])
new_reinfections[1:N_age, ] <- rbinom(outflow_S2[i,j], prob = norm_p_reinfection[i,j])
new_recoveries_two[1:N_age, ] <- rbinom(outflow_I2[i,j], prob = norm_p_gamma[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)

N_0 <- user(1000) # user-defined initial population size per patch

birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial patch population size
new_births[] <- rpois(birth_rate) #per day

births_protected[] <- rbinom(new_births[i], prob = seropoz_A4[i]) # protected by mAbs
births_not_protected[] <- new_births[i] - births_protected[i] #  NOT protected by mAbs

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
new_Sm[1,] <- Sm[1,j] - outflow_Sm[1,j] + births_protected[j]
new_Sm[2:N_age, ] <- Sm[i,j] - outflow_Sm[i,j]

new_S[1, ] <- S[1,j] - outflow_S[1,j] + new_waned_mAb[1,j] + births_not_protected[j] 
new_S[2:N_age, ] <- S[i,j] - outflow_S[i,j] + new_waned_mAb[i,j]

new_I[1:24, ] <-  I[i,j] - outflow_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j]
new_I[25, 13] <- I[i,j] - outflow_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j] + imported_cases
new_I[25, 1:12] <- I[i,j] - outflow_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j]
new_I[25, 14:N_patch] <- I[i,j] - outflow_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j]
new_I[26:N_age, ] <-  I[i,j] - outflow_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j]

new_R[1:N_age, ] <- R[i,j] - outflow_R[i,j] + new_recoveries[i,j] + new_recoveries_two[i,j]
new_S2[1:N_age, ] <- S2[i,j] - outflow_S2[i,j] + new_waned[i,j]
new_I2[1:N_age, ] <- I2[i,j] - outflow_I2[i,j] + new_reinfections[i,j]


## STEP 2 update with ageing

update(Sm[1, ]) <- if(tt %% 30 == 0) 0 else new_Sm[1,j]
update(Sm[2:48, ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j] else new_Sm[i,j]
update(Sm[N_age, ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j] + new_Sm[i,j] else new_Sm[i,j]

update(S[1, ]) <- if(tt %% 30 == 0) 0 else new_S[1,j]
update(S[2:48, ]) <- if(tt %% 30 == 0) new_S[i-1,j] else new_S[i,j]
update(S[N_age, ]) <- if(tt %% 30 == 0) new_S[i-1,j] + new_S[i,j] else new_S[i,j]

update(I[1, ]) <-  if(tt %% 30 == 0) 0 else new_I[1,j]
update(I[2:24, ]) <- if(tt %% 30 == 0) new_I[i-1,j] else new_I[i,j]
update(I[25, 13]) <- if(tt %% 30 == 0) new_I[i-1,j] else if(tt == imp_t[1] || tt == imp_t[2] || tt == imp_t[3] || tt == imp_t[4] || tt == imp_t[5]) 5 + new_I[i,j] else new_I[i,j]
update(I[25, 1:12]) <- if(tt %% 30 == 0) new_I[i-1,j] else new_I[i,j]
update(I[25, 14:N_patch]) <- if(tt %% 30 == 0) new_I[i-1,j] else new_I[i,j]
update(I[26:48, ]) <- if(tt %% 30 == 0) new_I[i-1,j] else new_I[i,j]
update(I[N_age, ]) <- if(tt %% 30 == 0) new_I[i-1,j] + new_I[i,j] else new_I[i,j]

update(R[1, ]) <- if(tt %% 30 == 0) 0 else new_R[1, j]
update(R[2:48, ]) <- if(tt %% 30 == 0) new_R[i - 1, j] else new_R[i, j]
update(R[N_age, ]) <- if(tt %% 30 == 0) new_R[i - 1, j] + new_R[i, j] else new_R[i, j]

update(S2[1, ]) <- if(tt %% 30 == 0) 0 else new_S2[1, j]
update(S2[2:48, ]) <- if(tt %% 30 == 0) new_S2[i - 1, j] else new_S2[i, j]
update(S2[N_age, ]) <- if(tt %% 30 == 0) new_S2[i - 1, j] + new_S2[i, j] else new_S2[i, j]

update(I2[1, ]) <- if(tt %% 30 == 0) 0 else new_I2[1, j]
update(I2[2:48, ]) <- if(tt %% 30 == 0) new_I2[i - 1, j] else new_I2[i, j]
update(I2[N_age, ]) <- if(tt %% 30 == 0) new_I2[i - 1, j] + new_I2[i, j] else new_I2[i, j]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size for use in FOI
N[1:N_age, ] <- Sm[i,j] + S[i,j] + I[i,j] + R[i,j] + S2[i,j] + I2[i,j]
## record adult seroprevalence for use in births_protected
seropoz_A4[] <- (sum(S2[N_age,i]) + sum(I[N_age,i]) + sum(I2[N_age,i]) + sum(R[N_age,i]))/ (sum(Sm[N_age,i]) + sum(S[N_age,i]) + sum(S2[N_age,i]) + sum(I[N_age,i]) + sum(I2[N_age,i]) + sum(R[N_age,i]))

##################################################################################################################################
# initial conditions
# equilibrium solution approximated to speed up balancing of demography
# no equilibrium solution for infection at this point
##################################################################################################################################

## initial states

S_ini_p[ , ] <- user()

initial(Sm[1:N_age, ]) <- 0
initial(S[1:N_age, ]) <- S_ini_p[i,j]
initial(I[1:N_age, ]) <- 0
initial(R[1:N_age, ]) <- 0
initial(S2[1:N_age, ]) <- 0
initial(I2[1:N_age, ]) <- 0

initial(tt) <- 1

##################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- c(N_age, N_patch)
dim(p_S) <- c(N_age, N_patch)
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- c(N_age, N_patch)
dim(p_I2) <- N_age

dim(Sm) <- c(N_age, N_patch)
dim(S) <- c(N_age, N_patch)
dim(I) <- c(N_age, N_patch)
dim(R) <- c(N_age, N_patch)
dim(S2) <- c(N_age, N_patch)
dim(I2) <- c(N_age, N_patch)

dim(S_ini_p) <- c(N_age, N_patch)
dim(norm_p_infection) <- c(N_age, N_patch)
dim(norm_p_infection_mAb) <- c(N_age, N_patch)
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- c(N_age, N_patch)
dim(norm_p_reinfection) <- c(N_age, N_patch)

dim(new_waned_mAb) <- c(N_age, N_patch)
dim(new_infections) <- c(N_age, N_patch)
dim(new_infections_mAb) <- c(N_age, N_patch)
dim(new_recoveries) <- c(N_age, N_patch)
dim(new_waned) <- c(N_age, N_patch)
dim(new_reinfections) <- c(N_age, N_patch)
dim(new_recoveries_two) <- c(N_age, N_patch)

dim(outflow_Sm) <- c(N_age, N_patch)
dim(outflow_S) <- c(N_age, N_patch)
dim(outflow_I) <- c(N_age, N_patch)
dim(outflow_R) <- c(N_age, N_patch)
dim(outflow_S2) <- c(N_age, N_patch)
dim(outflow_I2) <- c(N_age, N_patch)

dim(new_Sm) <- c(N_age, N_patch)
dim(new_S) <- c(N_age, N_patch)
dim(new_I) <- c(N_age, N_patch)
dim(new_R) <- c(N_age, N_patch)
dim(new_S2) <- c(N_age, N_patch)
dim(new_I2) <- c(N_age, N_patch)

dim(N) <- c(N_age, N_patch)
#dim(I_patch) <- N_patch
dim(imp_t) <- 5

dim(new_births) <- N_patch
dim(births_protected) <- N_patch
dim(births_not_protected) <- N_patch
dim(seropoz_A4) <- N_patch
dim(rate_infection) <- N_patch
dim(rate_infection_mAb) <- N_patch
dim(rate_reinfection) <- N_patch
dim(p_infection) <- N_patch
dim(p_infection_mAb) <- N_patch
dim(p_reinfection) <- N_patch
dim(external_I1) <- N_patch
dim(external_I2) <- N_patch