#############################
## Single Patch Model 2019 ##
#############################

# This simple stochastic model simulates MERS-CoV transmission in a single, homogenously mixed, 
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is aproximated as 360 days 
# It is written using Odin
# It can then be run from a user-edited R script "single_patch_2019_run.R"

N_age <- 27 #number of age classes


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
beta_mAb <- mAb_susc * beta # infection rate in the presence of maternally acquired Abs
beta_Ab <- Ab_susc * beta # reinfection rate (infection rate in the presence of Ab protection)
beta_I2 <- reduced_shed * beta # rate at which second infections infect naive susceptible dromedaries
beta_I2_mAb <- mAb_susc * reduced_shed * beta # rate at which second infections infect mAb protected calves
beta_I2_Ab <- Ab_susc * reduced_shed * beta # rate at which second infections infect Ab protected animals

# frequency dependent rate of infection
# when reduced_shed = 1 I and I2 are essentially the same compartment
rate_infection <- beta * (sum(I[1:N_age]) / sum(N[1:N_age]))  + beta_I2 * (sum(I2[1:N_age]) / sum(N[1:N_age])) 
rate_infection_mAb <- beta_mAb * (sum(I[1:N_age]) / sum(N[1:N_age])) + beta_I2_mAb * (sum(I2[1:N_age]) / sum(N[1:N_age]))
rate_reinfection <- beta_Ab * (sum(I[1:N_age]) / sum(N[1:N_age])) + beta_I2_Ab * (sum(I2[1:N_age]) / sum(N[1:N_age]))

#####################
## mortality rates ##
#####################
# user-defined age-dependent mortality rate
mu_6m <- user(0.0005) # death rate for 1st 6 months of life
mu_7_12m <- user(0.0005) # death rate for 7-12th months of life
mu_2nd_yr <- user(0.0005) # death rate for 2nd year of life
mu_3rd_yr <- user(0.00025) # death rate for 3rd year of life
mu_4th_yr <- user(0.00025) # death rate for 4th year of life
mu_adult_over_4 <- user(0.00025) # death rate in adulthood (>4 years)
# expand these across the age strata
mu[1:6] <- mu_6m
mu[7:12] <- mu_7_12m
mu[13:24] <- mu_2nd_yr
mu[25] <- mu_3rd_yr
mu[26] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

###########################
## recovery rate I --> R ## where R is non-infectious and completely immune to further infection
###########################
gamma <- user(0.05) # (/day) 

###################################################
## rate at which complete immunity wanes R --> S ## is it even legit to have complete immunity for any period?
###################################################
sigma <- user(0.0005) # (/day) to be taken from catalytic work eventually

###############################################################
## rate at which maternally-acquired immunity wanes Sm --> S ##
###############################################################
sigma_m <- user(0.006) # (/day) to be taken from catalytic work eventually

##############################################################################################################
# CONVERTING THESE RATES --> PROBABILITIES
# the above boil down to 7 rates which are converted to probabilities below
##############################################################################################################
#p_alpha <- 1 - exp(-alpha) # prob birth
p_infection <- 1 - exp(-rate_infection)
p_infection_mAb <- 1 - exp(-rate_infection_mAb)
p_reinfection <- 1 - exp(-rate_reinfection)
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
# R2 = recovered & completely immune for a 2nd+ time (allows those infected >once to be counted - R would not)

# probability of leaving each compartment for any reason (other than ageing which is dealt with later)
p_Sm[1:N_age] <- 1 - exp(- (sigma_m + rate_infection_mAb + mu[i])) # probability of leaving Sm 
p_S[1:N_age] <- 1 - exp(- (rate_infection + mu[i])) # probability of leaving S
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i])) # probability of leaving I
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i])) # probability of leaving R

p_S2[1:N_age] <- 1 - exp(- (rate_reinfection + mu[i])) # probability of leaving S2
p_I2[1:N_age] <- 1 - exp(- (gamma + mu[i])) # probability of leaving I2
p_R2[1:N_age] <- 1 - exp(- (sigma + mu[i])) # probability of leaving R2

# outflows (due to infection, recovery or death - ageing within a disease state is dealt with seperately)
outflow_Sm[1:N_age] <- rbinom(Sm[i], prob = p_Sm[i])
outflow_S[1:N_age] <- rbinom(S[i], prob = p_S[i])
outflow_I[1:N_age] <- rbinom(I[i], prob = p_I[i])
outflow_R[1:N_age] <- rbinom(R[i], prob = p_R[i])
outflow_S2[1:N_age] <- rbinom(S2[i], prob = p_S2[i])
outflow_I2[1:N_age] <- rbinom(I2[i], prob = p_I2[i])
outflow_R2[1:N_age] <- rbinom(R2[i], prob = p_R2[i])

# outflows due to ageing within a disease state
# first 2 years of life are classed into 24 month wide age classes, followed 2 1-year wide age classes
# then in class 27 when camels are +4 there is no further ageing - they remain there until they die

aged_Sm[1:24] <- if(tt %% 30 == 0) Sm[i] - outflow_Sm[i] else 0 # ageing happens every month
aged_Sm[25:26] <- if(tt %% 360 == 0) Sm[i] - outflow_Sm[i] else 0 # ageing happens every year
aged_S[1:24] <- if(tt %% 30 == 0) S[i] - outflow_S[i] else 0 # ageing happens every month
aged_S[25:26] <- if(tt %% 360 == 0) S[i] - outflow_S[i] else 0 # ageing happens every year
aged_I[1:24] <- if(tt %% 30 == 0) I[i] - outflow_I[i] else 0 # ageing happens every month
aged_I[25:26] <- if(tt %% 360 == 0) I[i] - outflow_I[i] else 0 # ageing happens every year
aged_R[1:24] <- if(tt %% 30 == 0) R[i] - outflow_R[i] else 0 # ageing happens every month
aged_R[25:26] <- if(tt %% 360 == 0) R[i] - outflow_R[i] else 0 # ageing happens every year
aged_S2[1:24] <- if(tt %% 30 == 0) S2[i] - outflow_S2[i] else 0 # ageing happens every month
aged_S2[25:26] <- if(tt %% 360 == 0) S2[i] - outflow_S2[i] else 0 # ageing happens every year
aged_I2[1:24] <- if(tt %% 30 == 0) I2[i] - outflow_I2[i] else 0 # ageing happens every month
aged_I2[25:26] <- if(tt %% 360 == 0) I2[i] - outflow_I2[i] else 0 # ageing happens every year
aged_R2[1:24] <- if(tt %% 30 == 0) R2[i] - outflow_R2[i] else 0 # ageing happens every month
aged_R2[25:26] <- if(tt %% 360 == 0) R2[i] - outflow_R2[i] else 0 # ageing happens every year

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
norm_p_sigma_m[1:N_age] <- p_sigma_m/(p_sigma_m + p_infection_mAb + p_mu[i])
norm_p_infection_mAb[1:N_age] <- p_infection_mAb/(p_sigma_m + p_infection_mAb + p_mu[i])
norm_p_infection[1:N_age] <- p_infection/(p_infection + p_mu[i])
norm_p_reinfection[1:N_age] <- p_reinfection/(p_reinfection + p_mu[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i])

# number of new infections, recoveries and newly susceptible
new_waned_mAb[1:N_age] <- rbinom(outflow_Sm[i], prob = norm_p_sigma_m[i])
new_infections_mAb[1:N_age] <- rbinom(outflow_Sm[i], prob = norm_p_infection_mAb[i])
new_infections[1:N_age] <- rbinom(outflow_S[i], prob = norm_p_infection[i])
new_recoveries[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_gamma[i])
new_waned[1:N_age] <- rbinom(outflow_R[i], prob = norm_p_sigma[i])
new_reinfections[1:N_age] <- rbinom(outflow_S2[i], prob = norm_p_reinfection[i])
new_recoveries_2[1:N_age] <- rbinom(outflow_I2[i], prob = norm_p_gamma[i]) 
new_waned_2[1:N_age] <- rbinom(outflow_R2[i], prob = norm_p_sigma[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)

birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial population size
#birth_rate <- N_0 * p_alpha * (1 + (delta *(cos(3 * cos(pi * tt / 360))))) # N_0 is the intial population size
new_births <- rpois(birth_rate) #per day

births_protected <- 0 #rbinom(new_births, prob = seropoz_A4) # of those born, these will be protected by mAbs
births_not_protected <- new_births #new_births - births_protected # of those born, these will NOT be protected by mAbs

#########################
## importation process ##
#########################

importation_rate <- user(0.01) # should be proportional to population size?
imported_cases <- rpois(importation_rate) #per day
imp_t <- user() # a user defined time at which cases are imported if you want a one off importation event

###############################################################################################################
# EQUATIONS for movement of individuals between disease states and age classes
# time-step = 1 day
###############################################################################################################

## tt %% 30 == 0 identifies the timesteps where ageing occurs in the first 24 age classes.
## when ageing occurs we need to allow ageing and change in disease states to happen in a single time-step
## -aged_I etc. in the else clause is actually redundant as it =0 except at ageing timepoints
## all importations (whether using rate or pulse) occur into age class 25 (2-3 years old)

# no need for tt %% here as no inflows
update(Sm[1]) <- Sm[1] - outflow_Sm[1] - aged_Sm[1] + births_protected
update(Sm[2:26]) <- Sm[i] - outflow_Sm[i] - aged_Sm[i] + aged_Sm[i-1]
update(Sm[N_age]) <- Sm[N_age] - outflow_Sm[i] + aged_Sm[(N_age - 1)]

#26 is separate because the width of 25 is 360 rather than 30
update(S[1]) <- if(tt %% 30 == 0) S[1] - outflow_S[1] - aged_S[1] + births_not_protected else S[1] - outflow_S[1] - aged_S[1] + new_waned_mAb[1] + births_not_protected 
update(S[2:25]) <- if(tt %% 30 == 0) S[i] - outflow_S[i] - aged_S[i] + aged_S[i-1] + new_waned_mAb[i-1] else S[i] - outflow_S[i] - aged_S[i] + aged_S[i-1] + new_waned_mAb[i]
update(S[26]) <- if(tt %% 360 == 0) S[26] - outflow_S[26] - aged_S[26] + aged_S[25] + new_waned_mAb[25] else S[26] - outflow_S[26] - aged_S[26] + aged_S[25] + new_waned_mAb[26]
update(S[N_age]) <- if(tt %% 360 == 0) S[N_age] - outflow_S[N_age] + aged_S[(N_age - 1)] + new_waned_mAb[N_age - 1] else S[N_age] - outflow_S[N_age] + aged_S[(N_age - 1)] + new_waned_mAb[N_age]

#26 also has imported cases going in either daily or on a set day (day set to not coincide with ageing)
update(I[1]) <-  if (tt %% 30 == 0) I[1] - outflow_I[1] - aged_I[1] else I[1] - outflow_I[1] - aged_I[1] + new_infections[1] + new_infections_mAb[1]
update(I[2:25]) <- if(tt %% 30 == 0) I[i] - outflow_I[i] - aged_I[i] + new_infections[i - 1] + new_infections_mAb[i - 1] + aged_I[i - 1] else I[i] - outflow_I[i] - aged_I[i] + new_infections[i] + new_infections_mAb[i] + aged_I[i - 1]
update(I[26]) <- if(tt %% 360 == 0) I[26] - outflow_I[26] - aged_I[26] + new_infections[25] + new_infections_mAb[25] + aged_I[25] + imported_cases else if(tt == imp_t) 1 + I[26] - outflow_I[26] - aged_I[26] + new_infections[26] + new_infections_mAb[26] + aged_I[25] + imported_cases else I[26] - outflow_I[26] - aged_I[26] + new_infections[26] + new_infections_mAb[26] + aged_I[25] + imported_cases
update(I[N_age]) <- if(tt %% 360 == 0) I[N_age] - outflow_I[N_age] + sum(new_infections[26:N_age]) + sum(new_infections_mAb[26:N_age]) + aged_I[(N_age - 1)] else I[N_age] - outflow_I[N_age] + new_infections[N_age] + new_infections_mAb[N_age] + aged_I[(N_age - 1)]

update(R[1]) <- if(tt %% 30 == 0) R[1] - outflow_R[1] - aged_R[1] else R[1] - outflow_R[1] - aged_R[1] + new_recoveries[1]
update(R[2:25]) <- if(tt %% 30 == 0) R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i - 1] + aged_R[i - 1] else R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[26]) <- if(tt %% 360 == 0) R[26] - outflow_R[26] - aged_R[26] + new_recoveries[25] + aged_R[25] else R[26] - outflow_R[26] - aged_R[26] + new_recoveries[26] + aged_R[25]
update(R[N_age]) <- if(tt %% 360 == 0) R[N_age] - outflow_R[N_age] + sum(new_recoveries[26:N_age]) + aged_R[(N_age - 1)] else R[N_age] - outflow_R[N_age] + new_recoveries[N_age] + aged_R[(N_age - 1)]

update(S2[1]) <- if(tt %% 30 == 0) S2[1] - outflow_S2[1] - aged_S2[1] else S2[1] - outflow_S2[1] - aged_S2[1] + new_waned[1] + new_waned_2[1]
update(S2[2:25]) <- if(tt %% 30 == 0) S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i - 1] + new_waned[i - 1] + new_waned_2[i - 1] else S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i - 1] + new_waned[i] + new_waned_2[i]
update(S2[26]) <- if(tt %% 360 == 0)S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i - 1] + new_waned[i - 1] + new_waned_2[i - 1] else S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i - 1] + new_waned[i] + new_waned_2[i]
update(S2[N_age]) <- if(tt %% 360 == 0) S2[N_age] - outflow_S2[N_age] + aged_S2[(N_age - 1)] + sum(new_waned[26:N_age]) + sum(new_waned_2[26:N_age]) else S2[N_age] - outflow_S2[N_age] + aged_S2[(N_age - 1)] + new_waned[N_age] + new_waned_2[N_age]

update(I2[1]) <- if(tt %% 30 == 0) I2[1] - outflow_I2[1] - aged_I2[1] else I2[1] - outflow_I2[1] - aged_I2[1] + new_reinfections[1]
update(I2[2:25]) <- if(tt %% 30 == 0) I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i - 1] + aged_I2[i - 1] else I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i] + aged_I2[i - 1]
update(I2[26]) <- if(tt %% 360 == 0) I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i - 1] + aged_I2[i - 1] else I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i] + aged_I2[i - 1]
update(I2[N_age]) <- if(tt %% 360 == 0) I2[N_age] - outflow_I2[N_age] + sum(new_reinfections[26:N_age]) + aged_I2[(N_age - 1)] else I2[N_age] - outflow_I2[N_age] + new_reinfections[N_age] + aged_I2[(N_age - 1)]

update(R2[1]) <- if(tt %% 30 == 0) R2[1] - outflow_R2[1] - aged_R2[1] else R2[1] - outflow_R2[1] - aged_R2[1] + new_recoveries_2[1]
update(R2[2:25]) <- if(tt %% 30 == 0) R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i - 1] + aged_R2[i - 1] else R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i] + aged_R2[i - 1]
update(R2[26]) <- if(tt %% 360 == 0) R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i - 1] + aged_R2[i - 1] else R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i] + aged_R2[i - 1]
update(R2[N_age]) <- if(tt %% 360 == 0) R2[N_age] - outflow_R2[N_age] + sum(new_recoveries_2[26:N_age]) + aged_R2[(N_age - 1)] else R2[N_age] - outflow_R2[N_age] + new_recoveries_2[N_age] + aged_R2[(N_age - 1)]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size and seroprevalence in each age group, and in adults >4

update(N[1:N_age]) <- Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i] + R2[i]

update(seroprevalence[1:N_age]) <- (I[i] + R[i] + S2[i] + I2[i] + R2[i]) / (Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i] + R2[i])

update(seropoz_A4) <- (sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]) + sum(R2[N_age]))/ (sum(Sm[N_age]) + sum(S[N_age]) + sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]) + sum(R2[N_age]))

##################################################################################################################################
# initial conditions
# equilibrium solution approximated to speed up balancing of demography
# no equilibrium solution for infection at this point
##################################################################################################################################

## initial states

initial(Sm[1:N_age]) <- 0
initial(S[1:N_age]) <- S_ini_p[i] # will be user-defined
initial(I[1:N_age]) <- 0 # will be user-defined
initial(R[1:N_age]) <- 0
initial(S2[1:N_age]) <- 0
initial(I2[1:N_age]) <- 0
initial(R2[1:N_age]) <- 0

initial(tt) <- 1
initial(N[1:N_age]) <- S_ini_p[i]
initial(seroprevalence[1:N_age]) <- 0
initial(seropoz_A4) <- 0

## initial population size for use in birthrate

N_0 <- user(1000) # user-defined

## setting initial conditions using the equilibrium solution for age distribution
                           
births_detr[1:360] <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * i / 360)))) # change to fixed N_0 = huge to avoid NaNs
#births_detr[1:360] <- 10000000 * p_alpha * (1 + (delta * (cos(3 * cos(pi * i / 360))))) # change to fixed N_0 = huge to avoid NaNs
births_det[1:360] <- rpois(births_detr[i]) 


## if we start the model with the equilibrium amount in each of the first 24 month-wide compartments,
## and no camels in the 3rd year of life (they would have just moved into the 4th year comp)
## then from here camels will start filling the yr 3 compartment every month and then every year this will
## empty into the 4th year and then the adult compartment. 
## Birthrate will be set to balance summed death rate of this age distribution - how?

a_max <- 32 ## estimated max number of years spent in the last age class before death (only 1% of the population above after 32 yrs)
            ## (calculated in "calculating_a_max.R using exp decay mod)

S_ini[1] <- 0
S_ini[2] <- (sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * mu[1]))
S_ini[3] <- (sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:2])))
S_ini[4] <- (sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:3])))
S_ini[5] <- (sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:4])))
S_ini[6] <- (sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:5])))
S_ini[7] <- (sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:6])))
S_ini[8] <- (sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:7])))
S_ini[9] <- (sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:8])))
S_ini[10] <- (sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:9])))
S_ini[11] <- (sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:10])))
S_ini[12] <- (sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:11])))
S_ini[13] <- sum(births_det[1:(360 - 330)]) * exp(- (30 * sum(mu[1:12])))
S_ini[14] <- (sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * sum(mu[1:13])))
S_ini[15] <- (sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:14])))
S_ini[16] <- (sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:15])))
S_ini[17] <- (sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:16])))
S_ini[18] <- (sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:17])))
S_ini[19] <- (sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:18])))
S_ini[20] <- (sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:19])))
S_ini[21] <- (sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:20])))
S_ini[22] <- (sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:21])))
S_ini[23] <- (sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:22])))
S_ini[24] <- (sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:23])))
S_ini[25] <- (sum(births_det[1:(360-330)]) * exp(-30*sum(mu[1:24])))
S_ini[26] <- sum(births_det[1:360]) *  exp(-(30*(sum(mu[1:24])) + 360*mu[25]))

yearly_influx <- sum(births_det[1:360]) * exp(-(30 * (sum(mu[1:24])) + 360*sum(mu[25:26]))) #annual inflow into the final age compartment
yr[1:a_max] <- i  # number of years of influx before max age expectancy reached
# camels remaining from each cohort to enter in the last 32 years influx that remain
cohort_remaining[1:a_max] <- yearly_influx * exp(-360*yr[i]*mu[N_age]) 
S_ini[N_age] <- sum(cohort_remaining[1:a_max])

#S_ini[N_age] <- (a_max * (sum(births_det[1:360]) * exp(- ((30 * sum(mu[1:24])) + 360 * sum(mu[25:26]))))) * exp(-(a_max*360*mu[27]))


# getting proportion of animals in each age-class, multiplying by N_0 and rounding to whole animals

S_ini_p[1:N_age] <- round((S_ini[i] / sum(S_ini[1:N_age])) * N_0)



################################################################################################################################

# OUTPUTS

################################################################################################################################

#########################################
## number of individuals in each state ##
#########################################

# total number of individuals still protected by maternal immunity
output(M) <- sum(Sm[1:N_age])

output(S_1) <- sum(S[1:N_age]) # susceptible individuals never infected
output(I_1) <- sum(I[1:N_age]) # individuals infectious for the 1st time
output(R_1) <- sum(R[1:N_age]) # individuals recovered from a 1st infection
output(S_2) <- sum(S2[1:N_age]) # susceptible individuals whose immunity has waned
output(I_2) <- sum(I2[1:N_age]) # individuals infectious for the 2nd+ time
output(R_2) <- sum(R2[1:N_age]) # individuals recovered from 2nd+ infections

output(Stot) <- sum(S[1:N_age]) + sum(S2[1:N_age]) + sum(Sm[1:N_age]) # total number of susceptible individuals
output(Itot) <- sum(I[1:N_age]) + sum(I2[1:N_age]) # total number of infectious individuals
output(Rtot) <- sum(R[1:N_age]) + sum(R2[1:N_age]) # total number of recovered individuals 

output(Ntot) <- sum(N[1:N_age]) # total number of individuals
output(N_J) <- sum(Sm[1:24]) + sum(S[1:24]) + sum(S2[1:24]) + sum(I[1:24]) + sum(I2[1:24]) + sum(R[1:24]) + sum(R2[1:24])
output(N_A) <- sum(Sm[25:N_age]) + sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age])

####################
## seroprevalence ##
####################
output(seropoz_A) <- 100 * (sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))/(sum(Sm[25:N_age]) + sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))
output(seropoz_J) <- 100 * (sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24])) / (sum(Sm[1:24]) + sum(S[1:24]) + sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24]))
output(seropz_tot) <- 100 * (sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))/(sum(Sm[1:N_age]) + sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))

############################
## age at first infection ##
############################
output(reinf_1) <- new_reinfections[1]
output(reinf_2) <- new_reinfections[2]
output(incidence_new_inf) <- sum(new_infections[1:N_age])
output(incidence_indig_inf) <- sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age])
output(total_incidence) <- (sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age]))

###########
## other ##
###########
output(inf) <- rate_infection
output(birthrate) <- birth_rate
output(births) <- new_births
output(importations) <- imported_cases
output(yy) <- yr[12]

################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- N_age
dim(p_S) <- N_age
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- N_age
dim(p_I2) <- N_age
dim(p_R2) <- N_age

dim(Sm) <- N_age
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S2) <- N_age
dim(I2) <- N_age
dim(R2) <- N_age

dim(S_ini) <- N_age
dim(S_ini_p) <- N_age
dim(yr) <- a_max
dim(cohort_remaining) <- a_max
dim(norm_p_infection) <- N_age
dim(norm_p_infection_mAb) <- N_age
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- N_age
dim(norm_p_reinfection) <- N_age

dim(new_waned_mAb) <- N_age
dim(new_infections) <- N_age
dim(new_infections_mAb) <- N_age
dim(new_recoveries) <- N_age
dim(new_waned) <- N_age
dim(new_reinfections) <- N_age
dim(new_recoveries_2) <- N_age
dim(new_waned_2) <- N_age

dim(aged_Sm) <- 26
dim(aged_S) <- 26
dim(aged_I) <- 26
dim(aged_R) <- 26
dim(aged_S2) <- 26
dim(aged_I2) <- 26
dim(aged_R2) <- 26

dim(outflow_Sm) <- N_age
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age
dim(outflow_S2) <- N_age
dim(outflow_I2) <- N_age
dim(outflow_R2) <- N_age
dim(births_det) <- 360
dim(births_detr) <- 360
dim(seroprevalence) <- N_age
dim(N) <- N_age
