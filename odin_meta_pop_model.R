N_age <- 27 ## number of age categories
n_patch <- 25

## RATES ##
alpha <- user(0.0005) # birth rate
beta <- user(0.3) # infection rate
mAb_susc <- user(0) # proportion of susceptibility experienced if mAbs present
Ab_susc <- user(0) # proportion of susceptibility experienced if previoulsy infected. default = 0
reduced_shed <- user(1) # proportion of shedding seen in reinfections default = 1 (no difference)
beta_mAb <- mAb_susc * beta # reinfection rate (infection rate in the presence of Ab protection)
beta_Ab <- Ab_susc * beta
beta_I2 <- reduced_shed * beta
beta_I2_mAb <- mAb_susc * reduced_shed * beta
beta_I2_Ab <- Ab_susc * reduced_shed * beta



## INCLUDING INFECTION FROM CONNECTED PATCHES ##

connectivity <- 0.05 ##conection strength between patches (0<x<1)

## FOI from first-time infections in other patches
FOI[1] <- ((sum(I[1:N_age,2])/N[2]) + (sum(I[1:N_age,6])/N[6]))
FOI[2] <- ((sum(I[1:N_age,1])/N[1]) + (sum(I[1:N_age,3])/N[3]) + (sum(I[1:N_age,7])/N[7]))
FOI[3] <- ((sum(I[1:N_age,2])/N[2]) + (sum(I[1:N_age,4])/N[4]) + (sum(I[1:N_age,8])/N[8]))
FOI[4] <- ((sum(I[1:N_age,3])/N[3]) + (sum(I[1:N_age,5])/N[5]) + (sum(I[1:N_age,9])/N[9]))
FOI[5] <- ((sum(I[1:N_age,4])/N[4]) + (sum(I[1:N_age,10])/N[10]))
FOI[6] <- ((sum(I[1:N_age,1])/N[1]) + (sum(I[1:N_age,7])/N[7]) + (sum(I[1:N_age,11])/N[11]))
FOI[7] <- ((sum(I[1:N_age,2])/N[2]) + (sum(I[1:N_age,6])/N[6]) + (sum(I[1:N_age,8])/N[8]) + (sum(I[1:N_age, 12])/N[12]))
FOI[8] <- ((sum(I[1:N_age,3])/N[3]) + (sum(I[1:N_age,7])/N[7]) + (sum(I[1:N_age,9])/N[9]) + (sum(I[1:N_age, 13])/N[13]))
FOI[9] <- ((sum(I[1:N_age,4])/N[4]) + (sum(I[1:N_age,8])/N[8]) + (sum(I[1:N_age,10])/N[10]) + (sum(I[1:N_age, 14])/N[14]))
FOI[10] <- ((sum(I[1:N_age,5])/N[5]) + (sum(I[1:N_age,9])/N[9]) + (sum(I[1:N_age,15])/N[15]))
FOI[11] <- ((sum(I[1:N_age,6])/N[6]) + (sum(I[1:N_age,12])/N[12]) + (sum(I[1:N_age,16])/N[16]))
FOI[12] <- ((sum(I[1:N_age,7])/N[7]) + (sum(I[1:N_age,11])/N[11]) + (sum(I[1:N_age,13])/N[13]) + (sum(I[1:N_age, 17])/N[17]))
FOI[13] <- ((sum(I[1:N_age,8])/N[8]) + (sum(I[1:N_age,12])/N[12]) + (sum(I[1:N_age,14])/N[14]) + (sum(I[1:N_age, 18])/N[18]))
FOI[14] <- ((sum(I[1:N_age,9])/N[9]) + (sum(I[1:N_age,13])/N[13]) + (sum(I[1:N_age,15])/N[15]) + (sum(I[1:N_age, 19])/N[19]))
FOI[15] <- ((sum(I[1:N_age,10])/N[10]) + (sum(I[1:N_age,14])/N[14]) + (sum(I[1:N_age,20])/N[20]))
FOI[16] <- ((sum(I[1:N_age,11])/N[11]) + (sum(I[1:N_age,17])/N[17]) + (sum(I[1:N_age,21])/N[21]))
FOI[17] <- ((sum(I[1:N_age,12])/N[12]) + (sum(I[1:N_age,16])/N[16]) + (sum(I[1:N_age,18])/N[18]) + (sum(I[1:N_age, 22])/N[22]))
FOI[18] <- ((sum(I[1:N_age,13])/N[13]) + (sum(I[1:N_age,17])/N[17]) + (sum(I[1:N_age,19])/N[19]) + (sum(I[1:N_age, 23])/N[23]))
FOI[19] <- ((sum(I[1:N_age,14])/N[14]) + (sum(I[1:N_age,18])/N[18]) + (sum(I[1:N_age,20])/N[20]) + (sum(I[1:N_age, 24])/N[24]))
FOI[20] <- ((sum(I[1:N_age,15])/N[15]) + (sum(I[1:N_age,19])/N[19]) + (sum(I[1:N_age,25])/N[25]))
FOI[21] <- ((sum(I[1:N_age,16])/N[16]) + (sum(I[1:N_age,22])/N[22]))
FOI[22] <- ((sum(I[1:N_age,17])/N[17]) + (sum(I[1:N_age,21])/N[21]) + (sum(I[1:N_age,23])/N[23]))
FOI[23] <- ((sum(I[1:N_age,18])/N[18]) + (sum(I[1:N_age,22])/N[22]) + (sum(I[1:N_age,24])/N[24]))
FOI[24] <- ((sum(I[1:N_age,19])/N[19]) + (sum(I[1:N_age,23])/N[23]) + (sum(I[1:N_age,25])/N[25]))
FOI[25] <- ((sum(I[1:N_age,20])/N[20]) + (sum(I[1:N_age,24])/N[24]))

## FOI from second-plus-time infections in other patches
FOI_2[1] <- ((sum(I2[1:N_age,2])/N[2]) + (sum(I2[1:N_age,6])/N[6]))
FOI_2[2] <- ((sum(I2[1:N_age,1])/N[1]) + (sum(I2[1:N_age,3])/N[3]) + (sum(I2[1:N_age,7])/N[7]))
FOI_2[3] <- ((sum(I2[1:N_age,2])/N[2]) + (sum(I2[1:N_age,4])/N[4]) + (sum(I2[1:N_age,8])/N[8]))
FOI_2[4] <- ((sum(I2[1:N_age,3])/N[3]) + (sum(I2[1:N_age,5])/N[5]) + (sum(I2[1:N_age,9])/N[9]))
FOI_2[5] <- ((sum(I2[1:N_age,1])/N[1]) + (sum(I2[1:N_age,3])/N[3]) + (sum(I2[1:N_age,7])/N[7]))
FOI_2[6] <- ((sum(I2[1:N_age,4])/N[4]) + (sum(I2[1:N_age,10])/N[10]))
FOI_2[7] <- ((sum(I2[1:N_age,1])/N[1]) + (sum(I2[1:N_age,7])/N[7]) + (sum(I2[1:N_age,11])/N[11]))
FOI_2[8] <- ((sum(I2[1:N_age,2])/N[2]) + (sum(I2[1:N_age,6])/N[6]) + (sum(I2[1:N_age,8])/N[8]) + (sum(I2[1:N_age, 12])/N[12]))
FOI_2[9] <- ((sum(I2[1:N_age,3])/N[3]) + (sum(I2[1:N_age,7])/N[7]) + (sum(I2[1:N_age,9])/N[9]) + (sum(I2[1:N_age, 13])/N[13]))
FOI_2[10] <- ((sum(I2[1:N_age,4])/N[4]) + (sum(I2[1:N_age,8])/N[8]) + (sum(I2[1:N_age,10])/N[10]) + (sum(I2[1:N_age, 14])/N[14]))
FOI_2[11] <- ((sum(I2[1:N_age,5])/N[5]) + (sum(I2[1:N_age,9])/N[9]) + (sum(I2[1:N_age,15])/N[15]))
FOI_2[12] <- ((sum(I2[1:N_age,6])/N[6]) + (sum(I2[1:N_age,12])/N[12]) + (sum(I2[1:N_age,16])/N[16]))
FOI_2[13] <- ((sum(I2[1:N_age,7])/N[7]) + (sum(I2[1:N_age,11])/N[11]) + (sum(I2[1:N_age,13])/N[13]) + (sum(I2[1:N_age, 17])/N[17]))
FOI_2[14] <- ((sum(I2[1:N_age,8])/N[8]) + (sum(I2[1:N_age,12])/N[12]) + (sum(I2[1:N_age,14])/N[14]) + (sum(I2[1:N_age, 18])/N[18]))
FOI_2[15] <- ((sum(I2[1:N_age,9])/N[9]) + (sum(I2[1:N_age,13])/N[13]) + (sum(I2[1:N_age,15])/N[15]) + (sum(I2[1:N_age, 19])/N[19]))
FOI_2[16] <- ((sum(I2[1:N_age,10])/N[10]) + (sum(I2[1:N_age,14])/N[14]) + (sum(I2[1:N_age,20])/N[20]))
FOI_2[17] <- ((sum(I2[1:N_age,11])/N[11]) + (sum(I2[1:N_age,17])/N[17]) + (sum(I2[1:N_age,21])/N[21]))
FOI_2[18] <- ((sum(I2[1:N_age,12])/N[12]) + (sum(I2[1:N_age,16])/N[16]) + (sum(I2[1:N_age,18])/N[18]) + (sum(I2[1:N_age, 22])/N[22]))
FOI_2[19] <- ((sum(I2[1:N_age,13])/N[13]) + (sum(I2[1:N_age,17])/N[17]) + (sum(I2[1:N_age,19])/N[19]) + (sum(I2[1:N_age, 23])/N[23]))
FOI_2[20] <- ((sum(I2[1:N_age,14])/N[14]) + (sum(I2[1:N_age,18])/N[18]) + (sum(I2[1:N_age,20])/N[20]) + (sum(I2[1:N_age, 24])/N[24]))
FOI_2[21] <- ((sum(I2[1:N_age,15])/N[15]) + (sum(I2[1:N_age,19])/N[19]) + (sum(I2[1:N_age,25])/N[25]))
FOI_2[22] <- ((sum(I2[1:N_age,16])/N[16]) + (sum(I2[1:N_age,22])/N[22]))
FOI_2[23] <- ((sum(I2[1:N_age,17])/N[17]) + (sum(I2[1:N_age,21])/N[21]) + (sum(I2[1:N_age,23])/N[23]))
FOI_2[24] <- ((sum(I2[1:N_age,18])/N[18]) + (sum(I2[1:N_age,22])/N[22]) + (sum(I2[1:N_age,24])/N[24]))
FOI_2[25] <- ((sum(I2[1:N_age,20])/N[20]) + (sum(I2[1:N_age,24])/N[24]))


rate_infection[1:n_patch] <- beta * (sum(I[1:N_age,i]) / N[i])  + 
  beta_I2 * (sum(I2[1:N_age, i]) / N[i]) + 
  (beta*connectivity*FOI[i]) +
  (beta_I2*connectivity*FOI_2[i])

rate_infection_mAb[1:n_patch] <- beta_mAb * (sum(I[1:N_age,i]) / N[i])  + 
  beta_I2_mAb * (sum(I2[1:N_age, i]) / N[i]) + 
  (beta_mAb*connectivity*FOI[i]) +
  (beta_I2_mAb*connectivity*FOI_2[i])

rate_reinfection[1:n_patch] <- beta_Ab * (sum(I[1:N_age,i]) / N[i])  + 
  beta_I2_Ab * (sum(I2[1:N_age, i]) / N[i]) + 
  (beta_Ab*connectivity*FOI[i]) +
  (beta_I2_Ab*connectivity*FOI_2[i])

# 
# rate_infection <- beta * (sum(I[1:N_age]) / N)  + beta_I2 * (sum(I2[1:N_age]) / N)
# rate_infection_mAb <- beta_mAb * (sum(I[1:N_age]) / N) + beta_I2_mAb * (sum(I2[1:N_age]) / N)
# rate_reinfection <- beta_Ab * (sum(I[1:N_age]) / N) + beta_I2_Ab * (sum(I2[1:N_age]) / N)

mu_6m <- user(0.005) # death rate for 1st 6 months of life, user-defined, default = 0.005
mu_2y <- user(0.001) # death rate for the rest of the 1st 2 yrs of life
mu_adult <- user(0.0005) # death rate in adulthood (2 yrs +)
mu[1:6] <- mu_6m
mu[7:24] <- mu_2y
mu[25:27] <- mu_adult

gamma <- user(0.05) # recovery rate
sigma <- user(0.0005) ## waning immunity

## converting rates to probabilities

p_alpha <- 1 - exp(-alpha)
p_infection[] <- 1 - exp(-rate_infection[i])
p_infection_mAb[] <- 1 - exp(-rate_infection_mAb[i])
p_reinfection[] <- 1 - exp(-rate_reinfection[i])
p_mu[] <- 1 - exp(-mu[i])
p_gamma <- 1 - exp(-gamma)
p_sigma <- 1 - exp(-sigma)

p_Sm[,] <- 1 - exp(- (rate_infection_mAb[j] + mu[i])) # probability of leaving Sm for those <6m protected by maternal Abs
p_S[,] <- 1 - exp(- (rate_infection[j] + mu[i])) # probability of leaving S (with the exception of through ageing)
p_I[] <- 1 - exp(- (gamma + mu[i])) # probability of leaving I (with the exception of through ageing)
p_R[] <- 1 - exp(- (sigma + mu[i])) # probability of leaving R (with the exception of through ageing)

p_S2[,] <- 1 - exp(- (rate_reinfection[j] + mu[i]))
p_I2[] <- 1 - exp(- (gamma + mu[i]))
p_R2[] <- 1 - exp(- (sigma + mu[i]))


# BIRTHS

## birth process: any new individual will enter 'S[1]'

delta <- user() #modulates the seasonality of births
pi <- 3.14159
N_0 <- user(1000) # initial population size
birth_rate <- N_0 * alpha * (1 + delta *(cos(3 * cos(pi * tt / 360))))
new_births <- rpois(birth_rate) 

# Determining how many of the newborns have maternal antibody protection 
# (as a function of adult >4 yr seroprevalence)

births_protected[] <- rbinom(new_births, prob = seropoz_A4[i]) # of those born, these will be protected by maternal antibodies
births_not_protected[] <- new_births - births_protected[i] # of those born, these will NOT be protected by maternal antibodies

# IMPORTATION OF INFECTION

## importation process
# importation_rate <- user(0.01)
# imported_cases <- rpois(importation_rate)

# STATE TRANSITIONS

# STEP 1 - HOW MANY INDIVIDUALS LEAVE EACH COMPARTMENT?

# outflows (due to infection, recovery or death - ageing is dealt with seperately)
outflow_Sm[,] <- rbinom(Sm[i,j], prob = p_Sm[i,j])
outflow_S[,] <- rbinom(S[i,j], prob = p_S[i,j])
outflow_I[,] <- rbinom(I[i,j], prob = p_I[i])
outflow_R[,] <- rbinom(R[i,j], prob = p_R[i])
outflow_S2[,] <- rbinom(S2[i,j], prob = p_S2[i,j])
outflow_I2[,] <- rbinom(I2[i,j], prob = p_I2[i])
outflow_R2[,] <- rbinom(R2[i,j], prob = p_R2[i])

# STEP 2 - OUT OF THOSE THAT HAVE LEFT, HOW MANY MOVED DISEASE STATE AND HOW MANY DIED?

#normalising the probabilities 
norm_p_infection_mAb[,] <- p_infection_mAb[j]/(p_infection_mAb[j] + p_mu[i])
norm_p_infection[,] <- p_infection[j]/(p_infection[j] + p_mu[i])
norm_p_reinfection[,] <- p_reinfection[j]/(p_reinfection[j] + p_mu[i])
norm_p_gamma[] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[] <- p_sigma/(p_sigma + p_mu[i])

new_infections[,] <- rbinom(outflow_S[i,j], prob = norm_p_infection[i,j])
new_infections_mAb[1:6,] <- rbinom(outflow_Sm[i,j], prob = norm_p_infection_mAb[i,j])
new_recoveries[,] <- rbinom(outflow_I[i,j], prob = norm_p_gamma[i])
new_waned[,] <- rbinom(outflow_R[i,j], prob = norm_p_sigma[i])
new_reinfections[,] <- rbinom(outflow_S2[i,j], prob = norm_p_reinfection[i,j])
new_recoveries_i2[,] <- rbinom(outflow_I2[i,j], prob = norm_p_gamma[i]) 
new_waned_i2[,] <- rbinom(outflow_R2[i,j], prob = norm_p_sigma[i])

# AGEING 

# number of individuals leaving each compartment through ageing alone (but remaining in S or I or R)
# individuals which have entered age compartment 27 remain there until they die

aged_Sm[1:6,] <- if(tt %% 30 == 0) Sm[i,j] - outflow_Sm[i,j] else 0

aged_S[1:24,] <- if(tt %% 30 == 0) S[i,j] - outflow_S[i,j] else 0

aged_S[25:26,] <- if(tt %% 360 == 0) S[i,j] - outflow_S[i,j] else 0 #13th comp has width 1 yr

aged_I[1:24,] <- if(tt %% 30 == 0) I[i,j] - outflow_I[i,j] else 0

aged_I[25:26,] <- if(tt %% 360 == 0) I[i,j] - outflow_I[i,j] else 0 #13th comp has width 1 yr

aged_R[1:24,] <- if(tt %% 30 == 0) R[i,j] - outflow_R[i,j] else 0

aged_R[25:26,] <- if(tt %% 360 == 0) R[i,j] - outflow_R[i,j] else 0 #13th comp has width 1 yr

aged_S2[1:24,] <- if(tt %% 30 == 0) S2[i,j] - outflow_S2[i,j] else 0

aged_S2[25:26,] <- if(tt %% 360 == 0) S2[i,j] - outflow_S2[i,j] else 0 #13th comp has width 1 yr

aged_I2[1:24,] <- if(tt %% 30 == 0) I2[i,j] - outflow_I2[i,j] else 0

aged_I2[25:26,] <- if(tt %% 360 == 0) I2[i,j] - outflow_I2[i,j] else 0 #13th comp has width 1 yr

aged_R2[1:24,] <- if(tt %% 30 == 0) R2[i,j] - outflow_R2[i,j] else 0

aged_R2[25:26,] <- if(tt %% 360 == 0) R2[i,j] - outflow_R2[i,j] else 0 #13th comp has width 1 yr


## equations for transitions between compartments
## S[1:6] = individuals protected by maternal antibodies (<6 months of age)
## S[7:N_age] = fully susceptible individuals
## I = infectious individuals
## R = recovered individuals

# imp_t <- user()

update(Sm[1,]) <- Sm[1,j] - outflow_Sm[1,j] - aged_Sm[1,j] + births_protected[j]
update(Sm[2:6,]) <- Sm[i,j] - outflow_Sm[i,j] - aged_Sm[i,j] + aged_Sm[i-1,j]


update(S[1,]) <- S[1,j] - outflow_S[1,j] - aged_S[1,j] + births_not_protected[j]
update(S[2:6,]) <- S[i,j] - outflow_S[i,j] - aged_S[i,j] + aged_S[i-1,j]
update(S[7,]) <- S[7,j] - outflow_S[7,j] - aged_S[7,j] + aged_S[6,j] + aged_Sm[6,j]
update(S[8:26,]) <- S[i,j] - outflow_S[i,j] - aged_S[i,j] + aged_S[i-1,j]
update(S[N_age,]) <- S[N_age,j] - outflow_S[N_age,j] + aged_S[(N_age - 1),j]

update(I[1,]) <-  if (tt %% 30 ==0) I[1,j] - outflow_I[1,j] - aged_I[1,j] else I[1,j] - outflow_I[1,j] - aged_I[1,j] + new_infections[1,j] + new_infections_mAb[1,j]
update(I[2:6,]) <- if(tt %% 30 == 0) I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i - 1,j] + new_infections_mAb[i - 1,j] + aged_I[i - 1,j] else I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i,j] + new_infections_mAb[i,j] + aged_I[i - 1,j]
update(I[7,]) <- if(tt %% 30 == 0) I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i - 1,j] + new_infections_mAb[i - 1,j] + aged_I[i - 1,j] else I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i,j] + aged_I[i - 1,j]
update(I[8:24,]) <- if(tt %% 30 == 0) I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i - 1,j] + aged_I[i - 1,j] else I[i,j] - outflow_I[i,j] - aged_I[i,j] + new_infections[i,j] + aged_I[i - 1,j]
update(I[25,]) <- if(tt %% 30 == 0) I[25,j] - outflow_I[25,j] - aged_I[25,j] + new_infections[24,j] + new_infections_mAb[24,j] + aged_I[24,j] else  I[25,j] - outflow_I[25,j] - aged_I[25,j] + new_infections[25,j] + aged_I[24,j]
update(I[26,]) <- if(tt %% 360 == 0) I[26,j] - outflow_I[26,j] - aged_I[26,j] + new_infections[25,j] + aged_I[25,j] else I[26,j] - outflow_I[26,j] - aged_I[26,j] + new_infections[26,j] + aged_I[25,j]
update(I[N_age,]) <- if(tt %% 360 == 0) I[N_age,j] - outflow_I[N_age,j] + sum(new_infections[26:N_age,j]) + aged_I[(N_age - 1),j] else I[N_age,j] - outflow_I[N_age,j] + new_infections[N_age,j] + aged_I[(N_age - 1),j]

update(R[1,]) <- if(tt %% 30 == 0) R[1,j] - outflow_R[1,j] - aged_R[1,j] else R[1,j] - outflow_R[1,j] - aged_R[1,j] + new_recoveries[1,j]
update(R[2:25,]) <- if(tt %% 30 == 0) R[i,j] - outflow_R[i,j] - aged_R[i,j] + new_recoveries[i - 1,j] + aged_R[i - 1,j] else R[i,j] - outflow_R[i,j] - aged_R[i,j] + new_recoveries[i,j] + aged_R[i - 1,j]
update(R[26,]) <- if(tt %% 360 == 0) R[26,j] - outflow_R[26,j] - aged_R[26,j] + new_recoveries[25,j] + aged_R[25,j] else R[26,j] - outflow_R[26,j] - aged_R[26,j] + new_recoveries[26,j] + aged_R[25,j]
update(R[N_age,]) <- if(tt %% 360 == 0) R[N_age,j] - outflow_R[N_age,j] + sum(new_recoveries[26:N_age,j]) + aged_R[(N_age - 1),j] else R[N_age,j] - outflow_R[N_age,j] + new_recoveries[N_age,j] + aged_R[(N_age - 1),j]

update(S2[1,]) <- if(tt %% 30 == 0) S2[1,j] - outflow_S2[1,j] - aged_S2[1,j] else S2[1,j] - outflow_S2[1,j] - aged_S2[1,j] + new_waned[1,j] + new_waned_i2[1,j]
update(S2[2:25,]) <- if(tt %% 30 == 0) S2[i,j] - outflow_S2[i,j] - aged_S2[i,j] + aged_S2[i-1,j] + new_waned[i - 1,j] + new_waned_i2[i - 1,j] else S2[i,j] - outflow_S2[i,j] - aged_S2[i,j] + aged_S2[i - 1,j] + new_waned[i,j] + new_waned_i2[i,j]
update(S2[26,]) <- if(tt %% 360 == 0)S2[i,j] - outflow_S2[26,j] - aged_S2[26,j] + aged_S2[25,j] + new_waned[26 - 1,j] + new_waned_i2[25,j] else S2[26,j] - outflow_S2[26,j] - aged_S2[26,j] + aged_S2[25,j] + new_waned[26,j] + new_waned_i2[26,j]
update(S2[N_age,]) <- if(tt %% 360 == 0) S2[N_age,j] - outflow_S2[N_age,j] + aged_S2[(N_age - 1),j] + sum(new_waned[26:N_age,j]) + sum(new_waned_i2[26:N_age,j]) else S2[N_age,j] - outflow_S2[N_age,j] + aged_S2[(N_age - 1),j] + new_waned[N_age,j] + new_waned_i2[N_age,j]

update(I2[1,]) <- if(tt %% 30 == 0) I2[1,j] - outflow_I2[1,j] - aged_I2[1,j] else I2[1,j] - outflow_I2[1,j] - aged_I2[1,j] + new_reinfections[1,j]
update(I2[2:25,]) <- if(tt %% 30 == 0) I2[i,j] - outflow_I2[i,j] - aged_I2[i,j] + new_reinfections[i - 1,j] + aged_I2[i - 1,j] else I2[i,j] - outflow_I2[i,j] - aged_I2[i,j] + new_reinfections[i,j] + aged_I2[i - 1,j]
update(I2[26,]) <- if(tt %% 360 == 0)I2[26,j] - outflow_I2[26,j] - aged_I2[26,j] + new_reinfections[25,j] + aged_I2[25,j] else I2[26,j] - outflow_I2[26,j] - aged_I2[26,j] + new_reinfections[26,j] + aged_I2[25,j]
update(I2[N_age,]) <- if(tt %% 360 == 0) I2[N_age,j] - outflow_I2[N_age,j] + sum(new_reinfections[26:N_age,j]) + aged_I2[(N_age - 1),j] else I2[N_age,j] - outflow_I2[N_age,j] + new_reinfections[N_age,j] + aged_I2[(N_age - 1),j]

update(R2[1,]) <- if(tt %% 30 == 0) R2[1,j] - outflow_R2[1,j] - aged_R2[1,j] else R2[1,j] - outflow_R2[1,j] - aged_R2[1,j] + new_recoveries_i2[1,j]
update(R2[2:25,]) <- if(tt %% 30 == 0) R2[i,j] - outflow_R2[i,j] - aged_R2[i,j] + new_recoveries_i2[i - 1,j] + aged_R2[i - 1,j] else R2[i,j] - outflow_R2[i,j] - aged_R2[i,j] + new_recoveries_i2[i,j] + aged_R2[i - 1,j]
update(R2[26,]) <- if(tt %% 360 == 0) R2[26,j] - outflow_R2[26,j] - aged_R2[26,j] + new_recoveries_i2[25,j] + aged_R2[25,j] else R2[26,j] - outflow_R2[26,j] - aged_R2[26,j] + new_recoveries_i2[26,j] + aged_R2[25,j]
update(R2[N_age,]) <- if(tt %% 360 == 0) R2[N_age,j] - outflow_R2[N_age,j] + sum(new_recoveries_i2[26:N_age,j]) + aged_R2[(N_age - 1),j] else R2[N_age,j] - outflow_R2[N_age,j] + new_recoveries_i2[N_age,j] + aged_R2[(N_age - 1),j]


update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work


## record total population size and seroprevalence in each age group, and in adults >4

update(N[]) <- (sum(Sm[1:6,i]) + sum(S[1:N_age,i]) + sum(I[1:N_age,i]) + sum(R[1:N_age,i]) + sum(S2[1:N_age,i]) + sum(I2[1:N_age,i]) + sum(R2[1:N_age,i]))
N_tot <- sum(N[])

#update(seroprevalence[1:6,]) <- (I[i,j] + R[i,j] + S2[i,j] + I2[i,j] + R2[i,j]) / (Sm[i,j] + S[i,j] + I[i,j] + R[i,j] + S2[i,j] + I2[i,j] + R2[i,j])
#update(seroprevalence[7:N_age,]) <- (I[i,j] + R[i,j] + S2[i,j] + I2[i,j] + R2[i,j]) / (S[i,j] + I[i,j] + R[i,j] + S2[i,j] + I2[i,j] + R2[i,j])
update(seropoz_A4[]) <- (sum(S2[N_age,i]) + sum(I[N_age,i]) + sum(I2[N_age,i]) + sum(R[N_age,i]) + sum(R2[N_age,i]))/ (sum(S[N_age,i]) + sum(S2[N_age,i]) + sum(I[N_age,i]) + sum(I2[N_age,i]) + sum(R[N_age,i]) + sum(R2[N_age,i]))

update(S_patch[]) <- sum(S[1:N_age,i]) + sum(S2[1:N_age,i]) + sum(Sm[1:6,i])
update(I_patch[]) <- sum(I[1:N_age,i]) + sum(I2[1:N_age,i])
update(R_patch[]) <- sum(R[1:N_age,i]) + sum(R2[1:N_age,i])

# INITIALISATION OF THE MODEL

## initial states

initial(Sm[1:6,]) <- 0
initial(S[1:26,]) <- S_ini_p[i,j] # will be user-defined
initial(S[N_age,]) <- 0
initial(I[1:N_age,]) <- I_ini[i,j] # will be user-defined
initial(R[1:N_age,]) <- 0
initial(S2[1:26,]) <- 0
initial(S2[N_age,]) <- S_ini_p[N_age,j]
initial(I2[1:N_age,]) <- 0
initial(R2[1:N_age,]) <- 0
initial(S_patch[]) <- sum(S_ini_p[1:N_age, i])
initial(I_patch[]) <- sum(I_ini[1:N_age,i])
initial(R_patch[]) <- 0
initial(tt) <- 1
#initial(seroprevalence[1:N_age,]) <- 0
initial(seropoz_A4[]) <- 0
initial(N[]) <- sum(S_ini_p[1:N_age,i]) + sum(I_ini[1:N_age,i])

# DEMOGRAPHIC EQUILIBRIUM SOLUTION 

## setting initial conditions using the equilibrium solution for age distribution
births_detr[1:360] <- 10000 * p_alpha * (1 + cos(3 * cos(pi * i / 360))) # change to fixed N_0 = 100 to avoid NaNs
births_det[1:360] <- rpois(births_detr[i]) 

## if we start the model with the equilibrium amount in each of the first month-wide compartments,
## and no camels in the 2nd year of life (they would have just moved into the adult compartment),
## then from here camels will start filling the yr 2 compartment every month and then eveyr year this will
## empty into the adult compartment. Birthrate will be set to balance summed death rate of this age distribution.

a_max <- 14 ## estimated max age of camels in years (as death rate of adults +2yrs is 1/5.5 years #7.5 years on av)

S_ini[1,] <- 0
S_ini[2,] <- round((sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * mu[1])), 0)
S_ini[3,] <- round((sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:2]))),0)
S_ini[4,] <- round((sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:3]))),0)
S_ini[5,] <- round((sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:4]))),0)
S_ini[6,] <- round((sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:5]))),0)
S_ini[7,] <- round((sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:6]))),0)
S_ini[8,] <- round((sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:7]))),0)
S_ini[9,] <- round((sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:8]))),0)
S_ini[10,] <- round((sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:9]))),0)
S_ini[11,] <- round((sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:10]))),0)
S_ini[12,] <- round((sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:11]))),0)
S_ini[13,] <- round(sum(births_det[1:(360 - 330)]) * exp(- (30 * sum(mu[1:12]))),0)
S_ini[14,] <- round((sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * sum(mu[1:13]))), 0)
S_ini[15,] <- round((sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:14]))),0)
S_ini[16,] <- round((sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:15]))),0)
S_ini[17,] <- round((sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:16]))),0)
S_ini[18,] <- round((sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:17]))),0)
S_ini[19,] <- round((sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:18]))),0)
S_ini[20,] <- round((sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:19]))),0)
S_ini[21,] <- round((sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:20]))),0)
S_ini[22,] <- round((sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:21]))),0)
S_ini[23,] <- round((sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:22]))),0)
S_ini[24,] <- round((sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:23]))),0)
S_ini[25,] <- 0
S_ini[26,] <- round(sum(births_det[1:360]) *  exp(-(30*(sum(mu[1:25])) + 360*mu[25])), 0)
S_ini[N_age,] <- round(a_max * (sum(births_det[1:360]) * exp(- ((30 * sum(mu[1:24])) + 360 * sum(mu[25:26])))),0)

I_ini[12,1] <-1
I_ini[12,2:n_patch] <- 0
I_ini[1:11,] <- 0
I_ini[13:14,] <- 0

S_ini_p[1:N_age,] <- round(((S_ini[i,j] / sum(S_ini[1:N_age,j])) * N_0), 0)

## useful outputs

#########################################
## number of individuals in each state ##
#########################################

# total number of individuals still protected by maternal immunity
output(M) <- sum(Sm[,])

output(S_1) <- sum(S[,]) # susceptible individuals never infected
output(I_1) <- sum(I[,]) # individuals infectious for the 1st time
output(R_1) <- sum(R[,]) # individuals recovered from a 1st infection
output(S_2) <- sum(S2[,]) # susceptible individuals whose immunity has waned
output(I_2) <- sum(I2[,]) # individuals infectious for the 2nd+ time
output(R_2) <- sum(R2[,]) # individuals recovered from 2nd+ infections

output(Stot) <- sum(S[,]) + sum(S2[,]) + sum(Sm[,]) # total number of susceptible individuals
output(Itot) <- sum(I[,]) + sum(I2[,]) # total number of infectious individuals
output(Rtot) <- sum(R[,]) + sum(R2[,]) # total number of recovered individuals 

output(Ntot) <- N_tot # total number of individuals
output(N_J) <- sum(Sm[1:6,]) + sum(S[1:13,]) + sum(S2[1:13,]) + sum(I[1:13,]) + sum(I2[1:13,]) + sum(R[1:13,]) + sum(R2[1:13,])
output(N_A) <- sum(S[25:N_age,]) + sum(S2[25:N_age,]) + sum(I[25:N_age,]) + sum(I2[25:N_age,]) + sum(R[25:N_age,]) + sum(R2[25:N_age,])


####################
## seroprevalence ##
####################

# output(seropoz_A) <- 100 * (sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))/(sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))
# output(seropoz_J) <- 100 * (sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24])) / (sum(Sm[1:24]) + sum(S[1:24]) + sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24]))
# output(seropz_tot) <- 100 * (sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))/(sum(Sm[1:24]) + sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))
# 
# ############################
# ## age at first infection ##
# ############################
# 
# output(reinf_1) <- new_reinfections[1]
# output(reinf_2) <- new_reinfections[2]
# output(incidence_new_inf) <- sum(new_infections[1:N_age])
# output(incidence_indig_inf) <- sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age])
# output(total_incidence) <- (sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age]))

###########
## other ##
###########
# output(inf) <- rate_infection
output(birthrate) <- birth_rate
output(births) <- new_births
# output(importations) <- imported_cases

## dim calls needed for arrays
dim(S_patch) <- n_patch
dim(I_patch) <- n_patch
dim(R_patch) <- n_patch
dim(Sm) <- c(6,n_patch)
dim(S) <- c(N_age, n_patch)
dim(I) <- c(N_age, n_patch)
dim(R) <- c(N_age, n_patch)
dim(S2) <- c(N_age, n_patch)
dim(I2) <- c(N_age, n_patch)
dim(R2) <- c(N_age, n_patch)
dim(N) <- n_patch
dim(S_ini) <- c(N_age, n_patch)
dim(S_ini_p) <- c(N_age, n_patch)
dim(I_ini) <- c(N_age, n_patch)
dim(new_infections) <- c(N_age, n_patch)
dim(new_infections_mAb) <- c(6, n_patch)
dim(new_recoveries) <- c(N_age, n_patch)
dim(new_waned) <- c(N_age, n_patch)
dim(new_reinfections) <- c(N_age, n_patch)
dim(new_recoveries_i2) <- c(N_age, n_patch)
dim(new_waned_i2) <- c(N_age, n_patch)
dim(aged_Sm) <- c(6, n_patch)
dim(aged_S) <- c(26, n_patch)
dim(aged_I) <- c(26, n_patch)
dim(aged_R) <- c(26, n_patch)
dim(outflow_Sm) <- c(6, n_patch)
dim(outflow_S) <- c(N_age, n_patch)
dim(outflow_I) <- c(N_age, n_patch)
dim(outflow_R) <- c(N_age, n_patch)
dim(aged_S2) <- c(26, n_patch)
dim(aged_I2) <- c(26, n_patch)
dim(aged_R2) <- c(26, n_patch)
dim(outflow_S2) <- c(N_age, n_patch)
dim(outflow_I2) <- c(N_age, n_patch)
dim(outflow_R2) <- c(N_age, n_patch)
dim(births_det) <- 360
dim(births_detr) <- 360
#dim(seroprevalence) <- c(N_age, n_patch)
dim(seropoz_A4) <- n_patch
# because mu varies with age
dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- c(6,n_patch)
dim(p_S) <- c(N_age, n_patch)
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- c(N_age, n_patch)
dim(p_I2) <- N_age
dim(p_R2) <- N_age
dim(norm_p_infection) <- c(N_age, n_patch)
dim(norm_p_infection_mAb) <- c(6, n_patch)
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_reinfection) <- c(N_age, n_patch)
dim(births_protected) <- n_patch
dim(births_not_protected) <- n_patch
dim(FOI) <- n_patch #the patch specific frequency term denoting which other patches infectious individuals contribute to the FOI in each patch
dim(FOI_2) <- n_patch # same as above but for 2nd+ times infected mindividuals in other patches
dim(rate_infection) <- n_patch
dim(rate_infection_mAb) <- n_patch
dim(rate_reinfection) <- n_patch
dim(p_infection) <- n_patch
dim(p_infection_mAb) <- n_patch
dim(p_reinfection) <- n_patch