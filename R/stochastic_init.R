
###############################################################
## stochastic initialisation at rough zoographic equilibrium ##
###############################################################
stoch_init <- function(alpha, delta, N_0, mu, N_age, n_r, n_c){

births_detr <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * seq(1:360) / 360))))
births_det <- vector(length = length(births_detr)) 

for(i in 1:(length(births_detr))){
  births_det[i] <- rpois(n = 1, lambda = births_detr[i])
}

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

# setting initial number of camels at demographic equilibrium
S_ini <- vector(length = N_age)
S_ini_p <- vector(length = N_age)
S_ini[1] <- 0

#for the first 48 month wide age classes
for(i in 1:(N_age-1)){
  S_ini[i] <- sum(births_det[ind1[i]:ind2[i]]) * exp(- (30 * sum(mu[1:(i - 1)])))
}
# special treatment for the final open-ended compartment 49
a_max <- 32 ## estimated max number of years spent in the last age class before death (only 1% of the population above after 32 yrs)
## (calculated in "calculating_a_max.R using exp decay mod)
yearly_influx <- sum(births_det[]) * exp(-(30 * (sum(mu[1:48])))) #annual inflow into the final age compartment
yr <- c(1:a_max)  # number of years of influx before max age expectancy reached
cohort_remaining <- vector(length = a_max)
# camels remaining from each cohort to enter in the last 32 years influx that remain
for(i in 1:a_max){
  cohort_remaining[i] <- yearly_influx * exp(- (360 * yr[i] * mu[N_age]))  
}

S_ini[N_age] <- sum(cohort_remaining[1:a_max])

# getting proportion of animals in each age-class, multiplying by N_0 and rounding to whole animals
for(i in 1:N_age){
  S_ini_p[i] <- round((S_ini[i] / sum(S_ini[1:N_age])) * N_0)
}
# expanding out for all patches
S_ini_a <- array(S_ini_p, dim = c(N_age, n_r, n_c)) 

return(S_ini_a)
}
