# MSIRS2I2 dust METAPOP 2021

###################
## compile model ##
###################

library(odin.dust)
msirs_meta_dust <- odin.dust::odin_dust("models/dustmeta.R")

#############################
## user defined variables: ##
#############################

N_age <- 49
N_patch <- 25

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

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
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

# input the time period that you wish to run the model for (in days)
time_period <- 9000
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day input that day here
imp_t <- 1510000  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

###############################################################
## stochastic initialisation at rough zoographic equilibrium ##
###############################################################

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
S_ini_m <- matrix(S_ini_p, nrow = length(S_ini_p), ncol = N_patch) 

###################
## run the model ##
###################

n_particles <- 2L

model <- msirs_meta_dust$new(pars = list(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = sigma_m,
                         Ab_susc = Ab_susc, mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu = mu, 
                         N_0 = N_0, importation_rate = importation_rate, imp_t = imp_t, delta = delta,
                         S_ini_p = S_ini_m), 
                         step = 1, 
                         n_particles = n_particles, 
                         n_threads = 1L, 
                         seed = 1L) 
n_steps <- 360

x <- array(NA, dim = c(model$info()$len, n_particles, n_steps))

for (t in seq_len(n_steps)) {
  x[ , , t] <- model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")

