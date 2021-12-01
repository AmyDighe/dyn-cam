# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model_vax <- odin::odin("models/odin_single_patch_model_vaccination_proportion.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################


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

# vaccination parameters

v_gamma <- gamma # rate of recovery from infection in vaccinated animals
v_sigma <- sigma # rate of waning complete infection induced immunity in vaccinated animals
v_sigma_m <- sigma_m # rate of waning of mAbs in vaccinated animals

v_mAb_susc <- mAb_susc # relative susceptibility of vaccinated animals with mAbs
v_susc <- 0.75 # relative susceptibility of vaccinated naive animals
v_Ab_susc <- 0.25  # relative susceptibility of vaccinated previously infected animals

v_shed <- 1/92 # infectiousness of vaccinated naive animals cf naive unvaccinated
v_reduced_shed <- 1/724 # infectiousness of vaccinated previously infected animals cf naive unvaccinated

# age dependent rates of vaccination
vax <- c(rep(0, 6), 0.8, rep(0, 42))

rho <- 0 #rate at which vaccine induced immunity wanes

###############
## run model ##
###############

last_imp <- imp_t[5]

#tic("single_patch_vax")
# include any user-defined parameters as arguments here
x <- sir_model_vax$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, Ab_susc = Ab_susc, 
                   mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                   mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, N_0 = N_0,
                   importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2,
                   v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_mAb_susc = v_mAb_susc, 
                   v_Ab_susc = v_Ab_susc, v_shed = v_shed, v_reduced_shed = v_reduced_shed, vax = vax)

#out <- as.data.frame(replicate(100, x$run(t)[, 403]))
out <- as.data.frame(x$run(t))

#toc()
plot(out$`vS[8]`/(out$`N_pop[8]`))
