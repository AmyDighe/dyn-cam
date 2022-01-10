sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################


# input a value for birth rate (per camel per day) 
alpha <- 0.000565 #0.0006 # 90% female * 50% reproductive age * 45.2% fecundity

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 
gamma <- 1/duration_infection

# input a value for the average rate of waning of mAbs
sigma_m <- 4.42/360 # from the catalytic work

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- c(90,1,1,1,90,30,90,1,90) # default = 
sigma_v <- 1/duration_immunity # default = 

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc_v <- c(0.75, 1, 0.75, 0.25, 1, 1, 0.75, 0.75, 1) # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed_v <- c(rep(1/92, 6),rep(0.1, 3)) # based on AUC from shedding in Alharbi 

# input values for the age dependent removal rate - balance birthrate
mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

# input an initial population size
N_0 <- 1000000

# input the time period that you wish to run the model for (in days)
time_period <- 19710 
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day input that day here
imp_t <- 269  + (360 * seq(0, 4, by = 1))

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

###############
## run model ##
###############
scenario <- names(beta_list)
last_imp <- imp_t[5]

for(j in 1:9){
  j = 7
  Ab_susc <- Ab_susc_v[j]
  sigma <- sigma_v[j]
  reduced_shed <- reduced_shed_v[j]
  beta_vector <- beta_list[[j]]
  mean_foi <- vector(length = length(beta_vector))
  
for(i in 1:(length(beta_vector))){
beta <- beta_vector[i]

x <- sir_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = sigma_m, Ab_susc = Ab_susc, 
               mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
               mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, N_0 = N_0,
               importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2)

if(beta<1){
  nruns <- 1000
} else {
  nruns <- 100
}
out <- as.data.frame(replicate(nruns, x$run(t)[, c(398, 401, 403, 405)]))
out_I1 <- out[,grep("I_1", colnames(out))]
out_I2 <- out[,grep("I_2", colnames(out))]
out_Itot <- out[,grep("Itot", colnames(out))]
out_N <- out[,grep("Ntot", colnames(out))]
idx_persist <- which(out_Itot[last_imp + (50*360),] > 0)
out_I1_persist <- out_I1[,idx_persist]
out_I2_persist <- out_I2[,idx_persist]
out_N_persist <- out_N[,idx_persist]

if(length(idx_persist) == 0){
  mean_foi[i] <- NA
} else if(length(idx_persist) > 1){
  I1_N <- colMeans(out_I1_persist[(last_imp + 25*360):(last_imp + 50*360),] / out_N_persist[(last_imp + 25*360):(last_imp + 50*360),])  
  I2_N <- colMeans(out_I2_persist[(last_imp + 25*360):(last_imp + 50*360),] / out_N_persist[(last_imp + 25*360):(last_imp + 50*360),])
  foi <- beta * I1_N + beta * reduced_shed * I2_N
  mean_foi[i] <- mean(foi)
} else if(length(idx_persist) == 1) {
  I1_N <- mean(out_I1_persist[(last_imp + 25*360):(last_imp + 50*360)] / out_N_persist[(last_imp + 25*360):(last_imp + 50*360)])
  I2_N <- mean(out_I2_persist[(last_imp + 25*360):(last_imp + 50*360)] / out_N_persist[(last_imp + 25*360):(last_imp + 50*360)])
  foi <- beta * I1_N + beta * reduced_shed * I2_N
  mean_foi[i] <- mean(foi)
  }

print(paste("i =", i, "foi= ", mean_foi[i], sep = " ")) }

saveRDS(file = paste("results/test_mean_foi_", scenario[j], ".rds", sep = ""), object = mean_foi)

}
