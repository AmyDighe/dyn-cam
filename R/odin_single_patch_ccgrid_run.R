# MSIRSIR

# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input a value for birth rate (per camel per day) 
alpha <- 0.000565 #0.0006 # 90% female * 50% reproductive age * 45.2% fecundity

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 
gamma <- 1/duration_infection

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 #1/92 # based on AUC from shedding in Alharbi 

# input values for the age dependent removal rate - balance birthrate

mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

# input the time period that you wish to run the model for (in days)
time_period <- 5671
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

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

persist_15 <- vector(length = dim(par_grid)[1])
period <- vector(length = dim(par_grid)[1])
ACF_example <- vector(mode = "list", length = dim(par_grid)[1])

for(i in 1:(dim(par_grid)[1])){
  # include any user-defined parameters as arguments here
  x <- sir_model(alpha = alpha, beta = par_grid$beta[i], gamma = gamma, 
                 sigma = par_grid$waning[i],sigma_m = par_grid$waning_mAb[i],
                 Ab_susc = Ab_susc, mAb_susc = mAb_susc, reduced_shed = reduced_shed, 
                 mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,mu_3rd_yr = mu_3rd_yr, 
                 mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, 
                 N_0 = par_grid$pop[i], importation_rate = importation_rate, 
                 imp_t = par_grid$import_time[i], delta = par_grid$seasonality[i], 
                 ind1 = ind1, ind2 = ind2) 
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  # calculate what proportion of runs establish endemicity
  persist_15[i] <- sum(out[5670,] > 0)
  # calculate periodicity of endemic runs using autocorrelation
  if(persist_15[i] > 0){
    idx_persist <- which(out[5670,] > 0)
    out_persist <- out[,idx_persist]
    period_v <- vector(length = persist_15[i])
    
    if(persist_15[i] > 1){
      for(j in 1:(persist_15[i])){
        ACF <- acf(out_persist[3000:(3000 + (360*5)), j], lag.max = 360*5, plot = FALSE)
        acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
        period_v[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period[i] <- mean(period_v)
    } else if(persist_15[i] == 1) {
      ACF <- acf(out_persist[3000:(3000 + (360*5))], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
      period_v <- acf_df[which.max(acf_df$acf), ]$lag 
      period[i] <- period_v
    }
    ACF_example[[i]] <- ACF
  } else {
    period[i] <- NA
    ACF_example[[i]] <- NA
  }
  
  # tracker
  print(i)
  
}

saveRDS(persist_15, file = "persist_15.rds")
saveRDS(period, file = "period.rds")
saveRDS(ACF_example, file = "ACF_example.rds")