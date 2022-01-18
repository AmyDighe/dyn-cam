# FINDING BIFURCATION POINT - PERIODICITY

sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input a value for birth rate (per camel per day) 
alpha <- 0.000565

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14
gamma <- 1/duration_infection

# input a value for the average rate of waning of mAbs
sigma_m <- 4.42/360 # from the catalytic work

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75 # default = 
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

# input the time period that you wish to run the model for (in days)
time_period <- 19710
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

estab <- vector(length = (dim(par_grid)[1]))
persist <- estab
period_av <- estab
period_sd <- estab
period_med <- estab

imp_t <- par_grid_bifurc$import_time[1] + (360 * seq(0,4, by = 1))
e_t <- imp_t[5] + (3 * 360)
m_t <- imp_t[5] + (25 * 360)
p_t <- imp_t[5] + (50 * 360)

for(i in 1:(dim(par_grid_bifurc)[1])){

  x <- sir_model$new(alpha = alpha, beta = par_grid_bifurc$beta[i], gamma = gamma, 
                     sigma = par_grid_bifurc$waning[i], sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                     mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4,
                     N_0 = par_grid_bifurc$pop[i], importation_rate = importation_rate, imp_t = imp_t, 
                     delta = par_grid_bifurc$seasonality[i], ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  estab[i] <- sum(out[e_t,] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 0){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    period <- vector(length = persist[i])
    
    if(persist[i] > 1){
      for(j in 1:(persist[i])){
        ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
        acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
        period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period_av[i] <- mean(period)
      period_sd[i] <- sd(period)
      period_med[i] <- median(period)
      
    } else if(persist[i] == 1) {
      ACF <- acf(out_persist[m_t:p_t], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
      period <- acf_df[which.max(acf_df$acf), ]$lag 
      
      period_av[i] <- period
      period_sd[i] <- NA
      period_med[i] <- period
    }
    
  } else {
    period_av[i] <- NA
    period_sd[i] <- NA
    period_med[i] <- NA
  }
  print(i)
}

saveRDS(file = "results/persist_estab_period_bifurc_sp.rds", object = data.frame(persist = persist[1:96], m_persist = m_persist[1:96], estab = estab[1:96],
                                                                                 period_av = period_av[1:96], period_med = period_med[1:96], period_sd = period_sd[1:96]))
                                                                          