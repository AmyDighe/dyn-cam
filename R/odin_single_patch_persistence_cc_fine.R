## FINER POP DIVISIONS TO LOOK AT CC

sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from shedding in Alharbi 

# input the time period that you wish to run the model for (in days)
time_period <- 19710
t <- seq(0:time_period)

estab <- vector(length = (dim(par_grid_fine_pop)[1]))
m_persist <- estab
persist <- estab
period_av <- estab
period_sd <- estab
period_med <- estab

for(i in 1:(dim(par_grid_fine_pop)[1])){
  
  imp_t <- par_grid_fine_pop$import_time[i] + (360 * seq(0,4, by = 1))
  e_t <- imp_t[5] + (3 * 360)
  m_t <- imp_t[5] + (25 * 360)
  p_t <- imp_t[5] + (50 * 360)
  
  x <- sir_model$new(alpha = alpha, beta = par_grid_fine_pop$beta[i], gamma = gamma, 
                     sigma = par_grid_fine_pop$waning[i], sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                     mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4,
                     N_0 = par_grid_fine_pop$pop[i], importation_rate = importation_rate, imp_t = imp_t, 
                     delta = par_grid_fine_pop$seasonality[i], ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  estab[i] <- sum(out[e_t,] > 0)
  m_persist[i] <- sum(out[m_t, ] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 1){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    period <- vector(length = persist[i])

      for(j in 1:(persist[i])){
        ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
        acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
        period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period_av[i] <- mean(period)
      period_sd[i] <- sd(period)
      period_med[i] <- median(period)
      
    } else {
    period_av[i] <- NA
    period_sd[i] <- NA
    period_med[i] <- NA
  }
  print(i)
}

saveRDS(file = "results/persist_estab_period_cc_fine_sp.rds", object = data.frame(persist = persist, m_persist = m_persist, estab = estab,
                                                                          period_av = period_av, period_sd = period_sd, period_med = period_med))
