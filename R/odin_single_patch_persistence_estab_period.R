sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input the time period that you wish to run the model for (in days)
time_period <- 360*35
t <- seq(0:time_period)

estab <- vector(length = (dim(par_grid_core)[1]))
m_persist <- estab
persist <- estab
period <- matrix(data = NA, nrow = 100, ncol = dim(par_grid_core)[1])

e_t <- 360*15
m_t <- 360*25
p_t <- 360*35

for(i in 1:(dim(par_grid_core)[1])){
  x <- sir_model$new(alpha = alpha, beta = par_grid_core$beta[i], gamma = gamma, 
                     sigma = par_grid_core$waning[i], sigma_m = sigma_m, 
                     Ab_susc = par_grid_core$susc[i], mAb_susc = mAb_susc, 
                     reduced_shed = par_grid_core$shedding[i], mu = mu,
                     N_0 = par_grid_core$pop[i], importation_rate = importation_rate, 
                     imp_t = imp_t, delta = par_grid_core$seasonality[i], 
                     ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)
  out <- as.data.frame(x$run(1:100)[,354])

  out <- as.data.frame(replicate(100, x$run(t)[, 354]))
  
  estab[i] <- sum(out[e_t,] > 0)
  m_persist[i] <- sum(out[m_t, ] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 50){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    
    for(j in 1:(persist[i])){
      ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[180:(360*5)], lag = ACF$lag[180:(360*5)])
      # write something in here about making sure acf is above the significance level
      period[j,i] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
  } else {
    period[,i] <- NA
  }
  print(paste(i, "persist=", persist[i], sep = " "))
}

 saveRDS(file = "results/persist_estab_period_sp_010222.rds", object = list(persist = data.frame(persist = persist, m_persist = m_persist, estab = estab),
                                                                  period = period))
 
 ############################################################# 
 ## example with annual, biennial and triennial periodicity ##
 #############################################################
 
 # input the time period that you wish to run the model for (in days)
 time_period <- 360*35
 t <- seq(0:time_period)
 
 estab <- vector(length = (dim(par_grid_core)[1]))
 m_persist <- estab
 persist <- estab
 period <- matrix(data = NA, nrow = 100, ncol = dim(par_grid_core)[1])
 
 e_t <- 360*15
 m_t <- 360*25
 p_t <- 360*35
 
 for(i in c(319, 315)){
   x <- sir_model$new(alpha = alpha, beta = par_grid_core$beta[i], gamma = gamma, 
                      sigma = par_grid_core$waning[i], sigma_m = sigma_m, 
                      Ab_susc = par_grid_core$susc[i], mAb_susc = mAb_susc, 
                      reduced_shed = par_grid_core$shedding[i], mu = mu,
                      N_0 = par_grid_core$pop[i], importation_rate = importation_rate, 
                      imp_t = imp_t, delta = par_grid_core$seasonality[i], 
                      ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)
   
   out <- as.data.frame(replicate(100, x$run(t)[, c(354, 369)]))
   saveRDS(out, file = paste("results/persistence/time_series_examples_period/out", i, ".rds", sep = ""))
   
   out_Itot <- out[,grep("Itot", colnames(out))]
   
   estab[i] <- sum(out_Itot[e_t,] > 0)
   m_persist[i] <- sum(out_Itot[m_t, ] > 0)
   persist[i] <- sum(out_Itot[p_t, ] > 0)
   
   if(persist[i] > 50){
     
     idx_persist <- which(out_Itot[p_t, ] > 0)
     out_persist <- out_Itot[,idx_persist]
     
     for(j in 1:(persist[i])){
       ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
       acf_df <- data.frame(acf = ACF$acf[180:(360*5)], lag = ACF$lag[180:(360*5)])
       # write something in here about making sure acf is above the significance level
       period[j,i] <- acf_df[which.max(acf_df$acf), ]$lag 
     }
   } else {
     period[,i] <- NA
   }
   print(paste(i, "persist=", persist[i], sep = " "))
   saveRDS(file = paste("results/persistence/time_series_examples_period/period", i, ".rds", sep = ""), object = period[,i])
 }
 

#####################################################
## looking into behavior when beta = 1 & delta = 1 ##
#####################################################

estab <- vector(length = 3)
m_persist <- estab
persist <- estab
period_av <- estab
period_sd <- estab
period_med <- estab

imp_t <- 151 + (360 * seq(0,4, by = 1))
e_t <- imp_t[5] + (3 * 360)
m_t <- imp_t[5] + (25 * 360)
p_t <- imp_t[5] + (50 * 360)

beta_v <- c(0.5, 1, 1.5)

for(i in 1:3){
  x <- sir_model$new(alpha = alpha, beta = beta_v[i], gamma = gamma, 
                     sigma = 1/90, sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu = mu,
                     N_0 = 50000, importation_rate = importation_rate, imp_t = imp_t, 
                     delta = 1, ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  matplot(y = out[(imp_t[5]):p_t,], type = "l",
          lty = 1, col = "#00000022")
  
  estab[i] <- sum(out[e_t,] > 0)
  m_persist[i] <- sum(out[m_t, ] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 0){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    period <- vector(length = persist[i])
    
    if(persist[i] > 1){
      for(j in 1:(persist[i])){
        ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
        acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
        period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period_av[i] <- mean(period)
      period_sd[i] <- sd(period)
      period_med[i] <- median(period)
      
    } else if(persist[i] == 1) {
      ACF <- acf(out_persist[m_t:p_t], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
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













## getting some examples out
results$id <- 1:360
results %>% filter(pop == 1000000,
                   beta == 0.25, seasonality == 0,
                   import_time == 151)

for(i in 1:(dim(par_grid)[1])){
  i = 109
  
  imp_t <- par_grid$import_time[i] + (360 * seq(0,4, by = 1))
  e_t <- imp_t[5] + (3*360)
  p_t <- imp_t[5] + (18*360)
  
  x <- sir_model$new(alpha = alpha, beta = par_grid$beta[i], gamma = gamma, 
                     sigma = par_grid$waning[i], sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu = mu,
                     N_0 = par_grid$pop[i], importation_rate = importation_rate, imp_t = imp_t, 
                     delta = par_grid$seasonality[i], ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  estab[i] <- sum(out[e_t,] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 0){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    period <- vector(length = persist[i])
    
    if(persist[i] > 1){
      for(j in 1:(persist[i])){
        ACF <- acf(out_persist[6840:p_t, j], lag.max = 360*5, plot = F)
        acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
        period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period_av[i] <- mean(period)
      period_sd[i] <- sd(period)
      period_med[i] <- median(period)
      
    } else if(persist[i] == 1) {
      ACF <- acf(out_persist[e_t:p_t], lag.max = 360*5, plot = TRUE)
      acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
      period <- acf_df[which.max(acf_df$acf), ]$lag 
      
      period_av[i] <- period
      period_sd[i] <- NA
    }
    
  } else {
    period_av[i] <- NA
    period_sd[i] <- NA
  }
  print(i)
}

out_persist$step <- t

matplot(y = out[m_t:p_t,1:5], type = "l",
        lty = 1, lwd = 1.5)
