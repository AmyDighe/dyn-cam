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

estab <- vector(length = (dim(par_grid)[1]))
m_persist <- estab
persist <- estab
period_av <- estab
period_sd <- estab
period_med <- estab

for(i in 1:(dim(par_grid)[1])){
  
  imp_t <- par_grid$import_time[i] + (360 * seq(0,4, by = 1))
  e_t <- imp_t[5] + (3 * 360)
  m_t <- imp_t[5] + (25 * 360)
  p_t <- imp_t[5] + (50 * 360)
  
  x <- sir_model$new(alpha = alpha, beta = par_grid$beta[i], gamma = gamma, 
                     sigma = par_grid$waning[i], sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu = mu,
                     N_0 = par_grid$pop[i], importation_rate = importation_rate, imp_t = imp_t, 
                     delta = par_grid$seasonality[i], ind1 = ind1, ind2 = ind2)

  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
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

 saveRDS(file = "results/persist_estab_period_sp.rds", object = data.frame(persist = persist, m_persist = m_persist, estab = estab,
                                                                  period_av = period_av, period_sd = period_sd, period_med = period_med))
d <- readRDS("results/persist_estab_period_sp.rds")
results <- cbind(d, par_grid)
results <- results %>% mutate(persist_estab = 100 * persist/estab)
#results$persist_estab[is.nan(results$persist_estab)]<-0

waning.labs <- c("\U03C3 = 30 days", "\U03C3 = 90 days")
names(waning.labs) <- unique(par_grid$waning)
seasonality.labs <- c("\U03B4 = 0.0", "\U03B4 = 0.5", "\U03B4 = 1.0")
names(seasonality.labs) <- unique(par_grid$seasonality)

# persistance by population size
ggplot(results)+
  geom_line(aes(x = log10(pop), y = persist, col = as.character(beta), lty = as.character(import_time)),
            size = 1)+
  facet_grid(waning~seasonality, labeller = labeller(waning = waning.labs, seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  scale_linetype_discrete(name = "import time\nwithin the year")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/persistence.png")

# persist/estab
ggplot(results)+
  geom_line(aes(x = log10(pop), y = persist_estab, col = as.character(beta), lty = as.character(import_time)),
            size = 1)+
  facet_grid(waning~seasonality, labeller = labeller(waning = waning.labs, seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  scale_linetype_discrete(name = "import time\nwithin the year")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/persist_estab.png")

# established
ggplot(results)+
  geom_line(aes(x = log10(pop), y = estab, col = as.character(beta), lty = as.character(import_time)),
            size = 1)+
  facet_grid(waning~seasonality, labeller = labeller(waning = waning.labs, seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  scale_linetype_discrete(name = "import time\nwithin the year")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/estab5.png")

# try estab with only 1 importation?
# run with more fine population classes

# periodicity
ggplot(results %>% filter(import_time == 151) %>% filter(persist>50))+
  geom_point(aes(x = log10(pop), y = period_av/360, col = as.character(beta)), size = 4)+
  scale_color_viridis_d(name = "beta", alpha = 0.5)+
  facet_grid(waning~seasonality, 
             labeller = labeller(waning = waning.labs, seasonality = seasonality.labs))+
  theme_minimal()+
  geom_errorbar(aes(x = log10(pop), 
                    ymin = (period_av - period_sd)/360, ymax = (period_av + period_sd)/360,
                    col = as.character(beta)))+
  geom_hline(yintercept = 1, col = "red", lty = 2)+
  geom_hline(yintercept = 2, col = "red", lty = 2)+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("period (years)")
ggsave(filename = "figs/mean_periodicity.png")

ggplot(results %>% filter(import_time == 151) %>% filter(persist>50))+
  geom_point(aes(x = log10(pop), y = period_med/360, col = as.character(beta)), size = 4)+
  scale_color_viridis_d(name = "beta", alpha = 0.5)+
  facet_grid(waning~seasonality, 
             labeller = labeller(waning = waning.labs, seasonality = seasonality.labs))+
  theme_minimal()+
  geom_hline(yintercept = 1, col = "red", lty = 2)+
  geom_hline(yintercept = 2, col = "red", lty = 2)+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("period (years)")
ggsave(filename = "figs/median_periodicity.png")


## looking into behavior when beta = 1 & delta = 1

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
        acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
        period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
      } 
      period_av[i] <- mean(period)
      period_sd[i] <- sd(period)
      period_med[i] <- median(period)
      
    } else if(persist[i] == 1) {
      ACF <- acf(out_persist[e_t:p_t], lag.max = 360*5, plot = TRUE)
      acf_df <- data.frame(acf = ACF$acf[100:(360*5)], lag = ACF$lag[100:(360*5)])
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
