# FINDING BIFURCATION POINT - PERIODICITY

sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75 # default = 

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from shedding in Alharbi 

# input the time period that you wish to run the model for (in days)
time_period <- 19710
t <- seq(0:time_period)

estab <- vector(length = (dim(par_grid_bifurc)[1]))
persist <- estab
period_av <- estab
period_sd <- estab
period_med <- estab

par_grid_fine_beta_extrapop <- expand.grid(waning = waning[2], 
                                  seasonality = seasonality[2:3], 
                                  beta = c(0.3, 0.35, 0.4, 0.45, 1.5), 
                                  pop = c(c(1,5) %o% 10^7),
                                  import_time = import_time[1])

par_grid_bifurc <- par_grid_fine_beta_extrapop

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

saveRDS(file = "results/persistence/sp_fine_beta_extrapop.rds", object = data.frame(persist = persist[1:20], estab = estab[1:20],
                                                                                 period_av = period_av[1:20], period_med = period_med[1:20], period_sd = period_sd[1:20]))

## alternative plot idea with proportions

estab <- vector(length = length(pop_og))
persist <- estab
period <- matrix(data = NA, nrow = 100, ncol = length(pop_og))

for(i in 1:(length(pop_og))){
  
  x <- sir_model$new(alpha = alpha, beta = 1, gamma = gamma, 
                     sigma = 1/90, sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                     mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4,
                     N_0 = pop_og[i], importation_rate = importation_rate, imp_t = imp_t, 
                     delta = 1, ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  estab[i] <- sum(out[e_t,] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 1){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]

      for(j in 1:(persist[i])){
        ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
        acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
        # write something in here about making sure acf is above the significance level
        period[j,i] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
      } else {
    period[,i] <- NA
    }
  print(i)
}

saveRDS(file = "alternative_periodicity_1.rds", object = period)

# bin values for stacked barchart

period_025 <- readRDS("results/persistence/alternative_periodicity_025.rds")
period_05<- readRDS("results/persistence/alternative_periodicity_05.rds")
period_1<- readRDS("results/persistence/alternative_periodicity_1.rds")

period <- rbind(period_025[,4:10], period_05[,4:10], period_1[,4:10]) #get rid of pop sizes where persistance is always<50

cuts <- apply(period, 2, cut, c(0, 350, 370, 710, 730, 1070, 1090, 1430, 1450, Inf), 
              labels= c("other", "1 yr", "other", "2 yrs", "other", "3 yrs", "other", "4 yrs", "other"))
cuts <- as.data.frame(cuts)
names(cuts) <- pop_og[4:10]
cuts$beta <- rep(c(0.25, 0.5, 1), each = 100)
cuts_long <- pivot_longer(cuts, cols = 1:7, names_to = "population", values_to = "periodicity")

cuts_long$population <- factor(cuts_long$population,levels = as.character(pop_og[4:10]))
cuts_long$periodicity <- factor(cuts_long$periodicity, levels = c(NA, "other", "4 yrs",
                                                                  "3 yrs", "2 yrs", "1 yr"),
                                exclude = NULL)
beta.labs <- c("\U03B2 = 0.25", "\U03B2 = 0.5", "\U03B2 = 1.0")
names(beta.labs) <- c(0.25, 0.5, 1)

ggplot(data = cuts_long, aes(fill = periodicity))+
  geom_bar(aes(x = population), position = "stack")+
  scale_fill_viridis_d(direction = -1)+
  theme_minimal()+
  theme(text= element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs))+
  ylab("percentage of runs")
ggsave(filename = "results/persistence/periodicity_alternative.png")


# corr <- acf(out_persist[m_t:p_t, j], lag.max=360*5,type="correlation",plot=TRUE,na.action=na.pass)
# significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(out_persist[m_t:p_t, j])))


## example of the ones with median 3 year periodicity
  x <- sir_model$new(alpha = alpha, beta = 0.25, gamma = gamma, 
                     sigma = 1/90, sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                     reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                     mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4,
                     N_0 = 50000, importation_rate = importation_rate, imp_t = imp_t, 
                     delta = 1, ind1 = ind1, ind2 = ind2)
  
  out <- as.data.frame(replicate(100, x$run(t)[, 403]))
  
  established <- sum(out[e_t,] > 0)
  persisted <- sum(out[p_t, ] > 0)
  
  if(persisted > 1){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    period <- vector(length = length(persisted))
    for(j in 1:(persisted)){
      ACF <- acf(out_persist[m_t:p_t, j], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
      # write something in here about making sure acf is above the significance level
      period[j] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
  } else {
    period <- NA
  }

plot(out_persist[,3], type = "l")


