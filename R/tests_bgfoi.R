## testing background foi

sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input the time period that you wish to run the model for (in days)
time_period <- 20000 
t <- seq(0:time_period)

imp_t <- 1

x <- sir_model$new(alpha = 0.000565, beta = 1, gamma = 1/14, sigma = 1/90, 
                   sigma_m = 4.42/360, Ab_susc = 0.75, mAb_susc = 0, 
                   reduced_shed = 1/90, mu = mu, N_0 = 1e+7,
                   importation_rate = 0, imp_t = imp_t, 
                   delta = 1, ind1 = ind1, ind2 = ind2, foi_bg_usr = 0.000015)

out <- as.data.frame(x$run(t))[, c(349, 352, 354, 356)]
plot(out$Itot, type = "l")

out <- as.data.frame(replicate(100, x$run(t)[, c(349, 352, 354, 356)]))
matplot(out[,grep("Itot", names(out))], type = "l")
out_Itot <- out[,grep("Itot", names(out))]
persist <- sum(out_Itot[19999,] > 0)

######################################################################################
## periodicity barplot - is it more stable with background? how is the persistence? ##
######################################################################################
time_period <- 14400 
t <- seq(0:time_period)

foi_bg <- c(0, 0.0001, 0.00001)
delta_v <- c(0,1)
pargrid_foibg_test <- expand.grid(foi_bg = foi_bg,
                                  delta_v = delta_v)
pargrid_foibg_test$imp_t_v <- c(151, rep(50000, 3),
                                151, rep(50000, 3))

pop_test <- pop_og[c(-1, -2, -3)]

e_t <- 13*360
m_t <- 25*360
p_t <- 40*360

for(j in 1:(dim(pargrid_foibg_test)[1])){
  estab <- vector(length = length(pop_test))
  persist <- estab
  period <- matrix(data = NA, nrow = 100, ncol = length(pop_test))
  imp_t <- pargrid_foibg_test$imp_t_v[j] + (360 * seq(0, 4, by = 1))
  
  for(i in 1:(length(pop_test))){
    
  x <- sir_model$new(alpha = alpha, beta = 0.25, gamma = gamma, sigma = 1/90, 
                   sigma_m = sigma_m, Ab_susc = 0.75, mAb_susc = mAb_susc, 
                   reduced_shed = 1/92, mu = mu, N_0 = pop_test[i],
                   importation_rate = 0, imp_t = imp_t, 
                   delta = pargrid_foibg_test$delta_v[j], ind1 = ind1, ind2 = ind2, 
                   foi_bg_usr = pargrid_foibg_test$foi_bg[j])  
  
  out <- as.data.frame(replicate(100, x$run(t)[, 354]))
  
  matplot(out[m_t:p_t,], type = "l", col = "black")
  abline(v = c(m_t, p_t))
  
  estab[i] <- sum(out[e_t,] > 0)
  persist[i] <- sum(out[p_t, ] > 0)
  
  if(persist[i] > 1){
    
    idx_persist <- which(out[p_t, ] > 0)
    out_persist <- out[,idx_persist]
    
    for(k in 1:(persist[i])){
      ACF <- acf(out_persist[m_t:p_t, k], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[150:(360*5)], lag = ACF$lag[150:(360*5)])
      # write something in here about making sure acf is above the significance level
      period[k,i] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
  } else {
    period[,i] <- NA
  }
  print(c(i,j))
  
  }
  
  saveRDS(file = paste("results/bplot", pargrid_foibg_test$foi_bg[j],
                       pargrid_foibg_test$delta_v[j], ".rds", sep = "_"), object = period)
}

period_0_0 <- readRDS("results/bplot_0_0_.rds")
period_0001_0 <- readRDS("results/bplot_0.001_0_.rds")
period_0005_0 <- readRDS("results/bplot_0.005_0_.rds")
period_001_0 <- readRDS("results/bplot_0.01_0_.rds")
period_0_1 <- readRDS("results/bplot_0_1_.rds")
period_0001_1 <- readRDS("results/bplot_0.001_1_.rds")
period_0005_1 <- readRDS("results/bplot_0.005_1_.rds")
period_001_1 <- readRDS("results/bplot_0.01_1_.rds")


period <- rbind(period_0_0, period_0001_0,
                period_0005_0, period_001_0,
                period_0_1, period_0001_1,
                period_0005_1, period_001_1) 


cuts <- apply(period, 2, cut, c(0, 350, 370, 710, 730, 1070, 1090, 1430, 1450, Inf), 
              labels= c("other", "1 yr", "other", "2 yrs", "other", "3 yrs", "other", "4 yrs", "other"))
cuts <- as.data.frame(cuts)
names(cuts) <- pop_test
cuts$foibg <- rep(c(0, 0.001, 0.005, 0.01, 0, 0.001, 0.005, 0.01), each = 100)
cuts$delta <- rep(c(0,0,0,0,1,1,1,1), each = 100)
cuts_long <- pivot_longer(cuts, cols = 1:7, names_to = "population", values_to = "periodicity")

cuts_long$population <- factor(cuts_long$population,levels = as.character(pop_test))
cuts_long$periodicity <- factor(cuts_long$periodicity, levels = c(NA, "other", "4 yrs",
                                                                  "3 yrs", "2 yrs", "1 yr"),
                                exclude = NULL)


ggplot(data = cuts_long, aes(fill = periodicity))+
  geom_bar(aes(x = population), position = "stack")+
  scale_fill_viridis_d(direction = -1)+
  theme_minimal()+
  theme(text= element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(delta~foibg)+
  ylab("percentage of runs")
ggsave(filename = "results/test_foibg_period.png")




## what about for the metapopulation model, how is the persistence cf the sp with bg foi?