sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input an initial population size
N_0 <- 1e+7

# input the time period that you wish to run the model for (in days)
time_period <- 12600 
t <- seq(0:time_period)

# if you want importation to occur on a specific day input that day here
imp_t <- 151  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

###############
## run model ##
###############
scenario <- names(beta_list)
last_imp <- imp_t[5]

for(j in 1:18){
  
  Ab_susc <- par_grid_R0$Ab_susc[j]
  sigma <- par_grid_R0$sigma[j]
  reduced_shed <- par_grid_R0$red_shed[j]
  beta_vector <- beta_list[[j]]
  mean_foi <- vector(length = length(beta_vector))
  
for(i in 1:(length(beta_vector))){
beta <- beta_vector[i]

x <- sir_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = sigma_m, Ab_susc = Ab_susc, 
               mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu = mu, N_0 = N_0,
               importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2,
               foi_bg_usr = 0)

if(beta<1){
  nruns <- 1000
} else {
  nruns <- 100
}
out <- as.data.frame(replicate(nruns, x$run(t)[, c(349, 352, 354, 356, 373)]))
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