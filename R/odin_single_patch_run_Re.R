sir_model <- odin::odin("models/odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input an initial population size
N_0 <- 1e+7

# if you want importation to occur on a specific day input that day here
imp_t <- 151  + (360 * seq(0, 4, by = 1))
last_imp <- imp_t[5]

# input the time period that you wish to run the model for (in days)
time_period <- 50*360 + last_imp 
t <- seq(0:time_period)

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

###############
## run model ##
###############
core <- c(2,3,5,8,9,11,14,15,17)

scenario_core <- scenario[core]
par_grid_R0_core <- par_grid_R0[core,]
#Se <- matrix(nrow = 23, ncol = 9)
#colnames(Se) <- scenario_core

beta_esti_lower_core <- beta_esti_lower[,core]

for(j in 5:9){
  Ab_susc <- par_grid_R0_core$Ab_susc[j]
  sigma <- par_grid_R0_core$sigma[j]
  reduced_shed <- par_grid_R0_core$red_shed[j]
  beta_vector <- beta_esti_lower_core[,j]
  mean_foi <- vector(length = length(beta_vector))
  
  for(i in 1:23){
    beta <- beta_vector[i]
    
    x <- sir_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma,
                       sigma_m = sigma_m, Ab_susc = Ab_susc, mAb_susc = mAb_susc,
                       reduced_shed = reduced_shed, mu = mu, N_0 = N_0,
                       importation_rate = importation_rate, imp_t = imp_t,
                       delta = delta, ind1 = ind1, ind2 = ind2, foi_bg_usr = 0)
    
    if(beta<0.15){
      nruns <- 100
    } else {
      nruns <- 10
    }

    out <- as.data.frame(replicate(nruns, x$run(t)[, c(354, 373)]))
    out_Itot <- out[,grep("Itot", colnames(out))]
    out_Se <- out[, grep("Se", colnames(out))]
    idx_persist <- which(out_Itot[last_imp + (50*360),] > 0)
    out_Se_persist <- out_Se[,idx_persist]
    
    if(length(idx_persist) == 0){
      Se[i,j] <- NA
    } else if(length(idx_persist) > 1){
      Se[i,j] <- mean(as.matrix(out_Se_persist))
    } else if(length(idx_persist) == 1) {
      Se[i,j] <- mean(out_Se_persist)
    }
    
    print(paste("i =", i, "Se= ", Se[i,j], sep = " ")) }
  
  saveRDS(file = "results/Se_lower.rds", object = Se)
  
}
