# MSIRS2I2 dust METAPOP 2021

###################
## compile model ##
###################

library(odin.dust)
msirs_meta_dust <- odin.dust::odin_dust("models/dustmeta.R")

source("R/stochastic_init.R") # function for starting model at zoographic equilibrium

#############################
## user defined variables: ##
#############################

n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of subpops

# correction for balancing ex vs int foi

correction_ex <- matrix(c(2, 3, 3, 3, 2,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          2, 3, 3, 3, 2), nrow = 5, ncol = 5)

connectivity <- c(0.01, 0.001, 0.1)

###################
## run the model ##
###################

###########################
###########################
#### test connectivity ####
###########################
###########################


n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta = par_grid_metapop$seasonality[1], 
                      N_0 = par_grid_metapop$pop[1],
                      mu, N_age, n_r = n_r, n_c = n_c)

storage.mode(S_ini_p) <- "double"

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, 
                                               beta = par_grid_metapop$beta[1],
                                               gamma = gamma, 
                                               sigma = par_grid_metapop$waning[1],
                                               sigma_m = sigma_m, 
                                               Ab_susc = par_grid_metapop$susc[1],
                                               mAb_susc = mAb_susc, 
                                               reduced_shed = par_grid_metapop$shedding[1], mu = mu, 
                                               N_0 = par_grid_metapop$pop[1], 
                                               importation_rate = importation_rate, 
                                               imp_t = imp_t, 
                                               delta =  par_grid_metapop$seasonality[1], 
                                               connectivity = connectivity[1],
                                               correction_ex = correction_ex,
                                               S_ini_p = S_ini_p,
                                               foi_bg_usr = foi_bg_usr), 
                                   step = 0, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Itot) # just extract Itot
#msirs_model$set_index(msirs_model$info()$index$Ntot) # just extract Ntot
steps <- seq(1, 3600, by = 10)
state <- msirs_model$simulate(steps)
traces <- t(drop(state[1,,]))
matplot(steps, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)



# PERSISTENCE ANALYSIS

steps <- seq(0, 12600, by = 10)

S_ini_list <- list(length = dim(par_grid_metapop)[1]) ## stochastic initialization
for(i in 1:(dim(par_grid_metapop)[1])){
  S_ini_p <- stoch_init(alpha, delta = par_grid_metapop$seasonality[i],
                        N_0 = par_grid_metapop$pop[i],
                        mu, N_age, n_r = n_r, n_c = n_c)
  storage.mode(S_ini_p) <- "double"
  S_ini_list[[i]] <- S_ini_p
}

e_t <- 13 * 360
m_t <- 25 * 360
p_t <- 35 * 360

for(j in 2:(length(connectivity))){
  estab <- vector(length = (dim(par_grid_metapop)[1]))
  m_persist <- estab
  persist <- estab

for(i in 1:(dim(par_grid_metapop)[1])){
  msirs_model$update_state(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                       alpha = alpha, 
                                       beta = par_grid_metapop$beta[i],
                                       gamma = gamma, 
                                       sigma = par_grid_metapop$waning[i],
                                       sigma_m = sigma_m, 
                                       Ab_susc = par_grid_metapop$susc[i],
                                       mAb_susc = mAb_susc, 
                                       reduced_shed = par_grid_metapop$shedding[i], mu = mu, 
                                       N_0 = par_grid_metapop$pop[i], 
                                       importation_rate = importation_rate, 
                                       imp_t = imp_t, 
                                       delta =  par_grid_metapop$seasonality[i], 
                                       connectivity = connectivity[j],
                                       correction_ex = correction_ex,
                                       S_ini_p = S_ini_list[[i]],
                                       foi_bg_usr = foi_bg_usr), step = 0)
  
  state <- msirs_model$simulate(steps)
  out <- as.data.frame(t(drop(state[1,,])))
   
  estab[i] <- sum(out[e_t/10,] > 0)
  m_persist[i] <- sum(out[m_t/10, ] > 0)
  persist[i] <- sum(out[p_t/10, ] > 0)
  print(paste("i = ", i, "persist = ",persist[i], sep = " "))

}

saveRDS(file = paste("results/persistence/dustmeta_persist_", connectivity[j], ".rds", sep = ""), 
        object = data.frame(persist = persist, m_persist = m_persist, estab = estab))
}

