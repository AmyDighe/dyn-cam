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

# input a value for the baseline effective contact rate
beta <- 0.5

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 90 # default = 
sigma <- 1/duration_immunity # default = 

# input a value between 0 and 1 for susceptibility of individuals with mAbs/Abs
## 0 would mean mAbs/Abs afford complete protection
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75 # default = 

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from shedding in Alharbi 

# if rather than a rate you want importation to occur on a specific day
imp_t <- 151  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# set level of connectivity between patches
connectivity <- 0.05

# correcting balance of internal and external foi
correction_ex <- matrix(data = c(2, 3, 3, 3, 2,
                                 3, 4, 4, 4, 3,
                                 3, 4, 4, 4, 3,
                                 3, 4, 4, 4, 3,
                                 2, 3, 3, 3, 2),
                        nrow = 5, ncol = 5, byrow = TRUE)
#initial pop per patch
N_0 <- 4e+05

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta, N_0, mu, N_age, n_r = n_r, n_c = n_c)
storage.mode(S_ini_p) <- "double"

###################
## run the model ##
###################
time_period <- 14400
n_particles <- 100L

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, beta = beta, gamma = gamma, 
                                               sigma = sigma, sigma_m = sigma_m, 
                                               Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                               reduced_shed = reduced_shed, mu = mu, N_0 = N_0, 
                                               importation_rate = importation_rate, imp_t = imp_t, 
                                               delta = delta, connectivity = connectivity, 
                                               correction_ex = correction_ex,
                                               S_ini_p = S_ini_p, foi_bg_usr = 0), 
                         step = 1, 
                         n_particles = n_particles, 
                         n_threads = 2L, 
                         seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Itot) # just extract I1
steps <- seq(1, 600, by = 5)
tic()
state <- msirs_model$simulate(steps)
toc()
# test plot of a single age group in a single patch
traces <- t(drop(state[1,,]))
matplot(steps, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)

# rearrange state to get the high dimensional array back
ix <- seq_along(msirs_model$info()$index$I)

dimensions <- msirs_model$info()$dim$I
test_x <- array(state[ix, , ], c(dimensions, dim(state)[-1]))

## output pulse heatmap gif to check connectivity is working

x_1p <- test_x[,,,1,] # just 1/100 of the runs
x_sum_age <- rowSums(x_1p) # sum over all ages in each patch
x_sum_age_df <- as.data.frame(apply(x_sum_age, 2, c)) # flatten array
x_sum_age_df$time <- rep(steps, each = 5) #indicate time
x_sum_age_df$nr <- rep(1:5, times = 120) #indicate row of grid
x_long <- tidyr::pivot_longer(x_sum_age_df, c(-time, -nr), names_to = "nc", 
                              names_prefix = "V", values_to = "I1")

a2 <- ggplot(data = x_long%>%filter(time>150))+
  geom_tile(aes(x = nc, y = nr, fill = I1))+
  transition_time(time)
animate(a2, nframes = length(unique(time)))

anim_save(filename = "results/dust_meta_intro.gif")

traces <- t(drop(state[1,,]))
matplot(steps, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)

## checking connectivity - what are some reasonable values?

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, beta = beta, gamma = gamma, 
                                               sigma = sigma, sigma_m = sigma_m, 
                                               Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                               reduced_shed = reduced_shed, mu = mu, N_0 = N_0, 
                                               importation_rate = importation_rate, imp_t = imp_t, 
                                               delta = delta, connectivity = 0.01, 
                                               S_ini_p = S_ini_p), 
                                   step = 1, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Itot) # just extract I1

msirs_model$update_state(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                     alpha = alpha, beta = beta, gamma = gamma, 
                                     sigma = sigma, sigma_m = sigma_m, 
                                     Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                     reduced_shed = reduced_shed, mu = mu, N_0 = N_0, 
                                     importation_rate = importation_rate, imp_t = imp_t, 
                                     delta = delta, connectivity = 0.1, 
                                     S_ini_p = S_ini_p), step = 1)

steps <- seq(1, 6000, by = 10)
tic()
state <- msirs_model$simulate(steps)
toc() #411.64
traces <- t(drop(state[1,,]))
matplot(steps, out, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)
