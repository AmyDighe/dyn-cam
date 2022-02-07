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

# set level of connectivity between patches
connectivity <- c(0.01, 0.1, 0.05, 0.001)

###################
## run the model ##
###################

###########################
###########################
#### test connectivity ####
###########################
###########################

par_grid_metapop_test <- expand_grid(beta = 5.5/14, 
                                     waning = 1/90,
                                     shedding = 0.25,
                                     seasonality = 1,
                                     susc = 0.75,
                                     pop = c(1e+03, 5e+03, 7.5e+03)/25,
                                     connectivity = c(0.1, 0.05, 0.01, 0.005, 0.001))

n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta = par_grid_metapop_test$seasonality[1], N_0 = 10,#par_grid_metapop_test$pop[1],
                      mu, N_age, n_r = n_r, n_c = n_c)

storage.mode(S_ini_p) <- "double"

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, 
                                               beta = par_grid_metapop_test$beta[1],
                                               gamma = gamma, 
                                               sigma = par_grid_metapop_test$waning[1],
                                               sigma_m = sigma_m, 
                                               Ab_susc = par_grid_metapop_test$susc[1],
                                               mAb_susc = mAb_susc, 
                                               reduced_shed = par_grid_metapop_test$shedding[1], mu = mu, 
                                               N_0 = 10,#par_grid_metapop_test$pop[1], 
                                               importation_rate = importation_rate, 
                                               imp_t = imp_t, 
                                               delta =  par_grid_metapop_test$seasonality[1], 
                                               connectivity = connectivity[1],
                                               correction_ex = correction_ex,
                                               S_ini_p = S_ini_p,
                                               foi_bg_usr = foi_bg_usr), 
                                   step = 0, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Ntot)
msirs_model$set_index(msirs_model$info()$index$Itot) # just extract Itot

# TEST
steps <- seq(1, 3600, by = 10)
state <- msirs_model$simulate(steps)
traces <- t(drop(state[1,,]))
matplot(steps, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)



# PERSISTENCE ANALYSIS
estab <- vector(length = (dim(par_grid_metapop_test)[1]))
m_persist <- estab
persist <- estab
steps <- seq(0, 12600, by = 10)

S_ini_list <- list(length = dim(par_grid_metapop_test)[1]) ## stochastic initialization
for(i in 1:(dim(par_grid_metapop_test)[1])){
  S_ini_list[[i]] <- stoch_init(alpha, delta = par_grid_metapop_test$seasonality[i],
                                N_0 = par_grid_metapop_test$pop[i],
                                mu, N_age, n_r = n_r, n_c = n_c)
}

e_t <- 13 * 360
m_t <- 25 * 360
p_t <- 35 * 360

for(i in 1:(dim(par_grid_metapop_test)[1])){
  
  msirs_model$update_state(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                       alpha = alpha, 
                                       beta = par_grid_metapop_test$beta[i],
                                       gamma = gamma, 
                                       sigma = par_grid_metapop_test$waning[i],
                                       sigma_m = sigma_m, 
                                       Ab_susc = par_grid_metapop_test$susc[i],
                                       mAb_susc = mAb_susc, 
                                       reduced_shed = par_grid_metapop_test$shedding[i], mu = mu, 
                                       N_0 = par_grid_metapop_test$pop[i], 
                                       importation_rate = importation_rate, 
                                       imp_t = imp_t, 
                                       delta =  par_grid_metapop_test$seasonality[i], 
                                       connectivity = par_grid_metapop_test$connectivity[i],
                                       correction_ex = correction_ex,
                                       S_ini_p = S_ini_p,
                                       foi_bg_usr = foi_bg_usr), step = 0)
  
  state <- msirs_model$simulate(steps)
  out <- as.data.frame(t(drop(state[1,,])))
   
  estab[i] <- sum(out[e_t/10,] > 0)
  m_persist[i] <- sum(out[m_t/10, ] > 0)
  persist[i] <- sum(out[p_t/10, ] > 0)
  print(i)

}

saveRDS(file = "results/persistence/dustmeta_persist_test.rds", 
        object = data.frame(persist = persist, m_persist = m_persist, estab = estab))

persist_test <- cbind(par_grid_metapop_test, persist, estab, p_e = persist/estab)



d <- readRDS("results/dustmeta_persist20.01.rds")
results <- cbind(d, par_grid_metapop)
results <- results %>% mutate(persist_estab = 100 * persist/estab)
#results$persist_estab[is.nan(results$persist_estab)]<-0

seasonality.labs <- c("\U03B4 = 0.0", "\U03B4 = 0.5", "\U03B4 = 1.0")
names(seasonality.labs) <- unique(par_grid_metapop$seasonality)

# persistance by population size
ggplot(results)+
  geom_line(aes(x = log10(25*pop), y = persist, col = as.character(beta)),
            size = 1)+
  facet_grid(~seasonality, labeller = labeller(seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/persistence_metapop.png")

# persist/estab
ggplot(results)+
  geom_line(aes(x = log10(25*pop), y = persist_estab, col = as.character(beta)),
            size = 1)+
  facet_grid(~seasonality, labeller = labeller(seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/persist_estab_metapop.png")

# established
ggplot(results)+
  geom_line(aes(x = log10(25*pop), y = estab, col = as.character(beta)),
            size = 1)+
  facet_grid(~seasonality, labeller = labeller(seasonality = seasonality.labs))+
  theme_minimal()+
  scale_color_viridis_d(name = "beta")+
  theme(text = element_text(size = 20))+
  xlab("log(population size)")+
  ylab("% of model runs in which endemic transmission persists")
ggsave(filename = "figs/estab_metapop.png")

## make a table comparing critical community size for sp and metapop
e <- readRDS("results/persist_estab_period_sp.rds")
results_sp <- cbind(e, par_grid)
results_sp <- results_sp %>% mutate(persist_estab = 100 * persist/estab,
                                    tag = "sp") %>% 
  filter(import_time ==151, waning == 1/90)%>%
  select(-period_av, -period_sd, -period_med)

results_mp <- results%>%
  mutate(pop = pop * 25,
         tag = "mp")

results_comp <- rbind(results_sp, results_mp)

results_comp <- results_comp %>%
  select(-waning, -import_time)

results_comp <- results_comp[order(results_comp$beta, results_comp$seasonality),]
results_comp_wide <- pivot_wider(results_comp, names_from = tag, 
                                  values_from = c(persist,m_persist, estab, persist_estab))
write.csv(file = "results/persistence_cf.csv", results_comp_wide)
