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

N_age <- 49 # number of age classes
n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of subpops

# input a value for birth rate (per camel per day) 
alpha <- 0.000565

# input a value for the baseline effective contact rate
beta <- 0.5

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14
gamma <- 1/duration_infection

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 60 # default = 
sigma <- 1/duration_immunity # default = 

# input a value for the average rate of waning of mAbs
sigma_m <- 4.42/365 # from the catalytic work

# input a value between 0 and 1 for susceptibility of individuals with mAbs/Abs
## 0 would mean mAbs/Abs afford complete protection
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from shedding in Alharbi 

# input values for the age dependent removal rate - balance birthrate
mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

mu <- vector(length = N_age)
mu[1:12] <- mu_1st_yr
mu[13:24] <- mu_2nd_yr
mu[25:36] <- mu_3rd_yr
mu[37:48] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day
imp_t <- 151  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# set level of connectivity between patches
connectivity <- 0.05

#initial pop per patch
N_0 <- 400000

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta, N_0, mu, N_age, n_r = n_r, n_c = n_c)


###################
## run the model ##
###################

n_particles <- 100L

msirs_model <- msirs_meta_dust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, beta = beta, gamma = gamma, 
                                               sigma = sigma, sigma_m = sigma_m, 
                                               Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                               reduced_shed = reduced_shed, mu = mu, N_0 = N_0, 
                                               importation_rate = importation_rate, imp_t = imp_t, 
                                               delta = delta, connectivity = connectivity, 
                                               S_ini_p = S_ini_p), 
                                   step = 1, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Itot) # just extract Itot

# TEST
steps <- seq(1, 600, by = 5)
state <- msirs_model$simulate(steps)
traces <- t(drop(state[1,,]))
matplot(steps, traces, type = "l", lty = 1, col = "#00000022",
        xlab = "Time", ylab = "Number infected (I)")
lines(steps, rowMeans(traces), col = "red", lwd = 2)

# PERSISTENCE ANALYSIS
estab <- vector(length = (dim(par_grid_metapop)[1]))
m_persist <- estab
persist <- estab
steps <- seq(1, 19710, by = 10)
connectivity <- 0.01

S_ini_list <- list(length = 10) ## stochastic initialization
for(i in 1:length(pop)){
  S_ini_list[[i]] <- stoch_init(alpha, delta, N_0 = unique(par_grid_metapop$pop)[i],
                                mu, N_age, n_r = n_r, n_c = n_c)
}
names(S_ini_list) <- as.character(unique(par_grid_metapop$pop))

for(i in 1:(dim(par_grid_metapop)[1])){
  imp_t <- par_grid_metapop$import_time[i] + (360 * seq(0,4, by = 1))
  e_t <- imp_t[5] + (3 * 360)
  m_t <- imp_t[5] + (25 * 360)
  p_t <- imp_t[5] + (50 * 360)
  
  msirs_model$update_state(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                       alpha = alpha, beta = par_grid_metapop$beta[i], gamma = gamma, 
                                       sigma = par_grid_metapop$waning[i], sigma_m = sigma_m, 
                                       Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
                                       reduced_shed = reduced_shed, mu = mu, N_0 = par_grid_metapop$pop[i], 
                                       importation_rate = importation_rate, imp_t = imp_t, 
                                       delta = par_grid_metapop$seasonality[i], 
                                       connectivity = connectivity, 
                                       S_ini_p = S_ini_list[[as.character(par_grid_metapop$pop[i])]]),
                           step = 1)
  
  state <- msirs_model$simulate(steps)
  out <- as.data.frame(t(drop(state[1,,])))
  
  estab[i] <- sum(out[e_t/10,] > 0)
  m_persist[i] <- sum(out[m_t/10, ] > 0)
  persist[i] <- sum(out[p_t/10, ] > 0)
  print(i)
}

saveRDS(file = paste("results/dustmeta_persist", connectivity, ".rds", sep = ""), 
        object = data.frame(persist = persist, m_persist = m_persist, estab = estab))

d <- readRDS("results/dustmeta_persist0.01.rds")
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
