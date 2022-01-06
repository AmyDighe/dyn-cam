# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model_vax <- odin::odin("models/odin_single_patch_model_vaccination_proportion.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################


# input a value for birth rate (per camel per day) 
alpha <- 0.000565

# input a value for the baseline effective contact rate
beta <- c(0.5, 0.25, 1)

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14
gamma <- 1/duration_infection

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 90 # default = 
sigma <- 1/duration_immunity # default = 

# input a value for the average rate of waning of mAbs
sigma_m <- 4.42/360 # from the catalytic work

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75 # default = 
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

# input an initial population size
N_0 <- 1000000

# input the time period that you wish to run the model for (in days)
time_period <- 9601
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day input that day here
imp_t <- 151  + (360 * seq(0, 4, by = 1))

# set a level of seasonality for births (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# index for summing births for each age class
ind1 <- rep(0,12)
ind2 <- rep(0,12)

for(y in 2:13){
  ind1[y-1] <- 360 - ((y - 1) * 30) + 1 
  ind2[y-1] <- 360 - ((y - 2) * 30)
}

# repeating for 4 years to cover all age classes
ind1 <- rep(ind1, 4)
ind2 <- rep(ind2, 4)

# vaccination parameters

v_gamma <- gamma # rate of recovery from infection in vaccinated animals
v_sigma <- sigma # rate of waning complete infection induced immunity in vaccinated animals
v_sigma_m <- sigma_m # rate of waning of mAbs in vaccinated animals

v_mAb_susc <- mAb_susc # relative susceptibility of vaccinated animals with mAbs
v_susc <- c(0.75, 1) # relative susceptibility of vaccinated naive animals
v_Ab_susc <- c(0.25, 0.5)  # relative susceptibility of vaccinated previously infected animals

v_shed <- c(1/92, 1/10, 1/200) # infectiousness of vaccinated naive animals cf naive unvaccinated
v_reduced_shed <- c(1/724, 1/100, 1/1448) # infectiousness of vaccinated previously infected animals cf naive unvaccinated

# age dependent rates of vaccination

age_idx_vax <- c(1:48)
rho <- c(1/(3*360),1/360, 1/(10*360)) # rate at which vaccine induced immunity wanes per day
coverage <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

###############
## run model ##
###############

last_imp <- imp_t[5] 
pars <- tidyr::expand_grid(beta, rho, nesting(v_susc, v_Ab_susc), nesting(v_shed, v_reduced_shed)) 

incidence <- matrix(nrow = length(coverage), ncol = dim(pars)[1])
output <- matrix(nrow = length(coverage), ncol = dim(pars)[1])

for(j in 1:(dim(pars)[1])){
  
for(i in 1:(length(coverage))){
  vaxp <- c(rep(0, 3), coverage[i], rep(0, 45))
  
x <- sir_model_vax$new(alpha = alpha, beta = pars$beta[j], gamma = gamma, sigma = sigma, sigma_m = sigma_m, Ab_susc = Ab_susc, 
                       mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                       mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, N_0 = N_0,
                       importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2,
                       v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_mAb_susc = v_mAb_susc, v_susc = pars$v_susc[j],
                       v_Ab_susc = pars$v_Ab_susc[j], v_shed = pars$v_shed[j], v_reduced_shed = pars$v_reduced_shed[j], 
                       vaxp = vaxp, rho = pars$rho[j])

out <- as.data.frame(x$run(t)[,c(1, 613, 595, 598, 601, 604)])
n_runs <- 100
out <- as.data.frame(replicate(n_runs, x$run(t)[,c(1, 613, 595, 598, 601, 604)]))

out_inc <- out[6000 : 9600, grep("total_incidence", colnames(out))]
incidence[i, j] <- mean(colSums(out_inc[ , ]))
# output the number of infectious camel days stratified by age group and the total 
out_I1 <- out[6000 : 9600, grep("^I_1", colnames(out))]
out_I2 <- out[6000 : 9600, grep("^I_2", colnames(out))]
out_vI1 <- out[6000 : 9600, grep("vI_1", colnames(out))]
out_vI2 <- out[6000 : 9600, grep("vI_2", colnames(out))]

output[i, j] <- mean(colSums(out_I1)) +
                   reduced_shed * mean(colSums(out_I2)) +
                   pars$v_shed[j] * mean(colSums(out_vI1)) +
                   pars$v_reduced_shed[j] * mean(colSums(out_vI2))
print(paste(i, j, sep = " "))
}
  saveRDS(output[,j], file = paste("results/vax_heatmap/output_beta_", pars$beta[j],
                                   "_rho_", signif(pars$rho[j], 2),
                                   "_susc_", pars$v_susc[j], "_shed_", pars$v_shed[j], ".rds", sep = ""))
  saveRDS(incidence[,j], file = paste("results/vax_heatmap/inc_beta_", pars$beta[j],
                                   "_rho_", signif(pars$rho[j], 2),
                                   "_susc_", pars$v_susc[j], "_shed_", pars$v_shed[j], ".rds", sep = ""))
}
saveRDS(output, file = "results/vax_heatmap/output_all.rds")
saveRDS(incidence, file = "results/vax_heatmap/inc_all.rds")

# get scenario names
scenarios <- vector(length = 6)
for(i in 1:6){
  scenarios[i] <- paste("susc = ", pars$v_susc[i], "/", pars$v_Ab_susc[i],
                        ", ", "shed = ", signif(pars$v_shed[i], 2), "/",
                        signif(pars$v_reduced_shed[i], 2), sep = "")
}


output <- readRDS("results/vax_heatmap/output_all.rds")
outputdf <- as.data.frame(output)
names(outputdf) <- scenarios
outputdf$coverage <- coverage
output_long <- outputdf %>% pivot_longer(cols = starts_with("\U03B2"), 
                                         names_to = "scenario", names_prefix = "V",
                                         values_to = "infectious_camel_days")

ggplot(data = output_long[grep("\U03B2 = 0.25", output_long$scenario),])+
  geom_tile(aes(y = coverage, x = scenario, fill = infectious_camel_days))


baseline <- matrix(rep(output[1,], 10), nrow = 10, ncol = 54, byrow = T) #w/ coverage = 0
output_diff <- 100 * (baseline - output[2:11, ]) / baseline # % reduction 

outputdiff_df <- as.data.frame(output_diff)
outputdiff_df$coverage <- coverage[2:11]
outputdiff_long <- outputdiff_df %>% pivot_longer(cols = starts_with("V"), 
                                         names_to = "scenario", names_prefix = "V",
                                         values_to = "percentage_reduction")
outputdiff_long$scenario <- rep(scenarios, 90)
final_df <- cbind(outputdiff_long, pars)

rho.labs <- c("\U03C1 = 3yrs", "\U03C1 = 1yr", "\U03C1 = 10yrs")
names(rho.labs) <- rho
beta.labs <- c("\U03B2 = 0.5", "\U03B2 = 0.25", "\U03B2 = 1.0")
names(beta.labs) <- beta


ggplot(final_df)+
  geom_tile(aes(y = coverage, x = reorder(scenario, percentage_reduction), fill = percentage_reduction))+
  facet_grid(rho~beta, labeller = labeller(rho = rho.labs, beta = beta.labs))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("scenario")+
  scale_fill_viridis_c(name = "% reduction in\ninfectious camel\ndays")+
  theme(plot.margin = unit(c(0,0,0,3), "cm"))+
  theme(text = element_text(size = 18))
ggsave(filename="results/vax_heatmap/heatmap.png")  

