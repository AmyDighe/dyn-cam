# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model_vax <- odin::odin("models/odin_single_patch_model_vaccination_proportion.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input a value for birth rate (per camel per day) 
alpha <- 0.000565

# input a value for the baseline effective contact rate
beta <- c(0.25, 0.5, 1)

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
v_susc <- 0.75 # relative susceptibility of vaccinated naive animals
v_Ab_susc <- 0.25  # relative susceptibility of vaccinated previously infected animals

v_shed <- 1/92 # infectiousness of vaccinated naive animals cf naive unvaccinated
v_reduced_shed <- 1/724 # infectiousness of vaccinated previously infected animals cf naive unvaccinated

# age dependent rates of vaccination

age_idx_vax <- c(1:48)
vaxp <- rep(0, 49)
rho <- c(1/360, 1/(3*360), 1/(10*360)) # rate at which vaccine induced immunity wanes per day
coverage <- 0.8

###############
## run model ##
###############

last_imp <- imp_t[5]
pars <- expand.grid(beta = beta, rho = rho)
incidence <- matrix(nrow = dim(pars)[1], ncol = length(age_idx_vax) + 1)
age_output <- array(dim = c(length(age_idx_vax) + 1, 49, dim(pars)[1]))

for(i in 4:(dim(pars)[1])){
  beta <- pars$beta[i]
  rho <- pars$rho[i]
  
  for(j in 1:(length(age_idx_vax)+1)){
    vaxp <- rep(0, 49)
    if(j < 49){
      vaxp[age_idx_vax[j]] <- coverage
    }
    # include any user-defined parameters as arguments here
    x <- sir_model_vax$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = sigma_m, Ab_susc = Ab_susc, 
                   mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,
                   mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, N_0 = N_0,
                   importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2,
                   v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_mAb_susc = v_mAb_susc, 
                   v_Ab_susc = v_Ab_susc, v_shed = v_shed, v_reduced_shed = v_reduced_shed, vaxp = vaxp, rho = rho)

    #out <- as.data.frame(x$run(t)[, c(103:(103+48), 250:(250+48), 397:(397+48), 544:(544+48), 613)])
    n_runs <- 100
    out <- as.data.frame(replicate(n_runs, x$run(t)[, c(103:(103+48), 250:(250+48), 397:(397+48), 544:(544+48), 613)])) #save just infectious comps and incidence
    # output total incidence
    out_inc <- out[6000 : 9600, grep("total_incidence", colnames(out))]
    incidence[i, j] <- mean(colSums(out_inc[ , ]))
    # output the number of infectious camel days stratified by age group and the total 
    out_I1 <- out[6000 : 9600, grep("^I\\[", colnames(out))]
    out_I2 <- out[6000 : 9600, grep("^I2", colnames(out))]
    out_vI1 <- out[6000 : 9600, grep("vI\\[", colnames(out))]
    out_vI2 <- out[6000 : 9600, grep("vI2", colnames(out))]

      output <- matrix(colSums(out_I1) +
        reduced_shed * colSums(out_I2) +
        v_shed * colSums(out_vI1) +
        v_reduced_shed * colSums(out_vI2), nrow = n_runs, ncol = 49, byrow = T)

  age_output[j, , i] <- colMeans(output)

print(paste(i,j, sep = " "))
  } #end of j loop
  saveRDS(incidence[i, ], file = paste("results/vax_age_opti/inc/", "incidence_2_", "beta_", pars$beta[i], "_rho_", signif(pars$rho[i], 2), ".rds", sep = ""))
  saveRDS(age_output[, ,i], file = paste("results/vax_age_opti/output/", "output_2_", "beta_", pars$beta[i], "_rho_", signif(pars$rho[i], 2), ".rds", sep = ""))
  } #end of i loop

# trying to vectorise this:
# output <- data.frame(nrow = 49*n_runs, ncol = 1)
# for(k in 1:(49*n_runs)){
#   output[,k] <- sum(out_I1[,k]) +
#     reduced_shed * sum(out_I2[,k]) +
#     v_shed * sum(out_vI1[,k]) +
#     v_reduced_shed * sum(out_vI2[,k])
# }
# names(output) <- names(out_vI2)
# 
# for(age in 1:49){
#   age_output[j, age, i] <- mean(unlist(output[,grep(paste("\\[", age, "\\]", sep = ""), colnames(output))]))
# }
scenario <- vector(length = 9)

for(i in 1:9){
  scenario[i] <- paste("\U03B2", " = ", pars$beta[i], ", ", " \U03C1", " = ", 1/(360 * pars$rho[i]), "yrs", sep = "")
}
scenario_v <- scenario
## Processing RESULTS

# for some reason reading multiple rds files isn't working - it can't find files
inc_1 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.25_rho_0.0028.rds")
inc_2 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.5_rho_0.0028.rds")
inc_3 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_1_rho_0.0028.rds")
inc_4 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.25_rho_0.00093.rds")
inc_5 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.5_rho_0.00093.rds")
inc_6 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_1_rho_0.00093.rds")
inc_7 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.25_rho_0.00028.rds")
inc_8 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_0.5_rho_0.00028.rds")
inc_9 <- readRDS("results/vax_age_opti/inc/incidence_2_beta_1_rho_0.00028.rds")

inc_df <- as.data.frame(cbind(inc_1, inc_2, inc_3, inc_4, inc_5, inc_6, inc_7, inc_8, inc_9))
inc_df$age_targeted <- 1:49
inc_df_long <- gather(inc_df, -age_targeted, key = "scenario", value = "incidence")
inc_df_long$scenario <- rep(scenario, each = 49)
matplot(y = inc_df, type = "l")

output_1 <- readRDS("results/vax_age_opti/output/output_2_beta_0.25_rho_0.0028.rds")
output_2 <- readRDS("results/vax_age_opti/output/output_2_beta_0.5_rho_0.0028.rds")
output_3 <- readRDS("results/vax_age_opti/output/output_2_beta_1_rho_0.0028.rds")
output_4 <- readRDS("results/vax_age_opti/output/output_2_beta_0.25_rho_0.00093.rds")
output_5 <- readRDS("results/vax_age_opti/output/output_2_beta_0.5_rho_0.00093.rds")
output_6 <- readRDS("results/vax_age_opti/output/output_2_beta_1_rho_0.00093.rds")
output_7 <- readRDS("results/vax_age_opti/output/output_2_beta_0.25_rho_0.00028.rds")
output_8 <- readRDS("results/vax_age_opti/output/output_2_beta_0.5_rho_0.00028.rds")
output_9 <- readRDS("results/vax_age_opti/output/output_2_beta_1_rho_0.00028.rds")

output_df <- as.data.frame(rbind(output_1, output_2, output_3, output_4, output_5, output_6, output_7, output_8, output_9))
output_df$total <- rowSums(output_df)
output_df$scenario <- rep(scenario, each = 49)
output_df$age_targeted <- as.numeric(rep(1:49, 9))
output_df_tot <- data.frame(total = output_df$total, 
                            scenario = output_df$scenario, 
                            age_targeted = output_df$age_targeted)

p1 <- ggplot(data = output_df_tot%>% filter(age_targeted<49))+
  geom_line(aes(x = age_targeted, y = total, col = scenario), lwd = 1)+
  theme_minimal()+
  ylab("total infectious camel days") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D")+
  geom_hline(yintercept = output_df_tot$total[output_df_tot$age_targeted == 49],
             lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))

ggsave(filename = "vax_age.png")

p2 <- ggplot(data = inc_df_long%>%filter(age_targeted<49))+
  geom_line(aes(x = age_targeted, y = incidence, col = scenario), lwd = 1)+
  theme_minimal()+
  ylab("incidence") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D")+
  geom_hline(yintercept = inc_df_long$incidence[inc_df_long$age_targeted == 49],
             lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))

ggsave(filename = "inc_vax_age.png")

## looking at whether infectious camel days increases for any age group
output_list <- list(output_1, output_2, output_3, output_4, 
                    output_5, output_6, output_7, output_8, output_9)

output_comp <- list(length = 9)
names(output_comp) <- scenario
for(i in 1:9){
  op <- output_list[[i]]
  output_comp[[i]] <- op[1:48, -(50:53)]- op[rep(49,48),-(50:53)]
}

output_cf <- do.call(rbind.data.frame, output_comp)
#output_cf$scenario <- rep(scenario, each = 48)
output_cf$age_targeted <- as.factor(rep(1:48, 9))
output_cf_long <- gather(output_cf, -c(scenario, age_targeted), key = "age_class", value = "infectiousness")
output_cf_long$age_class <- rep(1:49, each = (48*9))

ggplot(data = output_cf_long)+
  geom_line(aes(x = age_class, y = infectiousness, col = age_targeted))+
  scale_color_viridis_d()+
  facet_wrap(~scenario)+
  geom_hline(yintercept = 0, col = "red")+
  xlab("age (months)") +
  ylab("change in infectious camel days")+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))
ggsave(filename = "output_vax_age_by_age.png")

dummy_df <- data.frame(age_targeted = as.factor(1:48), Z = as.numeric(1:48))
for(i in 1:9){
ggplot(data = output_cf_long%>% filter(scenario == scenario_v[i]))+
  geom_line(aes(x = age_class, y = infectiousness))+
  scale_color_viridis_d()+
  facet_wrap(~age_targeted)+
  geom_vline(data = dummy_df, aes(xintercept = Z), lty = 2)+
  geom_hline(yintercept = 0, col = "red")+
  xlab("age (months)") +
  ylab("change in infectious camel days")+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
    ggtitle(scenario_v[i])

ggsave(filename = paste("output_vax_scenario", i, ".png", sep = ""))
}
