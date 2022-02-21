
sir_model_vax <- odin::odin("models/odin_single_patch_model_vaccination_proportion.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 90 # default = 
sigma <- 1/duration_immunity # default = 

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 0.75 # default = 

# input an initial population size
N_0 <- 1000000

# input the time period that you wish to run the model for (in days)
time_period <- 10800
t <- seq(0:time_period)

# set a level of seasonality for births (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# vaccination parameters

v_gamma <- gamma # rate of recovery from infection in vaccinated animals
v_sigma <- sigma # rate of waning complete infection induced immunity in vaccinated animals
v_sigma_m <- sigma_m # rate of waning of mAbs in vaccinated animals
v_mAb_susc <- mAb_susc # relative susceptibility of vaccinated animals with mAbs
coverage <- 0.8

###############
## run model ##
###############

incidence <- matrix(nrow = dim(pars)[1], ncol = 25)
incidence_l <- incidence
incidence_u <- incidence
age_output <- array(dim = c(25, 26, dim(pars)[1])) #dim(age targeted, age, scenario)
age_output_l <- age_output
age_output_u <- age_output

for(i in c(37:48)){
  beta <- pars$beta[i]
  reduced_shed <- pars$reduced_shed[i]
  v_shed <- pars$v_shed[i]
  v_reduced_shed <- pars$v_reduced_shed[i]
  
  for(j in 1:25){
    vaxp <- rep(0, 49)
    if(j <25){
      vaxp[2*j] <- coverage}
      
    # include any user-defined parameters as arguments here
    x <- sir_model_vax$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = sigma_m, Ab_susc = Ab_susc, 
                   mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu = mu, N_0 = N_0,
                   importation_rate = importation_rate, imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2,
                   v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_susc  = pars$v_susc[i], v_mAb_susc = v_mAb_susc, 
                   v_Ab_susc = pars$v_Ab_susc[i], v_shed = v_shed, v_reduced_shed = v_reduced_shed, vaxp = vaxp, rho = pars$rho[i],
                   foi_bg_usr = foi_bg_usr)

    if(beta < 0.5){
      n_runs <- 200
    } else {
      n_runs <- 100
    }
    
    out <- as.data.frame(replicate(n_runs, x$run(t)[, c(103, seq(103+1, 103+47, by = 2),103+48, 
                                                        250, seq(250+1, 250+47, by = 2), 250+48,
                                                        397, seq(397+1, 397+47, by = 2), 397+48,
                                                        544, seq(544+1, 544+47, by = 2), 544+48,
                                                        613)])) #save just infectious comps and incidence
    # output total incidence
    out_inc <- out[7200 : 10800, grep("total_incidence", colnames(out))]
    incidence[i, j] <- mean(colSums(out_inc[ , ]))
    incidence_l[i, j] <- quantile(colSums(out_inc[ , ]), 0.05)
    incidence_u[i, j] <- quantile(colSums(out_inc[ , ]), 0.95)
    # output the number of infectious camel days stratified by age group and the total 
    out_I1 <- out[7200 : 10800, grep("^I\\[", colnames(out))]
    out_I2 <- out[7200 : 10800, grep("^I2", colnames(out))]
    out_vI1 <- out[7200 : 10800, grep("vI\\[", colnames(out))]
    out_vI2 <- out[7200 : 10800, grep("vI2", colnames(out))]

      output <- matrix(colSums(out_I1) +
        reduced_shed * colSums(out_I2) +
        v_shed * colSums(out_vI1) +
        v_reduced_shed * colSums(out_vI2), nrow = n_runs, ncol = 26, byrow = T)

  age_output[j, , i] <- colMeans(output)
  age_output_l[j, , i] <- apply(output, 2, quantile, probs = 0.05, na.rm = TRUE)
  age_output_u[j, , i] <- apply(output, 2, quantile, probs = 0.95, na.rm = TRUE)

print(paste(i,j, incidence[i,j], sep = " "))
  } #end of j loop
  saveRDS(incidence[i, ], file = paste("results/vax_age_opti/inc/", "incidence_3_", i, ".rds", sep = ""))
  saveRDS(age_output[, ,i], file = paste("results/vax_age_opti/output/", "output_3_", i, ".rds", sep = ""))
  saveRDS(incidence_l[i, ], file = paste("results/vax_age_opti/inc/", "incidence_3l_", i, ".rds", sep = ""))
  saveRDS(age_output_l[, ,i], file = paste("results/vax_age_opti/output/", "output_3l_", i, ".rds", sep = ""))
  saveRDS(incidence_u[i, ], file = paste("results/vax_age_opti/inc/", "incidence_3u_", i, ".rds", sep = ""))
  saveRDS(age_output_u[, ,i], file = paste("results/vax_age_opti/output/", "output_3u_", i, ".rds", sep = ""))
  } #end of i loop

## Processing RESULTS
n_scenario <- 48
inc_mat <- matrix(NA, nrow = n_scenario, ncol = 25)
for(i in 1:12){
  inc_mat[i,] <- readRDS(paste("results/vax_age_opti/inc/", 
                               "incidence_3_", 
                               "beta_", 
                               pars$beta[i], 
                               "_rho_", 
                               signif(pars$rho[i], 2), 
                               ".rds", sep = ""))
  inc_mat_l[i,] <- readRDS(paste("results/vax_age_opti/inc/", 
                                 "incidence_3_", 
                                 "beta_", 
                                 pars$beta[i], 
                                 "_rho_", 
                                 signif(pars$rho[i], 2), 
                                 ".rds", sep = ""))
}

for(i in 13:48){
  inc_mat[i,] <- readRDS(paste("results/vax_age_opti/inc/", 
                               "incidence_3_", 
                               i,
                               ".rds", sep = ""))
}

scenario <- c(rep(1, 12), rep(2, 12), rep(3, 12), rep(4, 12))
beta_rho <- vector(length = 48)
for(i in 1:48){
beta_rho[i] <- paste("\U03B2", " = ", format(round(pars$beta[i],2), nsmall = 2), 
                     ", \U03C1", " = ", 1/(360*pars$rho[i]), " yrs", sep = "")}
inc_df <- as.data.frame(inc_mat)
inc_df$scenario <- scenario
inc_df$beta_rho <- beta_rho
inc_df_long <- pivot_longer(inc_df, c(-scenario, -beta_rho), names_to = "age_targeted", names_prefix = "V",
                            values_to = "incidence")
inc_df_long$age_targeted <- factor(2*as.numeric(inc_df_long$age_targeted), levels = c(seq(2,50, 2)))

inc_main <- inc_df_long%>%filter(scenario != 2)
beta_rho_ordered_main <- inc_main%>%filter(age_targeted == 48)
beta_rho_ordered_main <- beta_rho_ordered_main[order(beta_rho_ordered_main$incidence),]
inc_main$beta_rho <- factor(inc_main$beta_rho, levels = unique(beta_rho_ordered_main$beta_rho))

inc_s2 <- inc_df_long%>%filter(scenario==2)
beta_rho_ordered_s2 <- inc_s2%>%filter(age_targeted == 48)
beta_rho_ordered_s2 <- beta_rho_ordered_s2[order(beta_rho_ordered_s2$incidence),]
inc_s2$beta_rho <- factor(inc_s2$beta_rho, levels = unique(beta_rho_ordered_s2$beta_rho))

ggplot(data = inc_main%>%filter(age_targeted !=50, scenario ==1))+
  geom_line(aes(x = age_targeted, y = incidence/N_0, group = beta_rho,
                col = beta_rho), lwd = 1)+
  #geom_ribbon(aes(ymin = incidence_l, ymax = incidence_u))
  theme_minimal()+
  ylab("cumulative incidence over 10 years per animal") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D", name = "parameters")+
  geom_hline(yintercept = ((inc_df_long%>%filter(scenario == 1, age_targeted == 50))$incidence)/(N_0),
              lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))
ggsave(filename = "figs/inc_vax_age1.png")

ggplot(data = inc_s2%>%filter(age_targeted !=50))+
  geom_line(aes(x = age_targeted, y = incidence/N_0, group = beta_rho,
                col = beta_rho), lwd = 1)+
  #geom_ribbon(aes(ymin = incidence_l, ymax = incidence_u))
  theme_minimal()+
  ylab("cumulative incidence over 10 years per animal") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D", name = "parameters")+
  geom_hline(yintercept = ((inc_df_long%>%filter(scenario == 2, age_targeted == 50))$incidence)/(N_0),
             lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))
ggsave(filename = "figs/inc_vax_age2.png")

ggplot(data = inc_main%>%filter(age_targeted !=50, scenario ==3))+
  geom_line(aes(x = age_targeted, y = incidence/N_0, group = beta_rho,
                col = beta_rho), lwd = 1)+
  #geom_ribbon(aes(ymin = incidence_l, ymax = incidence_u))
  theme_minimal()+
  ylab("cumulative incidence over 10 years per animal") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D", name = "parameters")+
  geom_hline(yintercept = ((inc_df_long%>%filter(scenario == 3, age_targeted == 50))$incidence)/(N_0),
             lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))
ggsave(filename = "figs/inc_vax_age3.png")

ggplot(data = inc_main%>%filter(age_targeted !=50, scenario ==4))+
  geom_line(aes(x = age_targeted, y = incidence/N_0, group = beta_rho,
                col = beta_rho), lwd = 1)+
  #geom_ribbon(aes(ymin = incidence_l, ymax = incidence_u))
  theme_minimal()+
  ylab("cumulative incidence over 10 years per animal") +
  xlab("age targeted (months)") +
  scale_color_viridis_d(option = "D", name = "parameters")+
  geom_hline(yintercept = ((inc_df_long%>%filter(scenario == 4, age_targeted == 50))$incidence)/N_0,
             lty = 2, col = "grey", lwd = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_text(vjust= -1))+
  theme(axis.title.y = element_text(vjust = +3))+
  theme(plot.margin = unit(c(0,1,1,1), "cm"))

ggsave(filename = "figs/inc_vax_age4.png")



output_list <- vector(mode = "list", length = n_scenario)
for(i in 1:12){
  output_list[[i]] <- readRDS(paste("results/vax_age_opti/output/", 
                                "output_3_", "beta_", 
                                pars$beta[i], "_rho_", 
                                signif(pars$rho[i], 2),
                                ".rds", sep = ""))
}

for(i in 13:n_scenario){
  output_list[[i]] <- readRDS(paste("results/vax_age_opti/output/", 
                                    "output_3_", i,
                                    ".rds", sep = ""))
}


output_comp <- vector(mode = "list", length = n_scenario)
names(output_comp) <- 1:n_scenario
for(i in 1:n_scenario){
  op <- output_list[[i]]
  output_comp[[i]] <- op[1:24, ]- op[rep(25,24), ]
}

output_cf <- do.call(rbind.data.frame, output_comp)
output_cf$age_targeted <- as.factor(rep(seq(2,48,2), n_scenario))
output_cf$scenario <- as.factor(rep(1:n_scenario, each = 24))
output_cf_long <- pivot_longer(output_cf, c(-age_targeted, -scenario), 
                               names_to = "age_class", values_to = "infectiousness")
output_cf_long$age_class <- rep(c(1, seq(2, 48, 2), 49), times = (24*n_scenario))

scenario_v <- 1:12
dummy_df <- data.frame(age_targeted = as.factor(seq(2, 48, 2)), Z = as.numeric(seq(2, 48, 2)))
for(i in 1:n_scenario){
  p <- ggplot(data = output_cf_long%>% filter(scenario == scenario_v[i]))+
    geom_line(aes(x = age_class, y = infectiousness))+
    scale_color_viridis_d()+
    facet_wrap(~age_targeted)+
    geom_vline(data = dummy_df, aes(xintercept = Z), lty = 2)+
   # geom_hline(yintercept = 0, col = "red")+
    xlab("age (months)") +
    ylab("change in infectious camel days")+
    theme_minimal()+
    theme(text = element_text(size = 20))+
    theme(axis.title.x = element_text(vjust= -1))+
    theme(axis.title.y = element_text(vjust = +3))
  print(p)
}


