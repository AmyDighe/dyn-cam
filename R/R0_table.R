## estimating R0 from foi data 
mean_foi_MSIS1 <- readRDS("results/mean_foi_MSIS1.rds")
mean_foi_MSIS2 <- readRDS("results/mean_foi_MSIS2.rds")
mean_foi_MSIS3 <- readRDS("results/mean_foi_MSIS3.rds")
mean_foi_MSIS4 <- readRDS("results/mean_foi_MSIS4.rds")
mean_foi_MSIS5 <- readRDS("results/mean_foi_MSIS5.rds")
mean_foi_MSIRS1 <- readRDS("results/mean_foi_MSIRS1.rds")
mean_foi_MSIRS2 <- readRDS("results/mean_foi_MSIRS2.rds")
mean_foi_MSIRS3 <- readRDS("results/mean_foi_MSIRS3.rds")
mean_foi_comb <- readRDS("results/mean_foi_comb1.rds")
mean_foi_MSIRS5 <- readRDS("results/mean_foi_MSIRS5.rds")

## interpolating the beta for the catalytic foi 
## converting foi to R0 with CIs

foi_to_R0_simple <- function(beta_vector, foi_vector, foi_cat, duration_infection){
  interpol_beta <- approx(x = foi_vector, y = beta_vector,
                          xout = foi_cat, method = "linear")
  R0 <- interpol_beta$y * duration_infection
  return(list(R0 = R0,
              beta = interpol_beta$y))
}

R0_table <- function(beta_vector, foi_vector, duration_infection){
  R0_tab <- data.frame(study = names(foi_target_daily),
                         R0 = foi_to_R0_simple(beta_vector,
                                                foi_vector,
                                                foi_target_daily,
                                                duration_infection)$R0,
                          R0_lower = foi_to_R0_simple(beta_vector,
                                                      foi_vector,
                                                      foi_lower_daily,
                                                      duration_infection)$R0,
                          R0_upper = foi_to_R0_simple(beta_vector,
                                                      foi_vector,
                                                      foi_upper_daily,
                                                      duration_infection)$R0)
  return(R0_tab)
}


# combine CrI and central values into one column for comparison table

R0_col <- function(R0){
  R0[2:4] <- round(R0[2:4], digits = 2)
  output <- paste(R0$R0," (", R0$R0_lower, ", ", R0$R0_upper, ")", sep = "")
  return(output)
}



#### UPDATED

mean_foi_MSIRS_75_90_92 <- readRDS("results/test_mean_foi_MSIRS_75_90_92.rds")
mean_foi_MSIS_100_1_92 <- readRDS("results/test_mean_foi_MSIS_100_1_92.rds")
mean_foi_MSIRS_75_90_10 <- readRDS("results/test_mean_foi_MSIRS_75_90_10.rds")

beta_esti_MSIRS_75_90_92 <- foi_to_R0_simple(beta_list[[1]], mean_foi_MSIRS_75_90_92,
                                   foi_cat = foi_target_daily_ordered, 14)

R0_MSIRS_75_90_92 <- R0_table(beta_vector = beta_list[[1]], 
                              mean_foi_MSIRS_75_90_92,
                              14) # suggests 0.25, 0.5 and 1 is a good 3 to try
R0_MSIS_100_1_92 <- R0_table(beta_vector = beta_list[[2]], 
                              mean_foi_MSIS_100_1_92,
                              14) # suggests 0.25, 0.5 and 1 is a good 3 to try
R0_MSIRS_75_90_10 <- R0_table(beta_vector = beta_list[[7]], 
                              mean_foi_MSIRS_75_90_10,
                              14) # suggests 0.25, 0.5 and 1 is a good 3 to try

# plot beta against FOI for core results and show selection of 0.25, 0.50, 1.00
plot_dat <- data.frame(beta = beta_esti_MSIRS_75_90_92$beta, study = 1:23)

ggplot(plot_dat)+
  geom_point(aes(x = study, y = beta))+
  theme_minimal()+
  geom_hline(yintercept = c(0.25, 0.5, 1.0), col = "red", lty = 2)+
  theme(text = element_text(size = 20))

# create table of results so far for slides 11/01/22

R0_comp <- data.frame(study = R0_MSIRS_75_90_92$study,
                      region = REGION_1,
                      MSIRS_75_90_92 = R0_col(R0_MSIRS_75_90_92),
                      MSIS_100_1_92 = R0_col(R0_MSIS_100_1_92),
                      MSIRS_75_90_10 = R0_col(R0_MSIRS_75_90_10))

R0_comp_order <- R0_comp[order(R0_comp$region),]
write.csv(file = "results/R0/R0_table_core.csv", R0_comp_order)
