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
                                                duration_infection),
                          R0_lower = foi_to_R0_simple(beta_vector,
                                                      foi_vector,
                                                      foi_lower_daily,
                                                      duration_infection),
                          R0_upper = foi_to_R0_simple(beta_vector,
                                                      foi_vector,
                                                      foi_upper_daily,
                                                      duration_infection))
  return(R0_tab)
}

mean_foi_MSIRS_75_90_92 <- readRDS("results/mean_foi_MSIRS_75_90_92.rds")

beta_esti$beta <- foi_to_R0_simple(beta_list[[1]], mean_foi_MSIRS_75_90_92,
                              foi_cat = foi_target_daily_ordered, 14)

R0_MSIRS_75_90_92 <- R0_table(beta_vector = beta_list[[1]], 
                              mean_foi_MSIRS_75_90_92,
                              14) # suggests 0.25, 0.5 and 1 is a good 3 to try



R0_MSIS1 <- R0_table(beta_vector_MSIS1, mean_foi_MSIS1, 14)
R0_MSIS2 <- R0_table(beta_vector_MSIS2, mean_foi_MSIS2, 14)
R0_MSIS3 <- R0_table(beta_vector_MSIS3, mean_foi_MSIS3, 14)
R0_MSIS4 <- R0_table(beta_vector_MSIS4, mean_foi_MSIS4, 14)
R0_MSIS5 <- R0_table(beta_vector_MSIS5, mean_foi_MSIS5, 14)
R0_MSIRS1 <- R0_table(beta_vector_MSIRS1, mean_foi_MSIRS1, 14)
R0_MSIRS2 <- R0_table(beta_vector_MSIRS2, mean_foi_MSIRS2, 14)
R0_MSIRS3 <- R0_table(beta_vector_MSIRS3, mean_foi_MSIRS3, 14)
R0_comb <- R0_table(beta_vector_comb1, mean_foi_comb, 14)
R0_MSIRS5 <- R0_table(beta_vector_MSIRS5, mean_foi_MSIRS5, 14)

saveRDS(file = "results/R0/R0_MSIS1.rds", R0_MSIS1)
saveRDS(file = "results/R0/R0_MSIS2.rds", R0_MSIS2)
saveRDS(file = "results/R0/R0_MSIS3.rds", R0_MSIS3)
saveRDS(file = "results/R0/R0_MSIS4.rds", R0_MSIS4)
saveRDS(file = "results/R0/R0_MSIS5.rds", R0_MSIS5)
saveRDS(file = "results/R0/R0_MSIRS1.rds", R0_MSIRS1)
saveRDS(file = "results/R0/R0_MSIRS2.rds", R0_MSIRS2)
saveRDS(file = "results/R0/R0_MSIRS3.rds", R0_MSIRS3)
saveRDS(file = "results/R0/R0_MSIRS5.rds", R0_MSIRS5)
saveRDS(file = "results/R0/R0_comb.rds", R0_comb)

# combine CrI and central values into one column for comparison table

R0_col <- function(R0){
  R0[2:4] <- round(R0[2:4], digits = 2)
  output <- paste(R0$R0," (", R0$R0_lower, ", ", R0$R0_upper, ")", sep = "")
  return(output)
}

R0_comp <- data.frame(study = R0_MSIRS1$study,
                      region = REGION_1,
                      MSIS_1 = R0_col(R0_MSIS2),
                      MSIS_075 = R0_col(R0_MSIS3),
                      MSIS_025 = R0_col(R0_MSIS4),
                      MSIRS_90 = R0_col(R0_MSIRS1),
                      MSIRS_60 = R0_col(R0_MSIRS2),
                      MSIRS_30 = R0_col(R0_MSIRS3),
                      comb_075_60 = R0_col(R0_comb),
                      MSIS_075_01 = R0_col(R0_MSIS5),
                      MSIS_1_01 = R0_col(R0_MSIS1),
                      MSIRS_90_01 = R0_col(R0_MSIRS5))

R0_comp_order <- R0_comp[order(R0_comp$region),]
write.csv(file = "results/R0/R0_table.csv", R0_comp_order)
