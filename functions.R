
###################
### R0 analyses ###
###################

## interpolating the beta for the catalytic foi 
## converting foi to R0 with CIs

foi_to_R0_simple <- function(beta_vector, foi_vector, foi_cat, duration_infection){
  interpol_beta <- approx(x = foi_vector, y = beta_vector,
                          xout = foi_cat, method = "linear",
                          rule = 2)
  R0 <- interpol_beta$y * duration_infection
  return(list(R0 = R0,
              beta = interpol_beta$y))
}

R0_table <- function(beta_vector, foi_vector, foi_df, duration_infection){
  R0_tab <- data.frame(study = foi_df$study,
                       R0 = foi_to_R0_simple(beta_vector,
                                             foi_vector,
                                             foi_df$mean,
                                             duration_infection)$R0,
                       R0_lower = foi_to_R0_simple(beta_vector,
                                                   foi_vector,
                                                   foi_df$lower,
                                                   duration_infection)$R0,
                       R0_upper = foi_to_R0_simple(beta_vector,
                                                   foi_vector,
                                                   foi_df$upper,
                                                   duration_infection)$R0)
  return(R0_tab)
}


# combine CrI and central values into one column for comparison table

R0_col <- function(R0){
  R0[2:4] <- round(R0[2:4], digits = 2)
  output <- paste(R0$R0," (", R0$R0_lower, ", ", R0$R0_upper, ")", sep = "")
  return(output)
}




# extrapolation example

