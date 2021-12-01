# dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)
library(tictoc) # for benchmarking
# parameters

foi_fits <- readRDS("data/fit_4bbreal.rds")
foi_studies <- rownames(readRDS("data/SEROPOS.rds")) # from gen_real in cat-cam project
foi_fits_sum <- as.data.frame((rstan::summary(foi_fits))$summary)
foi_target_daily <- foi_fits_sum[grep("foi", rownames(foi_fits_sum)), "mean" ] /360
foi_upper_daily <- foi_fits_sum[grep("foi", rownames(foi_fits_sum)), "97.5%" ] /360
foi_lower_daily <- foi_fits_sum[grep("foi", rownames(foi_fits_sum)), "2.5%" ] /360
names(foi_target_daily) <- foi_studies # name with study & country
names(foi_upper_daily) <- foi_studies # name with study & country
names(foi_lower_daily) <- foi_studies # name with study & country
foi_target_daily_ordered <- sort(foi_target_daily)# order for ease with tuning
foi_upper_daily_ordered <- sort(foi_upper_daily)# order for ease with tuning
foi_lower_daily_ordered <- sort(foi_lower_daily)# order for ease with tuning

## need beta values to map to foi between 0.00011 and 0.0259

# from catalytic model 4bb
duration_Ab_cat <- 1/0.09
duration_mAB_cat <- 1/4.42

# for grid for MSISI
waning <- 1/c(duration_Ab_cat*360, 5 *360, 2*360)

# setting up grid for MSIRSIR single patch model
waning <- 1/c(30, 30*3, 30*6)
seasonality <- c(0, 0.5, 1)
beta <- c(0.25, 0.5, 1.4) # corresponding to foi of 0.35, 1, 3.5 under baseline assumptions
pop <- c(c(1,5) %o% 10^(3:7))
waning_mAb <- 1/c(duration_mAB_cat * 360, 360 * 0.5)
import_time <- c(89, 179, 269, 359)

par_grid <- expand.grid(waning = waning, 
                        seasonality = seasonality, 
                        beta = beta, 
                        pop = pop, 
                        waning_mAb = waning_mAb, 
                        import_time = import_time)

# grid for meta-population model
waning_meta <- 1/90
seasonality_meta <- c(0, 1)
waning_mAb_meta <- 1/(duration_mAB_cat * 360)
import_time_meta <- 271
connectivity <- c(0.01, 0.1, 1)

par_grid_metapop <- expand.grid(waning = waning_meta, 
                        seasonality = seasonality_meta, 
                        beta = beta, 
                        pop = pop, 
                        waning_mAb = waning_mAb_meta, 
                        import_time = import_time_meta,
                        connectivity = connectivity)

# beta vectors - make sure they span foi's from the cat model
beta_vector_MSIS1 <- seq(0.08, 0.45, by = 0.01) # MSIS 1 day comp, 100% susc, 1/10 redshed
beta_vector_MSIS2 <- seq(0.08, 0.65, by = 0.01) # MSIS 1 day comp, 100% susc, 1/92 redshed
beta_vector_MSIS3 <- seq(0.08, 0.7, by = 0.01) # MSIS 1 day comp, 75% susc, 1/92 redshed
beta_vector_MSIS4 <- seq(0.08, 1.1, by = 0.01) # MSIS 1 day comp, 25% susc, 1/92 redshed
beta_vector_MSIRS1 <- seq(0.08, 0.75, by = 0.01) # MSIRS 90 day comp, 100% susc, 1/92 redshed
beta_vector_MSIRS2 <- seq(0.08, 0.7, by = 0.01) # MSIRS 60 day comp, 100% susc, 1/92 redshed
beta_vector_MSIRS3 <- seq(0.08, 0.7, by = 0.01) # MSIRS 30 day comp, 100% susc, 1/92 redshed
beta_vector_comb1 <- seq(0.08, 0.75, by = 0.01) # comb 60 day comp, 75% susc, 1/92 redshed
beta_vector_MSIS5 <- seq(0.08, 0.55, by = 0.01) # MSIS5 1 day comp, 75% susc, 1/10 redshed
beta_vector_MSIRS5 <- seq(0.08, 0.6, by = 0.01) # MSIS5 90 day comp, 100% susc, 1/10 redshed
