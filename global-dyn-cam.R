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

# setting up grid for persistance analysis
waning <- 1/c(30, 30*3)
seasonality <- c(0, 0.5, 1)
beta <- c(0.25, 0.5, 1)
pop <- c(c(1,5) %o% 10^(3:6))
import_time <- c(151, 269)

par_grid <- expand.grid(waning = waning, 
                        seasonality = seasonality, 
                        beta = beta, 
                        pop = pop,
                        import_time = import_time)

# for finding point of bifurcation
par_grid_bifurc <- expand.grid(waning = waning[2], 
                        seasonality = seasonality[2:3], 
                        beta = c(0.2, 0.3, 0.35, 0.4, 0.45, 1.5), 
                        pop = pop,
                        import_time = import_time[1])

# for finding the critical community size
par_grid_cc_fine <- expand.grid(waning = waning[2], 
                               seasonality = seasonality, 
                               beta = c(0.25, 0.5, 1.0, 1.5), 
                               pop = c(c(1,2.5,5,7.5) %o% 10^(3:4), c(1,5) %o% 10^(5:6)),
                               import_time = import_time)

# grid for meta-population model
 connectivity <- c(0.001, 0.01, 0.1)

par_grid_metapop <- expand.grid(waning = waning[2],
                        seasonality = seasonality,
                        beta = beta,
                        pop = pop/25,
                        import_time = import_time[1])

beta_list <- list(MSIRS_75_90_92 = seq(0.08, 3.6, by = 0.05),
                  MSIS_100_1_92 = seq(0.08, 3.0, by = 0.05),
                 MSIS_75_1_92 = seq(0.08, 3.2, by = 0.05),
                 MSIS_25_1_92 = seq(0.08, 3.6, by = 0.05),
                 MSIRS_100_90_92 = seq(0.08, 3.5, by = 0.05),
                 MSIRS_100_30_92 = seq(0.08, 3.2, by = 0.05),
                 MSIRS_75_90_10 = seq(0.08, 1.8, by = 0.05),
                 MSIS_75_1_10 = seq(0.08, 1.05, by = 0.05),
                 MSIRS_100_90_10 = seq(0.08, 1.8, by = 0.05))


REGION_1 <- c("South_Asia",rep("Africa", 3), rep("Middle_East", 4), rep("Africa", 3), rep("Middle_East", 5),
              rep("South_Asia", 2), rep("Africa", 3), "Middle_East", "Africa")
