# dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)
# parameters

#foi_fits <- readRDS("data/fit_4bbreal.rds")

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

