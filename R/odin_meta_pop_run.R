# MSIRSIR METAPOP 2021

# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model <- odin::odin("odin_meta_pop_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input a value for birth rate (per camel per day) 
alpha <- 0.000565 #0.0006 # 90% female * 50% reproductive age * 45.2% fecundity

# input a value for the effective contact rate (per camel per day) 
beta <- 0.6 #default = 

# input a value for the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 
gamma <- 1/duration_infection

# input a value for the average duration of complete immunity following infection (in days) 
duration_immunity <- 90 # default = 
sigma <- 1/duration_immunity # default = 

# input a value for the average duration of partial immunity following infection (in days) 
duration_mAb_immunity <- 360 * (3/12) # default = 
sigma_m <- 1/duration_mAb_immunity # default = 

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 #1/92 # based on AUC from shedding in Alharbi 

# input values for the age dependent removal rate - balance birthrate

mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

# input an initial population size
N_0 <- 10000

# input the time period that you wish to run the model for (in days)
time_period <- 7352
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day input that day here
imp_t <- 271  

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

# metapop specific
connectivity <- 0.05

###############
## run model ##
###############

# include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma,
               sigma_m = sigma_m,
               Ab_susc = Ab_susc, mAb_susc = mAb_susc, reduced_shed = reduced_shed, 
               mu_1st_yr = mu_1st_yr, mu_2nd_yr = mu_2nd_yr,mu_3rd_yr = mu_3rd_yr, 
               mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, 
               N_0 = N_0, importation_rate = importation_rate, 
               imp_t = imp_t, delta = delta, ind1 = ind1, ind2 = ind2, 
               connectivity = connectivity)

# for a single run
out <- as.data.frame(x$run(t))
#saveRDS(file = "out_metapop.rds", out)
out <- readRDS("out_metapop.rds")
plot(out$Itot, type = "l")
foi <- (beta * mean((out$I_1/out$Ntot)[2000:(2000 + (360*7))])) + (beta * reduced_shed * mean((out$I_2/out$Ntot)[2000:(2000 + (360*7))]))
print(foi)

# subset to include total I per patch
out_I_patch <- out[2000:7352,c(2:27)]
I_patch_long <- gather(data = out_I_patch, "patch", "I", -tt)

p3 <- ggplot(data = I_patch_long)+
  geom_line(aes(x = tt, y = I, col = patch))

cowplot::save_plot(p3, filename = "metapop_ts.png")

# grid coordinate
grid_x <- rep(1:5, times = 5)
grid_y <- rep(1:5, each = 5)
I_patch_long$grid_x <- rep(grid_x, each = dim(out_I_patch)[1])
I_patch_long$grid_y <- rep(grid_y, each = dim(out_I_patch)[1])

a <- ggplot(data = I_patch_long)+
  geom_tile(aes(x = grid_x, y = grid_y, fill = I))+
  transition_time(tt)
animate(a, nframes = 1000)
anim_save(filename = "endemic_circulation.gif")

# looking at inititial spike
out_I_patch <- out[150:500,c(2:27)]
I_patch_long <- gather(data = out_I_patch, "patch", "I", -tt)

ggplot(data = I_patch_long)+
  geom_line(aes(x = tt, y = I, col = patch))

# grid coordinate
grid_x <- rep(1:5, times = 5)
grid_y <- rep(1:5, each = 5)
I_patch_long$grid_x <- rep(grid_x, each = dim(out_I_patch)[1])
I_patch_long$grid_y <- rep(grid_y, each = dim(out_I_patch)[1])

a2 <- ggplot(data = I_patch_long)+
  geom_tile(aes(x = grid_x, y = grid_y, fill = I))+
  transition_time(tt)

animate(a2, nframes = 250)
anim_save(filename = "initial_spike.gif")



