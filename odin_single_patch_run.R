# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model <- odin::odin("odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################


# input a value for birth rate (per camel per day) 

alpha <- 0.0005#0.0002926 # default = 

# input a value for the effective contact rate (per camel per day) 

beta <- 0.6 #default = 

# input a value for the average duration of the infectious period (in days) 

duration_infection <- 14 # default = 

gamma <- 1/duration_infection

# input a value for the average duration of complete immunity following infection (in days) 

duration_immunity <- 90 # default = 

sigma <- 1/duration_immunity # default = 

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all

Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input values for the age dependent death rate

mu_6m <- 0.0005 # death rate for 1st 6 months of life
mu_7_12m <- 0.0005 # death rate for 7-12th months of life
mu_2nd_yr <- 0.0005 # death rate for 2nd year of life
mu_3rd_yr <- 0.0005 # death rate for 3rd year of life
mu_4th_yr <- 0.0005 # death rate for 4th year of life
mu_adult_over_4 <- 0.0005 # death rate in adulthood (>4 years)

# input an initial population size
N_0 <- 100000

# input the time period that you wish to run the model for (in days)
time_period <- 36000 
t <- seq(0:time_period)

# set importation rate for introducing infectious individuals
importation_rate <- 0

# if rather than a rate you want importation to occur on a specific day input that day here
imp_t <- 15000  

# set a level of seasonality for births (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 


###############
## run model ##
###############

# include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, Ab_susc = Ab_susc, 
               mAb_susc = mAb_susc, mu_6m = mu_6m, mu_7_12m = mu_7_12m, mu_2nd_yr = mu_2nd_yr,
               mu_3rd_yr = mu_3rd_yr, mu_4th_yr = mu_4th_yr, mu_adult_over_4 = mu_adult_over_4, N_0 = N_0,
               importation_rate = importation_rate, imp_t = imp_t, delta = delta)

out <- as.data.frame(x$run(t))
plot(out$Ntot, type= "l")

out_sum <- out[, c(2, 257, 247:253)]
out_sum_long <- reshape2::melt(out_sum, id.var = "tt")
ggplot(data = out_sum_long) + geom_line(aes(x = tt, y = value, color = variable))+theme_minimal()

par(mfrow = c(5, 6))
par(mar = c(0,3,0,0), las=1)

plot(out$'N[27]'[1300:2700], type= "l", xaxt = 'n')
abline(v = 2000 - 1300, col = "red")

