# removing the check on indexes now that all naked "i" have been changed to start at 1 in odin
options(odin.no_check_naked_index = TRUE)
sir_model <- odin::odin("odin_single_patch_model.R", verbose = FALSE, skip_cache = TRUE)
library(ggplot2)

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

# input a value between 0 and 1 for susceptibility experienced by individuals with mAbs and Abs
## 0 would mean mAbs/Abs afford complete protection from MERS (default)
## 1 would mean mAbs/Abs afford no protection at all
Ab_susc <- 1 # default = 
mAb_susc <- 0 # default = 0

# input value for the proportion of baseline naive infectiousness
# seen in reinfected animals
reduced_shed <- 1/92 # based on AUC from Alharbi 
  
# input values for the age dependent removal rate - balance birthrate
out_sum <- out[, c(2, 455, 445:451)]
out_sum_long <- reshape2::melt(out_sum, id.var = "tt")
ggplot(data = out_sum_long) + geom_line(aes(x = tt, y = value, color = variable))+theme_minimal()
# 
# par(mfrow = c(5, 6))
# par(mar = c(0,3,0,0), las=1)
# 
# plot(out$'N[27]'[1300:2700], type= "l", xaxt = 'n')
# abline(v = 2000 - 1300, col = "red")

