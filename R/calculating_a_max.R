
# estimating a_max for the equilibrium solution

x = seq(1:36000)
plot(100*exp(-0.0004*x)/(exp(0)), type = "l", xlim = c(0,17500))
abline(h=5) #marks the point where only 5% of the population remains
abline(h=1) # marks the point where only 1% of the population remains
abline(h=50) #marks the point where 50% of the population remains
abline(v=7200)
abline(v=11520)
abline(v=1800) 
# after 32 years (11520 days) only 1% of the population remain - lets take 32 as the maximum survival in age class 27
# camels enter age class 27 when they are ~5 years old implying the maximum camel age in this model is 37

#after 5 years 50% of the population remains, suggesting the mean life expectancy is ~10 years