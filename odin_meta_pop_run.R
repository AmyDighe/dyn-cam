# simple run script #

sir_model <- odin::odin("odin_meta_pop_model.R", verbose = FALSE, skip_cache = TRUE)


 ## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 3/14

## input a value for recovery rate, gamma

infectious_period <- 14 ## duration of infection in days
gamma <- 1/infectious_period


## input a value for R --> S2, sigma
sigma <- 1/180


## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0.25
mAb_susc <- 0.25

## input values for the age dependent death rate

mu_6m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 1000

## input the time period that you wish to run the model for
time_period <- 4700 
t <- seq(0:time_period)


## input a strength of seasonality of births (1 = max, 0 = no seasonality)

delta = 1

## run model

x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, Ab_susc = Ab_susc, mAb_susc = mAb_susc, 
               mu_6m = mu_6m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0, delta = delta)

out <- as.data.frame(x$run(t))

## run multiple iterations of the model

x_res <- as.data.frame(replicate(100, x$run(0:4700)[,c(3:27, 28:52, 53:77, 4278:4302, 4335, 4337, 4339, 4342)]))

par(mgp = c(3, 0.5, 0), las = 1, mar = c(2, 2, 0.7, 0.5), oma = c(2,2,2,2))
layout(matrix(c(1:25, rep(26,15)), byrow = TRUE, ncol = 5, nrow = 8))

library(scales)
sir_col <- c("aquamarine2", "#cc0044", "#8c8cd9", "black", "pink", "gold", "dimgrey")

for(i in 1:25){
matplot(x_res[,c(seq(from=i,to = i+(104*99), by=104),seq(from=i + 25,to = i+25+(104*99), by=104),
                 seq(from=i +50,to = i+50+(104*99), by=104),seq(from=77,to = i+75+(104*99), by=104))], xlab = "Time (years)", ylab = "",
        type = "l", col = c(rep(alpha(sir_col[1], 0.1), 100),rep(alpha(sir_col[2], 0.1), 100),rep(alpha(sir_col[3], 0.1), 100),rep(alpha(sir_col[4], 0.1), 100)),
        lwd = 3, lty = 1, xaxt = "n", las = 1)
}

matplot(x_res[,c(seq(from = 101, to = 101+(104*99), by = 104), seq(from = 102, to = 102+(104*99), by = 104),
                 seq(from = 103, to = 103+(104*99), by = 104),seq(from = 104, to = 104+(104*99), by = 104))], xlab = "Time (years)", ylab = "",
type = "l", col = c(rep(alpha(sir_col[1], 0.1), 100),rep(alpha(sir_col[2], 0.1), 100),rep(alpha(sir_col[3], 0.1), 100),rep(alpha(sir_col[4], 0.1), 100)),
lwd = 3, lty = 1, xaxt = "n", las = 1)

colMax <- function(data) sapply(data, max, na.rm = TRUE)

infy <- c(1,104*(seq(from = 1, to = 99, by = 1)))
infy_1 <- 26:50
comby <- expand.grid(infy, infy_1)
inf_indexes <- rowSums(comby)
meta_p_infectious <- x_res[,inf_indexes]
meta_p_total_inf <- x_res[,seq(from = 102, to = 102+(104*99), by = 104)]
fades <- colSums(meta_p_total_inf[,]<1) >0 ##the number of runs that dipped to zero ie. faded out
no_starts <- (colMax(meta_p_total_inf[1:360,])) <8 ##the number of runs that did not see an epidemic (never more than 1 infectious individuals in the first year)

perc_of_epis_that_persist <- 100 - (100*(sum(fades) - sum(no_starts))/(100-(sum(no_starts))))

matplot(meta_p_total_inf[1:50,], type = "l", lty=1, lwd=2)

for(i in 1:25){
matplot(x_res[1:10,seq(from=i + 25,to = i+25+(104*99), by=104)], type = "l", lty =1, col = alpha(sir_col[2], 0.2), lwd = 2)
#legend("topleft", legend = c("S", "I", "R", "N"), col = sir_col[1:4], lwd = 2)  
}



