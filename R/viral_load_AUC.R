## Relative infectiousness
library(bayestestR)
library(dplyr)
# there is a big difference between the amount of viral RNA shed by
# naive camels when infected for the first time cf. reinfected camels

# https://www.nature.com/articles/s41598-019-52730-4#Sec20
# data extracted using https://plotdigitizer.com 

vl_data <- read.csv("viral_load.csv")

AUC_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative"))$x,
                                y = (vl_data%>%filter(group == "seronegative"))$y,
                                method = "trapezoid")
AUC_seropos <- area_under_curve(x = (vl_data%>%filter(group == "seropositive"))$x,
                                y = (vl_data%>%filter(group == "seropositive"))$y,
                                method = "trapezoid")

AUC_seroneg/AUC_seroneg # = 92


AUC_vax_seroneg
AUC_vax_seropos
AUC_vax