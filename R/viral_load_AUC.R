## Relative infectiousness
library(bayestestR)
library(dplyr)
library(ggplot2)
# there is a big difference between the amount of viral RNA shed by
# naive camels when infected for the first time cf. reinfected camels

# https://www.nature.com/articles/s41598-019-52730-4#Sec20
# data extracted using https://plotdigitizer.com 

vl_data <- read.csv("data/viral_load.csv")

ggplot(vl_data)+
  geom_line(aes(x = x, y = y, col = group))+
  theme_bw()

AUC_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative"))$x,
                                y = (vl_data%>%filter(group == "seronegative"))$y,
                                method = "trapezoid")
AUC_seropos <- area_under_curve(x = (vl_data%>%filter(group == "seropositive"))$x,
                                y = (vl_data%>%filter(group == "seropositive"))$y,
                                method = "trapezoid")

AUC_seropos/AUC_seroneg # 1/92


AUC_vax_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative_vax"))$x,
                                 y = (vl_data%>%filter(group == "seronegative_vax"))$y,
                                 method = "trapezoid")
AUC_vax_seropos<- area_under_curve(x = (vl_data%>%filter(group == "seropositive_vax"))$x,
                                   y = (vl_data%>%filter(group == "seropositive_vax"))$y,
                                   method = "trapezoid")

AUC_vax_seropos/AUC_seropos # 0.1267937
AUC_vax_seroneg/AUC_seroneg
