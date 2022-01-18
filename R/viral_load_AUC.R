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

AUC_vax_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative_vax"))$x,
                                 y = (vl_data%>%filter(group == "seronegative_vax"))$y,
                                 method = "trapezoid")
AUC_vax_seropos<- area_under_curve(x = (vl_data%>%filter(group == "seropositive_vax"))$x,
                                   y = (vl_data%>%filter(group == "seropositive_vax"))$y,
                                   method = "trapezoid")

redshed_inf <- AUC_seropos/AUC_seroneg # seropositivity reduces shedding by 92x 
v_redshed_if_seropos <- AUC_vax_seropos/AUC_seropos # vax seropositive animals reduces their shedding by around 8x
redshed_v_pos_cf_v_neg <- AUC_vax_seropos/AUC_vax_seroneg # vax reduces shedding 1000x more in animals that are already seropositive cf seronegative animals
v_redshed_if_seroneg <- AUC_vax_seroneg/AUC_seroneg #  vaccine does not reduce shedding in seronegative animals
redshed_if_vax_and_inf <- AUC_vax_seropos/AUC_seroneg # vaccination and seropositivity reduces shedding by ~724x compared to unvaccinated naive animals
redshed_vax_neg_cf_unvax_pos <- AUC_vax_seroneg/AUC_seropos #infection alone is 126x better at reducing shed than vaccination alone

## what about log AUC values?

logAUC_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative"))$x,
                                y = log((vl_data%>%filter(group == "seronegative"))$y),
                                method = "trapezoid")
logAUC_seropos <- area_under_curve(x = (vl_data%>%filter(group == "seropositive"))$x,
                                y = log((vl_data%>%filter(group == "seropositive"))$y),
                                method = "trapezoid")

logAUC_vax_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative_vax"))$x,
                                    y = log((vl_data%>%filter(group == "seronegative_vax"))$y),
                                    method = "trapezoid")
logAUC_vax_seropos<- area_under_curve(x = (vl_data%>%filter(group == "seropositive_vax"))$x,
                                   y = log((vl_data%>%filter(group == "seropositive_vax"))$y),
                                   method = "trapezoid")
logredshed_inf <- logAUC_seropos/logAUC_seroneg # seropositivity reduces shedding by 2x 
logv_redshed_if_seropos <- logAUC_vax_seropos/logAUC_seropos # vax seropositive animals reduces their shedding by around 1.5x
logredshed_v_pos_cf_v_neg <- logAUC_vax_seropos/logAUC_vax_seroneg # vax reduces shedding 3x more in animals that are already seropositive cf seronegative animals
logv_redshed_if_seroneg <- logAUC_vax_seroneg/logAUC_seroneg #  vaccine does not reduce shedding in seronegative animals
logredshed_if_vax_and_inf <- logAUC_vax_seropos/logAUC_seroneg # vaccination and seropositivity reduces shedding by ~3x compared to unvaccinated naive animals
logredshed_vax_neg_cf_unvax_pos <- logAUC_vax_seroneg/logAUC_seropos #infection alone is 2x better at reducing shed than vaccination alone
