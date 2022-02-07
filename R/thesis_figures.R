
#######################
## Chapter 4 figures ##
#######################

## methods - births
delta_v <- c(0, 0.5, 1)

time_period <- 1800
t <- seq(0:time_period)

birthrate_df <- matrix(nrow = 1801, ncol = 3, byrow = F)
births_df <- matrix(nrow = 1801, ncol = 3, byrow = F)

for(i in 1:3){
  x <- sir_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, 
                     sigma_m = sigma_m, Ab_susc = 0.75, mAb_susc = mAb_susc, 
                     reduced_shed = 1/92, mu = mu, N_0 = 1e+5,
                     importation_rate = importation_rate, 
                     imp_t = imp_t <- 151  + (360 * seq(0, 4, by = 1)), 
                     delta = delta_v[i], ind1 = ind1, ind2 = ind2,
                     foi_bg_usr = 0)
  out <- as.data.frame(x$run(t)[,c(369, 370)])
  
  births_df[,i] <- out$births
  birthrate_df[,i] <- out$birthrate
}

birth_df <- as.data.frame(cbind(births_df, birthrate_df, time = t))
names(birth_df) <- c("d0.births", "d05.births", "d1.births",
                     "d0.birthrate", "d05.birthrate", "d1.birthrate", "time")

birth_long <- birth_df %>% gather(v, value, 1:6) %>%
  separate(v, c("var", "col")) %>%
  arrange(time)%>%
  spread("col", "value")
len <- length(seq(0, 360, by = 30))
len2 <- length(seq(15, 345, by = 30))

ggplot(data = birth_long %>% filter(time>179 & time<541))+
  geom_point(aes(x = time-179, y = births/1e+2, col = var), size = 1.2)+
  geom_line(aes(x = time-179, y = birthrate/1e+2, col = var), size = 1.8)+
  scale_color_discrete(name = "strength of\nseasonality, \U03B4",
                       labels = c("0.0", "0.5", "1.0"),
                       type = c("#C0E9EF","#C8C7F7",  "#D0A4FF"))+
  theme_minimal()+
  scale_x_continuous(name = "time of year", breaks = c(seq(15, 345, by = 30), seq(0, 360, by = 30)),
                   labels = c("July", "Aug", "Sept", "Oct",
                              "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June",
                              rep(c(""), len)))+
  annotate("rect", xmin = 90, xmax = 270, ymin = 0, ymax = 150/1e+2,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 20))+
  ylab("Births per 1000 /day")+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(axis.title.y = element_text(vjust= 2))+
  theme(axis.ticks.x = element_line(color = c(rep(NA, len2), rep("black", len))))+
  theme(axis.line.x.bottom=element_line(color="black"))+
  theme(panel.grid.major.x = element_line(color = c(rep(NA, len2), rep("gray95", len))))+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figs/ch4_methods_births.png")

#### read in results ####
sp_coarse <- readRDS("results/persistence/sp_coarse.rds")
sp_fine_pop <- readRDS("results/persistence/persist_estab_period_cc_fine_sp.rds")
sp_fine_beta <- readRDS("results/persistence/sp_fine_beta.rds")
sp_fine_beta_extrapop <- readRDS("results/persistence/sp_fine_beta_extrapop.rds")
mp_coarse_0.01 <- readRDS("results/persistence/mp2_0.01.rds")


results_coarse_sp <- cbind(sp_coarse, par_grid)
results_coarse_sp <- results_coarse_sp %>% mutate(persist_estab = 100 * persist/estab,
                                                  model = "sp")

results_fine_pop_sp <- cbind(sp_fine_pop, par_grid_fine_pop)
results_fine_pop_sp <- results_fine_pop_sp %>% mutate(persist_estab = 100 * persist/estab,
                                                      model = "sp")

results_fine_beta_sp <- cbind(sp_fine_beta, par_grid_fine_beta)
results_fine_beta_sp <- results_fine_beta_sp %>% mutate(persist_estab = 100 * persist/estab)

results_fine_beta_sp_extrapop <- cbind(sp_fine_beta_extrapop, par_grid_fine_beta_extrapop)
results_fine_beta_sp_extrapop <- results_fine_beta_sp_extrapop %>% mutate(persist_estab = 100 * persist/estab)

results_coarse_mp <- cbind(mp_coarse_0.01, par_grid_metapop)
results_coarse_mp <- results_coarse_mp %>% mutate(persist_estab = 100 * persist/estab,
                                                  model = "mp", pop = 25*pop)

# combine coarse and fine beta for sp periodicity

results_sp <- rbind(results_coarse_sp[,-2], results_fine_beta_sp[, -2], 
                    results_fine_beta_sp_extrapop)


# combine coarse and fine pop for sp persistence

results_sp_persist <- rbind(results_coarse_sp, results_fine_pop_sp)

# combine coarse sp and mp for comparison

results_sp_mp <- rbind(results_coarse_mp, results_coarse_sp[,c(-4, -5, -6)])


## 1. R0 table


## 2. Persistence

sp <- readRDS("results/persistence/persist_estab_period_sp_010222.rds")

results_sp <- cbind(list.cbind(sp$persist), par_grid_core)

seasonality.labs <- c("\U03B4 = 0.5", "\U03B4 = 1.0")
names(seasonality.labs) <- unique(par_grid_core$seasonality)

shedding.labs <- c("1% shedding", "25% shedding", "50% shedding")
names(shedding.labs) <- unique(par_grid_core$shedding)

library(cowplot)
library(gridExtra)
library(lemon)

shed_1 <- ggplot(results_sp %>% filter(shedding == 0.01, pop < 1e+7))+
  geom_line(aes(x = log10(pop), y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = seasonality.labs,
                                                       shedding = shedding.labs),
             scale = "free")+
  scale_x_continuous(limits=c(2,6)) + 
  scale_y_continuous(limits=c(0,100))+
  theme_minimal()+
  scale_color_viridis_d(name = "R0")+
  theme(text = element_text(size = 20))+
  xlab("log10 population size")+
  ylab(" ")+
  theme(axis.title.x = element_blank())+
  theme(axis.line =element_line(color="black"))

shed_25 <- ggplot(results_sp %>% filter(shedding == 0.25, pop < 1e+7))+
  geom_line(aes(x = log10(pop), y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = seasonality.labs,
                                                       shedding = shedding.labs),
             scale = "free")+
  scale_x_continuous(limits=c(2,6)) + 
  scale_y_continuous(limits=c(0,100))+
  theme_minimal()+
  scale_color_viridis_d(name = "R0")+
  theme(text = element_text(size = 20))+
  xlab("log10 population size")+
  ylab("persistence (% of model runs)")+
  theme(axis.title.x = element_blank())+
  theme(strip.text.x = element_blank())+
  theme(axis.line =element_line(color="black"))#+
  #geom_line(data = persist_test, aes(y = log10(pop*25), x = persist, lty = as.factor(connectivity)), col = "red")

shed_50 <- ggplot(results_sp %>% filter(shedding == 0.50, pop < 1e+7))+
  geom_line(aes(x = log10(pop), y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = seasonality.labs,
                                                       shedding = shedding.labs),
             scale = "free")+
  scale_x_continuous(limits=c(2,6)) + 
  scale_y_continuous(limits=c(0,100))+
  theme_minimal()+
  scale_color_viridis_d(name = "R0")+
  theme(text = element_text(size = 20))+
  xlab("log10 population size")+
  ylab(" ")+
  theme(strip.text.x = element_blank())+
  theme(axis.line =element_line(color="black"))

sp_p <- grid.arrange(shed_1, shed_25, shed_50, ncol=1, nrow =3)


ggsave(filename = "figs/persistence_sp.png")

# p <- ggplot(results_sp %>% filter(shedding == 0.25, seasonality == 1, pop < 1e+04))+
#   geom_line(aes(x = log10(pop), y = persist, col = as.factor(beta*14)),
#             size = 1)+
#   theme_minimal()+
#   scale_color_viridis_d(name = "R0")+
#   theme(text = element_text(size = 20))+
#   xlab("log10 total population size")+
#   ylab("persistence (% of model runs)")+
#   theme(strip.text.x = element_blank())+
#   theme(axis.line =element_line(color="black"))+
#   geom_line(data = persist_test, aes(x = log10(pop*25), 
#                                      y = 100*p_e, lty = as.factor(connectivity)), 
#             col = "red", size = 1)
# 
# ggplot(persist_test)+
#   geom_line(aes(x = log10(pop*25), y = persist, col = as.factor(connectivity)))
# ggplot(persist_test)+
#   geom_line(aes(x = log10(pop*25), y = p_e, col = as.factor(connectivity)))

## critical community size table for sp

cc <- results_sp %>%
  arrange(beta, seasonality, pop) %>% 
  select(seasonality, beta, pop, persist_estab)

write.csv(file = "results/cc.csv", cc)



## 3. Periodicity
period_sp <- sp$period
cuts <- apply(period_sp, 2, cut, c(0, 350, 370, 710, 730, 1070, 1090, 1430, 1450, Inf), 
              labels= c("other", "1 yr", "other", "2 yrs", "other", "3 yrs", "other", "4 yrs", "other"))
cuts_long <- as.vector(t(cuts))
cuts_df <- cbind(period = cuts_long, par_grid_core)

cuts_df$pop <- factor(cuts_df$pop,levels = as.character(unique(par_grid_core$pop)))
cuts_df$period <- factor(cuts_df$period, levels = c(NA, "other", "4 yrs",
                                                                  "3 yrs", "2 yrs", "1 yr"))

cuts_df$R0 <- factor(c(rep(c("Low", "Medium", "High", "Very High"), times = 36000/4)),
                     levels = c("Low", "Medium", "High", "Very High"))



q1 <- ggplot(data = cuts_df%>%filter(seasonality == 1, !is.na(period)), 
       aes(x = pop, fill = period))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c("lightgrey", '#fde725', '#35b779', '#31688e', '#472d7b'))+
  theme_minimal()+
  theme(text= element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(shedding~R0, labeller = labeller(shedding = shedding.labs))+
  ylab("proportion of endemic simulations")+
  xlim(as.character(unique((par_grid_core %>% filter(pop>500))$pop)))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")+
  guides(fill = guide_legend(reverse=TRUE))

q2 <- ggplot(data = cuts_df%>%filter(seasonality == 0.5, !is.na(period)), 
             aes(x = pop, fill = period))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c("lightgrey", '#fde725', '#35b779', '#31688e', '#472d7b'))+
  theme_minimal()+
  theme(text= element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(shedding~R0, labeller = labeller(shedding = shedding.labs))+
  ylab("proportion of endemic simulations")+
  xlim(as.character(unique((par_grid_core %>% filter(pop>500))$pop)))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")+
  guides(fill = guide_legend(reverse=TRUE))

ggsave(filename = "results/persistence/periodicity_alternative.png")

# get examples of annual periodicity

out_annual_d1 <- readRDS("results/persistence/time_series_examples_period/out359.rds")
out_annual_d05 <- readRDS("results/persistence/time_series_examples_period/out355.rds")

### zoomed out
time_period <- 12600
out_annual_d1$t <- 0:time_period
out_a_long_d1 <- pivot_longer(out_annual_d1, cols = 1:200,
                              names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
out_annual_d05$t <- 0:time_period
out_a_long_d05 <- pivot_longer(out_annual_d05, cols = 1:200,
                               names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
out_a_long_d1$seasonality <- 1
out_a_long_d05$seasonality <- 0.5

out_a_long <- rbind(out_a_long_d05, out_a_long_d1)
birthrate <- data.frame(t = rep(c(0:time_period), 2), 
                        birthrate = c(out_annual_d1$birthrate.1, 
                                      out_annual_d05$birthrate.1),
                        seasonality = rep(c(1, 0.5), each = length(0:time_period)),
                        var = "jbirths")


ggplot(data = out_a_long %>% filter(t < p_t & t > (m_t + 2450), var == "Itot"))+#, run_no == "V1"))+
  geom_line(aes(x = t, y = values, 
                col = as.factor(seasonality),
                lty = var, group = interaction(run_no,seasonality,var)), alpha = 0.5)+
  geom_line(data = (birthrate %>% filter(t < p_t & t > (m_t + 2450))), 
            aes(x = t, y = 83.5*birthrate, col = as.factor(seasonality), lty = var, 
                group = interaction(seasonality,var)), size = 1)+
  scale_color_manual("strength of\n seasonality", values = c("#B2ABDA", '#472d7b'))+
  scale_linetype_manual("", values = c(1, 2), labels = c("infectious \nindividuals", "births"))+
  theme_classic()+
  scale_y_continuous(
    name = "Infectious individuals",
    sec.axis = sec_axis(~./83.5, name="births")
  )+
  scale_x_continuous(name = "time (years)", breaks = c(seq(p_t - (3*360), p_t, by = 360)),
                     labels = c(1:4))+
  annotate("rect", xmin = c(m_t+2450, p_t - (90+720), p_t - (90+360), p_t - 90), 
           xmax = c(p_t - 990, p_t - 630, p_t - 270, p_t), ymin = 0, ymax = 1000000,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 20))

## zoomed in 

len <- length(seq((p_t - 540), (p_t-180), by = 30))
len2 <- length(seq(p_t - 525, p_t-195, by = 30))

ggplot(data = out_a_long %>% filter(t < p_t & t > (p_t - 360), var == "Itot"))+#, run_no == "V1"))+
  geom_line(aes(x = t, y = values, group = interaction(run_no,seasonality,var), alpha = 0.5,
                col = as.factor(seasonality),
                lty = var), size = 1, alpha = 0.5)+
  geom_line(data = (birthrate %>% filter(t < p_t & t > (p_t - 360))), 
            aes(x = t, y = 83.5*birthrate, col = as.factor(seasonality), lty = var,
                group = interaction(seasonality, var)), size = 1)+
  theme_classic()+
  scale_color_manual("strength of\n seasonality", values = c("#B2ABDA", '#472d7b'))+
  scale_linetype_manual("", values = c(1, 2), labels = c("infectious \nindividuals", "births"))+
  scale_y_continuous(
    name = "Infectious individuals",
    sec.axis = sec_axis(~./83.5, name="births")
  )+
  scale_x_continuous(name = "time of year", breaks = c(seq(p_t - 345, p_t-15, by = 30), 
                                                       seq((p_t - 360), p_t, by = 30)),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "June",
                                "July", "Aug", "Sept", "Oct","Nov", "Dec",
                                rep(c(""), len)))+
  annotate("rect", xmin = p_t - 90, xmax = p_t, ymin = 0, ymax = 1000000,
           alpha = .05,fill = "blue")+
  annotate("rect", xmin = p_t - 360, xmax = p_t - 270, ymin = 0, ymax = 1000000,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 20))+
  ylab("Births per 1000 /day")+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(axis.title.y = element_text(vjust= 2))+
  theme(axis.ticks.x = element_line(color = c(rep(NA, len2), rep("black", len))))+
  theme(axis.line.x.bottom=element_line(color="black"))+
  theme(panel.grid.major.x = element_line(color = c(rep(NA, len2), rep("gray95", len))))+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position = "bottom")



# repeat the above with a smaller pop size (but still annual)

out_annual_d1 <- readRDS("results/persistence/time_series_examples_period/out319.rds")
out_annual_d05 <- readRDS("results/persistence/time_series_examples_period/out315.rds")

period_a <- readRDS("results/persistence/time_series_examples_period/period319.rds")

### zoomed out
time_period <- 12600
out_annual_d1$t <- 0:time_period
out_a_long_d1 <- pivot_longer(out_annual_d1, cols = 1:200,
                              names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
out_annual_d05$t <- 0:time_period
out_a_long_d05 <- pivot_longer(out_annual_d05, cols = 1:200,
                               names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
out_a_long_d1$seasonality <- 1
out_a_long_d05$seasonality <- 0.5

out_a_long <- rbind(out_a_long_d05, out_a_long_d1)
birthrate <- data.frame(t = rep(c(0:time_period), 2), 
                        birthrate = c(out_annual_d1$birthrate.1, 
                                      out_annual_d05$birthrate.1),
                        seasonality = rep(c(1, 0.5), each = length(0:time_period)),
                        var = "jbirths")

mean <- out_a_long %>% filter(t < p_t & t > (m_t + 2450), var == "Itot") %>%
  group_by(t,as.factor(seasonality)) %>% 
  summarise( mean_val = mean(values))

ggplot(data = out_a_long %>% filter(t < p_t & t > (m_t + 2450), var == "Itot"))+#, run_no == "V1"))+
  geom_line(aes(x = t, y = values, 
                col = as.factor(seasonality),
                lty = var, group = interaction(run_no,seasonality,var)), alpha = 0.5)+
  geom_line(data = (birthrate %>% filter(t < p_t & t > (m_t + 2450))), 
            aes(x = t, y = 83.5*birthrate, col = as.factor(seasonality), lty = var, 
                group = interaction(seasonality,var)), size = 1)+
  scale_color_manual("strength of\n seasonality", values = c("#B2ABDA", '#472d7b'))+
  scale_linetype_manual("", values = c(1, 2), labels = c("infectious \nindividuals", "births"))+
  theme_classic()+
  scale_y_continuous(
    name = "Infectious individuals",
    sec.axis = sec_axis(~./83.5, name="births")
  )+
  scale_x_continuous(name = "time (years)", breaks = c(seq(p_t - (3*360), p_t, by = 360)),
                     labels = c(1:4))+
  annotate("rect", xmin = c(m_t+2450, p_t - (90+720), p_t - (90+360), p_t - 90), 
           xmax = c(p_t - 990, p_t - 630, p_t - 270, p_t), ymin = 0, ymax = 2500,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 20))+
  geom_line(data = mean, aes(x = t, y = mean_val, group = mean$`as.factor(seasonality)`), size = 1.5,
            col = c(rep(c("#B2ABDA", "#380474"), each = 1149)), alpha = 1, lty = 1)

## zoomed in 

len <- length(seq((p_t - 540), (p_t-180), by = 30))
len2 <- length(seq(p_t - 525, p_t-195, by = 30))

ggplot(data = out_a_long %>% filter(t < p_t & t > (p_t - 360), var == "Itot"))+#, run_no == "V1"))+
  geom_line(aes(x = t, y = values, group = interaction(run_no,seasonality,var), alpha = 0.5,
                col = as.factor(seasonality),
                lty = var), size = 1, alpha = 0.5)+
  geom_line(data = (birthrate %>% filter(t < p_t & t > (p_t - 360))), 
            aes(x = t, y = 83.5*birthrate, col = as.factor(seasonality), lty = var,
                group = interaction(seasonality, var)), size = 1)+
  theme_classic()+
  scale_color_manual("strength of\n seasonality", values = c("#B2ABDA", '#472d7b'))+
  scale_linetype_manual("", values = c(1, 2), labels = c("infectious \nindividuals", "births"))+
  scale_y_continuous(
    name = "Infectious individuals",
    sec.axis = sec_axis(~./83.5, name="births")
  )+
  scale_x_continuous(name = "time of year", breaks = c(seq(p_t - 345, p_t-15, by = 30), 
                                                       seq((p_t - 360), p_t, by = 30)),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "June",
                                "July", "Aug", "Sept", "Oct","Nov", "Dec",
                                rep(c(""), len)))+
  annotate("rect", xmin = p_t - 90, xmax = p_t, ymin = 0, ymax = 2500,
           alpha = .05,fill = "blue")+
  annotate("rect", xmin = p_t - 360, xmax = p_t - 270, ymin = 0, ymax = 2500,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 20))+
  ylab("Births per 1000 /day")+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(axis.title.y = element_text(vjust= 2))+
  theme(axis.ticks.x = element_line(color = c(rep(NA, len2), rep("black", len))))+
  theme(axis.line.x.bottom=element_line(color="black"))+
  theme(panel.grid.major.x = element_line(color = c(rep(NA, len2), rep("gray95", len))))+
  theme(panel.grid.minor = element_blank())+
  geom_line(data = mean%>% filter(t < (p_t) & t > (p_t - 360)), aes(x = t, y = mean_val, 
                                                                        group = `as.factor(seasonality)`), 
            size = 1.5,
            col = c(rep(c("#B2ABDA", "#380474"), each = 359)), alpha = 1, lty = 1)



## example biennial, triennial

out_biennial <- readRDS("results/persistence/time_series_examples_period/out109.rds")
out_triennial <- readRDS("results/persistence/time_series_examples_period/out349.rds")
out_other <- readRDS("results/persistence/time_series_examples_period/out317.rds")

period_b <- readRDS("results/persistence/time_series_examples_period/period109.rds")
period_t <- readRDS("results/persistence/time_series_examples_period/period349.rds")
period_o <- readRDS("results/persistence/time_series_examples_period/period317.rds")

select_b <- which(period_b > 700 & period_b < 740)
out_biennial <- out_biennial[, grep("Itot", colnames(out_biennial))]
out_biennial <- out_biennial[, select_b]
out_biennial$t <- 0:time_period
out_b_long <- pivot_longer(out_biennial, 1:(length(select_b)), 
                           names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
###################
## every 2 years ##
###################

ggplot()+
  geom_line(data = out_b_long %>% filter(t>6840 & t<(6840+(10*360))), 
            aes(x = (t - 6840)/360, y =values, group = run_no), col = '#31688e', alpha = 0.1)+
  geom_line(data = out_b_long %>% filter(t>6840 & t<(6840+(10*360)), run_no == "4"),
            aes(x = (t - 6840)/360, y =values, group = run_no), size = 2, col = '#31688e')+
  xlab("time (years)")+
  ylab("infectious animals")+
  theme_classic()+
  theme(text = element_text(size = 20))

bi_acf <- acf(x = (out_b_long %>% filter(t>6840 & t<(6840+(10*360)), run_no == "4"))$values,
    lag.max = 5*360, plot = T)

###################
## every 3 years ##
###################

select_t <- which(period_t > 1070 & period_t < 1090)
out_triennial <- out_triennial[, grep("Itot", colnames(out_triennial))]
out_triennial <- out_triennial[, select_t]
out_triennial$t <- 0:time_period
out_t_long <- pivot_longer(out_triennial, 1:(length(select_t)), 
                           names_to = c("var", "run_no"), names_sep = "\\.", values_to = "values")
ggplot()+
  geom_line(data = out_t_long %>% filter(t>m_t & t<p_t), 
            aes(x = (t - m_t)/360, y =values, group = run_no), col = '#35b779', alpha = 0.2)+
  geom_line(data = out_t_long %>% filter(t>m_t & t<p_t, run_no == "36"),
            aes(x = (t - m_t)/360, y =values, group = run_no), size = 2, col = '#35b779')+
  xlab("time (years)")+
  ylab("infectious animals")+
  theme_classic()+
  theme(text = element_text(size = 20))

tri_acf <- acf(x = (out_t_long %>% filter(t>m_t & t< p_t, run_no == "36"))$values,
    lag.max = 5*360, plot = T)
  # maybe 27, 35, 36

###################
## every 4 years ##
###################

ggplot()+
  geom_line(data = out_t_long %>% filter(t>6840 & t<(6840+(10*360))), 
            aes(x = (t - 6840)/360, y =values, group = run_no), col = '#fde725', alpha = 0.2)+
  geom_line(data = out_t_long %>% filter(t>6840 & t<(6840+(10*360)), run_no == "83"),
            aes(x = (t - 6840)/360, y =values, group = run_no), size = 2, col = '#fde725')+
  xlab("time (years)")+
  ylab("infectious animals")+
  theme_classic()+
  theme(text = element_text(size = 20))

qua_acf <- acf(x = (out_t_long %>% filter(t>6840 & t<(6840+(10*360)), run_no == "83"))$values,
    lag.max = 5*360, plot = T)

#####################
## acf comparisons ##
#####################

ACF <- data.frame(lag = bi_acf$lag, biennial = bi_acf$acf, triennial = tri_acf$acf, 
                  quadrennial = qua_acf$acf)
ACF_long <- pivot_longer(ACF, 2:4, names_to = "period", values_to = "acf")
ACF_long$period <- factor(ACF_long$period, levels= c("biennial", "triennial", "quadrennial"))


ggplot(ACF_long)+
  geom_segment(aes(x=lag/360,xend=lag/360,y=0,yend=acf, col = period))+
  scale_colour_manual(values = c('#fde725', '#35b779', '#31688e'))+
  theme_minimal()+
  ylab("acf")+
  xlab("lag (years)")+
  geom_hline(yintercept = c(-0.05, +0.05), lty = 2, col = "red", size = 1)+
  theme(text = element_text(size = 20))


### what's going on with "other"
out_other$t <- 0:time_period
out_o_long <- pivot_longer(out_other, 1:100, 
                           names_to = "run_no", values_to = "Itot")
ggplot(data = out_o_long %>% filter((t < p_t & t > m_t)))+#, run_no == "V10"))+
  geom_line(aes(x = t, y = Itot, group = run_no))+
  theme_classic()










#######################
## Chapter 5 figures ##
#######################

## 1. age optimisation of targetted vaccination a and b


## 2. vaccine tpp


## 3. vaccine impact on persistence
