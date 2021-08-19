
# plot persistence analysis
persist <- readRDS("persist_15.rds")
period <- readRDS("period.rds")
ACF_examp <- readRDS("ACF_example.rds")

par_grid$persistence <- persist
par_grid$period <- period
par_grid$beta <- as.character(par_grid$beta)
par_grid$import_time<- as.character(par_grid$import_time)

p1 <- ggplot(data = par_grid%>%filter(waning_mAb == 1/(duration_mAB_cat * 360)), 
       aes(x = log10(pop), y = persistence))+
  geom_line(aes(col = beta, lty = import_time))+
  facet_grid(seasonality~waning, labeller = label_both)+
  theme_bw()

p2 <- ggplot(data = par_grid%>%filter(waning_mAb == 1/(0.5 * 360)), 
       aes(x = log10(pop), y = persistence))+
  geom_line(aes(col = beta, lty = import_time))+
  facet_grid(seasonality~waning, labeller = label_both)+
  theme_bw()

cowplot::save_plot(p1, filename = "persistence_1.png")
cowplot::save_plot(p2, filename = "persistence_2.png")

# table of periodicity

par_grid_endemic <- par_grid[which(!is.na(par_grid$period)),]

p4 <- ggplot(data = par_grid_endemic%>%filter(waning_mAb == 1/(duration_mAB_cat * 360))%>%
         filter(import_time == "271"), aes(x = log10(pop), y = period/360))+
  geom_point(aes(col = beta, group = beta, size = persistence))+
  facet_grid(seasonality~waning, labeller = label_both)+
  theme_bw() +
  geom_abline(slope = 0, intercept = 1, colour = "firebrick", lty = 2)+
  geom_abline(slope = 0, intercept = 2, colour = "firebrick", lty = 2)+
  ylim(0,2.5)

p5 <- ggplot(data = par_grid_endemic%>%filter(waning_mAb == 1/(0.5 * 360))%>%
         filter(import_time == "271"), aes(x = log10(pop), y = period/360))+
  geom_point(aes(col = beta, group = beta, size = persistence))+
  facet_grid(seasonality~waning, labeller = label_both)+
  theme_bw() +
  geom_abline(slope = 0, intercept = 1, colour = "firebrick", lty = 2)+
  geom_abline(slope = 0, intercept = 2, colour = "firebrick", lty = 2)+
  ylim(0,2.5)

cowplot::save_plot(p4, filename = "periodicity.png")
cowplot::save_plot(p5, filename = "periodicity2.png")
