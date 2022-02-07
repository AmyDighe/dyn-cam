## estimating R0 from foi data 

beta_esti <- matrix(nrow = 23, ncol = 18)
beta_esti_lower <- beta_esti
beta_esti_upper <- beta_esti
mean_foi <- vector(mode = "list", length = 18)
core <- c(2,3,5,8,9,11,14,15,17)

for(i in 1:18){
  mean_foi[[i]] <- readRDS(file = paste("results/mean_foi_", scenario[i], ".rds", sep = ""))
}

for(j in 1:18){
  
  beta_esti[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                     foi_vector = mean_foi[[j]],
                                     foi_cat = foi_df$mean, 14)$beta
  beta_esti_lower[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                     foi_vector = mean_foi[[j]],
                                     foi_cat = foi_df$lower, 14)$beta
  beta_esti_upper[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                     foi_vector = mean_foi[[j]],
                                     foi_cat = foi_df$upper, 14)$beta
}

## beta vs foi graphs

par(mfrow = c(3, 6))

for(i in 1:18){
  plot(x = beta_esti[,i], y = foi_df$mean, type = "p")
}

for(i in 1:18){
  plot(x = beta_esti_lower[,i], y = foi_df$lower, type = "p")
}

for(i in 1:18){
  plot(x = beta_esti_upper[,i], y = foi_df$upper, type = "p")
}

dat <- data.frame(foi_df, R0_1 = 14 * beta_esti[,2],
                  R0_25 = 14 * beta_esti[,14],
                  R0_50 = 14 * beta_esti[,8])

cuts <- data.frame(vals_1 = c(3.5, 7, 14, 21),
                   vals_25 = c(2.5, 4, 5.5, 7.5), 
                   vals_50 = c(1.75, 2.5, 3.5, 4.5),
                   y = rep(0.01, 4))
p1 <- ggplot(dat)+
  geom_point(aes(x = R0_1, y = mean, col = region), size = 3)+
  theme_minimal()+
  ylab("FOI")+
  xlab("R0 when shedding = 1%")+
  scale_color_viridis_d()+
  geom_vline(xintercept = cuts$vals_1, lty = 2, col = "firebrick")+
  geom_text(data = cuts, aes(x = vals_1 + 0.5, y = y, label = vals_1), 
            col = "firebrick")

p2 <- ggplot(dat)+
  geom_point(aes(x = R0_25, y = mean, col = region), size = 3)+
  theme_minimal()+
  ylab("FOI")+
  xlab("R0 when shedding = 25%")+
  scale_color_viridis_d()+
  geom_vline(xintercept = cuts$vals_25, lty = 2, col = "firebrick")+
  geom_text(data = cuts, aes(x = vals_25 + 0.1, y = y, label = vals_25), 
            col = "firebrick")

p3 <- ggplot(dat)+
  geom_point(aes(x = R0_50, y = mean, col = region), size = 3)+
  theme_minimal()+
  ylab("FOI")+
  xlab("R0 when shedding = 50%")+
  scale_color_viridis_d()+
  geom_vline(xintercept = cuts$vals_50, lty = 2, col = "firebrick")+
  geom_text(data = cuts, aes(x = vals_50 + 0.1, y = y, label = vals_50), 
            col = "firebrick")

prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)
prow

legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p_all <- plot_grid(prow, legend, rel_widths = c(3, .4))

ggsave(filename= "figs/foi_R0.png", p_all)


R0_tab <- vector(mode = "list", length = 18)
names(R0_tab) <- scenario
for(j in 1:18){
 R0_tab[[j]] <- R0_table(beta_vector = beta_list[[j]],
                         foi_vector = mean_foi[[j]],
                         foi_df = foi_df,
                         duration_infection = 14) 
}

R0_full <- bind_rows(R0_tab, .id = "column_label")
R0_appendix <- R0_full %>%
  mutate(R0_comp = paste(round(R0,1), " (", round(R0_lower,1), ", ",
                         round(R0_upper, 1), ")", sep = ""))
R0_appendix_wide <- R0_appendix %>% 
  dplyr::select(-R0, -R0_lower, -R0_upper)%>%
  pivot_wider(names_from = column_label, values_from = R0_comp)%>%
  mutate(region = REGION_1)

R0_appendix_order <- R0_appendix_wide[order(R0_appendix_wide$region),]
write.csv(file = "results/R0/R0_appendix.csv", R0_appendix_order)

# create table of all Re results

Se <- readRDS("results/Se.rds")
Se_l <- readRDS("results/Se_lower.rds")
Se_u <- readRDS("results/Se_upper.rds")

Re <- Se * R0_esti[,core]
Re_lower <- Se_l * R0_esti_lower[,core]
Re_upper <- Se_u * R0_esti_upper[,core]

Re_tab <- vector(mode = "list", length = 9)
names(Re_tab) <- scenario_core

for(j in 1:9){
 Re_tab[[j]] <-  data.frame(study = foi_df$study,
                            Re = R0_tab[[core[j]]]$R0 * Se[,j],
                            Re_lower = R0_tab[[core[j]]]$R0_lower * Se_l[,j],
                            Re_upper = R0_tab[[core[j]]]$R0_upper * Se_u[,j])
}

Re_full <- bind_rows(Re_tab, .id = "column_label")

Re_appendix <- Re_full %>%
  mutate(Re_comp = paste(round(Re,1), " (", round(Re_lower,1), ", ",
                         round(Re_upper, 1), ")", sep = ""))
Re_appendix_wide <- Re_appendix %>% 
  dplyr::select(-Re, -Re_lower, -Re_upper)%>%
  pivot_wider(names_from = column_label, values_from = Re_comp)%>%
  mutate(region = REGION_1)

Re_appendix_order <- Re_appendix_wide[order(Re_appendix_wide$region),]
write.csv(file = "results/R0/Re_appendix.csv", Re_appendix_order)


# create graphs for core R0 and Re results
Re_long <- pivot_longer(Re_full, 3:5, names_to = "quantity", 
                        values_to = "estimate")
R0_long <- pivot_longer(R0_full, 3:5, names_to = "quantity", 
                        values_to = "estimate")
Re_R0_long <- rbind(Re_long, R0_long)

saveRDS(file = "results/R0/Re_R0.rds", Re_R0_long)



