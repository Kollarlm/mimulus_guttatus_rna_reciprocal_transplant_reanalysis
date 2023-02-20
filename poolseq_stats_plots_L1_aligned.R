# Plotting output pf python scripts from Gould et al 2017 
# Leslie M Kollar

## Loading in libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggtext)

## Loading in data files

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned"
setwd(path)

# Files for pi
files=list.files(pattern = ".pi")
L1_aligned_dataset_pi = do.call(rbind, lapply(files,fread))


# Files for gstat
files=list.files(pattern = ".gstat")
L1_aligned_dataset_gstat = do.call(rbind, lapply(files,fread))

## Removed windows with average depth greater than 2 standard deviations

L1_SD <- L1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 179.8778
L1_filt <- (2*L1_SD)

L1_mean <- L1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=mean(win_depth, na.rm=TRUE)) # mean = 166.059

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

# gstat
L1_SD_g <- L1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=sd(avg.win.depth, na.rm=TRUE)) # standard deviation = 179.8778
L1_filt_g <- (2*L1_SD_g)

L1_mean_g <- L1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=mean(avg.win.depth, na.rm=TRUE)) # mean = 

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat %>% # removes data greater than 2 SD from mean
  filter(abs(scale(avg.win.depth)) < 2)

## Adding in column for the ratio of pi

L1_aligned_dataset_pi_filtered$piRatio <- (L1_aligned_dataset_pi_filtered$freqIA/L1_aligned_dataset_pi_filtered$freqCP)

## Take top and bottom 1% of piRatio
obs_L1 <-  nrow(L1_aligned_dataset_pi_filtered)

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% arrange(desc(piRatio))
L1_pi_top1 <- L1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% arrange((piRatio))
L1_pi_bottom1 <- L1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)

## Take top 1% of gstat

obs_L1_gstat <-  nrow(L1_aligned_dataset_gstat_filtered)

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat_filtered %>% arrange(desc(G_stat))
L1_pi_top1_gstat <- L1_aligned_dataset_gstat_filtered %>% 
  filter(row_number() < obs_L1_gstat * 0.01)

## Preparing the data to plot manhattan plot in ggplot

# Remvoing Chr so that the chromsomes are in order

L1_aligned_dataset_pi_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_pi_filtered$scaff))
L1_aligned_dataset_pi_filtered$scaff <- as.numeric(L1_aligned_dataset_pi_filtered$scaff)

L1_aligned_dataset_gstat_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_gstat_filtered$scaff))
L1_aligned_dataset_gstat_filtered$scaff <- as.numeric(L1_aligned_dataset_gstat_filtered$scaff)

data_cum_L1_pi <- L1_aligned_dataset_pi_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% 
  inner_join(data_cum_L1_pi, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)
  
data_cum_L1_gstat <- L1_aligned_dataset_gstat_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat_filtered %>% 
  inner_join(data_cum_L1_gstat, by ="scaff") %>% 
  mutate(bp_cum=win_end + bp_add)

## Scatter plot of piRatio

L1_pi_plot <- ggplot(L1_aligned_dataset_pi_filtered, aes(x=bp_cum, y= log10(piRatio))) +
  geom_point() +
  geom_hline(yintercept =log10(1.716585), color="red") +
  geom_hline(yintercept = log10(0.5020223), color="red") +
  annotate(geom = "rect", xmin=(119037830+850429), xmax=(119037830+7604769), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(119037830+1032334), xmax=(119037830+1246126), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67547331+13650670), xmax=(67547331+17847181), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(232897274+5329939), xmax=(232897274+7791197), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

  

# Scatter plot of read depth
L1_depth_plot <- ggplot(L1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = avg.win.depth)) +
  stat_smooth( method = "loess", span = 0.1, se = FALSE)+
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(232813426+5329939), xmax=(232813426+7791197), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)
  

## Scatter plot of gstat
L1_gstat_plot <- ggplot(L1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="red") +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(232813426+5329939), xmax=(232813426+7791197), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

## Arrange plots on one graph (Figure2)

L1_alinged_Fig2 <- ggarrange(L1_gstat_plot, L1_pi_plot,L1_depth_plot,nrow = 3, ncol = 1 )

tiff("./figures/L1_aligned_gouldFig2.png", width = 30, height = 20, units = "cm",res = 300)
print(L1_alinged_Fig2)
dev.off()

## Figure 3 from Gould et al 2017 (G stat values across chromosmes 5, 8, and 14)


L1_gstat_plot_8 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="red") +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)

L1_gstat_plot_5 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="red") +
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)

L1_gstat_plot_14 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="red") +
  annotate(geom = "rect", xmin=(232813426+5329939), xmax=(232813426+7791197), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

## Arrange plots on one graph (Figure3)

L1_alinged_Fig3 <- ggarrange(L1_gstat_plot_5, L1_gstat_plot_8, L1_gstat_plot_14,nrow = 3, ncol = 1 )

tiff("./figures/L1_aligned_gouldFig3.png", width = 30, height = 20, units = "cm",res = 300)
print(L1_alinged_Fig3)
dev.off()

## Distribution of G Statistic window values inside and outside the chromosome inversion (figure 4)

# Add column for within inversion or not within
L1_aligned_dataset_gstat_filtered_5_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum < 67527869+13650670) 
  L1_aligned_dataset_gstat_filtered_5_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_5_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum > 67527869+17847181)
L1_aligned_dataset_gstat_filtered_5_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_5_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_5_NOTinv2, L1_aligned_dataset_gstat_filtered_5_NOTinv1)
  
L1_aligned_dataset_gstat_filtered_5_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(between(bp_cum, (67527869+13650670),(67527869+17847181) ))
L1_aligned_dataset_gstat_filtered_5_inv$karyo <- "inv"
 
L1_aligned_dataset_gstat_filtered_5_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_5_inv, L1_aligned_dataset_gstat_filtered_5_NotINV_fig)

Gstat_chr5_density <- L1_aligned_dataset_gstat_filtered_5_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 5 Inversion")

################################################################################

L1_aligned_dataset_gstat_filtered_8_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+850429)) 

L1_aligned_dataset_gstat_filtered_8_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+7604769))
L1_aligned_dataset_gstat_filtered_8_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_8_NOTinv1, L1_aligned_dataset_gstat_filtered_8_NOTinv2)

L1_aligned_dataset_gstat_filtered_8_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+850429),(118986430+7604769)))
L1_aligned_dataset_gstat_filtered_8_inv$karyo <- "inv"

L1_aligned_dataset_gstat_filtered_8_inv_wosmall1 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))

L1_aligned_dataset_gstat_filtered_8_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_8_inv_wosmall1, L1_aligned_dataset_gstat_filtered_8_inv_wosmall2, L1_aligned_dataset_gstat_filtered_8_NotINV_fig)

Gstat_chr8_density <- L1_aligned_dataset_gstat_filtered_8_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 8 Inversion")

################################################################################

L1_aligned_dataset_gstat_filtered_8small_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8small_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8small_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))
L1_aligned_dataset_gstat_filtered_8small_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8small_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_8_NOTinv1, L1_aligned_dataset_gstat_filtered_8_NOTinv2)

L1_aligned_dataset_gstat_filtered_8small_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+1032334), (118986430+1246126)))
L1_aligned_dataset_gstat_filtered_8small_inv$karyo <- "inv"

L1_aligned_dataset_gstat_filtered_8small_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_8small_inv, L1_aligned_dataset_gstat_filtered_8small_NotINV_fig)

Gstat_chr8small_density <- L1_aligned_dataset_gstat_filtered_8small_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 8small Inversion")
############################################################################################################################

L1_aligned_dataset_gstat_filtered_14_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum < (232813426+5329939)) 

L1_aligned_dataset_gstat_filtered_14_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_14_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum > (232813426+7791197))
L1_aligned_dataset_gstat_filtered_14_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_14_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_14_NOTinv2, L1_aligned_dataset_gstat_filtered_14_NOTinv1)

L1_aligned_dataset_gstat_filtered_14_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(between(bp_cum, (232813426+5329939), (232813426+7791197)))
L1_aligned_dataset_gstat_filtered_14_inv$karyo <- "inv"

L1_aligned_dataset_gstat_filtered_14_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_14_inv, L1_aligned_dataset_gstat_filtered_14_NotINV_fig)

Gstat_chr14_density <- L1_aligned_dataset_gstat_filtered_14_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 14 Inversion")

figure4_density_gstat <- ggarrange(Gstat_chr5_density, Gstat_chr8_density, Gstat_chr8small_density, Gstat_chr14_density, nrow = 4, ncol = 1)


tiff("./figures/figure4_density_gstat.png", width = 20, height = 30, units = "cm",res = 300)
print(figure4_density_gstat)
dev.off()

#######################################################
## Checking for significance of elevated G stat in non inverted, inverted, and small inv in chr8
L1_aligned_dataset_gstat_filtered_8_inv_wosmall1 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))

L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined <- rbind(L1_aligned_dataset_gstat_filtered_8_inv_wosmall2, L1_aligned_dataset_gstat_filtered_8_inv_wosmall2)

L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined$karyo <- "invLARGE"


L1_aligned_dataset_gstat_filtered_8small_inv_anova <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+1032334), (118986430+1246126)))
L1_aligned_dataset_gstat_filtered_8small_inv_anova$karyo <- "invSMALL"

L1_aligned_dataset_gstat_filtered_8ALL_inv_anova <- rbind(L1_aligned_dataset_gstat_filtered_8small_inv_anova,L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined)

L1_aligned_dataset_gstat_filtered_8small_anova <- rbind(L1_aligned_dataset_gstat_filtered_8ALL_inv_anova, L1_aligned_dataset_gstat_filtered_8_NotINV_fig)

library(car)
summary(aov(G_stat ~ karyo,L1_aligned_dataset_gstat_filtered_8small_anova))
boxplot(G_stat ~ karyo,L1_aligned_dataset_gstat_filtered_8small_anova)
