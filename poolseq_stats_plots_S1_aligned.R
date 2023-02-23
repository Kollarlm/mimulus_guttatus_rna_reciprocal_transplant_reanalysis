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

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/S1_aligned"
setwd(path)

# Files for pi
files=list.files(pattern = ".pi")
S1_aligned_dataset_pi = do.call(rbind, lapply(files,fread))


# Files for gstat
files=list.files(pattern = ".gstat")
S1_aligned_dataset_gstat = do.call(rbind, lapply(files,fread))

# Files for Tajimas D
files=list.files(pattern = ".tajD")
S1_aligned_dataset_tajD = do.call(rbind, lapply(files,fread))

# Removed NAs for tajimas D
S1_aligned_dataset_nona_tajD <- na.omit(S1_aligned_dataset_tajD)
colnames(S1_aligned_dataset_nona_tajD) <- c("gene", "scaffold", "gene_start", "gene_end","gene_len", "IA_gene_cov", "IA_S", "TajD_IA", "CP_gene_cov", "CP_S", "TajD_CP", "tot_gene_cov", "tot_S", "TajD_tot")

## Removed windows with average depth greater than 2 standard deviations

S1_SD <- S1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 179.8778
S1_filt <- (2*S1_SD)

S1_mean <- S1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=mean(win_depth, na.rm=TRUE)) # mean = 166.059

S1_aligned_dataset_pi_filtered <- S1_aligned_dataset_pi %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

# gstat
S1_SD_g <- S1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=sd(avg.win.depth, na.rm=TRUE)) # standard deviation = 179.8778
S1_filt_g <- (2*S1_SD_g)

S1_mean_g <- S1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=mean(avg.win.depth, na.rm=TRUE)) # mean = 

S1_aligned_dataset_gstat_filtered <- S1_aligned_dataset_gstat %>% # removes data greater than 2 SD from mean
  filter(abs(scale(avg.win.depth)) < 2)

## Adding in column for the ratio of pi

S1_aligned_dataset_pi_filtered$piRatio <- (S1_aligned_dataset_pi_filtered$freqIA/S1_aligned_dataset_pi_filtered$freqCP)

## Take top and bottom 1% of piRatio
obs_S1 <-  nrow(S1_aligned_dataset_pi_filtered)

S1_aligned_dataset_pi_filtered <- S1_aligned_dataset_pi_filtered %>% arrange(desc(piRatio))
S1_pi_top1 <- S1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_S1 * 0.01)

S1_aligned_dataset_pi_filtered <- S1_aligned_dataset_pi_filtered %>% arrange((piRatio))
S1_pi_bottom1 <- S1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_S1 * 0.01)

## Take top 1% of gstat

obs_S1_gstat <-  nrow(S1_aligned_dataset_gstat_filtered)

S1_aligned_dataset_gstat_filtered <- S1_aligned_dataset_gstat_filtered %>% arrange(desc(G_stat))
S1_pi_top1_gstat <- S1_aligned_dataset_gstat_filtered %>% 
  filter(row_number() < obs_S1_gstat * 0.01)

######## Preparing the data to plot manhattan plot in ggplot ########

# Remvoing Chr so that the chromsomes are in order

S1_aligned_dataset_pi_filtered$scaff <- gsub("chr","", as.character(S1_aligned_dataset_pi_filtered$scaff))
S1_aligned_dataset_pi_filtered$scaff <- as.numeric(S1_aligned_dataset_pi_filtered$scaff)

S1_aligned_dataset_gstat_filtered$scaff <- gsub("chr","", as.character(S1_aligned_dataset_gstat_filtered$scaff))
S1_aligned_dataset_gstat_filtered$scaff <- as.numeric(S1_aligned_dataset_gstat_filtered$scaff)

data_cum_S1_pi <- S1_aligned_dataset_pi_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

S1_aligned_dataset_pi_filtered <- S1_aligned_dataset_pi_filtered %>% 
  inner_join(data_cum_S1_pi, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)
  
data_cum_S1_gstat <- S1_aligned_dataset_gstat_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

S1_aligned_dataset_gstat_filtered <- S1_aligned_dataset_gstat_filtered %>% 
  inner_join(data_cum_S1_gstat, by ="scaff") %>% 
  mutate(bp_cum=win_end + bp_add)

## Scatter plot of piRatio

S1_pi_plot <- ggplot(S1_aligned_dataset_pi_filtered, aes(x=bp_cum, y= log10(piRatio))) +
  geom_point() +
  geom_hline(yintercept =log10(3.510945), color="red") +
  geom_hline(yintercept = log10(0.7889974), color="red") +
  annotate(geom = "rect", xmin=(120007578+858736), xmax=(120007578+6465310), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(120007578+5998986), xmax=(120007578+6260288), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67165139+14125309), xmax=(67165139+18136144), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(233542315+5892556), xmax=(233542315+8677481), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

  

# Scatter plot of read depth
S1_depth_plot <- ggplot(S1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = avg.win.depth)) +
  stat_smooth( method = "loess", span = 0.1, se = FALSE)+
  annotate(geom = "rect", xmin=(120007578+858736), xmax=(120007578+6465310), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(120007578+5998986), xmax=(120007578+6260288), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67141889+14125309), xmax=(67141889+18136144), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(233542315+5892556), xmax=(233542315+8677481), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)
  

## Scatter plot of gstat
S1_gstat_plot <- ggplot(S1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.71607
), color="red") +
  annotate(geom = "rect", xmin=(120007578+858736), xmax=(120007578+6465310), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(120007578+5998986), xmax=(120007578+6260288), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67141889+14125309), xmax=(67141889+18136144), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(233542315+5892556), xmax=(233542315+8677481), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

## Arrange plots on one graph (Figure2)

S1_alinged_Fig2 <- ggarrange(S1_gstat_plot, S1_pi_plot,S1_depth_plot,nrow = 3, ncol = 1 )

tiff("./figures/S1_aligned_gouldFig2.png", width = 30, height = 20, units = "cm",res = 300)
print(S1_alinged_Fig2)
dev.off()

## Figure 3 from Gould et al 2017 (G stat values across chromosmes 5, 8, and 14)


S1_gstat_plot_8 <- S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.71607), color="red") +
  annotate(geom = "rect", xmin=(120007578+858736), xmax=(120007578+6465310), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(120007578+5998986), xmax=(120007578+6260288), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)

S1_gstat_plot_5 <- S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.71607), color="red") +
  annotate(geom = "rect", xmin=(67141889+14125309), xmax=(67141889+18136144), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)


S1_gstat_plot_14 <- S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="red") +
  annotate(geom = "rect", xmin=(233542315+5892556), xmax=(233542315+8677481), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

## Arrange plots on one graph (Figure3)

S1_alinged_Fig3 <- ggarrange(S1_gstat_plot_5, S1_gstat_plot_8, S1_gstat_plot_14,nrow = 3, ncol = 1 )

tiff("./figures/S1_aligned_gouldFig3.png", width = 30, height = 20, units = "cm",res = 300)
print(S1_alinged_Fig3)
dev.off()

## Distribution of G Statistic window values inside and outside the chromosome inversion (figure 4)

# Add column for within inversion or not within
S1_aligned_dataset_gstat_filtered_5_NOTinv1 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum < 67169025+14125309) 
  S1_aligned_dataset_gstat_filtered_5_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_5_NOTinv2 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum > 67169025+18136144)
S1_aligned_dataset_gstat_filtered_5_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_5_NotINV_fig <- rbind(S1_aligned_dataset_gstat_filtered_5_NOTinv2, S1_aligned_dataset_gstat_filtered_5_NOTinv1)
  
S1_aligned_dataset_gstat_filtered_5_inv <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(between(bp_cum, (67169025+14125309),(67169025+18136144) ))
S1_aligned_dataset_gstat_filtered_5_inv$karyo <- "inv"
 
S1_aligned_dataset_gstat_filtered_5_figure4 <- rbind(S1_aligned_dataset_gstat_filtered_5_inv, S1_aligned_dataset_gstat_filtered_5_NotINV_fig)

Gstat_chr5_density <- S1_aligned_dataset_gstat_filtered_5_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 5 Inversion")

################################################################################

S1_aligned_dataset_gstat_filtered_8_NOTinv1 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (120007578+858736)) 

S1_aligned_dataset_gstat_filtered_8_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_8_NOTinv2 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (120007578+6465310))
S1_aligned_dataset_gstat_filtered_8_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_8_NotINV_fig <- rbind(S1_aligned_dataset_gstat_filtered_8_NOTinv1, S1_aligned_dataset_gstat_filtered_8_NOTinv2)

S1_aligned_dataset_gstat_filtered_8_inv <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (120007578+858736),(120007578+6465310)))
S1_aligned_dataset_gstat_filtered_8_inv$karyo <- "inv"

S1_aligned_dataset_gstat_filtered_8_inv_wosmalS1 <-  S1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (120007578+5998986)) 

S1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  S1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (120007578+6260288))

S1_aligned_dataset_gstat_filtered_8_figure4 <- rbind(S1_aligned_dataset_gstat_filtered_8_inv_wosmalS1, S1_aligned_dataset_gstat_filtered_8_inv_wosmall2, S1_aligned_dataset_gstat_filtered_8_NotINV_fig)

Gstat_chr8_density <- S1_aligned_dataset_gstat_filtered_8_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 8 Inversion")

################################################################################

S1_aligned_dataset_gstat_filtered_8small_NOTinv1 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (120007578+5998986)) 

S1_aligned_dataset_gstat_filtered_8small_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_8small_NOTinv2 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (120007578+6260288))
S1_aligned_dataset_gstat_filtered_8small_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_8small_NotINV_fig <- rbind(S1_aligned_dataset_gstat_filtered_8_NOTinv1, S1_aligned_dataset_gstat_filtered_8_NOTinv2)

S1_aligned_dataset_gstat_filtered_8small_inv <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (120007578+5998986), (120007578+6260288)))
S1_aligned_dataset_gstat_filtered_8small_inv$karyo <- "inv"

S1_aligned_dataset_gstat_filtered_8small_figure4 <- rbind(S1_aligned_dataset_gstat_filtered_8small_inv, S1_aligned_dataset_gstat_filtered_8small_NotINV_fig)

Gstat_chr8small_density <- S1_aligned_dataset_gstat_filtered_8small_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
  labs("Chromosome 8small Inversion")
############################################################################################################################

S1_aligned_dataset_gstat_filtered_14_NOTinv1 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum < (233542315+5892556)) 

S1_aligned_dataset_gstat_filtered_14_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_14_NOTinv2 <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum > (233542315+8677481))
S1_aligned_dataset_gstat_filtered_14_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_gstat_filtered_14_NotINV_fig <- rbind(S1_aligned_dataset_gstat_filtered_14_NOTinv2, S1_aligned_dataset_gstat_filtered_14_NOTinv1)

S1_aligned_dataset_gstat_filtered_14_inv <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(between(bp_cum, (233542315+5892556), (233542315+8677481)))
S1_aligned_dataset_gstat_filtered_14_inv$karyo <- "inv"

S1_aligned_dataset_gstat_filtered_14_figure4 <- rbind(S1_aligned_dataset_gstat_filtered_14_inv, S1_aligned_dataset_gstat_filtered_14_NotINV_fig)

Gstat_chr14_density <- S1_aligned_dataset_gstat_filtered_14_figure4 %>% 
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
S1_aligned_dataset_gstat_filtered_8_inv_wosmalS1 <-  S1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (120007578+5998986)) 

S1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  S1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (120007578+6260288))

S1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined <- rbind(S1_aligned_dataset_gstat_filtered_8_inv_wosmall2, S1_aligned_dataset_gstat_filtered_8_inv_wosmall2)

S1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined$karyo <- "invLARGE"


S1_aligned_dataset_gstat_filtered_8small_inv_anova <-  S1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (120007578+5998986), (120007578+6260288)))
S1_aligned_dataset_gstat_filtered_8small_inv_anova$karyo <- "invSMALL"

S1_aligned_dataset_gstat_filtered_8ALL_inv_anova <- rbind(S1_aligned_dataset_gstat_filtered_8small_inv_anova,S1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined)

S1_aligned_dataset_gstat_filtered_8small_anova <- rbind(S1_aligned_dataset_gstat_filtered_8ALL_inv_anova, S1_aligned_dataset_gstat_filtered_8_NotINV_fig)

library(car)
summary(aov(G_stat ~ karyo,S1_aligned_dataset_gstat_filtered_8small_anova))
boxplot.8.gstat <- boxplot(G_stat ~ karyo,S1_aligned_dataset_gstat_filtered_8small_anova)

png("./figures/boxplot.8.gstat.png", width = 20, height = 30, units = "cm",res = 300)
print(boxplot.8.gstat)
dev.off()

######### barplots of Tajimas D ##########

S1_aligned_dataset_nona_tajD$scaffold <- gsub("chr","", as.character(S1_aligned_dataset_nona_tajD$scaffold))
S1_aligned_dataset_nona_tajD$scaffold <- as.numeric(S1_aligned_dataset_nona_tajD$scaffold)

data_cum_S1_D <- S1_aligned_dataset_nona_tajD %>% 
  group_by(scaffold) %>% 
  summarise(max_bp=max(gene_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaffold, bp_add)

S1_aligned_dataset_nona_tajD <- S1_aligned_dataset_nona_tajD %>% 
  inner_join(data_cum_S1_D, by ="scaffold") %>% 
  mutate(bp_cum = gene_end + bp_add)

# Adding inverted vs noninverted to the chromosomes of interest.
S1_aligned_dataset_nona_tajD_5_NOTinv1 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum < 66634174+14125309) 
S1_aligned_dataset_nona_tajD_5_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_5_NOTinv2 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum > 66634174+18136144)
S1_aligned_dataset_nona_tajD_5_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_5_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_5_NOTinv2, S1_aligned_dataset_nona_tajD_5_NOTinv1)

S1_aligned_dataset_nona_tajD_5_inv <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(between(bp_cum, (66634174+14125309),(66634174+18136144) ))
S1_aligned_dataset_nona_tajD_5_inv$karyo <- "inv 5"

S1_aligned_dataset_nona_tajD_5 <- rbind(S1_aligned_dataset_nona_tajD_5_inv, S1_aligned_dataset_nona_tajD_5_NotINV_fig)

################################################################################

S1_aligned_dataset_nona_tajD_8_NOTinv1 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118479470+858736)) 
S1_aligned_dataset_nona_tajD_8_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_8_NOTinv2 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118479470+6465310))
S1_aligned_dataset_nona_tajD_8_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_8_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_8_NOTinv1, S1_aligned_dataset_nona_tajD_8_NOTinv2)

S1_aligned_dataset_nona_tajD_8_inv <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118479470+858736),(118479470+6465310)))
S1_aligned_dataset_nona_tajD_8_inv$karyo <- "inv 8"

#removing small inversion from 8 inversion
S1_aligned_dataset_nona_tajD_8_inv_wosmalS1 <-  S1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118479470+5998986)) 

S1_aligned_dataset_nona_tajD_8_inv_wosmall2 <-  S1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118479470+6260288))

S1_aligned_dataset_nona_tajD_8 <- rbind(S1_aligned_dataset_nona_tajD_8_inv_wosmalS1, S1_aligned_dataset_nona_tajD_8_inv_wosmall2, S1_aligned_dataset_nona_tajD_8_NotINV_fig)

################################################################################

S1_aligned_dataset_nona_tajD_8small_inv <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118479470+5998986), (118479470+6260288)))
S1_aligned_dataset_nona_tajD_8small_inv$karyo <- "inv 8 SMALL"

############################################################################################################################

S1_aligned_dataset_nona_tajD_14_NOTinv1 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum < (232192457+5892556)) 

S1_aligned_dataset_nona_tajD_14_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_14_NOTinv2 <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum > (232192457+8677481))
S1_aligned_dataset_nona_tajD_14_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_14_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_14_NOTinv2, S1_aligned_dataset_nona_tajD_14_NOTinv1)

S1_aligned_dataset_nona_tajD_14_inv <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(between(bp_cum, (232192457+5892556), (232192457+8677481)))
S1_aligned_dataset_nona_tajD_14_inv$karyo <- "inv 14"

S1_aligned_dataset_nona_tajD_14 <- rbind(S1_aligned_dataset_nona_tajD_14_inv, S1_aligned_dataset_nona_tajD_14_NotINV_fig)

#######################################################

# TajD dataset
### Inverted genes didnt pass?
tajd <- rbind(S1_aligned_dataset_nona_tajD_5, S1_aligned_dataset_nona_tajD_8,S1_aligned_dataset_nona_tajD_8small_inv,S1_aligned_dataset_nona_tajD_14)

S1_aligned_dataset_nona_tajD_noninv <-  S1_aligned_dataset_nona_tajD %>% 
  filter(scaffold %in% c("1","2","3","4","6","7","9","10","11","12","13"))

S1_aligned_dataset_nona_tajD_noninv$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_final <- rbind(S1_aligned_dataset_nona_tajD_noninv,tajd)

data1 <- S1_aligned_dataset_nona_tajD_final[,c(2,8,17)]
colnames(data1) <- c("chromosome", "TajD", "Karyo")
data1$tajtype <- "IA"


data2 <- S1_aligned_dataset_nona_tajD_final[,c(2,11,17)]
colnames(data2) <- c("chromosome", "TajD", "Karyo")
data2$tajtype <- "CP"

data3 <- S1_aligned_dataset_nona_tajD_final[,c(2,14,17)]
colnames(data3) <- c("chromosome", "TajD", "Karyo")
data3$tajtype <- "total"

loath <- rbind(data1, data2, data3)

tajD_pot_S1 <- loath %>% 
  ggplot(aes(fill=tajtype, y= TajD, x= Karyo)) +
  geom_boxplot()

png("./figures/tajD_pot_S1.png", width = 20, height = 30, units = "cm",res = 300)
print(tajD_pot_S1)
dev.off()

write.csv(S1_pi_top1_gstat, file = "S1_pi_top1_gstat.csv")
write.csv(S1_pi_top1, file="S1_pi_top1.csv")
write.csv(S1_pi_top1, file="S1_pi_bottom1.csv")