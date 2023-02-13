# Plotting output pf python scripts from Gould et al 2017 
# Leslie M Kollar

## Loading in libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpubr)

## Loading in data files

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/S1_aligned"
setwd(path)

# Files for pi
files=list.files(pattern = "matching.pi")
S1_aligned_dataset_pi = do.call(rbind, lapply(files,fread))

# Files for gstat
files=list.files(pattern = "matching.gstat")
S1_aligned_dataset_gstat = do.call(rbind, lapply(files,fread))

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
obs_S1 <-  nrow(S1_aligned_dataset_gstat_filtered)

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

## Scatter plot of piRatio

S1_pi_plot <- ggplot(S1_aligned_dataset_pi_filtered, aes(x=win.start, y = log10(piRatio))) +
  geom_point() +
  geom_hline(yintercept =log10(5.008270), color="red") +
  geom_hline(yintercept = log10(0.6356174), color="red")

## Scatter plot of read depth
S1_depth_plot <- ggplot(S1_aligned_dataset_gstat_filtered, aes(x=win_start, y = avg.win.depth)) +
  stat_smooth( method = "loess", span = 0.1, se = FALSE)
  

## Scatter plot of gstat
S1_gstat_plot <- ggplot(S1_aligned_dataset_gstat_filtered, aes(x=win_start, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(33.67083), color="red") 

## Arrange plots on one graph

S1_alinged_Fig1 <- ggarrange(S1_gstat_plot, S1_pi_plot,S1_depth_plot,nrow = 3, ncol = 1 )

tiff("./S1_aligned_gouldFig1.png", width = 30, height = 20, units = "cm",res = 300)
print(S1_alinged_Fig1)
dev.off()




