## # Plotting output for F python scripts from Gould et al 2017 
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
files=list.files(pattern = ".fst")
L1_aligned_dataset_fst = do.call(rbind, lapply(files,fread))

# Removing windowns with 2 SD difference from mean

L1_SD <- L1_aligned_dataset_fst %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 227.9529
L1_filt <- (2*L1_SD)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

#Taking the top 1% of fst

obs_L1 <-  nrow(L1_aligned_dataset_fst_filtered)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst_filtered %>% arrange(desc(win_Fst))
L1_fst_top1 <- L1_aligned_dataset_fst_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)

#Changing chromosome so they are in order

L1_aligned_dataset_fst_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_fst_filtered$scaff))
L1_aligned_dataset_fst_filtered$scaff <- as.numeric(L1_aligned_dataset_fst_filtered$scaff)

data_cum_L1_fst <- L1_aligned_dataset_fst_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst_filtered %>% 
  inner_join(data_cum_L1_fst, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)

# Scatter plot for FST

L1_fst_plot <- L1_aligned_dataset_fst_filtered %>% 
  ggplot(aes(x=bp_cum, y= win_Fst)) +
  geom_point() +
  geom_hline(yintercept =0.01738773, color="red") +
  annotate(geom = "rect", xmin=(118943830+850429), xmax=(	 118943830+7604769), ymin=-Inf, ymax=Inf,
           color = "#2C77BF", 
           fill = "#2C77BF", alpha = 0.5)+
  annotate(geom = "rect", xmin=(118943830+1032334), xmax=(118943830+1246126), ymin=-Inf, ymax=Inf,
           color = "green", 
           fill = "green", alpha = 0.5)+
  annotate(geom = "rect", xmin=(67485786+13650670), xmax=(67485786+17847181), ymin=-Inf, ymax=Inf,
           color = "orange", 
           fill = "orange", alpha = 0.5)+
  annotate(geom = "rect", xmin=(233048607+5329939), xmax=(233048607+7791197), ymin=-Inf, ymax=Inf,
           color = "yellow", 
           fill = "yellow", alpha = 0.5)

  tiff("./figures/L1_fst_plot.png", width = 30, height = 20, units = "cm",res = 300)
  print(L1_fst_plot)
  dev.off()
