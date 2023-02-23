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

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/S1_aligned"
setwd(path)

# Files for pi
files=list.files(pattern = ".fst")
S1_aligned_dataset_fst = do.call(rbind, lapply(files,fread))

# Removing windowns with 2 SD difference from mean

S1_SD <- S1_aligned_dataset_fst %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 227.9529
S1_filt <- (2*S1_SD)

S1_aligned_dataset_fst_filtered <- S1_aligned_dataset_fst %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

#Taking the top 1% of fst

obs_S1 <-  nrow(S1_aligned_dataset_fst_filtered)

S1_aligned_dataset_fst_filtered <- S1_aligned_dataset_fst_filtered %>% arrange(desc(win_Fst))
S1_fst_top1 <- S1_aligned_dataset_fst_filtered %>% 
  filter(row_number() < obs_S1 * 0.01)

#Changing chromosome so they are in order

S1_aligned_dataset_fst_filtered$scaff <- gsub("chr","", as.character(S1_aligned_dataset_fst_filtered$scaff))
S1_aligned_dataset_fst_filtered$scaff <- as.numeric(S1_aligned_dataset_fst_filtered$scaff)

data_cum_S1_fst <- S1_aligned_dataset_fst_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

S1_aligned_dataset_fst_filtered <- S1_aligned_dataset_fst_filtered %>% 
  inner_join(data_cum_S1_fst, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)

