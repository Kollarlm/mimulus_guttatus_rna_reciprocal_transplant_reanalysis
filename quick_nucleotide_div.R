library(ggplot2)
library(tidyverse)
library(tidyr)

## These files were generated with bwa -> 
path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis"
setwd(path)

pi.S1 <- read.table("./data/S1/S1_aligned_S1_gould_WGS_sort_duplicates_bcf.windowed.pi", header=TRUE)
pi.L1 <- read.table("./data/L1/L1_aligned_L1_gould_WGS_sort_duplicates_bcf.windowed.pi", header=TRUE)


hist(pi.S1$PI,br=20)
hist(pi.L1$PI,br=20)


boxplot(pi.all$PI,ylab="diversity")

pi.chr8 <- subset(pi.all, pi.all$CHROM == "chr8")

plot(pi.chr8$BIN_END,pi.chr8$PI,xlab="position",ylab="diversity")

S1.pi.plot <- pi.all %>% 
  filter(CHROM == "chr8") %>% 
  ggplot(aes(x=BIN_END, y= PI)) +
  geom_point()+
  geom_vline(xintercept = 858736, color="red") +
  geom_vline(xintercept = 6465310, color="red") +
  geom_vline(xintercept = 5998986, color= "blue")+
  geom_vline(xintercept = 6260288, color= "blue")

tiff("./figures/S1.pi.plot.png", width = 30, height = 20, units = "cm",res = 300)
print(S1.pi.plot)
dev.off()

L1.pi.plot <- pi.L1 %>% 
  filter(CHROM == "chr8") %>% 
  ggplot(aes(x=BIN_END, y= PI)) +
  geom_point()+
  geom_vline(xintercept = 850429, color="red") +
  geom_vline(xintercept = 7604769, color="red") +
  geom_vline(xintercept = 1032334, color= "blue")+
  geom_vline(xintercept = 1246126, color= "blue")

tiff("./figures/L1.pi.plot.png", width = 30, height = 20, units = "cm",res = 300)
print(L1.pi.plot)
dev.off()
