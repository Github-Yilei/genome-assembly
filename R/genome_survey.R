
#############
# Description: this script is used to estimate the genome size by the data from jellyfish histo
# Author: YiLei
# Date: 2021/06/01
############
```
file <- jellyfish_17k.histo

kmer_count <- read.delim(file, head = F, sep = " ", stringsAsFactors = F)
names(kmer_count) <- c('bin', 'frequencies')

kmer_count$sum <- kmer_count[,1] * kmer_count[,2]

kmer_count$freq <- kmer_count$sum/sum(kmer_count$sum)*100

Knum <- sum(kmer_count$sum)

Kdepth <- which(kmer_count[,3] == max(kmer_count[c(1:1000), 3]))

library(ggplot2)

windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),
             TMN=windowsFont("Times New Roman"),
             ARL=windowsFont("Arial"))

default_theme <- theme(panel.background = element_blank(), 
      panel.grid.major.y = element_line(colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      text = element_text(family = "TMN"))

ggplot(data = kmer_count) + geom_point(aes(x = bin, y = freq))+
  xlim(0, 100) + geom_vline(xintercept = Kdepth, color = "red") +
  labs(x = "depth",
       y = "Frequency(%)", 
       title = "Kmer distribution") + 
  default_theme

# genome size
## The first error peak can be removed 
genome_size <-  Knum/Kdepth/1000000
print(genome_size)
