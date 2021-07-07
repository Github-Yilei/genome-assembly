
#############
# Description: this script is used to estimate the genome size by the data from jellyfish histo
# Author: YiLei
# Date: 2021/06/01
############
```
file <- jellyfish_17k.histo

kmer_count <- read.delim(file, head = F, sep = " ", stringsAsFactors = F)

# fixing data type and colnames.
names(kmer_count) <- c('depth', 'frequencies')
kmer_count$depth <- as.numeric(kmer_count$depth)
kmer_count$frequencies <- as.numeric(kmer_count$frequencies)

kmer_count$sum <- kmer_count$depth * kmer_count$frequencies
kmer_count$freq <- kmer_count$sum/sum(kmer_count$sum)*100

library(ggplot2)
windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),
             TMN=windowsFont("Times New Roman"),
             ARL=windowsFont("Arial"))

default_theme <- theme(panel.background = element_blank(), 
      panel.grid.major.y = element_line(colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      text = element_text(family = "TMN"))

ggplot(data = kmer_count) + geom_point(aes(x = depth, y = freq))+
  xlim(3, 150) +
  ylim(0, 1.5) + 
  geom_vline(xintercept = Kdepth, color = "red") +
  labs(x = "depth",
       y = "Frequency(%)", 
       title = "Kmer distribution") + 
  default_theme


# genome size
## The 1st or 2ed point might be error peaks that can be removed .
Knum <- sum(kmer_count$sum[3:nrow(kmer_count)])
Kdepth <- which(kmer_count[,3] == max(kmer_count[c(3:1000), 3]))
genome_size <-  Knum/Kdepth/1000000
print(genome_size)


