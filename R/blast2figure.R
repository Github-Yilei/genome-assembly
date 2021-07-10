library(ggplot2)
file <- "C:\\Users\\Admin\\Desktop\\"
df1 <- read.delim(paste0(file, 'parsed_mapped_blastn.txt'), header = FALSE)
names(df1) <- c("qseqid", "sseqid", "pident", "nident", "qlen",
               "slen", "evalue", "bitscore")

df2 <- read.delim(paste0(file, 'parsed_unmapped_blastn.txt'), header = FALSE)
names(df2) <- c("qseqid", "sseqid", "pident", "nident", "qlen",
               "slen", "evalue", "bitscore")

ggplot() + geom_density(data = df1, aes(x = nident, fill = "mapped"))+
   geom_density(data = df2, aes(x = nident*10, fill = "unmapped")) +  
   scale_x_continuous(sec.axis = sec_axis(~ ./10, name = "nident")) + 
  theme_bw()
