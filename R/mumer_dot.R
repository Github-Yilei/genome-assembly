file <- "C:\\Users\\Admin\\Desktop\\"
df1 <- read.delim(paste0(file, "tri_chr_filtered.coords"), header = F, 
                       sep = "\t",
                       stringsAsFactors = FALSE)

colnames(df1) <- c("ref_start", "ref_end", "qry_start", "qry_end", "ref_len", "qry_len", 
                  "identiy", "ref_tag","qry_tag" )

df <- subset(df1,ref_tag == "sampe.unique.REduced.paired_only.counts_AAGCTT.9g1" & 
               qry_tag == "sampe.unique.REduced.paired_only.counts_AAGCTT.9g1_RagTag_1")


df$structure <- ifelse(df[,3] < df[,4], "normal", "inversion")


ggplot(data = df) + geom_segment(aes(x = ref_start/10000, y = qry_start/10000,
                                     xend = ref_end/10000, yend = qry_end/10000, 
                                     col = structure), size = 1) + 
  scale_color_manual(values = c("blue", "red")) + 
  theme_bw() + 
  xlab(label = "Canu chr1") +
  ylab(label = "Hap1 chr1")
