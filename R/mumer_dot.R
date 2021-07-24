file <- "C:\\Users\\Admin\\Desktop\\"
df1 <- read.delim(paste0(file, "tri_chr_filtered.coords"), header = F, 
                       sep = "\t",
                       stringsAsFactors = FALSE)
colnames(df1) <- c("ref_start", "ref_end", "qry_start", "qry_end", "ref_len", "qry_len", 
                  "identiy", "ref_tag","qry_tag" )

df <- subset(df1,ref_tag == "sampe.unique.REduced.paired_only.counts_AAGCTT.9g1" & 
               qry_tag == "sampe.unique.REduced.paired_only.counts_AAGCTT.9g1_RagTag_2")
x_range <- range(c(df$ref_start, df$ref_end))
y_range <- range(c(df$qry_start, df$qry_end))
plot.new()

plot.window(xlim = x_range,
            ylim = y_range)

for( i in 1:nrow(df)){
  
  if (df[i,3] < df[i,4]){
    lines(x = df[i,1:2], y = df[i,3:4], col = "red")
  } else{
    lines(x = df[i,1:2], y = df[i,3:4], col = "blue")
  }
  
}

box()
axis(1, at = seq(0, x_range[2], 10000), labels = seq(0, x_range[2], 10000) / 10000)
axis(2, at = seq(0, y_range[2], 10000), labels = seq(0, y_range[2], 10000) / 10000)
