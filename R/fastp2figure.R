file <- "C:\\Users\\Admin\\Desktop\\parse_fastp_result.txt"
df <- read.delim(file, header = T, sep = ',')


df_spl <- tidyr::separate(data = df, 
                          col = X, 
                          into = c("sample", "class"), 
                          sep = ";")
