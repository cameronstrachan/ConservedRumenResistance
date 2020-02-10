library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Must supply blast output files as argument", call.=FALSE)
} 

files <- c()

for (x in 1:length(args)){
  files <- c(files, args[x])
}

for (file in files){

  input_file <- paste("~/master/rumen2/dataflow/02-blast-out/", file, sep = "")
  file_prefix <- as.character(strsplit(file, ".txt")[1])
  output_file <- paste("~/master/rumen2/dataflow/03-analysis/", file_prefix, '.csv', sep = '')

  df <- read.delim(input_file, header=FALSE)
  colnames(df)[1:13] <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq")

  df <- df %>%
    mutate(per_cov = (length / qlen) * 100) %>%
    filter(per_cov > 98) %>%
    filter(pident > 99)

  df_freq <- as.data.frame(table(df$qseqid)) %>%
    filter(Freq > 0) %>%
    arrange(desc(Freq))
  
  colnames(df_freq) <- c('gene', 'frequency')
  
  write.csv(df_freq, output_file)
  
}