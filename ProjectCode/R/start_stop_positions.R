library(tidyverse)

df_card <- read.delim("~/master/rumen2/dataflow/02-blast-out/stewart2019_mags_genes_300_cp_pathogen_mapped_card.txt", header=FALSE)
colnames(df_card)[1:13] <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq")

df_card <- df_card %>%
  filter(pident > 95) %>%
  separate(sseqid, into = c("gene", "organism"), sep = "\\[")

df_card$gene <- gsub("_", "", df_card$gene)
df_card$organism <- gsub("\\]_", "", df_card$organism)
df_card$organism <- gsub("_", " ", df_card$organism)

df_card <- df_card %>%
  separate(gene, into = c('rm', "accession", "ARD", "gene"), sep = "\\|")

#

files <- list.files('~/master/rumen2/dataflow/03-analysis/', pattern = '^start_stop')
list_dfs <- list()
i <- 1

for (file in files){
  df <- read.csv(paste('~/master/rumen2/dataflow/03-analysis/', file, sep = ''), header=TRUE)
  df$file <- file
  list_dfs[[i]] <- df
  i <- i + 1
}

df <- bind_rows(list_dfs)

df$file <- gsub("start_stop_stewart2019_mags_genes_300_cp_", "", df$file)
df$file <- gsub(".txt", "", df$file)
df$file <- gsub("_", " ", df$file)
df$X <- NULL

df_start_stop <- df_card %>%
  select(qseqid, gene) %>%
  inner_join(df)

colnames(df_start_stop) <- c('gene_name', 'ard', 'pathogen_genome_id', 'pident',  'pathogen_start', 'pathogen_end', 'rumen_start', 'rumen_end', 'rumen_direction', 'pathogen')

df_start_stop$pathogen_start_dir <- NA
df_start_stop$pathogen_end_dir <- NA
df_start_stop$rumen_start_dir <- NA
df_start_stop$rumen_end_dir <- NA

for (x in 1:nrow(df_start_stop)){
  df_start_stop[x, 'pathogen_start_dir'] <- min(df_start_stop[x, 'pathogen_start'], df_start_stop[x, 'pathogen_end'])
  df_start_stop[x, 'pathogen_end_dir'] <- max(df_start_stop[x, 'pathogen_start'], df_start_stop[x, 'pathogen_end'])
  df_start_stop[x, 'rumen_start_dir'] <- min(df_start_stop[x, 'rumen_start'], df_start_stop[x, 'rumen_end'])
  df_start_stop[x, 'rumen_end_dir'] <- max(df_start_stop[x, 'rumen_start'], df_start_stop[x, 'rumen_end'])
}


df_start_stop_clean <- df_start_stop %>%
  select(gene_name, ard, pathogen, pathogen_genome_id, pident, pathogen_start_dir, pathogen_end_dir, rumen_start_dir, rumen_end_dir)

write.csv(df_start_stop_clean, '~/master/rumen2/dataflow/03-analysis/compiled_start_stop.txt')
