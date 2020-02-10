library(tidyverse)
library(stringi)

df_card <- read.delim("~/master/rumen2/dataflow/02-blast-out/stewart2019_mags_genes_300_cp_pathogen_mapped_card.txt", header=FALSE)
colnames(df_card)[1:13] <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq")

df_card <- df_card %>%
  filter(pident > 70) %>%
  separate(sseqid, into = c("gene", "organism"), sep = "\\[")

df_card$gene <- gsub("_", "", df_card$gene)
df_card$organism <- gsub("\\]_", "", df_card$organism)
df_card$organism <- gsub("_", " ", df_card$organism)

df_card <- df_card %>%
  separate(gene, into = c('rm', "accession", "ARD", "gene"), sep = "\\|")



for (i in 1:nrow(df_card)){
  df_card[i,"file"] <- stri_reverse(str_split_fixed(stri_reverse(df_card[i,"qseqid"]), "_", 3)[[3]])
}


df_meta <- read.csv("~/master/rumen2/dataflow/00-meta/stewart2019_ftp_meta.csv") %>%
  select(scientific_name, file)

df_meta$file <- gsub(".fa.gz", "", df_meta$file)

df_final <- inner_join(df_card, df_meta) %>%
  select(gene, scientific_name) %>%
  group_by(gene, scientific_name) %>%
  mutate(freq = length(gene)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(gene)

colnames(df_final) <- c("ARD", "Bin Assignment", "Frequency")
write.csv(df_final, "~/master/rumen2/newfigures/prevalence_gene_rumen.csv")

check <- df_card %>%
  filter(pident > 98)

length(unique(df_card$gene))
length(unique(check$gene))

df_cluster <- inner_join(df_card, df_meta) %>%
  group_by(file) %>%
  mutate(num_genes_mag = length(unique(gene))) %>%
  ungroup() %>%
  filter(num_genes_mag > 1)
  
