library(tidyverse)
library(stringi)

df <- read.delim("~/master/rumen/dataflow/02-blast/V1_V2_pathogens_rumen.txt", header=FALSE)
colnames(df) <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length", "sseq")

df <- df %>%
  select(-sseq) %>%
  filter(length > 200) %>%
  filter(pident > 60) %>%
  distinct()


for (i in 1:nrow(df)){
  df[i,"orf_id"] <- stri_reverse(str_split_fixed(stri_reverse(df[i,"sseqid"]), "_", 2)[[1]])
  df[i,"genome"] <- stri_reverse(str_split_fixed(stri_reverse(df[i,"sseqid"]), "_", 2)[[2]])
}

df <- df %>%
  separate(genome, into = c("genome", "contig"), sep = "_") %>%
  group_by(genome) %>%
  mutate(norfs_genome = length(unique(orf_id))) %>%
  ungroup()

df1 <- df %>%
  filter(norfs_genome >1) %>%
  rename(Accession = genome)

df2 <- read.csv("~/master/rumen/dataflow/00-meta/tableS2.csv")

df2$Accession <- as.character(df2$Accession)

df3 <- read_csv("~/master/rumen/dataflow/02-headers/pathogens_rumen.csv")

colnames(df3) <- c("sseqid", "header", "rm")

df <- inner_join(df1, df2) %>%
  select(qseqid, sseqid, pident, norfs_genome, Accession, Figure, Description.of.blast.hit..NCBI.,Host, Note) %>%
  distinct() %>%
  inner_join(df3) %>%
  select(-rm)

rm(list=setdiff(ls(), "df"))

df <- df %>%
  separate(header, into = c("rm", "start", "stop", "strand"), sep = "#")

df <- df %>%
 select(-rm, -qseqid, -pident, -Figure) %>%
 distinct()

write.csv(df, "~/master/rumen/dataflow/00-meta/genomes_with_ant6_duplication.csv")