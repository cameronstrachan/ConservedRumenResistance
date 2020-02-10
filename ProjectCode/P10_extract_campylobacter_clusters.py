# python libraries
import os, sys
import subprocess
import pandas as pd

# custom libraries
system = str(input('\n' + 'Local or Server (L or S):'))

if system == 'S':
    sys.path.insert(0, '/home/strachan/master/')
else:
    sys.path.insert(0, '/Users/cameronstrachan/master/')

from modules import seq_core_lin as sc
from modules import seq_gen_lin as sg

# after running cluster frequency comparison (Rmd), tables were exported
# that exported clusters with high levels of recombination from
# campylobacter coli and jejni

# extract the clusters

file_obj = sc.Fasta('campylobacter_coli.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_coli_extracted_clusters.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.extract_regions(df = 'dataflow/03-analysis/cluster_positions_ccoli.csv', col_start = "cluster_start", col_end = "cluster_end", col_contig = "pathogen_genome_id")

file_obj = sc.Fasta('campylobacter_jejuni.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_jejuni_extracted_clusters.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.extract_regions(df = 'dataflow/03-analysis/cluster_positions_cjejuni.csv', col_start = "cluster_start", col_end = "cluster_end", col_contig = "pathogen_genome_id")

# dereplicate

file_obj = sc.Fasta('campylobacter_coli_extracted_clusters.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_coli_extracted_clusters_derep.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.dereplicate()


file_obj = sc.Fasta('campylobacter_jejuni_extracted_clusters.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_jejuni_extracted_clusters_derep.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.dereplicate()

# run prodigal

file_obj = sc.Fasta('campylobacter_coli_extracted_clusters_derep.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_coli_extracted_clusters_derep.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

file_obj = sc.Fasta('campylobacter_jejuni_extracted_clusters_derep.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('campylobacter_jejuni_extracted_clusters_derep.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# but sequences on one line

file_obj = sc.Fasta('campylobacter_coli_extracted_clusters_derep.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('campylobacter_coli_extracted_clusters_derep_oneline.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.saveonelinefasta()


file_obj = sc.Fasta('campylobacter_jejuni_extracted_clusters_derep.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('campylobacter_jejuni_extracted_clusters_derep_oneline.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.saveonelinefasta()

# dereplicate prot seqs

file_obj = sc.Fasta('campylobacter_coli_extracted_clusters_derep_oneline.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('campylobacter_coli_extracted_clusters_derep_derep.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.dereplicate()


file_obj = sc.Fasta('campylobacter_jejuni_extracted_clusters_derep_oneline.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('campylobacter_jejuni_extracted_clusters_derep_derep.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.dereplicate()

# annotate genes against CARD
db_file = "card_db.fasta"
blastdb = 'dataflow/02-blast-db/'
blastout = 'dataflow/02-blast-out/'

file_obj = sc.Fasta('campylobacter_coli_extracted_clusters_derep_derep.fasta', 'dataflow/01-prot/')
file_obj.setOutputLocation(blastout)
outputfilename = "campylobacter_coli_extracted_clusters_derep_derep" + '_card.txt'
file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastp', db=db_file, dblocation=blastdb, max_target_seqs=1, evalue=1e-3, num_threads = 60)

file_obj = sc.Fasta('campylobacter_jejuni_extracted_clusters_derep_derep.fasta', 'dataflow/01-prot/')
file_obj.setOutputLocation(blastout)
outputfilename = "campylobacter_jejuni_extracted_clusters_derep_derep" + '_card.txt'
file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastp', db=db_file, dblocation=blastdb, max_target_seqs=1, evalue=1e-3, num_threads = 60)
