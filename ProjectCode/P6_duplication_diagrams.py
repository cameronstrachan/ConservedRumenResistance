import os, sys
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

# custom libraries
system = str(input('\n' + 'Local or Server (L or S):'))

if system == 'S':
    sys.path.insert(0, '/home/strachan/master/')
else:
    sys.path.insert(0, '/Users/cameronstrachan/master/')

from modules import seq_core_lin as sc
from modules import seq_gen_lin as sg

# This extracted contigs for the gene duplication analysis (Supp. Fig. 3)

genes_df = pd.read_csv('dataflow/00-meta/genomes_with_ant6_duplication.csv', low_memory=False)
genes = genes_df['Accession'].tolist()

file_obj = sc.Fasta('fig1_fig3_ncbi_nucl_hits.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('pathogens_duplicates.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.subsetfasta(seqlist = genes, headertag='_duplicate')

# From the above file, the regions with 50kB of AadE-Ia or AadE-Ib were trimmed
# out and extracted in genious, then the ORFs were predicted

file = "duplicate_gene_diagrams_trimmed.fasta"

file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
file_obj.setOutputName(file)
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# All ORFs were then blasted against each other to do the synteny analysis (sup fig 3)

file = "duplicate_gene_diagrams_trimmed.fasta"
blastdbdir = 'dataflow/02-blast-db/'

file_obj = sc.Fasta(file, 'dataflow/01-prot/')
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='prot')

indir = 'dataflow/01-prot/'
blastdir = 'dataflow/02-blast/'
file = "duplicate_gene_diagrams_trimmed.fasta"

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdir)
outputfilename = "duplicate_gene_diagrams_trimmed.txt"
file_obj.setOutputName(outputfilename)

blastdb = "duplicate_gene_diagrams_trimmed.fasta"

file_obj.runblast(blast='blastp', db=blastdb, dblocation=blastdbdir, max_target_seqs=5000, evalue=1e-3, num_threads = 40, max_hsps = 1)
