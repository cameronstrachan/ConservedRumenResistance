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

# The ORFs from the rumen genomes were made into a blast database

file = "rumen_genomes.fasta"
blastdbdir = 'dataflow/02-blast-db/'

file_obj = sc.Fasta(file, 'dataflow/01-prot/')
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='prot')

# The 3 resistance gene ORFs originally found in the 3 genomes from bacteroidetes (ANT6, Sat4 and Aph3)
# were put into a file. These were then blasted against the

file = "3_resistance_genes.fasta"
indir = 'dataflow/01-prot/'


file_obj = sc.Fasta(file, indir)
file_obj.setOutputLocation(blastdir)

outputfilename = "3_resistance_genes.txt"
blastdb = "rumen_genomes.fasta"

file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastp', db=blastdb, dblocation=blastdbdir, max_target_seqs=10000, evalue=1e-4, num_threads = 60, max_hsps = 1)

# The csv file used is the blast output, subsetted for the hits for ANT6.

genes_df = pd.read_csv('dataflow/00-meta/ANT6_rumen.csv', low_memory=False)
genes = genes_df['sseqid'].tolist()

file_obj = sc.Fasta('rumen_genomes.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_rumen.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.subsetfasta(seqlist = genes, headertag='RUMEN')

# ANT6 ORFS from ncbi were trimmed above 250 aa and below 250 aa

file_obj = sc.Fasta('ANT6_ncbi.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_ncbi_250.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.lengthcutoff(replaceheaders = False, length = 250, direction = 'above')

file_obj = sc.Fasta('ANT6_ncbi_250.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_ncbi_250_350.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.lengthcutoff(replaceheaders = False, length = 350, direction = 'below')

# ANT6 ORFS from the rumen were trimmed above 250 aa and below 250 aa

file_obj = sc.Fasta('ANT6_rumen.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_rumen_250.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.lengthcutoff(replaceheaders = False, length = 250, direction = 'above')

file_obj = sc.Fasta('ANT6_rumen_250.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_rumen_250_350.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.lengthcutoff(replaceheaders = False, length = 350, direction = 'below')

# The trimmed ORFS were then combined into a single file (copy/paste: ANT6_ncbi_rumen_250_350.fasta) and then the headers were renamed.

file_obj = sc.Fasta('ANT6_ncbi_rumen_250_350.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('ANT6_ncbi_rumen_250_350_rename.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.saveonelinefasta()

# Muscle and fast free were the run with the following parameters. 

#../bin/muscle -in dataflow/01-prot/ANT6_ncbi_rumen_250_350_rename.fasta -out dataflow/03-alignments/ANT6_ncbi_rumen_250_350.afa -maxiters 3 -diags -sv -distance1 kbit20_3
#../bin/FastTree dataflow/03-alignments/ANT6_ncbi_rumen_250_350.afa > dataflow/03-trees/ANT6_ncbi_rumen_250_350.afa.newick
