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

# Genomes that were downloaded (Figure 1 and Figure 2) were combined into fig1_fig3_ncbi_nucl_hits.fasta
# Note the file name was made before fig 3 was moved to fig 2. These were then combined
# with the rumen genomes (pathogens_rumen.fasta) and made into a blast database.

file = "fig1_fig3_ncbi_nucl_hits.fasta"

file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
file_obj.setOutputName(file)
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

seqs_concatn = ['rumen_genomes.fasta', 'fig1_fig3_ncbi_nucl_hits.fasta']

sg.concat(inputfolder='dataflow/01-prot/', outputpath='dataflow/01-prot/pathogens_rumen.fasta', filenames=seqs_concatn)

file = "pathogens_rumen.fasta"
blastdbdir = 'dataflow/02-blast-db/'

file_obj = sc.Fasta(file, 'dataflow/01-prot/')
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='prot')

# The two version of ANT6 (v1_v2_4309680.fasta) were then blasted against the pathogen and rumen genomes.

indir = 'dataflow/01-prot/'
blastdir = 'dataflow/02-blast/'
file = "v1_v2_4309680.fasta"

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdir)
outputfilename = "V1_V2_pathogens_rumen.txt"
file_obj.setOutputName(outputfilename)

blastdb = "pathogens_rumen.fasta"

file_obj.runblast(blast='blastp', db=blastdb, dblocation=blastdbdir, max_target_seqs=2000, evalue=1e-3, num_threads = 40, max_hsps = 1)

# A header mapping file was created for the pathogen and rumen genomes.

headerfile = 'dataflow/02-headers/'
file = "pathogens_rumen.fasta"

file_obj = sc.Fasta(file, 'dataflow/01-prot/')
file_obj.setOutputName(file)
file_obj.setOutputLocation(headerfile)
headers = file_obj.fasta2headermap()
df = pd.DataFrame.from_dict(headers, orient="index")
df['file'] = file
df.to_csv(headerfile + file.split('.fa')[0] + '.csv')

# Genomes were selected that showed had a ANT5 duplication ((>200 amino acids and >60% identity to the versions from B. thetaiotoamicron (4309680))
# See R1_find_duplications
# This extracted ORFs that were used for the Figure 3 tree (iqtree)

genes_df = pd.read_csv('dataflow/00-meta/genomes_with_ant6_duplication.csv', low_memory=False)
genes = genes_df['sseqid'].tolist()

file_obj = sc.Fasta('pathogens_rumen.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('pathogens_duplicates.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.subsetfasta(seqlist = genes, headertag='_duplicate')
