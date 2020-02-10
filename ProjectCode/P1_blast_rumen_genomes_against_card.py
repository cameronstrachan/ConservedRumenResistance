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
from modules import seq_scrape as ss

# Hungate1000 collection and Stewart et al. 2018 genomes downloaded and concatenated into rumen_genomes.fasta
file = "rumen_genomes.fasta"

# Run prodigal on rumen_genomes.fasta
file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
file_obj.setOutputName(file)
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# CARD database was downloaded as card_db.fasta
# make blast database from from card_db.fasta

file = "card_db.fasta"
indir = 'dataflow/01-prot/'
blastdbdir = 'dataflow/02-blast-db/'
blastdir = 'dataflow/02-blast/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='prot')

# blast rumen genomes against CARD database

file = "rumen_genomes.fasta"

file_obj = sc.Fasta(file, indir)
file_obj.setOutputLocation(blastdir)

outputfilename = "rumen_genomes_card.txt"
blastdb = "card_db.fasta"

file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastp', db=blastdb, dblocation=blastdbdir, max_target_seqs=1, evalue=1e-3, num_threads = 60)
