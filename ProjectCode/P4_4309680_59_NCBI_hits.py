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

# from phylogenetic analysis, decided to repeat the analysis carried out in figure 1, using the contig 4309680-submission.assembly_59, which
# contained a second ANT6 version. also included other contigs of interest from the tree. these were subsetted into a files called subclade_island.fasta

genes = ['4309680-submission.assembly_59', '3964017-submission.assembly_7', '3643350-assembly_6', '3394949-submission.assembly_17', 'RUG117_52']

file_obj = sc.Fasta('rumen_genomes.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('subclade_island.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.subsetfasta(seqlist = genes, headertag='none')

# 4309680-submission.assembly_59 was then blasted against the
# NCBI nucleotide collection (nr/nt) using web-based blastn and the full-length sequence for each of the top 50 hits
# was downloaded and concatenated into island2_pathogens.fasta. With the contigs selected above into rumen_genomes_island2_pathogens.fasta.


files = ['island2_pathogens.fasta', 'subclade_island.fasta']
sg.concat(inputfolder='dataflow/01-nucl/', outputpath='dataflow/01-nucl/rumen_genomes_island2_pathogens.fasta', filenames=files)

# a blast database was then made with the contigs of interest, including 4309680-submission.assembly_59

file = "subclade_island.fasta"
indir = 'dataflow/01-nucl/'
blastdir = 'dataflow/02-blast/'
blastdbdir = 'dataflow/02-blast-db/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='nucl')

file = "rumen_genomes_island2_pathogens.fasta"

# all of the pathogen and selected contigs were then blasted against the selected contigs.
# increase the HSPs in the blast so that regions that span different parts of the islands of interest can be plotted.


file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdir)
outputfilename = "second_island_single_gene_mapping.txt"
file_obj.setOutputName(outputfilename)

blastdb = "subclade_island.fasta"

file_obj.runblast(blast='blastn', db=blastdb, dblocation=blastdbdir, max_target_seqs=10, evalue=1e-3, num_threads = 60, max_hsps = 5)

# then called the ORFS from the selected contigs.

file_obj = sc.Fasta('subclade_island.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('subclade_island.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# created a header mapping file for the selected contigs.

headerfile = 'dataflow/02-headers/'
file = 'subclade_island.fasta'

file_obj = sc.Fasta('subclade_island.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('subclade_island.fasta')
file_obj.setOutputLocation(headerfile)
headers = file_obj.fasta2headermap()
df = pd.DataFrame.from_dict(headers, orient="index")
df['file'] = 'subclade_island.fasta'
df.to_csv(headerfile + file.split('.fa')[0] + '.csv')

# ran prodigal on all selected contigs and the pathogens.

file_obj = sc.Fasta('rumen_genomes_island2_pathogens.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('rumen_genomes_island2_pathogens.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# made a header mapping file for the selected ontigs and pathogens.

headerfile = 'dataflow/02-headers/'
file = 'rumen_genomes_island2_pathogens.fasta'

file_obj = sc.Fasta('rumen_genomes_island2_pathogens.fasta', 'dataflow/01-prot/')
file_obj.setOutputName('rumen_genomes_island2_pathogens.fasta')
file_obj.setOutputLocation(headerfile)
headers = file_obj.fasta2headermap()
df = pd.DataFrame.from_dict(headers, orient="index")
df['file'] = 'rumen_genomes_island2_pathogens.fasta'
df.to_csv(headerfile + file.split('.fa')[0] + '.csv')
