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

# 3 contigs were selected from rumen_genomes.fasta, based on containing highly conserved homologues in CARD
# subset the selected contigs into a file called rumen_genomes_resistance_genes.fasta

contigs = ['4309689-submission.assembly_79', '4309680-submission.assembly_52', 'RUG782_1']

file_obj = sc.Fasta('rumen_genomes.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('rumen_genomes_resistance_genes.fasta')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.subsetfasta(seqlist = contigs, headertag='resistance_genes')

# run prodigal on rumen_genomes_resistance_genes.fasta

file_obj = sc.Fasta('rumen_genomes_resistance_genes.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('rumen_genomes_resistance_genes.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

# the 3 contigs (Hungate collection 4309689_79 and 43809680_52, MAG RUG782_1) were then blasted against the
# NCBI nucleotide collection (nr/nt) using web-based blastn and the full-length sequence for each of the top 50 hits
# was downloaded and concatenated into resistance_island_blast_hits_concatenated.fasta
# A meta data file with the start and stop of each aligned region, from the online blast, for each unique sequence
# was unsed to extract specific regions into a file called resistance_island_blast_hits_concatenated_extractedCONTIGs.fasta

file = 'resistance_island_blast_hits_concatenated.fasta'

file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
outputfilename = file.split(".f")[0] + '_extractedCONTIGs' + '.fasta'
file_obj.setOutputName(outputfilename)
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.extractORFs_gff3(gff3_table_loc = 'dataflow/00-meta/resistance_blast_hit_cotigs.csv')

# a file was created to later map headers (full header to header truncated after the first white space)

indir = 'dataflow/01-nucl/'
headerfile = 'dataflow/02-headers/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(headerfile)
headers = file_obj.fasta2headermap()
df = pd.DataFrame.from_dict(headers, orient="index")
df['file'] = file
df.to_csv(headerfile + file.split('.fa')[0] + '.csv')

# The concatenated rumen sourced contigs (Hungate collection 4309689_79 and 43809680_52, MAG RUG782_1) were used to
# create a a blast database.

file = "rumen_genomes_resistance_genes.fasta"
indir = 'dataflow/01-nucl/'
blastdbdir = 'dataflow/02-blast-db/'
blastdir = 'dataflow/02-blast/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='nucl')

# The specific regions that were extracted from the online blast hits (resistance_island_blast_hits_concatenated_extractedCONTIGs.fasta)
# were concatenated with the 3 rumen sourced contigs (Hungate collection 4309689_79 and 43809680_52, MAG RUG782_1) and blasted
# against the 3 rumen contigs (rumen_genomes_resistance_genes.fasta).

file = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"

file_obj = sc.Fasta(file, indir)
file_obj.setOutputLocation(blastdir)

outputfilename = "resistance_island_mapping.txt"
blastdb = "rumen_genomes_resistance_genes.fasta"

file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastn', db=blastdb, dblocation=blastdbdir, max_target_seqs=10, evalue=1e-3, num_threads = 60, max_hsps = 5)

# all of the specific regions and contigs were then blasted against eachother, to give a table
# that can be used to figure out which sequenced are unique.


file = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"
indir = 'dataflow/01-nucl/'
blastdbdir = 'dataflow/02-blast-db/'
blastdir = 'dataflow/02-blast/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputName(file)
file_obj.setOutputLocation(blastdbdir)
file_obj.runmakeblastdb(dbtype='nucl')

file = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"

file_obj = sc.Fasta(file, indir)
file_obj.setOutputLocation(blastdir)

outputfilename = "resistance_island_mapping_allvall.txt"
blastdb = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"

file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastn', db=blastdb, dblocation=blastdbdir, max_target_seqs=50, evalue=1e-3, num_threads = 60, max_hsps = 1)

# a table was created to list the unique regions and contigs, which were subsetted into a file,
# translated into ORFs and blasted online to get some first annotations.

genes_df = pd.read_csv('dataflow/00-meta/resistance_blast_hit_cotigs_unique.csv', low_memory=False)
genes = genes_df['qseqid'].tolist()

file_obj = sc.Fasta("resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta", 'dataflow/01-nucl/')
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.setOutputName("resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen_unique.fasta")
file_obj.subsetfasta(seqlist = genes, headertag='unique')

file_obj = sc.Fasta('resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen_unique.fasta', 'dataflow/01-nucl/')
file_obj.setOutputName('resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen_unique.fasta')
file_obj.setOutputLocation('dataflow/01-prot/')
file_obj.runprodigal()

file_obj = sc.Fasta('resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen_unique.fasta', 'dataflow/01-prot/')
file_obj.setOutputLocation('dataflow/03-blast-tables/')
file_obj.runonlineblast()


# extract all of the islands from all the rumen data for blasting against the 3 islands of interest (Hungate collection 4309689_79 and 43809680_52, MAG RUG782_1).
# increase the HSPs in the blast so that regions that span different parts of the islands of interest can be plotted.

file = "rumen_genomes.fasta"
indir = 'dataflow/01-nucl/'
blastdbdir = 'dataflow/02-blast-db/'



file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
outputfilename = file.split(".f")[0] + '_extractedCONTIGs_all_rumen' + '.fasta'
file_obj.setOutputName(outputfilename)
file_obj.setOutputLocation('dataflow/01-nucl/')
file_obj.extractORFs_gff3(gff3_table_loc = 'dataflow/00-meta/resistance_blast_hit_cotigs_all_rumen.csv')

files = ["resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta", "rumen_genomes_extractedCONTIGs_all_rumen.fasta"]

sg.concat(inputfolder='dataflow/01-nucl/', outputpath='dataflow/01-nucl/rumen_genomes_extractedCONTIGs_all.fasta', filenames=files)


file = "rumen_genomes_extractedCONTIGs_all.fasta"
indir = 'dataflow/01-nucl/'
blastdir = 'dataflow/02-blast/'

file_obj = sc.Fasta(file, indir)
file_obj.setOutputLocation(blastdir)

outputfilename = "resistance_island_mapping2.txt"
blastdb = "rumen_genomes_resistance_genes.fasta"

file_obj.setOutputName(outputfilename)
file_obj.runblast(blast='blastn', db=blastdb, dblocation=blastdbdir, max_target_seqs=10, evalue=1e-3, num_threads = 60, max_hsps = 5)


# all prots against all prots from island

# file = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"
#
# file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
# file_obj.setOutputName(file)
# file_obj.setOutputLocation('dataflow/01-prot/')
# file_obj.runprodigal()
#
# file = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"
# indir = 'dataflow/01-nucl/'
# blastdbdir = 'dataflow/02-blast-db/'
#
# file_obj = sc.Fasta(file, 'dataflow/01-prot/')
# file_obj.setOutputName(file)
# file_obj.setOutputLocation(blastdbdir)
# file_obj.runmakeblastdb(dbtype='prot')
#
# blastdir = 'dataflow/02-blast/'
#
# file_obj = sc.Fasta(file, 'dataflow/01-prot/')
# file_obj.setOutputLocation(blastdir)
#
# outputfilename = "resistance_island_all_v_all_prot.txt"
# blastdb = "resistance_island_blast_hits_concatenated_extractedCONTIGs_3rumen.fasta"
#
# file_obj.setOutputName(outputfilename)
# #file_obj.runblast(blast='blastp', db=blastdb, dblocation=blastdbdir, max_target_seqs=100, evalue=1e-3, num_threads = 60, max_hsps = 1)
#
# headerfile = 'dataflow/02-headers/'
#
# file_obj = sc.Fasta(file, 'dataflow/01-prot/')
# file_obj.setOutputName(file)
# file_obj.setOutputLocation(headerfile)
# headers = file_obj.fasta2headermap()
# df = pd.DataFrame.from_dict(headers, orient="index")
# df['file'] = file
# df.to_csv(headerfile + file.split('.fa')[0] + '.csv')
