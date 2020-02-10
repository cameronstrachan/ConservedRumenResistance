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

blastdb = 'dataflow/02-blast-db/'
blastout = 'dataflow/02-blast-out/'
blastxml = 'dataflow/02-blast-xml/'
blastin = 'dataflow/01-nucl/'
prot_dir = 'dataflow/01-prot/'
analysis_folder = 'dataflow/03-analysis/'
dbs = ['campylobacter_coli.fasta', 'campylobacter_jejuni.fasta', 'staphylococcus_aureus.fasta', 'clostridioides_difficile.fasta', 'acinetobacter_baumannii.fasta', 'streptococcus_pneumoniae.fasta']
output_files_blast = []

# step 1 - trim rumen genes file those those above 300bp
file_obj = sc.Fasta('stewart2019_metagenomes_genes.fasta', blastin)
file_obj.setOutputName('stewart2019_metagenomes_genes_1000.fasta')
file_obj.setOutputLocation(blastin)
file_obj.lengthcutoff(replaceheaders = False, length = 1000, direction = 'above')

input_file = 'stewart2019_metagenomes_genes_1000.fasta'

# step 2 - loop through the different pathogen genome database files, blastn
for db_file in dbs:
    file_obj = sc.Fasta(input_file, blastin)
    outname = input_file.split('.fa')[0] + '_' + db_file.split('.fa')[0] + '.txt'
    file_obj.setOutputName(outname)
    file_obj.setOutputLocation(blastout)
    #file_obj.runblast(blast='blastn', db=db_file, dblocation=blastdb, max_target_seqs=5000, evalue=1e-100, num_threads = 60, max_hsps = 1)
    output_files_blast.append(outname)

# step 3 - summarize frequency of hits
command = 'Rscript summarize_hit_frequency.R '
output_files_freq = []

for file in output_files_blast:
    command = command + file + ' '
    outname = file.split('.tx')[0] + '.csv'
    output_files_freq.append(outname)

#os.system(command)

# step 4 - extract genes from the frequency tables for annotation and convert to ORFS
genes = []

for file in output_files_freq:
    csv_file = analysis_folder + file
    df = pd.read_csv(csv_file, low_memory=False)
    gene_ids = df['gene'].tolist()
    genes = genes + gene_ids

genes_unique = list(set(genes))

file_obj = sc.Fasta(input_file, blastin)
outname = input_file.split('.fa')[0] + '_pathogen_mapped.fasta'
file_obj.setOutputName(outname)
file_obj.setOutputLocation(blastin)
file_obj.subsetfasta(seqlist = genes_unique, headertag='none')

file_obj = sc.Fasta(outname, blastin)
file_obj.setOutputName(outname)
file_obj.setOutputLocation(prot_dir)
file_obj.translateORFs()

# step 5 - annotate the ORFs
#file_obj = sc.Fasta(outname, prot_dir)
#file_obj.setOutputLocation(blastxml)
#file_obj.runonlineblast(numhits=1)

#blastfiles = [outname]
#xmlfiles = [f for f in os.listdir(blastxml) if f.endswith(".xml")]
#sg.blastxmltotable(xmlinputfolder=blastxml, blastinputfolder=blastin,outputpath='dataflow/03-analysis/compiled_annotations.txt', xmlfilenames=xmlfiles, blastfilename=blastfiles)

# step 6 - make header map file to extract locations of genes

file_obj = sc.Fasta(input_file, blastin)
headers = file_obj.fasta2headermap()

analysis_folder = 'dataflow/03-analysis/'
df = pd.DataFrame.from_dict(headers, orient="index")
df['file'] = input_file

df.index.name = 'id'
df.reset_index(inplace=True)

df.columns = ['id', 'full_header', 'file']

df['id_unnumbered'] = 'NA'

df['id_unnumbered'] = df.apply(lambda x: str(x['id']).rsplit('_', 1)[0], axis=1)

df = df[df['id'].isin(genes_unique)]

df.to_csv(analysis_folder + input_file.split('.fa')[0] + '_mapped_headers.csv')

# step 7 - annotate the mapped gene products to the CARD database

file_obj = sc.Fasta(outname, prot_dir)
file_obj.setOutputLocation(blastout)

outputfilename = outname.split('.fa')[0] + '_card.txt'
db_file = "card_db.fasta"

file_obj.setOutputName(outputfilename)

file_obj.runblast(blast='blastp', db=db_file, dblocation=blastdb, max_target_seqs=1, evalue=1e-3, num_threads = 60)
