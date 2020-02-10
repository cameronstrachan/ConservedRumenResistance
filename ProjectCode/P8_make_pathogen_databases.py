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

# concatenate all of the downloaded pathogen assemblies and then make blast DBs for each

dirs = ['staphylococcus_aureus', 'campylobacter_jejuni', 'campylobacter_coli', 'clostridioides_difficile', 'acinetobacter_baumannii', 'streptococcus_pneumoniae']
head_dir = 'dataflow/01-nucl/'

for dir in dirs:
    path_dir = head_dir + dir + '/'
    unzip_command = 'gunzip ' + path_dir + '*.gz'
    os.system(unzip_command)
    lis = [f for f in os.listdir(path_dir) if f.endswith(".fna")]
    output_file = head_dir + dir + '.fasta'
    sg.concat(inputfolder=path_dir, outputpath=output_file, filenames=lis)

files = ['staphylococcus_aureus', 'campylobacter_jejuni', 'campylobacter_coli', 'clostridioides_difficile', 'acinetobacter_baumannii', 'streptococcus_pneumoniae']

for file in files:
    file_obj = sc.Fasta(file, 'dataflow/01-nucl/')
    file_obj.setOutputName(file)
    file_obj.setOutputLocation('dataflow/02-blast-db/')
    file_obj.runmakeblastdb(dbtype='nucl')
