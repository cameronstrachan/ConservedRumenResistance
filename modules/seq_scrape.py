#!/usr/bin/env python3

from Bio import Entrez
import subprocess

'''
This is simply some functions involved with scraping
genomic data
'''

def ncbigenomescrape(searchterm, searchterm2='assembly[title]', location='none'):

	Entrez.email = "strachc@gmail.com"

	searchterm = str(searchterm)
	query = searchterm + '[orgn]' + ' AND ' + searchterm2

	handle = Entrez.esearch(db="nucleotide", term=query, retmax=5000)
	searchrecord = Entrez.read(handle)
	gen_ids = searchrecord["IdList"]

	num_genomes = len(gen_ids)

	print('Number of genomes: %1d' % num_genomes)

	if num_genomes == 0:
		print('No genomes downloaded')
		return('No genomes downloaded')
		cont = 'n'
#	elif num_genomes > 2:
#		print('Too many genomes available')
#		return('No genomes downloaded')
	else:
		for gen_id in gen_ids:
			handle = Entrez.efetch(db="nucleotide", id=gen_id, rettype="gb", retmode="text")
			genbankfile = handle.read()
			outputfile = open(location + str(searchterm) + '_' + gen_id + '_genome' + '.gb', 'w')
			outputfile.write(genbankfile)
			outputfile.close()
		return('Genomes downloaded')
		print('Done')


def srafastqdownlaod(accession, outputdir='dataflow/01-fastq'):

	command = '../bin/fastq-dump.2.9.2 ' + '--outdir ' + outputdir + ' --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip  ' + accession
	process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()

# Testing program
def main():
	pass

if __name__ == '__main__':
    main()
