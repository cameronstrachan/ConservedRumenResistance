# ConservedRumenResistance

## Overview

These are the scripts used for carrying out the analysis of highly conserved antibiotic 
resistance genes, described in the manuscript *The role of the rumen in generating 
specific antibiotic resistances* (Strachan et al. 2020). This is intended to allow
those interested to look up specific parameters and flags used and replicate any aspects
of the analysis.

# Description

Much of the basic fasta manipulation, ORF calling and blasting was done using python 
driver scripts, which call call functions from seq_core_lin.py and seq_gen_lin.py 
(modules). In these modules, you can find all the default parameters used, if not set
in the driver scripts. The driver scripts create an object each time, which is simply 
defined by the file input and file output location (this is always called file_obj in the 
scripts). The functions then act on the specified input file and create the specified 
output file. This was intended to make it easier to follow the steps carried out. 
  

## Flow of analysis

The driver scripts are named P1 - P10 and go in order of the analysis carried out in the
paper. The file types flow through the numbered folders in dataflow.

Fig 1 - P1/P2  
Supp Fig 1 - P3  
Fig 2 - P4  
Fig 3 - P5  
Supp Fig 3- P6  
Fig 4 - P7/P8/P9  
Fig 5 - P7/P8/P9  
Supp Fig 4 - P10  


## Intermediate tables and figures

Certain intermediate tables and figures were generated using the R script found in the
folder R or with the Rmarkdown notebooks found in markdown_notebooks. 

## Data downloading

Manual downloading of the data is described in the manuscript, but when are script was
used it is found in the folder data_downloading.