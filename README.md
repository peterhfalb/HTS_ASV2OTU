# Pipeline for post-processing DADA2 ITS2 ASV output to Denoised OTUs

### Zazzy Metabarcoding Pipeline created by Luis Morgado, University of Oslo
### Updated by Eivind Kverme Ronold, University of Oslo, September 2024
### Updated by Peter Falb, University of Minnesota for compatability at UM


Inside this folder is all the scripts and data necessary to run the pipeline on an ASV table. ASV table included is from the Coryolis project, Kennedy lab. ASV-table used is called: Coryolus_combined_sequences_taxa.txt

**NOTE:** This is a pipeline designed for an HPC cluster located and operated in Norway. Your local architecture and software availability may differ. Check all scripts for software and environments that are loaded, such as Anaconda environments used for calling certain perl/python  scripts. You will also have to set up your own local ITSx environment (unless your HPC cluster has it installed) and download and install mumu locally. 

**NOTE 2:** Original pipeline handles sequences directly from multiplexed libraries from the sequencing centre and demultiplexes, cleans and runs DADA2. This is a truncated version to adapt for post-processing an already generated ASV table.

## Step 1 - Script: 01_Make_sample_fasta_from_ASV_table.R

This is a preliminary script to take the ASV table and turn each sample into a separate fasta file. This is required in order to run the clustering as
it is set up in the pipeline. Make sure that the structure of the ASV table is correct. Each sample needs to be a column and there needs to be a column (or row names) containing each sequence. Also important to doublecheck the map-file that is generated Script sets new names on samples, S[0-9][0-9][0-9]. This naming is required for the pipeline to recognise samples. 

## Step 2 - Script: 02_parallell_ITSx_fungi.txt

First set up a virtual environment on your computing cluster. For UM-MSI, ssh into your account from terminal and run the following commands to set up the environment:

```bash
# Load conda (adjust module name to match your cluster)
module load conda

# Create environment and install ITSx from bioconda
conda create -n itsx_env -c bioconda -c conda-forge itsx

# Activate it
conda activate itsx_env

# Load VSEARCH as a module on top of the conda environment
module load vsearch

# Verify both work
ITSx --help
vsearch --version
```


This is formulated as a shell script for queueing in an HPC cluster (slurm) 
Adapt the queueing parameters to whatever system you are running
Script requires two dependencies, VSEARCH and ITSx. 
Adapt the activation of these to your own environment
**Method:** Sequences are split into batches using VSEARCH  within each sample for parallell processing with ITSx. This reduces computation time. Still, expect a runtime of many hours for data of 100 or more separate samples. This is the process that takes the longest time in the pipeline. In essence, ITSx identifies and removes the highly conserved 5.8S part of the ITS2 barcodes. This makes for better clustering later, as the 5.8S regions inflates the similarity between sequences. 

## STEP 3 - Script: 03_VSEARCH_cluster.txt

This script does many things. Most importantly, clusters the ASVs into OTUs at a custom threshold (97% default). 
Also formulated as a slurm script, adapt as needed. Supports multithreading for speed.
This script depends on two background scripts, rename.pl and uc2otutab_original.py. 
Make sure that the paths to background scripts are set up properly. These are perl and python scripts to organise sample names and generate an OTU table from the clustered centroids
Do not change anything in the background scripts.
Very important to follow the S[0-9][0-9][0-9] naming convention unless you want to do a lot of custom regex expressions in all the dependencies.
This script also includes the VSEARCH denovo chimera-check. This usually catches chimeras that leak through DADA2. It is basically identical to UCHIME-denovo from USEARCH.
Outputs both a centroid and an OTU-table. Remember to set the name for your project.
Usually completes in a few minutes on an HPC cluster with default settings.

## STEP 4 - Script: 04_mumu_curation_VSEARCH.txt

First, you need to clone mumu on your system. On terminal navigate to the directory where you want the package. Then clone the repository:
```bash
cd ~
mkdir packages # ONLY if this folder doesn't exist yet.
cd packages

ml gcc/13.1.0-5z64cho
git clone https://github.com/frederic-mahe/mumu.git
cd ./mumu/
make
make check
#make install  # SKIP THIS unless you have sudo permissions on your device. You will just reference the binary in the SLURM call.
```
*mumu* is a post-clustering clean up algorithm that is supposed to find rare variantssequences of common sequences that leak through the clustering steps.
It works in two steps. 
*First*, all sequences are blasted against each other. Sequences with high similarity are checked for patterns in the OTU table. If a rare sequence with high similarity to a common sequence also show a very similar occurrence pattern, is is merged with this "parent" sequence. Mostly used to clean up singleton and doubleton sequences. 
*mums* is a unix version of lulu, with some optimisation to run faster. In the script you are required to set the path to your folder with mumu installed. Also set the name of the input files to exactly match the output from VSEARCH. This script outputs a new centroid and a new OTU table. VSEARCH outputs are not deleted in case you want to compare. 

## STEP 5 - Script: 05_Taxonomic_annotation_SINTAX_ITS2.txt

SINTAX is a taxonomic assignment method developed by R. Edgar for USEARCH.
It is also implemented in VSEARCH so possible to load from there is you do not have a USEARCH license. As the script is set up, there is no need to change any 
filenames. It works directly on the mumu output files. Database needs to have a specific structure to the fasta headers, included in the "database" folder is the most
recent UNITE ITS2 database, formatted for optimalization with SINTAX. This database has had all sequences with an "unidentified" annotation anywhere in the header removed.
In addition there is an R script that will attach the taxonomy to the OTU table from mums.
If you want to run that, make sure to set the correct name for the mumu OTU-table that is called. 


### Fin.
