#!/bin/bash
#SBATCH --job-name=sintax_INSERTJOBNAMEHERE
#SBATCH --time=01:00:00
#SBATCH --nodes=1  # specify one node
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=4 
#SBATCH --mem-per-cpu=4G 

#SBATCH -p msismall

#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=INSERTEMAILHERE

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error


module purge
module load usearch/11.0_64bit

usearch -sintax Centroid_mumu_curated.fas -db database/UNITE_ITS_04.04.2024.ITS2_ALL_EUK_no_unidentified_for_sintax.fasta -tabbedout Taxonomy_sintax_ITS2.txt -strand plus -sintax_cutoff 0.8

module purge
module load R/4.2.2-gcc-8.2.0-vp7tyde

R < Combine_OTU_with_tax_ITS2_VSEARCH.R --save
