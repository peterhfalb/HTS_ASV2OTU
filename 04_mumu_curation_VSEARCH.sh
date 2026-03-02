#!/bin/bash
#SBATCH --job-name=mumu_INSERTJOBNAMEHERE
#SBATCH --time=00:30:00
#SBATCH --nodes=1  # specify one node
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=5G 

#SBATCH -p msismall

#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=INSERTEMAILHERE

# mumu curation script. Make sure to set proper names on input and output files
# source: https://github.com/frederic-mahe/mumu

module purge
ml blast-plus/2.13.0-gcc-8.2.0-vo4mr4d

cd /YOURWORKINGDIRECTORY/

# SETTING VARIABLES - this NEEDS to be the same as script 2
PROJ="COOL_PROJECT_NAME"

# Use the centroid file from the VSEARCH output
cp ${PROJ}.centroids OTU_centroids

# Fix the headers to only contain the sequence ID and nothing else
sed -i 's/;.*//' OTU_centroids

# Use the centroid file to make a blast database for internal matching
makeblastdb -in OTU_centroids -parse_seqids -dbtype nucl

# Blast all sequences against eachother to create the match.list
blastn -db OTU_centroids \
       -outfmt '6 qseqid sseqid pident' \
       -out match_list.txt \
       -qcov_hsp_perc 80 \
       -perc_identity 84 \
       -num_threads 16 \
       -query OTU_centroids

# Clean up the amplicon names before mumu
cp ${PROJ}.otutable pre_mumu_table.txt

awk 'BEGIN{FS=OFS="\t"} {sub(/;.*/, "", $1); print}' pre_mumu_table.txt > mumu_table.txt

rm pre_mumu_table.txt

# Run mumu curation
~/packages/mumu/mumu \
        --otu_table mumu_table.txt \
        --match_list match_list.txt \
        --log mumu.log \
        --new_otu_table ${PROJ}_mumu_curated.txt

# Generate a new centroid file from the mumu_curated table for taxonomic annotation
grep -A 1 -Ff <(awk 'NR>1 {print $1}' ${PROJ}_mumu_curated.txt) OTU_centroids | sed '/--/d' > Centroid_mumu_curated.fas