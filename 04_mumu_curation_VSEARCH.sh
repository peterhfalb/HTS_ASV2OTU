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


# Use the centroid file with this extension from the SWARM output as it contains the swarm centroids
cp [INSERT PROJECT NAME].centroids OTU_centroids


# Fix the headers to only contain the sequence ID from dada2 and nothing else
sed -i 's/;.*//' OTU_centroids


# Use the centroid file to make a blast database for internal matching
makeblastdb -in OTU_centroids -parse_seqids -dbtype nucl


# Blast all sequences against eachother to create the match.list which is one of the parameters for mumu-curation 
blastn -db OTU_centroids -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -num_threads 16 -query OTU_centroids

# Clean up the amplicon names before mumu 

cp [INSERT PROJECT NAME].otutable pre_mumu_table.txt

awk 'BEGIN{FS=OFS="\t"} {sub(/;.*/, "", $1); print}' pre_mumu_table.txt > mumu_table.txt

rm pre_mumu_table.txt

# This runs the mumu curation, output file will be matched with the data from SWARM output
~/packages/mumu/mumu \
		--otu_table mumu_table.txt \
		--match_list match_list.txt \
		--log mumu.log \
		--new_otu_table [INSERT PROJECT NAME]_mumu_curated.txt


# Generate a new centroid file from the mumu_curated table for taxonomic annotation
grep -A 1 -Ff <(awk '{print $1}' [INSERT PROJECT NAME]_mumu_curated.txt) OTU_centroids | sed '/--/d' > Centroid_mumu_curated.fas

