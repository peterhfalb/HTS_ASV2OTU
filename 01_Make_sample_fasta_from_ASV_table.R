# Script for turning an ASV table into sample-fasta files and a centroid for running Zazzy post-processing of ASV-tables

# Edited by Peter Falb on March 2, 2026 for input of phyloseq data format, and generalization of the script to be more flexible for different datasets.

library(tidyverse)
library(phyloseq)

# set wd
setwd("/Users/peterfalb/Documents/Research/Kennedy/FAB2_2022/ASVtoOTU/ITS/Clean_annotated_scripts")

# load the phyloseq data
ps <- readRDS("/Users/peterfalb/Documents/Research/Kennedy/FAB2_2022/data/phyloseqDiversity/FAB2.phylo.ITS.rare.div.rds")

# Set number of samples (+ 1 to catch the ASV column)
n_samples <- nsamples(ps) + 1

# Extract ASV table
asv_table <- data.frame(otu_table(ps), check.names = FALSE)

# Extract Sequence_ID from tax_table and add as first column named "ASV"
sequence_ids <- data.frame(tax_table(ps))[, "Sequence_ID", drop = FALSE]

asv_table <- cbind(ASV = sequence_ids$Sequence_ID, asv_table)

# Add a column for the total reads per ASV across all samples, make sure to catch all the sample columns
asv_table <- asv_table %>% 
  mutate(Total = rowSums(asv_table[, 2:n_samples])) %>% 
  arrange(desc(Total)) # Arrange be descending order of read counts

# Print the names just to check for sample columns.
names(asv_table)

# Extract the actual sequences in a vector
Seqs <-asv_table$ASV

# Set up a mapping file to make sure the names used through the bioinformatic clustering can be re-mapped to real sample names
sample_names <- names(asv_table[,2:n_samples]) # List of the sample names
# Make a column with names S001 - S090 to match the sample names
# This sample naming is very important as the Zazzy post-processing has many scripts that call samples by a S[0-9][0-9][0-9] regex. 
Map_file <- data.frame(
  "BioInf_Sample" = paste0("S", str_pad(sample_names, width = 3, side = "left", pad = "0")),
  "Sample_name" = sample_names)

# Check the sample names:
Map_file$BioInf_Sample

# Save the map file IMPORTANT :)
write.csv(Map_file, "Map_file.csv", row.names = F, quote = F)

# Rename the asv_table
# Remove the taxonomy column which won't be used
asv_table <- asv_table[,2:n_samples]
#Rename the columns
colnames(asv_table) <- Map_file$BioInf_Sample

# Extract the sequences, use fasta-format header style
Sequences <- data.frame("SeqID" = paste0(">", rownames(asv_table)),
                        "Sequence" = Seqs)
# Make all the seqids same length, this just makes all the tables look nicer
### NOT NECESSARY HERE ALSO NOT APPLICABLE
#Sequences$SeqID[1:9] <- str_replace(Sequences$SeqID[1:9], ">Seq", ">Seq000")
#Sequences$SeqID[10:99] <- str_replace(Sequences$SeqID[10:99], ">Seq", ">Seq00")
#Sequences$SeqID[100:999] <- str_replace(Sequences$SeqID[100:999], ">Seq", ">Seq0")

# Turn into a centroid file, this will be the pre-processing centroid
Centroid.vec <- data.frame(as.vector(t(cbind(Sequences$SeqID, Sequences$Sequence)))) 

# Write out the centroid
write.table(Centroid.vec, "Centroid.fasta", col.names = F, row.names = F, quote = F)

# Define a function to write a fasta file for each sample
# This is required for the current setup of Zazzy. 
make.fasta <- function(x, n){
  S <- Sequences %>% add_column("Reads" = x[,n])
  SS <- S %>% filter(Reads > 0)
  SSS <- SS %>% mutate(SeqID = str_replace(SeqID, SeqID, paste0(SeqID, ";size=", Reads))) %>% 
    arrange(desc(Reads)) %>% 
    select(c("SeqID", "Sequence"))
  
  # Write fasta
  Sample.vec <- data.frame(as.vector(t(cbind(SSS$SeqID, SSS$Sequence))))
  
  return(Sample.vec)
}

# Number of repetitions
n_reps <- ncol(asv_table)

# Initialize an empty list to store the combined columns
combined_list <- list()

# Loop to create the combined columns and name them S[number]
for (i in 1:n_reps) {
  combined <- make.fasta(asv_table, i)
  combined_list[[paste0(colnames(asv_table)[i])]] <- combined
}

# Write each sample file to a separate fasta
for (name in names(combined_list)) {
  file_name <- paste0(name, ".fas")
  write.table(combined_list[[name]], file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# There should now be a .fas file for each sample in the working folder. Each header should contain information about how many times the specific sequence occur in that specific sample. This can now be run directly through the Zazzy post-processing of ASVs

