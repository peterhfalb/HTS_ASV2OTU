# SLURM-based pipeline for post-processing DADA2 ASV tables to denoised OTUs

### Pipeline Overview and Authors

Full credit for the development of this pipeline goes to Eivind Ronold and Luis Morgado. I (PF) have updated and adapted all the scripts to work on University of Minnesota systems, as well as implemented various bug fixes and speed/code improvements. Individual scripts were consolidated into one pipeline script which compiles and installs all necessary packages for running the pipeline. The pipeline was also adapted to handle all types of sequencing data currently used in the Kennedy Lab (16S-V4, ITS1, ITS2, 18S-V4, 18S-AMF). 

The original Zazzy pipeline created by Luis Morgado was designed to handle sequences directly from multiplexed libraries from sequencing centers, then clean an run DADA2. This is a truncated version adapted to take pre-DADA2 processed ASV files outputted from Trevor Gould's MSI pipeline. 

**Zazzy Metabarcoding Pipeline created by Luis Morgado, University of Oslo**

**Updated by Eivind Kverme Ronold, University of Oslo, September 2024**

**Updated by Peter Falb, University of Minnesota for compatability at UM**

## Setup Instructions:

Login in to MSI via terminal on your mac (usage might be different on Windows, see https://msi.umn.edu/connecting/connecting-to-hpc-resources):

```bash
ssh -Y yourUMNusername@agate.msi.umn.edu #Fill in with your UMN username
```

Next, move into a personal directory where you would like the pipeline files to be installed (this is a permanent, single location you will use anytime you run the pipeline):

```bash
#For example:
cd Packages/KennedyLabPipelines/
```

Once in your preferred directory, clone the github repository:
```bash
git clone https://github.com/peterhfalb/HTS_ASV2OTU.git
```

To set up the pipeline script and packages move into the cloned repository and run the setup script:
```bash
# move into the cloned repository
cd HTS_ASV2OTU/ 

# run this line to update the pipeline script and
# answer the prompt to enter your email for SLURM notifications
bash setup.sh 
```
This script does 3 things, *FIRST* it updates the path to reflect where the pipeline code is located on the cluster, *SECOND* it prompts you for your email, which will be used to send SLURM notifications when you run the script, *THIRD* it installs all the packages needed to run the pipeline, and *FOURTH* it checks to make sure the taxonomy database files are correctly located within the Kennedy Lab shared directory. Package installation happens now because it is the most likely step in the pipeline where errors are going to occur. I have tried to set it up so that dependencies are properly handled via the MSI infrastructure, but if anything goes wrong, please screenshot the error and contact me (Peter Falb; falb0011@umn.edu). When you run the actual pipeline script, it will check again to make sure all packages are installed, and attempt to install them if they aren't.

**NOTE: It may take 10-15 minutes to install the necessary packages**

After proper installation, the ASV to OTU pipeline should now be ready to go.

## Usage Instructions

This script takes one input file - the ASV table outputted from Trevor Gould's DADA2 pipeline. This should be a table that has the ASVs as column one, sample abundances in the next columns, then taxonomy and taxonomy bootstrap values as final columns. To run the ASV to OTU pipeline, you will want to create a new project directory on the cluster, and then upload your ASV table. Then you will run the pipeline script.

**Step 1: Login to the cluster, and prepare files**
```bash
# login to the cluster (fill in your UMN username)
ssh -Y yourMSIusername@agate.msi.umn.edu 

# navigate to where you want your project to live
cd yourPreferredDirectory/ 

# create a directory for your project
mkdir yourProjectName_ASVtoOTU 

# navigate into the new directory
cd yourProjectName_ASVtoOTU 

# print and copy this path to the file directory
pwd 

# leave the ssh window
exit 

# copy your ASV table from your local computer to your new project directory on the cluster
# (fill in file paths and your UMN username)
scp /path/to/your/ASVTABLE.txt yourMSIusername@agate.msi.umn.edu:/path/to/your/yourProjectName_ASVtoOTU/ 

```

**Step 2: Run the pipeline**
To run the script, you need to ssh login to the cluster and go to the directory where the pipeline script is housed. Then, you will have to run the script with a SLURM submission, specifying 4-5 things:
1. Path to project directory
2. Path to ASV table
3. Project Name (this will determine how your output files are named)
4. Primer Set / Barcoding Region
5. Whether to run ITSx (*only applicable for ITS1 or ITS2 regions*) -- See note below

*ITSx is a program which removes the highly conserved regions flanking the ITS variable region. ITS primers often pick up these conserved regions, which can artificially inflate the sequence similarity (clustering) between different ITS ASVs. The pipeline will automatically run ITSx for ITS primer sets, but if you prefer to NOT run this step, add the flag --skip_itsx*

**Navigate to directory of script then run the command**
```bash
cd path/to/pipelineCode/directory/
sbatch ASVtoOTU_msiSLURM.sh <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx]

# PRIMER SET OPTIONS:
#   ITS1       — fungal ITS1 region (UNITE database, ITSx optional)
#   ITS2       — fungal ITS2 region (UNITE database, ITSx optional)
#   16S-V4     — bacterial 16S V4 region (SILVA database, no ITSx)
#   18S-V4     — microeukaryote 18S V4 region (PR2 database, no ITSx)
#   18S-AMF    — arbuscular mycorrhizal fungi 18S (MaarjAM database, no ITSx)

# EXAMPLES:
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 ITS2
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 ITS2 --skip-itsx
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 16S-V4
```

**A note about SBATCH/SLURM scripts if you are unfamiliar:** when you run the sbatch command, it submits the script as a 'SLURM' submission to the computing cluster. This means the 'job' will get in a queue to eventually run. Depending on what time of day/week you submit it, it could take anywhere from 2 seconds to 30 minutes to initiate (usually towards the lower end in my experience). Once it starts, you will get an email saying it started. Next, you will either get an email that the script COMPLETED or FAILED. If it FAILED, it will specify an Exit Code number, usually 2. If it failed, contact me (Peter F) and I can help troubleshoot. Every time you run the pipeline, the job will output a .out and .err file within the directory where the pipeline is stored. You don't really need to worry about these files, UNLESS something fails, then both will be helpful for understanding what error was thrown/what went wrong.

## Output Files:

The ASV to OTU pipeline will output a LOT of files (see below for descriptions of all), but you will likely be most interested in two:

1. *ProjectName_OTU_with_taxonomy_PrimerSet.txt* -  This is your new data table, functionally equivalent to the input ASV table, but with OTU-level abundances with new taxonomic assignments
2. *pipeline_run.log* - This is a log file which has logged a. everything the pipeline did and b. what each of the steps removed. You will want to look at the "Pipeline Summary - Quality Control" section which tells you how many ASVs/OTUs were removed in each step of the process. 

The *pipeline_run.log* also summarizes other files you may be interested in.

## Summary of main pipeline steps:

### Step 0 - ITSx (only relevant for ITS datasets)
 In essence, ITSx identifies and removes the highly conserved 5.8S part of the ITS2 barcodes. This makes for better clustering later, as the 5.8S regions inflates the similarity between sequences.

### Step 1 - VSEARCH Clustering

This step uses the VSEARCH algorithm to cluster ASVs at 97% similarity (weighted by sequence abundance; *see VSEARCH documentation*), then does an additional chimera check to find and remove any chimeras that might have slipped through DADA2. The number removed here should be quite small (check log file).

### Step 2 - MUMU/LULU OTU Curation

*mumu* is a post-clustering clean up algorithm that is supposed to find rare variants sequences of common sequences that leak through the clustering steps.
It works in two steps. All sequences are blasted against each other. Sequences with high similarity are checked for patterns in the OTU table. If a rare sequence with high similarity to a common sequence also show a very similar occurrence pattern, is is merged with this "parent" sequence. These sequences are assumed to be additional sequencing errors that have made it through all previous cleaning steps. Mostly used to clean up singleton and doubleton sequences. *mumu* is a unix version of lulu, with some optimisation to run faster.

### Step 3 - DADA2 RDP Taxonomy Assignment

After the OTUs have been created and curated, taxonomy is assigned to the centroid of each OTU using the exact same approach used by Trevor Gould in his DADA2 pipeline. This uses the *assignTaxonomy* function in DADA2, which uses the RDP Naive Bayesian classifier (doi: 10.1128/AEM.00062-07) to classify sequences. Taxonomy databases are unique to each primer (ITS1 and ITS2 have one database, UNITE sh) and are the same databases used by Trevor Gould to assign taxonomy. I have updated each taxonomy database to their most recent version, and they are located in the Kennedy lab shared taxonomy folder on the MSI computing cluster.

## Citations:

If you use this, please give credit to Luis Morgado for developing the original Zazzy Metabarcoding Pipeline (https://github.com/ekronold/Zazzy_metabarcoding_pipeline), and Eivind Kverme Ronold for making the initial updates for adaptation to DADA2 output ASV files.

**Additionally** please cite the primary packages and databases used in the pipeline:

#### Packages:

ITSx: Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for use in environmental sequencing. Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco, Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez, Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov, K. Martin Eriksson, R. Henrik Nilsson *Methods in Ecology and Evolution, 4: 914-919, 2013* (DOI: 10.1111/2041-210X.12073)

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Frédéric Mahé. (2026) mumu: post-clustering curation tool for metabarcoding data v.1.1.2 (https://github.com/frederic-mahe/mumu)

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188. (https://doi.org/10.1038/s41467-017-01312-x)

#### Databases

González-Miguéns, R., Gàlvez-Morante, À., Skamnelou, M., Antó, M., Casacuberta, E., Richter, D.J., Lara, E. et al. 2025. A novel taxonomic database for eukaryotic mitochondrial cytochrome oxidase subunit I gene (eKOI),with a focus on protists diversity. Database (Oxford). 2025:baaf057.

Chuvochina M, Gerken J, Frentrup M, Sandikci Y, Goldmann R, Freese HM, Göker M, Sikorski J, Yarza P, Quast C, Peplies J, Glöckner FO, Reimer LC (2026) SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Nucleic Acids Research, gkaf1247.

Abarenkov K, Nilsson RH, Larsson K-H, Taylor AFS, May TW, Frøslev TG, Pawlowska J, Lindahl B, Põldmaa K, Truong C, Vu D, Hosoya T, Niskanen T, Piirmann T, Ivanov F, Zirk A, Peterson M, Cheeke TE, Ishigami Y, Jansson AT, Jeppesen TS, Kristiansson E, Mikryukov V, Miller JT, Oono R, Ossandon FJ, Paupério J, Saar I, Schigel D, Suija A, Tedersoo L, Kõljalg U. 2024. The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Research, https://doi.org/10.1093/nar/gkad1039

Öpik, M., Vanatoa, A., Vanatoa, E., Moora, M., Davison, J., Kalwij, J.M., Reier, Ü., Zobel, M. 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). New Phytologist 188: 223-241.
