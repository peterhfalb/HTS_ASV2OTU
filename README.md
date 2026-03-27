# SLURM pipeline to convert DADA2 ASV tables to denoised OTUs

### Pipeline Overview and Authors

This pipeline script was written specifically for use by the Kennedy Lab at the University of Minnesota. It takes DADA2 ASV tables outputs from Trevor Gould's pipeline, and converts them to denoised OTUs. The goal here is to further reduce sequencing noise and errors coming off the sequencing platform, and will likely be most helpful for cleaning non-host associated bulk communities with extremely high sequence diversity. OTUs are clustered at 97% sequence identity, which very approximately represents the level of species. The pipeline is specifically designed for the Minnesota Agate computing cluster, and requires user access to the Kennedy lab shared folders for use. If you would like to run this on a different computing cluster, or don't have access to the Kennedy Lab folders, please contact me (Peter F; falb0011@umn.edu). 

Full credit for the development of this pipeline goes to Eivind Ronold and Luis Morgado. I (PF) have updated and adapted all the scripts to work on University of Minnesota systems, as well as implemented various bug fixes and speed/code improvements. Individual scripts were consolidated into one pipeline script which compiles and installs all necessary packages for running the pipeline. The pipeline was also adapted to handle all types of sequencing data currently used in the Kennedy Lab (16S-V4, ITS1, ITS2, 18S-V4, 18S-AMF). 

The original Zazzy pipeline created by Luis Morgado was designed to handle sequences directly from multiplexed libraries from sequencing centers, then clean and run DADA2. This is a truncated version adapted to take pre-DADA2 processed ASV files outputted from Trevor Gould's MSI pipeline. 

**Zazzy Metabarcoding Pipeline created by Luis Morgado, University of Oslo**

**Updated by Eivind Kverme Ronold, University of Oslo, September 2024**

**Updated by Peter Falb, University of Minnesota, March 2026**

## Installation Instructions:

Login in to MSI via terminal on your mac (usage might be different on Windows, see https://msi.umn.edu/connecting/connecting-to-hpc-resources):

```bash
# login to MSI via SSH secure shell
ssh -Y yourUMNusername@agate.msi.umn.edu # change 'yourUMNusername' to your UMN username
```

Next, move into a personal directory where you would like the pipeline files to be installed, or create a new directory to house your pipelines. This is a permanent, single location you will use anytime you run the pipeline, and is personal to you (it won't work for other users).

```bash
# [optional]: make a directory for pipelines/packages (if desired/needed)
mkdir pipelines # mkdir is "make directory"

# move into the directory where you want the pipeline installed
cd pipelines/ # cd is "change directory"
```

Once in your preferred directory, clone the HTS_ASV2OTU github repository:
```bash
git clone https://github.com/peterhfalb/HTS_ASV2OTU.git
```

To set up and install the pipeline, open the cloned repository and run the setup/installation script:
```bash
# move into the cloned repository
cd HTS_ASV2OTU/ 

# run this line of code to install and compile the pipeline package
bash setup.sh
# then type into the prompt asking for your email, to setup your email for SLURM notifications
# pipeline installation could take 10-15 minutes, but is only necessary the first time you run the pipeline
```
This script does 3 things, *FIRST* it updates the path to reflect where the pipeline code is located on the cluster, *SECOND* it prompts you for your email, which will be used to send SLURM notifications when you run the script, *THIRD* it installs all the packages needed to run the pipeline, and *FOURTH* it checks to make sure the taxonomy database files are correctly located within the Kennedy Lab shared directory. Package installation happens now because it is the most likely step in the pipeline where errors are going to occur. I have tried to set it up so that dependencies are properly handled via the MSI infrastructure, but if anything goes wrong, please screenshot the error and contact me (Peter Falb; falb0011@umn.edu). When you run the actual pipeline script, it will check again to make sure all packages are installed, and attempt to install them if they aren't.

**NOTE: It may take 10-15 minutes to install the necessary packages**

After proper installation, the ASV to OTU pipeline should now be ready to go.

## Usage Instructions:

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

# print and copy the path to the file directory
pwd # "print working directory" - copy what this prints

# leave the ssh window
exit 

# copy your ASV table from your local computer to your new project directory on the cluster
# (fill in file paths and your UMN username)
scp /path/to/your/ASVTABLE.txt yourMSIusername@agate.msi.umn.edu:/path/to/your/yourProjectName_ASVtoOTU/ 
# scp is 'secure copy', and allows you to copy files from your local device directly to the cluster

```

**Step 2: Run the pipeline**

To run the script, you need to ssh login to the cluster and go to the directory where the pipeline script is housed. Then, you will have to run the script with a SLURM submission, specifying 4-6 things:
1. Path to project directory
2. Path to ASV table
3. Project Name (this will determine how your output files are named)
4. Primer Set / Barcoding Region
5. Whether to run ITSx (*only applicable for ITS1 or ITS2 regions; SEE NOTE BELOW*)
6. Manually pick reference database (if you would like to specify a non-default reference database for your given group)

Expect a runtime of 10-45 minutes. It should not go much longer than that, but the SLURM script requests 2 hours of time on the cluster just in case.

*ITSx is a program which removes the highly conserved regions flanking the ITS variable region. ITS primers often pick up these conserved regions, which can artificially inflate the sequence similarity (clustering) between different ITS ASVs. The pipeline will automatically run ITSx for ITS primer sets, but if you prefer to NOT run this step, add the flag --skip_itsx*

*NOTE: recent testing has shown that ITSx seems to remove the synmock ASVs used by the Kennedy lab for the positive control. If using an ITS dataset with a synmock community, you may want to skip the ITSx step*

*Additionally, please avoid using taxonomic ranking names (e.g. Class, Order, etc.) as sample names in the input ASV table, as this will screw up the column filtering code*

**Navigate to directory of script then run the command**
```bash
# login to the cluster (fill in your UMN username)
ssh -Y yourMSIusername@agate.msi.umn.edu 

# navigate to the directory where you cloned and set up the pipeline scripts (HTS_ASV2OTU folder)
cd pipelines/HTS_ASV2OTU/ # you can also use cd .. to go up a directory

# submit the slurm job with the following commands:
sbatch ASVtoOTU_msiSLURM.sh <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx] [--db <database>]

# PRIMER SET OPTIONS:
#   ITS1       — fungal ITS1 region (UNITE database, ITSx optional)
#   ITS2       — fungal ITS2 region (UNITE database, ITSx optional)
#   16S-V4     — bacterial 16S V4 region (SILVA database, no ITSx)
#   18S-V4     — microeukaryote 18S V4 region (PR2 database, no ITSx)
#   18S-AMF    — arbuscular mycorrhizal fungi 18S (MaarjAM database, no ITSx)

# DATABASE OPTIONS (if manually chosen, using flag --db):
#   SILVA        - bacteria SSU 16S rRNA sequences
#   UNITE        - fungal ITS1 and ITS2 regions
#   PR2          - full-length SSU 18S sequences, covers whole eukaryote tree but with focus on protists
#   Maarjam      - AMF 18S SSU sequences
#   EukaryomeITS - ITS1 and ITS2 sequences with good coverage across the eukaryote tree
#   EukaryomeSSU - 18S SSU sequences with good coverage across the eukaryote tree (especially for AMF?)

# EXAMPLES:
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 ITS2
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 ITS2 --skip-itsx
#   sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/table.tsv FAB2 16S-V4
```

Replace the arguments above in <> with your filepaths, project name and primer set (exclude the <>), with a single space between each argument. Run the flag --skip-itsx at the end to skip the ITSx step.

**A note about SBATCH/SLURM scripts if you are unfamiliar:** when you run the sbatch command, it submits the script as a 'SLURM' submission to the computing cluster. This means the 'job' will get in a queue to eventually run. Depending on what time of day/week you submit it, it could take anywhere from 2 seconds to 30 minutes to initiate (usually towards the lower end in my experience). Once it starts, you will get an email saying it started. Next, you will either get an email that the script COMPLETED or FAILED. If it FAILED, it will specify an Exit Code number, usually 2. If it failed, contact me (Peter F) and I can help troubleshoot. Every time you run the pipeline, the job will output a .out and .err file within the directory where the pipeline is stored. You don't really need to worry about these files, UNLESS something fails, then both will be helpful for understanding what error was thrown/what went wrong.

**A note about reference databases:** Generally the default will probably work well, but for 18S AMF and microeukaryote data, the EukaryomeSSU dataset *may* get better taxonomic annotation. It might be worthwhile to test out using this database for taxonomy assignment. If you were feeling spicy, you could also test out using the EukaryomeITS database for taxonomy annotation of fungal ITS data, but at this point I can't comment on how it compares to UNITE taxonomy.

After your job has started, you can watch its progress in real time using the following command:
```bash
# first ensure you are still in the pipeline script folder on MSI
# then run:
tail -f pipeline_*.out

```

This will show you a running output of what the pipeline is doing. If you have run the pipeline before, it will also show you log outputs from previous runs. It's not necessary to do this step, but if you want to track the pipeline in real time, this is an easy way to do it.

After the pipeline is done, or if you want to exit the live view of whats happening, hit Control+C & Enter on Mac (not sure what the equivalent is for Windows).

## Output Files:

The ASV to OTU pipeline will output a LOT of files (see below for descriptions of all), but you will likely be most interested in two:

1. *ProjectName_OTU_with_taxonomy_PrimerSet.txt* -  This is your new data table, functionally equivalent to the input ASV table, but with OTU-level abundances with new taxonomic assignments. (ProjectName is your inputted project name, and PrimerSet will be your primer set)
2. *pipeline_run.log* - This is a log file which has logged a. everything the pipeline did and b. what each of the steps removed. You will want to look at the "Pipeline Summary - Quality Control" section which tells you how many ASVs/OTUs were removed in each step of the process. 

The *pipeline_run.log* also summarizes other files you may be interested in.

## Summary of main pipeline steps:

### Step 0 - ITSx (only relevant for ITS datasets)
 In essence, ITSx identifies and removes the highly conserved 5.8S part of the ITS2 barcodes. This makes for better clustering later, as the 5.8S regions inflates the similarity between sequences. It also discards sequences that are identified as not ITS, which may accidentally discard synthetic community members (synmock). If you run this and find that your synmock does not contain the expected high abundance synthetic community members, it may be because ITSx has removed them. Try re-running the pipeline without ITSx.

### Step 1 - VSEARCH Clustering

This step uses the VSEARCH algorithm to cluster ASVs at 97% similarity (weighted by sequence abundance; *see VSEARCH documentation*), then does an additional chimera check to find and remove any chimeras that might have slipped through DADA2. The number removed here should be quite small (check log file), but may vary wildly by primer set. 

**NOTE:** *Initial testing seems to suggest that up to 20% of the OTUs in 18S-V4 primer datasets are detected to be chimeric. My assumption with the 18S primer is that it is trying to pick up so much sequence diversity across the eukaryotic tree that it picks up a lot of errors a long the way. Nevertheless, I (PF) would like to do more testing to fully understand whats going on here.*

### Step 2 - MUMU/LULU OTU Curation

*mumu* is a post-clustering clean up algorithm that is supposed to find rare variants sequences of common sequences that leak through the clustering steps.
It works in two steps. All sequences are blasted against each other. Sequences with high similarity are checked for patterns in the OTU table. If a rare sequence with high similarity to a common sequence also shows a very similar occurrence pattern, it is merged with this "parent" sequence. These sequences are assumed to be additional sequencing errors that have made it through all previous cleaning steps. Mostly used to clean up singleton and doubleton sequences. *mumu* is a unix version of the LULU package (Frøslev et al., 2017), with some optimisation to run faster.

**NOTE:** *Similar to above, I have found in my preliminary testing of 18S-V4 datasets that more than 50% of OTUs are detected to have a parent OTU by mumu. Similar assumptions to above for why this is happening, but it is worth looking into this, and whether mumu parameters should be tweaked for 18S datasets.*

### Step 3 - DADA2 RDP Taxonomy Assignment

After the OTUs have been created and curated, taxonomy is assigned to the centroid of each OTU using the exact same approach used by Trevor Gould in his DADA2 pipeline. This uses the *assignTaxonomy* function in DADA2, which uses the RDP Naive Bayesian classifier (doi: 10.1128/AEM.00062-07) to classify sequences. Taxonomy databases are unique to each primer (ITS1 and ITS2 have one database: UNITE sh) and are the same databases used by Trevor Gould to assign taxonomy. I have updated each taxonomy database to their most recent version, and they are located in the Kennedy lab shared taxonomy folder on the MSI computing cluster.

This will likely take the longest of any step in the process -- somewhere between 5-45 minutes depending on how many OTUs there are, and the size of the reference database (PR2 takes the longest). 

## Citations

If you use this pipeline, please give credit to Luis Morgado for developing the original Zazzy Metabarcoding Pipeline (https://github.com/ekronold/Zazzy_metabarcoding_pipeline), and Eivind Kverme Ronold for making the initial updates for adaptation to DADA2 output ASV files.

*Additionally, please cite the primary packages and databases used in the pipeline:*


### Packages:

ITSx: Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for use in environmental sequencing. Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco, Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez, Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov, K. Martin Eriksson, R. Henrik Nilsson *Methods in Ecology and Evolution, 4: 914-919, 2013* (https://doi.org/10.1111/2041-210X.12073)

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Frédéric Mahé. (2026) mumu: post-clustering curation tool for metabarcoding data v.1.1.2 (https://github.com/frederic-mahe/mumu)

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188. (https://doi.org/10.1038/s41467-017-01312-x)

Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869


### Databases:

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C. et al. 2013. The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Res. 41:D597–604.

Chuvochina M, Gerken J, Frentrup M, Sandikci Y, Goldmann R, Freese HM, Göker M, Sikorski J, Yarza P, Quast C, Peplies J, Glöckner FO, Reimer LC (2026) SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Nucleic Acids Research, gkaf1247.

Abarenkov K, Nilsson RH, Larsson K-H, Taylor AFS, May TW, Frøslev TG, Pawlowska J, Lindahl B, Põldmaa K, Truong C, Vu D, Hosoya T, Niskanen T, Piirmann T, Ivanov F, Zirk A, Peterson M, Cheeke TE, Ishigami Y, Jansson AT, Jeppesen TS, Kristiansson E, Mikryukov V, Miller JT, Oono R, Ossandon FJ, Paupério J, Saar I, Schigel D, Suija A, Tedersoo L, Kõljalg U. 2024. The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Research, https://doi.org/10.1093/nar/gkad1039

Öpik, M., Vanatoa, A., Vanatoa, E., Moora, M., Davison, J., Kalwij, J.M., Reier, Ü., Zobel, M. 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). New Phytologist 188: 223-241.

Leho Tedersoo, Mahdieh S Hosseyni Moghaddam, Vladimir Mikryukov, Ali Hakimzadeh, Mohammad Bahram, R Henrik Nilsson, Iryna Yatsiuk, Stefan Geisen, Arne Schwelm, Kasia Piwosz, Marko Prous, Sirje Sildever, Dominika Chmolowska, Sonja Rueckert, Pavel Skaloud, Peeter Laas, Marco Tines, Jae-Ho Jung, Ji Hye Choi, Saad Alkahtani, Sten Anslan , EUKARYOME: the rRNA gene reference database for identification of all eukaryotes, Database, Volume 2024, 2024, baae043, https://doi.org/10.1093/database/baae043