#!/bin/bash
#SBATCH --job-name=ITSx_INSERTJOBNAMEHERE
#SBATCH --time=00:30:00
#SBATCH --nodes=1  # specify one node
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=40 
#SBATCH --mem-per-cpu=5G 

#SBATCH -p msismall

#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=INSERTEMAILHERE

set -e  # stop immediately on any error

module purge
module load conda

source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh

conda activate itsx_env

module load vsearch

cd YOURWORKINGDIRECTORY

# SETTING VARIABLES
VSEARCH=$(which vsearch)
SPLITS=40
mkdir -p uncut_fasta  # -p prevents error if directory already exists

for FILE in `ls S[0-9][0-9][0-9].fas` ; do

    TMP_FASTA1=$(mktemp)  # moved inside loop

    cat "${FILE}" > uncut_fasta/"${FILE}"

    faSplit sequence "${FILE}" $SPLITS splitted

    for f in `ls splitted*` ; do
        ITSx -i $f --complement F -t F --preserve T -o $f.out &
    done

    wait

    cat splitted*ITS2.fasta >> "${TMP_FASTA1}"

    rm splitted*

    "${VSEARCH}" --derep_fulllength "${TMP_FASTA1}" \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --minuniquesize 1 \
         --relabel_sha1 \
         --output "${FILE}" > /dev/null

    rm "${TMP_FASTA1}"
done