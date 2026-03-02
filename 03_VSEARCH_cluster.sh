#!/bin/bash
#SBATCH --job-name=VSEARCH_INSERTJOBNAMEHERE
#SBATCH --time=2:00:00
#SBATCH --nodes=1  # specify one node
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=10 
#SBATCH --mem-per-cpu=5G 

#SBATCH -p msismall

#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=INSERTEMAILHERE

module purge

cd /YOURWORKINGDIRECTORY/ #make sure to change this according to your system

#SETTING VARIABLES
VSEARCH=vsearch
THREADS=10
CLUSTER=0.97
LOW_ABUND_SAMPLE=2
PROJ="COOL_PROJECT_NAME"
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
TMP_FASTA3=$(mktemp)
TMP_FASTA4=$(mktemp)
TMP_FASTA5=$(mktemp)
TMP_FASTA7=$(mktemp)

#rename fasta headers in single files
ml perl/5.36.0-gcc-8.2.0-uuoo54n

cat S[0-9][0-9][0-9].fas >> "${TMP_FASTA7}"
perl back_ground_scripts/rename.pl
cat renamed*.fas > "${TMP_FASTA1}"
rm renamed*.fas

ml purge

ml vsearch/2.3.4
# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA7}" \
             --sizein \
             --sizeout \
             --minseqlength 10 \
             --fasta_width 0 \
             --relabel_sha1 \
	     --minuniquesize ${LOW_ABUND_SAMPLE} \
	     --sizein --sizeout \
             --threads ${THREADS} \
             --output "${TMP_FASTA2}" > /dev/null

# Sorting
"${VSEARCH}" --fasta_width 0 \
             --sortbysize "${TMP_FASTA2}" \
             --minseqlength 10 \
             --threads ${THREADS} \
             --output "${TMP_FASTA3}" --sizein --sizeout

#Clustering
"${VSEARCH}" --cluster_size "${TMP_FASTA3}" \
	     --id ${CLUSTER} \
	     --sizein --sizeout \
             --minseqlength 10 \
             --qmask none \
             --threads ${THREADS} \
	     --centroids "${TMP_FASTA4}"

# Chimera checking
"${VSEARCH}" --uchime_denovo "${TMP_FASTA4}" \
            --sizein --sizeout \
            --qmask none \
            --threads ${THREADS} \
            --nonchimeras "${TMP_FASTA5}"

# Sorting & filtering
"${VSEARCH}" --fasta_width 0 \
             --minseqlength 10 \
             --sortbysize "${TMP_FASTA5}" \
             --threads ${THREADS} \
             --output $PROJ.centroids --sizein --sizeout

#sed 's/;size=.*//' -i $PROJ.centroids

#mapping of raw reads against OTU representatives
"${VSEARCH}" --usearch_global "${TMP_FASTA1}" \
	     --db $PROJ.centroids  \
             --minseqlength 10 \
	     --strand plus \
	     --id ${CLUSTER} \
	     --maxaccepts 0 \
             --qmask none \
             --threads ${THREADS} \
	     --uc $PROJ.uc

ml purge

#preparing for making table
sed -i 's/	S0/	barcodelabel=S0/g' $PROJ.uc
sed -i 's/	S1/	barcodelabel=S1/g' $PROJ.uc
sed -i 's/	S2/	barcodelabel=S2/g' $PROJ.uc
sed -i 's/	S3/	barcodelabel=S3/g' $PROJ.uc
sed -i 's/	S4/	barcodelabel=S4/g' $PROJ.uc
sed -i 's/	S5/	barcodelabel=S5/g' $PROJ.uc
sed -i 's/	S6/	barcodelabel=S6/g' $PROJ.uc
sed -i 's/	S7/	barcodelabel=S7/g' $PROJ.uc
sed -i 's/	S8/	barcodelabel=S8/g' $PROJ.uc
sed -i 's/	S9/	barcodelabel=S9/g' $PROJ.uc

#make table
ml python/3.10.9_anaconda2023.03_libmamba
python back_ground_scripts/uc2otutab_original.py $PROJ.uc > $PROJ.otutable

rm -f $PROJ.uc "${TMP_FASTA1}" "${TMP_FASTA2}" "${TMP_FASTA3}" "${TMP_FASTA4}" "${TMP_FASTA5}" "${TMP_FASTA7}"
