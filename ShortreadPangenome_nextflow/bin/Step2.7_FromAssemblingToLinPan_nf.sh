#!/bin/bash

SRA=$1
REF=$2
LINPAN=$3
CONTIG=$4
GENOME=$5
N=$6
APP=$7
MASURCA=$8
GENOME_SIZE=$9
PREFIX=${10}
MASURCA_CONFIG_UNMAPPED2=${11}

# Please note that this script was based on Hu et al. (2020) in Legume Genomics: Methods and Protocols (and its scripts) AND Yao et al. (2015) in Genome Biology
module load biocontainers
module load biopython 
module load boost
module load megahit

cd ${SRA}/
cd mapped2/

python3 ${APP}/splitUP.py unmapped2_merged.fastq

mv unmapped2_merged.fastq_R1.fastq unmapped2_R1.fastq
mv unmapped2_merged.fastq_R2.fastq unmapped2_R2.fastq
mv unmapped2_merged.fastq_unpaired.fastq unmapped2_unpaired.fastq

# Assemblyng unmapped reads independently 
## MaSuRCA assembly using unmapped reads
mkdir -p ./masurca_unmapped2

cd ./masurca_unmapped2
cp ~/masurca_config_example ./${MASURCA_CONFIG_UNMAPPED2}
FRAG_MEAN=`awk 'NR==8 {print $6}' ../ihist2.picard.txt`
FRAG_SD=`awk 'NR==8 {print $7}' ../ihist2.picard.txt`
FWD_READS=${SRA}/mapped2/unmapped2_R1.fastq
REV_READS=${SRA}/mapped2/unmapped2_R2.fastq
SGL_READS=${SRA}/mapped2/unmapped2_unpaired.fastq

sed -i "s#PE= pe #PE= pe ${FRAG_MEAN} ${FRAG_SD} ${FWD_READS} ${REV_READS}#" ${MASURCA_CONFIG_UNMAPPED2}
sed -i "s#PE= se #PE= se ${FRAG_MEAN} ${FRAG_SD} ${SGL_READS}#" ${MASURCA_CONFIG_UNMAPPED2}

echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca ${MASURCA_CONFIG_UNMAPPED2} # configure the configuration file beforehand following the updated MaSuRCA manual
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../

## MegaHit assembly using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped2_R1.fastq -2 unmapped2_R2.fastq -r unmapped2_unpaired.fastq -o ./megahit_unmapped2 --out-prefix unmapped2_megahit 
echo "Megahit assembly of unmapped reads finished"

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd masurca_${g}/
  cd CA/
  cat primary.genome.scf.fasta alternative.genome.scf.fasta > ../${g}_masurca.fasta # a limitation in that we cannot distinguish primary and alternate contigs (of an individual or among individuals) since we pooled all the reads among samples 
  cd ../
  sortbyname.sh in=${g}_masurca.fasta out=${g}_masurca_sorted.fasta length descending overwrite=T
  reformat.sh minlength=200 in=${g}_masurca_sorted.fasta out=${g}_masurca_200.fa
  cd ../
  
  cd megahit_${g}/
  cp ${g}_megahit.contigs.fa ${g}_megahit.fasta
  sortbyname.sh in=${g}_megahit.fasta out=${g}_megahit_sorted.fasta length descending overwrite=T
  reformat.sh minlength=200 in=${g}_megahit_sorted.fasta out=${g}_megahit_200.fa
  cd ../
done
