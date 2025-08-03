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

# Please note that this script was based on Hu et al. (2020) in Legume Genomics: Methods and Protocols (and its scripts) AND Yao et al. (2015) in Genome Biology

cd ${SRA}/

# Extracting mapped and unmapped reads (on the linear representative genome) separately
module purge
module load biocontainers
module load samtools
module load bamtools
module load biopython

cd ./mapped

samtools view -@ ${N} -b -F 4 -f 64 mapped.sam > mapped_f.bam # mapped forward reads
samtools view -@ ${N} -b -F 4 -f 128 mapped.sam > mapped_r.bam # mapped reverse reads
samtools merge -@ ${N} mapped_merged.bam mapped_f.bam mapped_r.bam
samtools sort -@ ${N} -n -o mapped_sorted.bam mapped_merged.bam
bamtools convert -in mapped_sorted.bam -out mapped_merged.fastq -format fastq

python ${APP}/splitUP.py mapped_merged.fastq # splitUP.py script from Hu et al. (2020)

# ###
# # If the mapped_merged.fastq file is too large, it will cause an out-of-memory issue. In that case, split the file by 4n lines first and then concatenate them later like below.
# 
# cp mapped_merged.fastq mapped_merged_cp.fastq # copy the file for redundancy
# split -l 100000000 --numeric-suffixes=1 --additional-suffix=.fastq --verbose mapped_merged.fastq mapped_merged_
# 
# for i in {01..52}; do
# echo "printing tail ${i}"
# tail -4 mapped_merged_${i}.fastq # manually check if all the files are split at the end of a fastq-formatted-sequence line
# done
# 
# for i in {01..52}; do
# python ~/splitUP.py mapped_merged_${i}.fastq
# done
# 
# #For the unpaired, concatenate them and rerun splitUP.py again to avoid artefacts due to split command.
# 
# for i in {01..52}; do
# cat mapped_merged_${i}.fastq_unpaired.fastq >> mapped_merged_53.fastq # Do not forget later that this '53' file is a concatenated file of all '1'-'52' unpaired files.
# done 
# 
# python ~/splitUP.py mapped_merged_53.fastq
# mv mapped_merged_53.fastq_unpaired.fastq mapped_merged.fastq_unpaired.fastq
# 
# for i in {01..53}; do
# cat mapped_merged_${i}.fastq_R1.fastq >> mapped_merged.fastq_R1.fastq
# cat mapped_merged_${i}.fastq_R2.fastq >> mapped_merged.fastq_R2.fastq
# done

# ###

samtools view -@ ${N} -b -f 68 mapped.sam > unmapped_f.bam # unmapped forward reads
samtools view -@ ${N} -b -f 132 mapped.sam > unmapped_r.bam # unmapped reverse reads
samtools merge -@ ${N} unmapped_merged.bam unmapped_f.bam unmapped_r.bam
samtools sort -@ ${N} -n -o unmapped_sorted.bam unmapped_merged.bam
bamtools convert -in unmapped_sorted.bam -out unmapped_merged.fastq -format fastq

python ${APP}/splitUP.py unmapped_merged.fastq

mv mapped_merged.fastq_R1.fastq mapped_R1.fastq
mv mapped_merged.fastq_R2.fastq mapped_R2.fastq
mv mapped_merged.fastq_unpaired.fastq mapped_unpaired.fastq
mv unmapped_merged.fastq_R1.fastq unmapped_R1.fastq
mv unmapped_merged.fastq_R2.fastq unmapped_R2.fastq
mv unmapped_merged.fastq_unpaired.fastq unmapped_unpaired.fastq
