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
cd ./mapped

# Assembling mapped and unmapped reads independently, redundantly (MaSuRCA and MegaHit) modifying the approach of Yao et al. (2015)
module purge
module load biocontainers
module load megahit

## MegaHit assembly using mapped reads
echo "Megahit assembly of mapped reads started"
megahit -t ${N} -1 mapped_R1.fastq -2 mapped_R2.fastq -r mapped_unpaired.fastq -o ./megahit_mapped --out-prefix mapped_megahit 
echo "Megahit assembly of mapped reads finished"

## MegaHit assembly using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped_R1.fastq -2 unmapped_R2.fastq -r unmapped_unpaired.fastq -o ./megahit_unmapped --out-prefix unmapped_megahit 
echo "Megahit assembly of unmapped reads finished"

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in mapped unmapped; do
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

## Quast assessment - skip
#module load biocontainers
#module load quast

#echo "Masurca mapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_masurca_mapped --eukaryote ./masurca_mapped/mapped_masurca_200.fa
#echo "Masurca mapped assembly Quast finished"
 
#echo "Masurca unmapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_masurca_unmapped --eukaryote ./masurca_unmapped/unmapped_masurca_200.fa
#echo "Masurca mapped assembly Quast finished"

#echo "Megahit mapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_megahit_mapped --eukaryote ./megahit_mapped/mapped_megahit_200.fa
#echo "Megahit mapped assembly Quast finished"
 
#echo "Megahit unmapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_megahit_unmapped --eukaryote ./megahit_unmapped/unmapped_megahit_200.fa
#echo "Megahit unmapped assembly Quast finished"
