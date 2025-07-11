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

# Meta-assembling mapped and unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu

mkdir -p ./canu_mapped/
mkdir -p ./canu_unmapped/

## Canu assmebly using mapped contigs
cd ./canu_mapped
echo "Canu assembly of mapped reads started"
canu -p final_mapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_mapped/mapped_masurca_200.fa ${SRA}/mapped/megahit_mapped/mapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_mapped.contigs.fasta final_mapped.unassembled.fasta > final_mapped.fa
echo "Canu assembly of mapped reads finishted"
cd ../

## Canu assmebly using unmapped contigs
cd ./canu_unmapped
echo "Canu assembly of unmapped reads started"
canu -p final_unmapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_unmapped/unmapped_masurca_200.fa ${SRA}/mapped/megahit_unmapped/unmapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_unmapped.contigs.fasta final_unmapped.unassembled.fasta > final_unmapped.fa
echo "Canu assembly of unmapped reads finishted"
cd ../
