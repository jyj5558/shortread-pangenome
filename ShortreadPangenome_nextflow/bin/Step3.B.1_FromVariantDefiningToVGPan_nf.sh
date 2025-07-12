#!/bin/bash

SRA=$1
REF=$2
LINPAN=$3
CONTIG=$4
VGPAN=$5
PREFIX=$6
N=$7
DEPOT=$8
APP=$9
EMAIL=$10

# Constructing a VG pangenome
module purge
module load biocontainers
module load bwa
module load picard
module load boost
module load bcftools
module load htslib
module load nf-core

mkdir -p ${VGPAN}/augmented/called/snp/

## Calling SNPs (and small SVs < 50 bp) 
### Indexing reference
cd ${LINPAN}

PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.fa --OUTPUT ${PREFIX}_panref2_sorted.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.fa

### Mapping reads
cd ${VGPAN}/augmented
mkdir -p ./mapped/
cd ./mapped/

for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${LINPAN}/${PREFIX}_panref2_sorted.fa  ${SRA}/Nremoved/${i}_*1.fq_Ns_removed ${SRA}/Nremoved/${i}_*2.fq_Ns_removed > ${i}.sam
  PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT ${i}.sam --OUTPUT ${i}_sorted.bam --SORT_ORDER coordinate
  PicardCommandLine MarkDuplicates --INPUT ${i}_sorted.bam --OUTPUT ${i}_sorted.marked.bam --METRICS_FILE ${i}_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}_sorted.marked.bam
done

### Running GATK4 HaplotypeCaller
module load gatk4
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${LINPAN}/${PREFIX}_panref2_sorted.fa -I ${id}_sorted.marked.bam  -O ${VGPAN}/augmented/called/snp/${id}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
done
