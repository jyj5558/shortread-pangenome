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

# Read mapping using Bowtie2 
module purge
module load biocontainers
module load bowtie2
module load bbmap
module load samtools
module load r
 
bowtie2-build ${REF}/${GENOME}.fna ${REF}/${GENOME}
 
mkdir -p ./mapped
cd ./mapped
 
bowtie2 -p ${N} -I 0 -X 1000 -x ${REF}/${GENOME} -1 ../Nremoved/merged_Nremoved_1.fq -2 ../Nremoved/merged_Nremoved_2.fq --end-to-end --sensitive -S mapped.sam 2> bowtie.log
 
## Summary stats on bam files
samtools sort -@ ${N} -o mapping_sorted.bam mapped.sam

## Average insert size of the mapped paired-end reads for the MaSuRCA assembly
reformat.sh in=mapping_sorted.bam ihist=ihist.txt
${APP}/jdk-21.0.1/bin/java -Duser.country=US -Duser.language=en -jar ${APP}/picard.jar CollectInsertSizeMetrics -I mapping_sorted.bam -O ihist.picard.txt -H ihist.picard.pdf -M 0.05 #--DEVIATIONS # picard is better than reformat.sh
 
### Mapping rate, depth, and breadth
echo "mapping rate is" > ./mapping_rate.txt
samtools flagstat -@ ${N} ./mapping_sorted.bam >> ./mapping_rate.txt 
echo "depth is" > ./mapping_depth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> ./mapping_depth.txt 
echo "breadth is" > ./mapping_breadth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./mapping_breadth.txt


# Extracting mapped and unmapped reads (on the linear representative genome) separately
module purge
module load biocontainers
module load samtools
module load bamtools
module load biopython

samtools view -@ ${N} -b -F 4 -f 64 mapped.sam > mapped_f.bam # mapped forward reads
samtools view -@ ${N} -b -F 4 -f 128 mapped.sam > mapped_r.bam # mapped reverse reads
samtools merge -@ ${N} mapped_merged.bam mapped_f.bam mapped_r.bam
samtools sort -@ ${N} -n -o mapped_sorted.bam mapped_merged.bam
bamtools convert -in mapped_sorted.bam -out mapped_merged.fastq -format fastq

samtools view -@ ${N} -b -f 68 mapped.sam > unmapped_f.bam # unmapped forward reads
samtools view -@ ${N} -b -f 132 mapped.sam > unmapped_r.bam # unmapped reverse reads
samtools merge -@ ${N} unmapped_merged.bam unmapped_f.bam unmapped_r.bam
samtools sort -@ ${N} -n -o unmapped_sorted.bam unmapped_merged.bam
bamtools convert -in unmapped_sorted.bam -out unmapped_merged.fastq -format fastq
