#!/bin/bash
#SBATCH --job-name=
#SBATCH -A 
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

###-----Define frequently used variables-----###
N=128
PREFIX=short-read
BASE=/scratch/negishi/jeon96/swallow # in case of using simulated reads, BASE=/scratch/negishi/jeon96/swallow/sim
DEPOT=/depot/fnrdewoody
APP=/home/jeon96/app
REF=/scratch/negishi/jeon96/swallow/original/ref
SRA=/scratch/negishi/jeon96/swallow/original/sra
LINPAN=/scratch/negishi/jeon96/swallow/linpan
MCPAN=/scratch/negishi/jeon96/swallow/mcpan
VGPAN=/scratch/negishi/jeon96/swallow/vgpan
ORGPAN=/scratch/negishi/jeon96/swallow/original/pan
CONTIG=/scratch/negishi/jeon96/swallow/contig
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic
GENOME_SIZE=1.1g
CLEANED_SRA=/scratch/negishi/jeon96/swallow/original/sra/cleaned
SIM=/scratch/negishi/jeon96/swallow/sim
RefInd=bHirRus1_LinPan
POPGEN=/scratch/negishi/jeon96/swallow/popgen
BENCHMARK=/scratch/negishi/jeon96/swallow/benchmarking

# Please note that this script was not used in its entirety for the paper (only the "N removal" step was used since the previous sequence simulating step guaranteed quality and non-contamination). This is an example script to be used when working with sequences from real samples (i.e., for simulated sequences, this may not be used.). 

# Data Quality Control - may be skipped if using simulated reads of high quality without adapter sequences
cd ${SRA}

module purge
module load biocontainers
module load kraken2
module load biopython
module load cutadapt
module load fastqc
module load trim-galore

## Initial quality check using Fastqc
mkdir -p ./qc/
cd ./qc/

for g in `ls -lt ../raw/ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  fastqc ../raw/${g}R1_001.fastq --extract --quiet
  fastqc ../raw/${g}R2_001.fastq --extract --quiet
  
  cat ${g}1_fastqc/summary.txt ${g}2_fastqc/summary.txt > ${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}fastqc_summary.txt)
  echo "raw"
  echo "$FAIL"
  rm -r ${g}?_fastqc* 
done
cd ../

## Adapter & low quality reads removal using Trim_galore
mkdir -p ./cleaned/
cd ./cleaned/

for g in `ls -lt ../raw/ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ./ --paired ../raw/${g}R1_001.fastq ../raw/${g}R2_001.fastq
  cat ${g}1.fastq_trimming_report.txt ${g}2.fastq_trimming_report.txt > ${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
done
cd ../

## Contamination removal using Kraken 
mkdir -p ./filtered/
cd ./filtered/

DB=/scratch/bell/dewoody/LEPC/JY-test/db/contam_lib # make a "contam_lib" beforehand with potential contaminant sequences
KRAKEN2_NUM_THREADS=64
export KRAKEN2_DB_PATH="/scratch/bell/dewoody/LEPC/JY-test/db/contam_lib/"

mkdir -p ./filtered/
cd ./filtered/

for g in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  kraken2 --db ${DB} --quick --threads ${KRAKEN2_NUM_THREADS} --output ${g}cseq --paired --unclassified-out ${g}filtered#.fq ${SRA}/cleaned/${g}R1_001_val_1.fq.gz ${SRA}/cleaned/${g}R2_001_val_2.fq.gz --report ${i}_kraken_report.txt --use-names
done
cd ../

## N removal using "Remove_reads_with_Ns.py" script from Hu et al. (2020)
mkdir -p ./Nremoved/
cd ./Nremoved/

for i in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`
  do
  python ~/Remove_reads_with_Ns.py ../filtered/${i}_*_filtered_1.fq ../filtered/${i}_*_filtered_2.fq
done

for i in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`
  do
  cat ./${i}*_1.fq* >> merged_Nremoved_1.fq # pooling forward reads (metagenome-like assembly strategy; Hu et al. (2020), Yao et al. (2015))
  cat ./${i}*_2.fq* >> merged_Nremoved_2.fq # pooling backward reads
done
cd ../
