#!/bin/bash

SRA=$1
DB=$2
N=$3
 
# Please note that this script was not used in its entirety for the paper (only the "N removal" step was used since the previous sequence simulating step guaranteed quality and non-contamination). This is an example script to be used when working with sequences from real samples (i.e., for simulated sequences, this may not be used.). 
# Make a "contam_lib" (= $DB) with potential contaminant sequences to be used for Kraken2 beforehand 

# Data Quality Control - may be skipped if using simulated reads of high quality without adapter sequences
cd ${SRA}

module purge
module load biocontainers
module load biopython
module load cutadapt
module load fastqc
module load trim-galore
module load bbmap
module load gcc/12.2.0

export PATH=$HOME/app/kraken2-2.1.5:$PATH

## Initial quality check using Fastqc
mkdir -p ./qc/
cd ./qc/

for g in `ls -lt ../raw/ | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq` # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
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
 
for g in `ls -lt ../raw/ | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq` # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  do
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ./ --paired ../raw/${g}R1_001.fastq ../raw/${g}R2_001.fastq
  cat ${g}1.fastq_trimming_report.txt ${g}2.fastq_trimming_report.txt > ${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
done
# in cleaned file make sure fq are gzipped: 
for f in *.fq; do pigz -p 32 "$f"; done

cd ../

## Contamination removal using Kraken 
mkdir -p ./filtered/
cd ./filtered/

 
for g in `ls -lt ${SRA}/cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  kraken2 --db ${DB} \
  --threads ${N} \
  --unclassified-out ./${g}filt_R1.fq \
  --classified-out ./${g}cont_R1.fq \
  --output ./${g}_R1_out \
  --report ./${g}R1_kraken_report.txt \
  --confidence 0.3 \
  --minimum-hit-groups 2 \
  --use-names ${SRA}/cleaned/${g}R1_001_val_1.fq.gz
 
  kraken2 --db ${DB} \
  --threads ${N} \
  --unclassified-out ./${g}filt_R2.fq \
  --classified-out ./${g}cont_R2.fq \
  --output ./${g}_R2_out \
  --report ./${g}R2_kraken_report.txt \
  --confidence 0.3 \
  --minimum-hit-groups 2 \
  --use-names ${SRA}/cleaned/${g}R1_001_val_1.fq.gz
 
# # did not use
# #   k2 classify --db /scratch/negishi/allen715/kraken2-2.1.5/pluspf_08gb --threads ${N} --unclassified-out ./${g}filt_R1.fq --classified-out ./${g}cont_R1.fq --output ./${g}_R1_out -report ./${g}R1_kraken_report.txt --confidence 0.3 --minimum-hit-groups 2 --use-names ${SRA}/cleaned/${g}R1_001_val_1.fq.gz 
# #   k2 classify --db /scratch/negishi/allen715/kraken2-2.1.5/pluspf_08gb --threads ${N} --unclassified-out ./${g}filt_R2.fq --classified-out ./${g}cont_R2.fq --output ./${g}_R2_out -report ./${g}R2_kraken_report.txt --confidence 0.3 --minimum-hit-groups 2 --use-names ${SRA}/cleaned/${g}R2_001_val_2.fq.gz -
done

for g in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  repair.sh in1=${g}R1.fq in2=${g}R2.fq out1=${g}fixed_R1.fq out2=${g}fixed_R2.fq outs=${g}fixed_s.fq
done

cd ${SRA}

## N removal using "Remove_reads_with_Ns.py" script from Hu et al. (2020)
mkdir -p ./Nremoved/
cd ./Nremoved/

for i in $(ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1); do
    r1=$(ls ../filtered/${i}_*_filt_fixed_R1.fq | head -n 1)
    r2=$(ls ../filtered/${i}_*_filt_fixed_R2.fq | head -n 1)
    echo "Running on $r1 and $r2"
    python ~/Remove_reads_with_Ns.py "$r1" "$r2"
done

for i in $(ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1); do
    echo "$i"
done

mv ../filtered/*fq_Ns_removed .

# merge all R1 files
cat *_R1.fq_Ns_removed >> merged_Nremoved_1.fq

# merge all R2 files
cat *_R2.fq_Ns_removed >> merged_Nremoved_2.fq

cd ../

