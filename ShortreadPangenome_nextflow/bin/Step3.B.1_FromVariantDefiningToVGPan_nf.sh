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
cd ${VGPAN}/augmented

## Running sarek with simulated samples to augment the linear pangenome
### Preparing an input sheet
ls -lrt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1 | sort > sample; \
num_samples=$(wc -l < sample)
printf '1\n%.0s' $(seq 1 ${num_samples}) > lane ;  \
ls -1 ${SRA}/Nremoved/${i}_*1.fq_Ns_removed | sort > fastq_1;  \
ls -1 ${SRA}/Nremoved/${i}_*2.fq_Ns_removed | sort > fastq_2;  \
echo "patient,sample,lane,fastq_1,fastq_2" > header; \
paste sample sample lane fastq_1 fastq_2 | tr "\t" "," > content; \
cat header content > samplesheet.csv; \
rm header content sample lane fastq_1 fastq_2 primary alternate primary_fastq_1 primary_fastq_2 alternate_fastq_1 alternate_fastq_2 # generate samplesheet.csv as an input to sarek

# ### Running sarek pipeline
# nextflow run nf-core/sarek -r 3.3.2 -profile singularity -work-dir ${VGPAN}/augmented/ -resume -params-file ${VGPAN}/augmented/nf-params.json \
# --input ${VGPAN}/augmented/mapped.csv --outdir ${VGPAN}/called/snp/ --tools haplotypecaller --skip_tools baserecalibrator \
# --trim_fastq true --genome custom --fasta ${LINPAN}/${PREFIX}_panref2_sorted.fa --fasta_fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai --igenomes_ignore true \
# --email ${EMAIL}

module load gatk4
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  mkdir -p ${VGPAN}/called/snp/variant_calling/${i}/
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${LINPAN}/${PREFIX}_panref2_sorted.fa -I ${VGPAN}/aligned/$(i}.${PREFIX}_surject.sorted.marked.bam  -O $#{VGPAN}/called/snp/variant_calling/${i}/${i}.haplotypecaller.vcf.gz --tmp-dir .
done
