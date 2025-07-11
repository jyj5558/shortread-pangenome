#!/bin/bash
#SBATCH --job-name=orgpan_re
#SBATCH -A fnrdewoody
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

###-----Define frequently used variables-----###
N=128 # number of cores
PREFIX=short-swallow-out
BASE=/scratch/negishi/jeon96/swallow # in case of using simulated reads, BASE=/scratch/negishi/jeon96/swallow/sim
APP=/home/jeon96/app
#PREFIX=barnswallow # please define PREFIX as whatever you which
ORGPAN=/scratch/negishi/jeon96/swallow/original/original/pan/out
REF=/scratch/negishi/jeon96/swallow/original/original/ref
CLEANED_SRA=/scratch/negishi/jeon96/swallow/original/original/sra/cleaned
DEPOT=/depot/fnrdewoody
SRA=/scratch/negishi/jeon96/swallow/original/original/sra
LINPAN=/scratch/negishi/jeon96/swallow/linpan/linpan
MCPAN=/scratch/negishi/jeon96/swallow/mcpan/mcpan
VGPAN=/scratch/negishi/jeon96/swallow/vgpan/vgpan
RefInd=refp 
OUT=/scratch/negishi/allen715/pangenome/forNA2/out
REFDIR=/scratch/negishi/allen715/pangenome/forNA2/ref/
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed
GENOME_SIZE=1.1g
POPGEN=/scratch/negishi/jeon96/swallow/popgen
BENCHMARK=/scratch/negishi/jeon96/swallow/benchmarking


# Mapping test data reads

cd ${ORGPAN}/
mkdir -p ./aligned_re/tmp
cd ./aligned_re/

${APP}/vg convert -t ${N} ../${PREFIX}_mod_chopped_new_L.xg -p > ../${PREFIX}_org_new.pg

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
 
## Aligning the individuals
  #${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ${ORGPAN}/${PREFIX}_mod_chopped.xg -g ${ORGPAN}/${PREFIX}_mod_chopped_pruned.gcsa > ${i}_org_aln.gam
  cp ../aligned/${i}_org_aln.gam .
 
## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
  #${APP}/vg filter -t ${N} ${i}_org_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${ORGPAN}/${PREFIX}_mod_chopped.xg > ${i}_org_aln.filtered.gam 
 
done
 
#cat *_org_aln.filtered.gam > combined_org_aln.filtered.gam # -> for population genetic analyses

cd ${ORGPAN}/aligned_re/
 
## Augmenting the graph with all variation from the GAM -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg convert -t ${N} ${ORGPAN}/${PREFIX}_mod_chopped_new_L.xg -p > ${PREFIX}_org.pg
#${APP}/vg augment -t ${N} ${ORGPAN}/${PREFIX}_org.pg ${ORGPAN}/aligned/combined_org_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_orgSV.gam > ${PREFIX}_orgSV.pg 
 
## Indexing the augmented graph -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg mod -t ${N} -X 256 ${PREFIX}_orgSV.pg > ${PREFIX}_orgSV_chopped.pg
#${APP}/vg index -t ${N} -x ${PREFIX}_orgSV_chopped.xg ${PREFIX}_orgSV_chopped.pg
#${APP}/vg prune -t ${N} -k 45 ${PREFIX}_orgSV_chopped.pg > ${PREFIX}_orgSV_chopped_pruned.pg
#${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_orgSV_chopped_pruned.gcsa ${PREFIX}_orgSV_chopped_pruned.pg
 
## Indexing the augmented graph with -L option -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}_orgSV_chopped.pg -x ${PREFIX}_orgSV_chopped_new_L.xg #-L: preserve alt paths in the xg
#${APP}/vg convert -t ${N} ${PREFIX}_orgSV_chopped_new_L.xg -p > ${PREFIX}_orgSV_new.pg


# Variant calling (for SVs) for each individual # use highmem queue
cd ${ORGPAN}/aligned_re/
mkdir -p ../called_re/sv/vg/
 
## Computing the snarls
#${APP}/vg snarls -t ${N} ${PREFIX}_orgSV_new.pg > ${PREFIX}_orgSV.snarls # -> for population genetic analyses
${APP}/vg snarls -t ${N} ../${PREFIX}_new.pg > ${PREFIX}_org.snarls # -> for variant benchmarking
 
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
 
## Realigning the individuals -> for population genetic analyses
  #${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}_orgSV_chopped.xg -g ./${PREFIX}_orgSV_chopped_pruned.gcsa > ${i}_orgSV_aln.gam # -> for population genetic analyses

## Filtering secondary and ambiguous read mappings out of the gam -> for population genetic analyses
  #${APP}/vg filter -t ${N} ${i}_orgSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}_orgSV_chopped.xg > ${i}_orgSV_aln.filtered.gam # -> for population genetic analyses
 
## Computing the support
  #${APP}/vg pack -t ${N} -x ${PREFIX}_orgSV_new.pg -g ${i}_orgSV_aln.filtered.gam -Q 5 -o ${i}_orgSV.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for population genetic analyses
  ${APP}/vg pack -t ${N} -x ../${PREFIX}_new.pg -g ${i}_org_aln.gam -Q 5 -o ${i}_org.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for variant benchmarking
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
cd ${ORGPAN}/aligned_re/
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${APP}/vg call -t ${N} ../${PREFIX}_org_new.pg -r ${PREFIX}_org.snarls -k ${i}_org.pack -s ${i} -a -A -c 50 -C 100000 > ../called_re/sv/vg/${i}_orgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called
done

## Combining separately called SV vcf files
cd ${ORGPAN}/called_re/sv
 
module purge
module load biocontainers
module load bcftools
module load vcftools

## Compressing and indexing each vcf file first #before running this, need to change vcf file headers
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  sed 's/refp#0#//g' ${i}_orgSV.vcf > ${i}_orgSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_orgSV2.vcf -Oz -o ${i}_orgSV.sorted.vcf.gz
  bcftools index ${i}_orgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_orgSV2.vcf
done

bcftools merge -m all -Oz -o ${PREFIX}_orgSV.merged.vcf.gz --threads ${N} *_orgSV.sorted.vcf.gz 


# Variant calling based on surjected bamfiles for each individual
cd ${ORGPAN}/aligned_re/
mkdir -p ../called_re/snp

module purge
module load biocontainers
module load samtools
module load bwa
module load picard
module load nf-core
module load bbmap

## Creating a list of reference paths
${APP}/vg paths -x ${ORGPAN}/${PREFIX}.full.gbz -S ${RefInd} -L > ../${RefInd}.${PREFIX}-pg_paths.txt # use .full graph just to make path lists (used updated version of vg for this - v.1.53)

## Projecting each sample's reads to the RefInd
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  #${APP}/vg filter -t ${N} ${i}_org_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../${PREFIX}_mod_chopped.xg > ${i}_org_interleaved_aln.filtered.gam # -> for population genetic analyses
  ${APP}/vg filter -t ${N} ${i}_org_aln.gam -i -x ../${PREFIX}_mod_chopped.xg > ${i}_org_interleaved_aln.gam # -> for variant benchmarking
  #${APP}/vg surject -x ../${PREFIX}_mod_chopped.xg ${i}_org_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam # -> for population genetic analyses
  ${APP}/vg surject -x ../${PREFIX}_mod_chopped.xg ${i}_org_interleaved_aln.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam # -> for variant benchmarking
done

### Preprocessing the reference genome and bam files
module load gatk4
module load picard
module load samtools
module load bcftools
module load boost

#cd ${REF}
#sed 's/>/>refp#0#/g' GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed.fasta > ref.renamed.sorted.fa
#samtools faidx ref.renamed.sorted.fa
#PicardCommandLine CreateSequenceDictionary --REFERENCE ref.renamed.sorted.fa --OUTPUT ref.renamed.sorted.dict
#bwa index -a bwtsw ref.renamed.sorted.fa

cd ${ORGPAN}/aligned_re/

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### Indexing bam files
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}-pg_surject.sorted.bam ${i}.${PREFIX}-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}-pg_surject.sorted.bam   
done

### Marking duplicates
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam
done

### GATK HaplotypeCaller
for i in `cat ${BASE}/original/original/SRR_Acc_List_new.txt`; do
  echo "processing $i"
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${REF}/ref.renamed.sorted.fa -I ${i}.${PREFIX}-pg_surject.sorted.marked.bam  -O ${ORGPAN}/called_re/snp/${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
  echo "processed $i"
done # use sorted file instead of just renamed 

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
cd ${REF}
PicardCommandLine CreateSequenceDictionary --REFERENCE ref.renamed.fa --OUTPUT ref.renamed.dict
bwa index -a bwtsw ref.renamed.fa

cd ${ORGPAN}/aligned_re/
mkdir -p ../called_re/sv/delly
 
### Running delly2
cd ../called_re/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samples+="${ORGPAN}/aligned_re/${i}.${PREFIX}-pg_surject.sorted.marked.bam "
done 

${DEPOT}/apps/delly/src/delly call -g ${REF}/ref.renamed.sorted.fa ${samples} > orgpan_delly.vcf # used ref.renamed.sorted.fa not ref.renamed.fa 
bcftools view -Oz -o  orgpan_delly.bcf --threads ${N} orgpan_delly.vcf
bcftools index orgpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o orgpan_delly.filtered.bcf orgpan_delly.bcf
bcftools view -Ov -o orgpan_delly.filtered.vcf orgpan_delly.filtered.bcf 