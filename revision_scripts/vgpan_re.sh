#!/bin/bash
#SBATCH --job-name=vgpan_re
#SBATCH -A fnrdewoody
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
REF=/scratch/negishi/jeon96/swallow/original/original/ref
SRA=/scratch/negishi/jeon96/swallow/original/original/sra
LINPAN=/scratch/negishi/jeon96/swallow/linpan/linpan
MCPAN=/scratch/negishi/jeon96/swallow/mcpan/mcpan
VGPAN=/scratch/negishi/jeon96/swallow/vgpan/vgpan
ORGPAN=/scratch/negishi/jeon96/swallow/original/original/pan
CONTIG=/scratch/negishi/jeon96/swallow/contig
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic
GENOME_SIZE=1.1g
CLEANED_SRA=/scratch/negishi/jeon96/swallow/original/original/sra/cleaned
SIM=/scratch/negishi/jeon96/swallow/sim
RefInd=bHirRus1_LinPan
POPGEN=/scratch/negishi/jeon96/swallow/popgen
BENCHMARK=/scratch/negishi/jeon96/swallow/benchmarking

cd ${VGPAN}/
mkdir -p ./augmented_re/
cd ./augmented_re/
mkdir -p ./tmp/

cat ../augmented/*aln.gam > combined_aln.gam 

## Augmenting the graph with all variation from the GAM of the above step
#${APP}/vg convert -t ${N} ${PREFIX}2_mod_chopped_new_L.xg -p > ${PREFIX}2.pg
#${APP}/vg augment -t ${N} ${PREFIX}2.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_aug.gam > ${PREFIX}2_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 # -> for population genetic analyses 
${APP}/vg augment -t ${N} ../augmented/${PREFIX}2.pg combined_aln.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_aug.gam > ${PREFIX}2_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 # -> for variant benchmarking 

## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_aug.pg > ${PREFIX}2_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}2_aug_chopped.xg ${PREFIX}2_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_aug_chopped.pg > ${PREFIX}2_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_aug_chopped_pruned.gcsa ${PREFIX}2_aug_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}2_aug_chopped.pg -x ${PREFIX}2_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}2_aug_chopped_new_L.xg -p > ${PREFIX}2_aug_new.pg


cd ${VGPAN}/
mkdir -p ./aligned_re/tmp
cd ./aligned_re/

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../augmented_re/${PREFIX}2_aug_chopped.xg -g ../augmented_re/${PREFIX}2_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

#cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam # -> for population genetic analyses
 
## Augmenting the graph with all variation from the GAM -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg augment -t ${N} ../augmented/${PREFIX}2_aug_new.pg combined_aug_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_augSV.gam > ${PREFIX}2_augSV.pg 

# Indexing the augmented graph -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_augSV.pg > ${PREFIX}2_augSV_chopped.pg
#${APP}/vg index -t ${N} -x ${PREFIX}2_augSV_chopped.xg ${PREFIX}2_augSV_chopped.pg
#${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_augSV_chopped.pg > ${PREFIX}2_augSV_chopped_pruned.pg
#${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_augSV_chopped_pruned.gcsa ${PREFIX}2_augSV_chopped_pruned.pg

## Indexing the augmented graph with -L option -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}2_augSV_chopped.pg -x ${PREFIX}2_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
#${APP}/vg convert -t ${N} ${PREFIX}2_augSV_chopped_new_L.xg -p > ${PREFIX}2_augSV_new.pg


# Variant calling (for SVs) for each individual using vg
cd ${VGPAN}/aligned_re/
mkdir -p ../called_re/sv/vg/

## Computing the snarls
#${APP}/vg snarls -t ${N} ${PREFIX}2_augSV_new.pg > ${PREFIX}2_augSV.snarls # -> for population genetic analyses
${APP}/vg snarls -t ${N} ../augmented_re/${PREFIX}2_aug_new.pg > ${PREFIX}2_aug.snarls # -> for variant benchmarking

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do

## Realigning the individuals -> for population genetic analyses
  #${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}2_augSV_chopped.xg -g ./${PREFIX}2_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam # -> for population genetic analyses

## Filtering secondary and ambiguous read mappings out of the gam -> for population genetic analyses
  #${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}2_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam # -> for population genetic analyses

## Computing the support
  #${APP}/vg pack -t ${N} -x ${PREFIX}2_augSV_new.pg -g ${i}_augSV_aln.filtered.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for population genetic analyses
  ${APP}/vg pack -t ${N} -x ../augmented_re/${PREFIX}2_aug_new.pg -g ${i}_aug_aln.gam -Q 5 -o ${i}_aug.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for variant benchmarking
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  #${APP}/vg call -t ${N} ${PREFIX}2_augSV_new.pg -r ${PREFIX}2_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}_vgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called # -> for population genetic analyses
  ${APP}/vg call -t ${N} ../augmented_re/${PREFIX}2_aug_new.pg -r ${PREFIX}2_aug.snarls -k ${i}_aug.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called_re2/sv/vg/${i}_vgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called # -> for variant benchmarking
done


# Variant calling based on surjected bamfiles for each individual 
cd ${VGPAN}/aligned_re/
mkdir -p ../called_re/snp
#${APP}/vg gbwt -p -E -x ../augmented/${PREFIX}2-pg.vg --gbz-format -g ../augmented/${PREFIX}2-pg.gbz # make a gbz file 

module purge
module load biocontainers
module load samtools
module load picard
module load bcftools
module load boost
module load gatk4

## Projecting each sample's reads to the reference paths -> filtering suppressed for comparisons with linear assemblies but interleaving
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam # -> for population genetic analyses
  ${APP}/vg filter -t ${N} ../aligned_re/${i}_aug_aln.gam -i -x ../augmented_re/${PREFIX}2_aug_chopped.xg > ${i}_aug_interleaved_aln.gam # -> for variant benchmarking
  #${APP}/vg surject -x ../augmented/${PREFIX}2_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}2-pg_surject.bam # -> for population genetic analyses
  ${APP}/vg surject -x ../augmented_re/${PREFIX}2_aug_chopped.xg ${i}_aug_interleaved_aln.gam --threads ${N} --prune-low-cplx --interleaved -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}2-pg_surject.bam # -> for variant benchmarking
done

## Calling SNPs (and small SVs < 50 bp) 
### Indexing bam files
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}2-pg_surject.sorted.bam ${i}.${PREFIX}2-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}2-pg_surject.sorted.bam   
done

### Marking duplicates
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}2-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}2-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam
done

### Running GATK4 HaplotypeCaller
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do 
  echo "processing $i"
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${LINPAN}/${PREFIX}_panref2_sorted.fa -I ${VGPAN}/called_re/snp/preprocessing/markduplicates/${i}/${i}.md.cram  -O ${VGPAN}/called_re/snp/variant_calling/${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
  echo "processed $i"
done

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
cd ${VGPAN}/aligned_re/
mkdir -p ../called_re/sv/delly

### Running delly2
cd ../called_re/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samples+="${VGPAN}/aligned_re/${i}.${PREFIX}2-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > vgpan_delly.vcf 
bcftools view -Ob -o vgpan_delly.bcf vgpan_delly.vcf
bcftools index vgpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o vgpan_delly.filtered.bcf vgpan_delly.bcf
bcftools view -Ov -o vgpan_delly.filtered.vcf vgpan_delly.filtered.bcf