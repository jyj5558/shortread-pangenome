#!/bin/bash
#SBATCH --job-name=mcpan_re
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

mv /scratch/negishi/jeon96/swallow/linpan/temp/swallow/linpan /scratch/negishi/jeon96/swallow/linpan/

cd ${MCPAN}/${PREFIX}-pg/

cat *aln.gam > combined_aln.gam 

## Augmenting the graph with all variation from the GAM of the above step
#${APP}/vg convert -t ${N} ${PREFIX}_mod_chopped_new_L.xg -p > ${PREFIX}.pg
#${APP}/vg augment -t ${N} ${PREFIX}.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_aug.gam > ${PREFIX}_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 # -> for population genetic analyses 
${APP}/vg augment -t ${N} ${PREFIX}.pg combined_aln.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_aug.gam > ${PREFIX}_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 # -> for variant benchmarking 

## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_aug.pg > ${PREFIX}_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}_aug_chopped.xg ${PREFIX}_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_aug_chopped.pg > ${PREFIX}_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_aug_chopped_pruned.gcsa ${PREFIX}_aug_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}_aug_chopped.pg -x ${PREFIX}_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}_aug_chopped_new_L.xg -p > ${PREFIX}_aug_new.pg


mkdir -p ./aligned_re/tmp
cd ./aligned_re/

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../${PREFIX}_aug_chopped.xg -g ../${PREFIX}_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

#cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam # -> for population genetic analyses 

## Augmenting the graph with all variation from the GAM -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg augment -t ${N} ../${PREFIX}_aug_new.pg combined_aug_aln.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_augSV.gam > ${PREFIX}_augSV.pg 

# Indexing the augmented graph -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg mod -t ${N} -X 256 ${PREFIX}_augSV.pg > ${PREFIX}_augSV_chopped.pg
#${APP}/vg index -t ${N} -x ${PREFIX}_augSV_chopped.xg ${PREFIX}_augSV_chopped.pg
#${APP}/vg prune -t ${N} -k 45 ${PREFIX}_augSV_chopped.pg > ${PREFIX}_augSV_chopped_pruned.pg
#${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_augSV_chopped_pruned.gcsa ${PREFIX}_augSV_chopped_pruned.pg

## Indexing the augmented graph with -L option -> skipped for variant benchmarking in revision (but not skipped for population genetic analyses)
#${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}_augSV_chopped.pg -x ${PREFIX}_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
#${APP}/vg convert -t ${N} ${PREFIX}_augSV_chopped_new_L.xg -p > ${PREFIX}_augSV_new.pg


# Variant calling (for SVs) for each individual using vg
cd ${MCPAN}/${PREFIX}-pg/aligned_re/
mkdir -p ../called_re/sv/vg/

## Computing the snarls
#${APP}/vg snarls -t ${N} ${PREFIX}_augSV_new.pg > ${PREFIX}_augSV.snarls # -> for population genetic analyses
${APP}/vg snarls -t ${N} ../${PREFIX}_aug_new.pg > ../${PREFIX}_aug.snarls # -> for variant benchmarking

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do

## Realigning the individuals -> for population genetic analyses 
  #${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}_augSV_chopped.xg -g ./${PREFIX}_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam # -> for population genetic analyses
  
## Filtering secondary and ambiguous read mappings out of the gam -> for population genetic analyses
  #${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam # -> for population genetic analyses

## Computing the support
  #${APP}/vg pack -t ${N} -x ${PREFIX}_augSV_new.pg -g ${i}_augSV_aln.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for population genetic analyses
  ${APP}/vg pack -t ${N} -x ${PREFIX}_aug_new.pg -g ${i}_aug_aln.gam -Q 5 -o ${i}_aug.pack #-Q 5: ignore mapping and base qualitiy < 5 # -> for variant benchmarking
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  #${APP}/vg call -t ${N} ${PREFIX}_augSV_new.pg -r ${PREFIX}_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}_mcSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called # -> for population genetic analyses
  ${APP}/vg call -t ${N} ../${PREFIX}_aug_new.pg -r ../${PREFIX}_aug.snarls -k ${i}_aug.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called_re/sv/vg/${i}_mcSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called # -> for variant benchmarking
done


# Variant calling based on surjected bamfiles for each individual (If Sarek pipeline does not work, just manually call SNPs using GATK)
cd ${MCPAN}/${PREFIX}-pg/aligned_re/
mkdir -p ../called_re/snp

module purge
module load biocontainers
module load picard
module load bwa
module load boost
module load samtools
module load bbmap
module load gatk4

## Creating a list of reference paths
${APP}/vg paths -x ${MCPAN}/${PREFIX}-pg/${PREFIX}-pg.full.gbz -S ${RefInd} -L > ../${RefInd}.${PREFIX}-pg_paths.txt # use .full graph just to make comprehensive path lists

## Projecting each sample's reads to the RefInd
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam # -> for population genetic analyses
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -i -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_interleaved_aln.gam # -> for variant benchmarking
  #${APP}/vg surject -x ../${PREFIX}_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam # -> for population genetic analyses
  ${APP}/vg surject -x ../${PREFIX}_aug_chopped.xg ${i}_aug_interleaved_aln.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam # -> for variant benchmarking
done

## Calling SNPs (and small SVs < 50 bp) 
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

### Running GATK4 HaplotypeCaller
PicardCommandLine CreateSequenceDictionary --REFERENCE ${LINPAN}/${PREFIX}_panref2_resorted.renamed.fa --OUTPUT ${LINPAN}/${PREFIX}_panref2_resorted.renamed.dict

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do 
  echo "processing $i"
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R ${LINPAN}/${PREFIX}_panref2_resorted.renamed.fa -I ${i}.${PREFIX}-pg_surject.sorted.marked.bam  -O ${MCPAN}/${PREFIX}-pg/called_re/snp/${i}.haplotypecaller.vcf.gz --tmp-dir . --native-pair-hmm-threads ${N}
  echo "processed $i"
done

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.renamed.fa --OUTPUT ${PREFIX}_panref2_sorted.renamed.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.renamed.fa

mkdir -p ${MCPAN}/${PREFIX}-pg/called_re/sv/delly

### Running delly2
cd ${MCPAN}/${PREFIX}-pg/called_re/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samples+="${MCPAN}/${PREFIX}-pg/aligned_re/${i}.${PREFIX}-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa ${samples} > mcpan_delly.vcf 
bcftools view -Ob -o mcpan_delly.bcf mcpan_delly.vcf
bcftools index mcpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o mcpan_delly.filtered.bcf mcpan_delly.bcf
bcftools view -Ov -o mcpan_delly.filtered.vcf mcpan_delly.filtered.bcf