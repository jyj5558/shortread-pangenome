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


# Minigraph-Cactus Pangenome
# Mapping test data reads - when not working with simulated reads, use the pangenome files based on ${PREFIX}-pg.vg, not ${PREFIX}_aug.pg (i.e., without augmentation)
cd ${MCPAN}/${PREFIX}-pg/
mkdir -p ./aligned/tmp
cd ./aligned/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../${PREFIX}_aug_chopped.xg -g ../${PREFIX}_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) 
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam
#cat *_aug_aln.gam > combined_aug_aln.gam

## Augmenting the graph with all variation from the GAM
${APP}/vg augment -t ${N} ../${PREFIX}_aug_new.pg combined_aug_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_augSV.gam > ${PREFIX}_augSV.pg 
#${APP}/vg augment -t ${N} ../${PREFIX}_aug_new.pg combined_aug_aln.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_augSV.gam > ${PREFIX}_augSV.pg 

# Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_augSV.pg > ${PREFIX}_augSV_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}_augSV_chopped.xg ${PREFIX}_augSV_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_augSV_chopped.pg > ${PREFIX}_augSV_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_augSV_chopped_pruned.gcsa ${PREFIX}_augSV_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}_augSV_chopped.pg -x ${PREFIX}_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}_augSV_chopped_new_L.xg -p > ${PREFIX}_augSV_new.pg


# Variant calling (for SVs) for each individual using vg
cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/sv/vg/

## Computing the snarls
${APP}/vg snarls -t ${N} ${PREFIX}_augSV_new.pg > ${PREFIX}_augSV.snarls

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}_augSV_chopped.xg -g ./${PREFIX}_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam 
  ${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam

## Computing the support
  ${APP}/vg pack -t ${N} -x ${PREFIX}_augSV_new.pg -g ${i}_augSV_aln.filtered.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
  #${APP}/vg pack -t ${N} -x ${PREFIX}_augSV_new.pg -g ${i}_augSV_aln.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg call -t ${N} ${PREFIX}_augSV_new.pg -r ${PREFIX}_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}_mcSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called
done


# Filtering SV vcf files from vg
cd ${MCPAN}/${PREFIX}-pg/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools

## Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  sed 's/bHirRus1_LinPan#0#//g' ${i}_mcSV.vcf > ${i}_mcSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_mcSV2.vcf -Oz -o ${i}_mcSV.sorted.vcf.gz
  bcftools index ${i}_mcSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_mcSV2.vcf
done

## Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 

## Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSV.merged.vcf.gz > vcf-stats.txt

## In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.       1st Qu.    Median     Mean     3rd Qu.   Max.
# -26.684     9.542     21.751     35.119    43.435 23775.500
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.050    0.583    1.570    2.000 6173.960
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.000    0.120    0.320    0.583    1.000    1.545    2.500    4.320    6173.960

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0000  0.0800   0.3031  0.6000  1.0000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4012   1.3013   1.4679   1.4516   1.6741   1.9699

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.2819  0.2897  0.2963  0.3031  0.3061  0.4271

quit()

## Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --out ${PREFIX}_mcSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Annotating SVs in vcf file
bcftools sort ${PREFIX}_mcSV.filtered.recode.vcf -Oz -o ${PREFIX}_mcSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}_mcSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}_mcSV.filtered.popped.vcf # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree.

singularity run ${APP}/vcflib-1.09.simg
N=128
PREFIX=short-read # needs to define again inside Apptainer
vcfwave -I 1000 -t ${N} ${PREFIX}_mcSV.filtered.popped.vcf > ${PREFIX}_mcSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
truvari anno svinfo -m 50 -o ${PREFIX}_mcSV.decomposed.annotated.vcf ${PREFIX}_mcSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python3 ${APP}/SVextractor_fromTruvariAnno.py 50 ${PREFIX}_mcSV.decomposed.annotated.vcf ${PREFIX}_mcSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)


# Variant calling (for SNPs using Sarek and for SVs using delly2) based on surjected bamfiles for each individual (If Sarek pipeline does not work, just manually call SNPs using GATK)
cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/snp

module purge
module load biocontainers
module load picard
module load bwa
module load boost
module load nf-core
module load samtools
module load bbmap

## Creating a list of reference paths
${APP}/vg paths -x ${MCPAN}/${PREFIX}-pg/${PREFIX}-pg.full.gbz -S ${RefInd} -L > ../${RefInd}.${PREFIX}-pg_paths.txt # use .full graph just to make comprehensive path lists

## Projecting each sample's reads to the RefInd 
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -i -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_interleaved_aln.gam
  ${APP}/vg surject -x ../${PREFIX}_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam
  #${APP}/vg surject -x ../${PREFIX}_aug_chopped.xg ${i}_aug_interleaved_aln.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam
done

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### Indexing bam files
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}-pg_surject.sorted.bam ${i}.${PREFIX}-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}-pg_surject.sorted.bam   
done

### Preprocessing input files for sarek
#cd $CLUSTER_SCRATCH
#mkdir -p nf-core
#mkdir -p /scratch/negishi/jeon96/singularity/cache
#cd nf-core  
#nf-core download sarek # download sarek pipeline beforehand

cd ${MCPAN}/${PREFIX}-pg/aligned/

cat ${BASE}/original/SRR_Acc_List.txt > sample; \
#printf '1\n%.0s' {1..25} > lane ;  \
ls -1 ${MCPAN}/${PREFIX}-pg/aligned/*surject.sorted.bam > bam; \
ls -1 ${MCPAN}/${PREFIX}-pg/aligned/*surject.sorted.bam.bai > bai; \
echo "patient,sample,bam,bai" > header; \
paste sample sample bam bai | tr "\t" "," > content; \
cat header content > mapped.csv; \
rm header content sample bam bai # generate mapped.csv as an input to sarek

sed 's/>/>bHirRus1_LinPan#0#/g' ${LINPAN}/${PREFIX}_panref2_sorted.fa > ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa # accommodate "reference contig" names of the surjected bam files
samtools faidx ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa
cat ${PREFIX}_panref2_sorted.renamed.fa.fai | cut -f1 | grep "unmap_tig" | sort > unmap_tig.sorted
cat ${PREFIX}_panref2_sorted.renamed.fa.fai | cut -f1 | grep -v "unmap_tig" | sort > ref_tig.sorted
cat unmap_tig.sorted ref_tig.sorted > bamcontig.list
sortbyname.sh in=${PREFIX}_panref2_sorted.renamed.fa out=${PREFIX}_panref2_resorted.renamed.fa list=bamcontig.list overwrite=T # match the order of contigs of the linear pangenome reference with the order in the bam files for sarek
samtools faidx ${PREFIX}_panref2_resorted.renamed.fa 

# Generateing nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#    "input": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/aligned\/mapped.csv",
#    "step": "markduplicates",
#    "outdir": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/called\/snp\/",
#    "tools": "haplotypecaller",
#    "skip_tools": "baserecalibrator",
#    "trim_fastq": true,
#    "genome": "custom",
#    "fasta": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa",
#    "fasta_fai": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa.fai",
#    "igenomes_ignore": true,
#    "email": "jeon96@purdue.edu",
#    "multiqc_title": "short-read-mcpan_swallow"
#}

### Running sarek pipeline
nextflow run nf-core/sarek -r 3.3.2 -profile singularity -work-dir ${MCPAN}/${PREFIX}-pg/aligned/ -params-file ${MCPAN}/${PREFIX}-pg/aligned/nf-params.json #-resume 
#--step markduplicates --input ${MCPAN}/${PREFIX}-pg/aligned/mapped.csv --outdir ${MCPAN}/${PREFIX}-pg/called/snp --tools haplotypecaller --skip_tools baserecalibrator # used when resuming from "markduplicates" step

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
#cd ${LINPAN}
#PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.renamed.fa --OUTPUT ${PREFIX}_panref2_sorted.renamed.dict
#bwa index -a bwtsw ${PREFIX}_panref2_sorted.renamed.fa

cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/sv/delly

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam
done

### Running delly2
cd ../called/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${MCPAN}/${PREFIX}-pg/aligned/${i}.${PREFIX}-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa ${samples} > mcpan_delly.vcf 
bcftools view -Ob -o mcpan_delly.bcf mcpan_delly.vcf
bcftools index mcpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o mcpan_delly.filtered.bcf mcpan_delly.bcf
bcftools view -Ov -o mcpan_delly.filtered.vcf mcpan_delly.filtered.bcf


# Filtering vcf files
cd ${MCPAN}/${PREFIX}-pg/called/snp

module purge
module load biocontainers
module load bcftools
module load vcftools

## Filtering SNP vcf files
### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_mcSNP.sorted.vcf.gz
  bcftools index ${i}_mcSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o ${PREFIX}_mcSNP.merged.vcf.gz --threads ${N} *_mcSNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSNP.merged.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     73.64     121.80     157.47    194.96 148296.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.418    6.500 2741.100
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.143    3.000    4.000    4.278    5.000    5.500    6.000    7.000    8.000    2741.100

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8434  0.9600  0.9600
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.962   5.308   5.681   5.702   6.238   6.938

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8071  0.8277  0.8378  0.8434  0.8550  0.9499

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --out ${PREFIX}_mcSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' ${PREFIX}_mcSNP.filtered.recode.vcf > ${PREFIX}_mcSNP.filtered.renamed.vcf

## Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o mcpan_delly.filtered.vcf.gz mcpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats mcpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,2000)
summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     840.0     868.4    1200.0 10000.0
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02161  0.04000  0.24000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.004124  0.006186  0.009072  0.021608  0.012371  0.277938

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --out mcpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' mcpan_delly.filtered2.recode.vcf > mcpan_delly.filtered2.renamed.vcf
bcftools sort mcpan_delly.filtered2.renamed.vcf -Ov -o mcpan_delly.filtered2.sorted.vcf    



# VG Pangenome
# Mapping test data reads - when not working with simulated reads, use the pangenome files based on ${PREFIX}2-pg.vg, not ${PREFIX}2_aug.pg (i.e., without augmentation)
cd ${VGPAN}/
mkdir -p ./aligned/tmp
cd ./aligned/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../augmented/${PREFIX}2_aug_chopped.xg -g ../augmented/${PREFIX}2_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now) 
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

#cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam
cat *_aug_aln.gam > combined_aug_aln.gam

## Augmenting the graph with all variation from the GAM
${APP}/vg augment -t ${N} ../augmented/${PREFIX}2_aug_new.pg combined_aug_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_augSV.gam > ${PREFIX}2_augSV.pg 
#${APP}/vg augment -t ${N} ../augmented/${PREFIX}2_aug_new.pg combined_aug_aln.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_augSV.gam > ${PREFIX}2_augSV.pg 

# Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_augSV.pg > ${PREFIX}2_augSV_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}2_augSV_chopped.xg ${PREFIX}2_augSV_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_augSV_chopped.pg > ${PREFIX}2_augSV_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_augSV_chopped_pruned.gcsa ${PREFIX}2_augSV_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}2_augSV_chopped.pg -x ${PREFIX}2_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}2_augSV_chopped_new_L.xg -p > ${PREFIX}2_augSV_new.pg


# Variant calling (for SVs) for each individual using vg
cd ${VGPAN}/aligned/
mkdir -p ../called/sv/vg/

## Computing the snarls
${APP}/vg snarls -t ${N} ${PREFIX}2_augSV_new.pg > ${PREFIX}2_augSV.snarls

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}2_augSV_chopped.xg -g ./${PREFIX}2_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam 
  ${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}2_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam

## Computing the support
  ${APP}/vg pack -t ${N} -x ${PREFIX}2_augSV_new.pg -g ${i}_augSV_aln.filtered.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
  #${APP}/vg pack -t ${N} -x ${PREFIX}2_augSV_new.pg -g ${i}_augSV_aln.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg call -t ${N} ${PREFIX}2_augSV_new.pg -r ${PREFIX}2_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}_vgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called
done


# Filtering SV vcf files from vg
cd ${VGPAN}/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools

## Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ${i}_vgSV.vcf -Oz -o ${i}_vgSV.sorted.vcf.gz
  bcftools index ${i}_vgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_vgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

## Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}2_vgSV.merged.vcf.gz --threads ${N} *_vgSV.sorted.vcf.gz 

## Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}2_vgSV.merged.vcf.gz > vcf-stats.txt

## In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.       1st Qu.    Median     Mean     3rd Qu.   Max.
# -19.510     9.542     17.835     31.086    36.919 6701.480
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.000    0.600    1.662    2.040 2739.560
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.000    0.040    0.250    0.600    1.000    1.680    2.800    4.877    2739.560

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0400  0.2800   0.4412  0.9600  1.0000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4893   1.4052   1.6702   1.6089   1.8430   2.2308

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.3727  0.4326  0.4374  0.4412  0.4463  0.5403

quit()

## Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --out ${PREFIX}2_vgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Annotating SVs in vcf file
bcftools sort ${PREFIX}2_vgSV.filtered.recode.vcf -Oz -o ${PREFIX}2_vgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}2_vgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}2_vgSV.filtered.popped.vcf # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
PREFIX=short-read # needs to define again inside Apptainer
vcfwave -I 1000 -t ${N} ${PREFIX}2_vgSV.filtered.popped.vcf > ${PREFIX}2_vgSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
truvari anno svinfo -m 50 -o ${PREFIX}2_vgSV.decomposed.annotated.vcf ${PREFIX}2_vgSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python3 ${APP}/SVextractor_fromTruvariAnno.py 50 ${PREFIX}2_vgSV.decomposed.annotated.vcf ${PREFIX}2_vgSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)


# Variant calling (for SNPs using Sarek and for SVs using delly2) based on surjected bamfiles for each individual 
cd ${VGPAN}/aligned/
mkdir -p ../called/snp
#${APP}/vg gbwt -p -E -x ../augmented/${PREFIX}2-pg.vg --gbz-format -g ../augmented/${PREFIX}2-pg.gbz # make a gbz file 

module purge
module load biocontainers
module load nf-core
module load samtools
module load picard
module load bcftools
module load boost

## Projecting each sample's reads to the reference paths 
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam
  #${APP}/vg filter -t ${N} ${i}_aug_aln.gam -i -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_interleaved_aln.gam
  ${APP}/vg surject -x ../augmented/${PREFIX}2_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}2-pg_surject.bam
  #${APP}/vg surject -x ../augmented/${PREFIX}2_aug_chopped.xg ${i}_aug_interleaved_aln.gam --threads ${N} --prune-low-cplx --interleaved -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}2-pg_surject.bam
done

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### Indexing bam files
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}2-pg_surject.sorted.bam ${i}.${PREFIX}2-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}2-pg_surject.sorted.bam   
done

### Preprocessing input files for sarek
cat ${BASE}/original/SRR_Acc_List.txt > sample; \
#printf '1\n%.0s' {1..25} > lane ;  \
ls -1 ${VGPAN}/aligned/*surject.sorted.bam > bam; \
ls -1 ${VGPAN}/aligned/*surject.sorted.bam.bai > bai; \
echo "patient,sample,bam,bai" > header; \
paste sample sample bam bai | tr "\t" "," > content; \
cat header content > mapped.csv; \
rm header content sample bam bai # generate mapped.csv as an input to sarek

# Generating nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#    "input": "\/scratch\/negishi\/jeon96\/swallow\/vgpan\/aligned\/mapped.csv",
#    "step": "markduplicates",
#    "outdir": "\/scratch\/negishi\/jeon96\/swallow\/vgpan\/called\/snp\/",
#    "tools": "haplotypecaller",
#    "skip_tools": "baserecalibrator",
#    "trim_fastq": true,
#    "genome": "custom",
#    "fasta": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa",
#    "fasta_fai": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa.fai",
#    "igenomes_ignore": true,
#    "email": "jeon96@purdue.edu",
#    "multiqc_title": "short-read-vgpan_swallow"
#}

### Running sarek pipeline
nextflow run nf-core/sarek -r 3.3.2 -profile singularity -work-dir ${VGPAN}/aligned/ -params-file ${VGPAN}/aligned/nf-params.json #-resume 
#--step markduplicates --input ${VGPAN}/aligned/mapped.csv --outdir ${VGPAN}/called/snp --tools haplotypecaller --skip_tools baserecalibrator  

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
cd ${VGPAN}

cd ${VGPAN}/aligned/
mkdir -p ../called/sv/delly

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}2-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}2-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam
done

### Running delly2
cd ../called/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${VGPAN}/aligned/${i}.${PREFIX}2-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > vgpan_delly.vcf 
bcftools view -Ob -o vgpan_delly.bcf vgpan_delly.vcf
bcftools index vgpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o vgpan_delly.filtered.bcf vgpan_delly.bcf
bcftools view -Ov -o vgpan_delly.filtered.vcf vgpan_delly.filtered.bcf


# Filtering vcf files
cd ${VGPAN}/called/snp/

module purge
module load biocontainers
module load bcftools
module load vcftools

## Filtering SNP vcf files
### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_vgSNP.sorted.vcf.gz
  bcftools index ${i}_vgSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_vgSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o ${PREFIX}2_vgSNP.merged.vcf.gz --threads ${N} *_vgSNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}2_vgSNP.merged.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     75.64     121.84     159.29    194.96 138011.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.452    6.500 2599.000
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.250    3.000    4.000    4.333    5.000    5.500    6.000    7.000    8.000    2599.000

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8418  0.9600  0.9600
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.996   5.357   5.743   5.765   6.306   7.005

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8059  0.8259  0.8364  0.8418  0.8537  0.9488

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --out ${PREFIX}2_vgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o vgpan_delly.filtered.vcf.gz vgpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats vgpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,2000)
summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     780.0     809.2    1200.0 10000.0
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02266  0.04000  0.24000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.004234  0.006543  0.009623  0.022664  0.015397  0.283680

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --out vgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all



# NCBI linear representative genome
# Variant calling (for SNPs using Sarek and for SVs using delly2) by Andrew N. Black
cd ${REF}

## Indexing the lienar representative genome
samtools faidx GCF_015227805.2_bHirRus1.pri.v3_genomic.fna

## Creating dictionary
samtools dict GCF_015227805.2_bHirRus1.pri.v3_genomic.fna > GCF_015227805.2_bHirRus1.pri.v3_genomic.dict  

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
cd /scratch/negishi/blackan/forAB/cleaned_sra

### Preprocessing input files for sarek
ls -1 *val_1.fq.gz | cut -c 1-11 | uniq > sample; 
printf '1\n%.0s' {1..25} > lane;  
ls -1 *val_1.fq.gz > fastq_1;  
ls -1 *val_2.fq.gz > fastq_2; 
echo "patient,sample,lane,fastq_1,fastq_2" > header; 
paste sample sample lane fastq_1 fastq_2 | tr "\t" "," > tmp; 
cat header tmp > samplesheet.csv; 
rm header sample lane fastq_1 fastq_2

# Generating nf-params-pan.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#   "input": "\/scratch\/negishi\/blackan\/forAB\/cleaned_sra\/samplesheet.csv",
#   "outdir": "\/scratch\/negishi\/blackan\/forAB\/RESULTS\/pan\/",
#   "tools": "haplotypecaller",
#   "skip_tools": "baserecalibrator",
#   "trim_fastq": true,
#   "save_mapped": true,
#   "save_output_as_bam": true,
#   "genome": "custom",
#   "fasta": "\/scratch\/negishi\/blackan\/forAB\/ref\/GCF_015227805.2_bHirRus1.pri.v3_genomic.fna",
#   "fasta_fai": "\/scratch\/negishi\/blackan\/forAB\/ref\/GCF_015227805.2_bHirRus1.pri.v3_genomic.fna.fai",
#   "save_reference": true,
#   "igenomes_ignore": true,
#   "multiqc_title": "ref_swallow"
#}

### Running sarek pipeline
module load nextflow

cd /scratch/negishi/blackan/forAB/cleaned_sra
nextflow run $CLUSTER_SCRATCH/nf-core/nf-core-sarek-3.3.2/workflow/ \
-profile singularity -work-dir /scratch/negishi/blackan/forAB/RESULTS/pan/ -resume -params-file ./nf-params-pan.json 

### Compressing and indexing each vcf file first (By Andrew Black until here)
for i in `ls -1 *vcf.gz`; do   bcftools index $i -c -t ; done
bcftools merge *vcf.gz > ncbi_raw_snp.vcf

### Changing name and sorting  (By JongYoon Jeon from here)
cd ${REF}/called/

bcftools sort -Oz -o ncbi_SNP.sorted.vcf.gz ncbi_raw_snp.vcf

## Calling SVs using delly2
module load biocontainers
module load picard
module load bwa
module load boost

### Preprocessing the linear representative genome
cd ${REF}/
PicardCommandLine CreateSequenceDictionary --REFERENCE ${GENOME}.fna --OUTPUT ${GENOME}.dict
bwa index -a bwtsw ${GENOME}.fna

### Mapping for delly2
mkdir -p ./mapped/tmp/
cd ./mapped/
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
    bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${REF}/${GENOME}.fna  ${CLEANED_SRA}/${i}_1_val_1.fq ${CLEANED_SRA}/${i}_2_val_2.fq > ${i}.sam
    PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT ${i}.sam --OUTPUT ${i}_sorted.bam --SORT_ORDER coordinate
    PicardCommandLine MarkDuplicates --INPUT ${i}_sorted.bam --OUTPUT ${i}_sorted.marked.bam --METRICS_FILE ${i}_metrics.txt
    PicardCommandLine BuildBamIndex --INPUT ${i}_sorted.marked.bam
done

### Running delly2
cd ../
mkdir -p ./called/
cd ./called/
export OMP_NUM_THREADS=25 #equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
    samples+="${REF}/mapped/${i}_sorted.marked.bam "
done

${DEPOT}/apps/delly/src/delly call -g ${REF}/${GENOME}.fna ${samples} > ncbi_delly.vcf 
bcftools view -Ob -o ncbi_delly.bcf ncbi_delly.vcf
bcftools index ncbi_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o ncbi_delly.filtered.bcf ncbi_delly.bcf
bcftools view -Ov -o ncbi_delly.filtered.vcf ncbi_delly.filtered.bcf


# Filtering vcf files
## Filtering SNP vcf files
module purge
module load biocontainers
module load bcftools
module load vcftools

### Filtering with population-level parameters
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --missing-indv 
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --missing-site
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --depth
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --site-mean-depth
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ncbi_SNP.sorted.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     78.28     129.64      166.88    210.02 209853.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    4.000    5.368    6.082    7.000 5288.500
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.500    3.286    4.000    5.000    5.368   6.000    7.000    7.875    9.000    5288.500

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8425  0.9600  0.9600
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 3.180   5.885   6.360   6.398   7.079   7.834

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8096  0.8273  0.8369  0.8425  0.8530  0.9452

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --out ncbi_SNP.sorted.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Filtering SV vcf file from delly2
bcftools convert --thread ${N} -Oz -o ncbi_delly.filtered.vcf.gz ncbi_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --depth
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ncbi_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,2500)
summary(var_qual$qual)
#  Min.    1st Qu.  Median      Mean   3rd Qu. Max.
# 300.0     360.0     540.0     740.2    1015.5 10000.0
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu. Median    Mean   3rd Qu.  Max.
# 0.0000  0.0000  0.0000   0.0164  0.0400  0.2400
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#    Min.    1st Qu.   Median     Mean     3rd Qu.    Max.
# 0.001764  0.003919  0.006663  0.016390  0.009994  0.234764

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --out ncbi_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all



# Linear pangenome
# Variant calling (for SNPs using Sarek and for SVs using delly2) by Andrew N. Black
cd ${LINPAN}

## Indexing linear pangenome
samtools faidx short_panref2_sorted.fa

## Creating dictionary
samtools dict short_panref2_sorted.fa > short_panref2_sorted.dict 

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
cd /scratch/negishi/blackan/forAB/cleaned_sra

### Preprocessing input files for sarek
ls -1 *val_1.fq.gz | cut -c 1-11 | uniq > sample; 
printf '1\n%.0s' {1..25} > lane;  
ls -1 *val_1.fq.gz > fastq_1;  
ls -1 *val_2.fq.gz > fastq_2; 
echo "patient,sample,lane,fastq_1,fastq_2" > header; 
paste sample sample lane fastq_1 fastq_2 | tr "\t" "," > tmp; 
cat header tmp > samplesheet.csv; 
rm header sample lane fastq_1 fastq_2

# Generating nf-params-pan.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#   "input": "\/scratch\/negishi\/blackan\/forAB\/cleaned_sra\/samplesheet.csv",
#   "outdir": "\/scratch\/negishi\/blackan\/forAB\/RESULTS\/pan\/",
#   "tools": "haplotypecaller",
#   "skip_tools": "baserecalibrator",
#   "trim_fastq": true,
#   "save_mapped": true,
#   "save_output_as_bam": true,
#   "genome": "custom",
#   "fasta": "\/scratch\/negishi\/blackan\/forAB\/ref\/short_panref2_sorted.fa",
#   "fasta_fai": "\/scratch\/negishi\/blackan\/forAB\/ref\/short_panref2_sorted.fa.fai",
#   "save_reference": true,
#   "igenomes_ignore": true,
#   "multiqc_title": "ref_swallow"
#}

### Running sarek pipeline
module load nextflow

cd /scratch/negishi/blackan/forAB/cleaned_sra
nextflow run $CLUSTER_SCRATCH/nf-core/nf-core-sarek-3.3.2/workflow/ \
-profile singularity -work-dir /scratch/negishi/blackan/forAB/RESULTS/pan/ -resume -params-file ./nf-params-pan.json 

### Compressing and indexing each vcf file first (By Andrew Black until here)
for i in `ls -1 *vcf.gz`; do   bcftools index $i -c -t ; done
bcftools merge *vcf.gz > linpan_raw_snp.vcf

### Changing names and sorting (By JongYoon Jeon from here)
cd ${LINPAN}/called/

bcftools sort -Oz -o linpan_SNP.sorted.vcf.gz linpan_raw_snp.vcf

## Calling SVs using delly2
module load biocontainers
module load picard
module load bwa
module load boost

### Preprocessing the linear pangenome
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.renamed.fa --OUTPUT ${PREFIX}_panref2_sorted.renamed.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.renamed.fa

### Mapping for delly2
mkdir -p ./mapped/tmp/
cd ./mapped/
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${LINPAN}/${PREFIX}_panref2_sorted.fa  ${CLEANED_SRA}/${i}_1_val_1.fq ${CLEANED_SRA}/${i}_2_val_2.fq > ${i}.sam
  PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT ${i}.sam --OUTPUT ${i}_sorted.bam --SORT_ORDER coordinate
  PicardCommandLine MarkDuplicates --INPUT ${i}_sorted.bam --OUTPUT ${i}_sorted.marked.bam --METRICS_FILE ${i}_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}_sorted.marked.bam
done

### Running delly2
cd ../
mkdir -p ./called/
cd ./called/
export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${LINPAN}/mapped/${i}_sorted.marked.bam "
done

${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > linpan_delly.vcf 
bcftools view -Ob -o linpan_delly.bcf linpan_delly.vcf
bcftools index linpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o linpan_delly.filtered.bcf linpan_delly.bcf
bcftools view -Ov -o linpan_delly.filtered.vcf linpan_delly.filtered.bcf


# Filtering vcf files
## Filtering SNP vcf files
module purge
module load biocontainers
module load bcftools
module load vcftools

### Filtering with population-level parameters
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --missing-indv 
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --missing-site
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --depth
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --site-mean-depth
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats linpan_SNP.sorted.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     78.28     129.64      166.79    209.64 205949.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    4.000    5.333    6.076    7.000 5272.000
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.500    3.250    4.000    5.000    5.333   6.000    7.000    7.857    9.000    5272.000

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8427  0.9600  0.9600
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 3.178   5.878   6.361   6.394   7.065   7.833

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8098  0.8274  0.8370  0.8427  0.8531  0.9453

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                     
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --out linpan_SNP.sorted.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Filtering SV vcf files from delly2
bcftools convert --thread ${N} -Oz -o linpan_delly.filtered.vcf.gz linpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf linpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf linpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf linpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf linpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf linpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats linpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,2500)
summary(var_qual$qual)
# Min.   1st Qu. Median   Mean 3rd Qu. Max.
# 300     360     540     739    1009 10000
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.    1st Qu.  Median     Mean     3rd Qu.  Max.
# 0.00000  0.00000  0.00000   0.01652  0.04000  0.24000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#    Min.    1st Qu.   Median     Mean     3rd Qu.    Max.
# 0.002148  0.003905  0.006833  0.016525  0.009957  0.236236

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results 
vcftools --gzvcf linpan_delly.filtered.vcf.gz --out linpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all



# Original pangenome
# Variant calling (for SNPs using Sarek and for SVs using vg and delly2) by Natalie M. Allen
cd ${ORGPAN}

###-----Define frequently used variables-----### for original pangenome
N=128 # number of cores
BASE=/scratch/negishi/allen715/pangenome/forNA2 # e.g., /scratch/negishi/jeon96/swallow
APP=/home/allen715/vgtoolkit_1.53 # e.g., /home/jeon96/app
PREFIX=barnswallow # please define PREFIX as whatever you which
ORGPAN=/scratch/negishi/allen715/pangenome/forNA2/out # e.g., /scratch/negishi/jeon96/swallow/mcpan
REF=/scratch/negishi/allen715/pangenome/forNA2/ref # e.g., scratch/negishi/jeon96/swallow/original/ref 
CLEANED_SRA=/scratch/negishi/allen715/pangenome/forNA2/cleaned_sra # e.g., ${BASE}/original/sra/cleaned
DEPOT=/depot/fnrdewoody
RefInd=refp 
REF=refp
OUT=/scratch/negishi/allen715/pangenome/forNA2/out
REFDIR=/scratch/negishi/allen715/pangenome/forNA2/ref/
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed

# Mapping test data reads

cd ${ORGPAN}/
mkdir -p ./aligned/tmp
cd ./aligned/
 
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
 
## Aligning the individuals
${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ${ORGPAN}/${PREFIX}_mod_chopped.xg -g ${ORGPAN}/${PREFIX}_mod_chopped_pruned.gcsa > ${i}_org_aln.gam
 
## Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now)
${APP}/vg filter -t ${N} ${i}_org_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${ORGPAN}/${PREFIX}_mod_chopped.xg > ${i}_org_aln.filtered.gam
 
done
 
cat *_org_aln.filtered.gam > combined_org_aln.filtered.gam

cd ${ORGPAN}/aligned/
 
## Augmenting the graph with all variation from the GAM
${APP}/vg convert -t ${N} ${ORGPAN}/${PREFIX}_mod_chopped_new_L.xg -p > ${PREFIX}_org.pg
${APP}/vg augment -t ${N} ${ORGPAN}/${PREFIX}_org.pg ${ORGPAN}/aligned/combined_org_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_orgSV.gam > ${PREFIX}_orgSV.pg 
 
## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_orgSV.pg > ${PREFIX}_orgSV_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}_orgSV_chopped.xg ${PREFIX}_orgSV_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_orgSV_chopped.pg > ${PREFIX}_orgSV_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_orgSV_chopped_pruned.gcsa ${PREFIX}_orgSV_chopped_pruned.pg
 
## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}_orgSV_chopped.pg -x ${PREFIX}_orgSV_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}_orgSV_chopped_new_L.xg -p > ${PREFIX}_orgSV_new.pg


# Variant calling (for SVs) for each individual # use highmem queue
cd ${ORGPAN}/aligned/
mkdir -p ../called/sv
 
## Computing the snarls
${APP}/vg snarls -t ${N} ${PREFIX}_orgSV_new.pg > ${PREFIX}_orgSV.snarls
 
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
 
## Aligning the individuals
${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}_orgSV_chopped.xg -g ./${PREFIX}_orgSV_chopped_pruned.gcsa > ${i}_orgSV_aln.gam

## Filtering secondary and ambiguous read mappings out of the gam
${APP}/vg filter -t ${N} ${i}_orgSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}_orgSV_chopped.xg > ${i}_orgSV_aln.filtered.gam
 
## Computing the support
${APP}/vg pack -t ${N} -x ${PREFIX}_orgSV_new.pg -g ${i}_orgSV_aln.filtered.gam -Q 5 -o ${i}_orgSV.pack #-Q 5: ignore mapping and base qualitiy < 5
done

## Calling variants; run this step using highmem queue (otherwise it can't be finished)
cd ${ORGPAN}/aligned/
for i in `cat ${BASE}/SRR_Acc_List_new.txt`; do
  ${APP}/vg call -t ${N} ${PREFIX}_orgSV_new.pg -r ${PREFIX}_orgSV.snarls -k ${i}_orgSV.pack -s ${i} -a -A -c 50 -C 100000 > ../called/sv/${i}_orgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called
done

## Combining separately called SV vcf files
cd ${ORGPAN}/called/sv
 
module purge
module load biocontainers
module load bcftools
module load vcftools

## Compressing and indexing each vcf file first #before running this, need to change vcf file headers
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
  sed 's/refp#0#//g' ${i}_orgSV.vcf > ${i}_orgSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_orgSV2.vcf -Oz -o ${i}_orgSV.sorted.vcf.gz
  bcftools index ${i}_orgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_orgSV2.vcf
done

bcftools merge -m all -Oz -o ${PREFIX}_orgSV.merged.vcf.gz --threads ${N} *_orgSV.sorted.vcf.gz 


# Variant calling (for SNPs using Sarek and for SVs using delly2) based on surjected bamfiles for each individual
cd ${ORGPAN}/aligned/
mkdir -p ../called/snp

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
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
  ${APP}/vg filter -t ${N} ${i}_org_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../${PREFIX}_mod_chopped.xg > ${i}_org_interleaved_aln.filtered.gam
  ${APP}/vg surject -x ../${PREFIX}_mod_chopped.xg ${i}_org_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam
done

## Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### Indexing bam files
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}-pg_surject.sorted.bam ${i}.${PREFIX}-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}-pg_surject.sorted.bam   
done

### Preprocessing the reference genome and bam files
module purge
module load biocontainers
module load gatk4
module load picard
module load samtools
module load bcftools
module load boost

cd ${REF}
sed 's/>/>refp#0#/g' GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed.fasta > ref.renamed.sorted.fa
samtools faidx ref.renamed.sorted.fa
PicardCommandLine CreateSequenceDictionary --REFERENCE ref.renamed.sorted.fa --OUTPUT ref.renamed.sorted.dict
bwa index -a bwtsw ref.renamed.sorted.fa

cd ${ORGPAN}/aligned/

### GATK HaplotypeCaller
 
for i in `cat ${BASE}/SRR_Acc_List_new.txt`; do
  mkdir -p /scratch/negishi/allen715/pangenome/forNA2/out/called/snp/variant_calling/${i}/
  gatk  --java-options "-Xmx100G -XX:ParallelGCThreads=128" HaplotypeCaller -R ${REF}/ref.renamed.sorted.fa -I /scratch/negishi/allen715/pangenome/forNA2/out/called/snp/preprocessing/markduplicates/${i}/${i}.md.cram  -O /scratch/negishi/allen715/pangenome/forNA2/out/called/snp/variant_calling/${i}/${i}.haplotypecaller.vcf.gz --tmp-dir .
done # use sorted file instead of just renamed 

## Calling SVs using delly2
### Preprocessing the reference genome and bam files
cd ${REF}
PicardCommandLine CreateSequenceDictionary --REFERENCE ref.renamed.fa --OUTPUT ref.renamed.dict
bwa index -a bwtsw ref.renamed.fa

cd ${ORGPAN}/aligned/
mkdir -p ../called/sv/delly
 
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam
done
 
### Running delly2
cd ../called/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/SRR_Acc_List.txt`; do
  samples+="${ORGPAN}/aligned/${i}.${PREFIX}-pg_surject.sorted.marked.bam "
done 

${DEPOT}/apps/delly/src/delly call -g ${REF}/ref.renamed.sorted.fa ${samples} > orgpan_delly.vcf # used ref.renamed.sorted.fa not ref.renamed.fa 
bcftools view -Oz -o  orgpan_delly.bcf --threads ${N} orgpan_delly.vcf
bcftools index orgpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o orgpan_delly.filtered.bcf orgpan_delly.bcf
bcftools view -Ov -o orgpan_delly.filtered.vcf orgpan_delly.filtered.bcf 


# Filtering vcf files (By JongYoon Jeon from here)
cd ${ORGPAN}/called/snp/

module purge
module load biocontainers
module load bcftools
module load vcftools

## Filtering SNP vcf files
### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_orgSNP.sorted.vcf.gz
  bcftools index ${i}_orgSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining separately called SNP vcf files
bcftools merge -m all -Oz -o barnswallow_orgSNP.merged.vcf.gz --threads ${N} *_orgSNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSNP.merged.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     76.84     121.84     157.44    194.96 169562.00
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.383    6.500 3814.730
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.333    3.000    4.000    4.333    5.000    5.500    6.000    7.000    8.000    3814.730

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min. 1st Qu.  Median  Mean 3rd Qu.  Max.
# 0.000  0.840  0.920   0.843  0.960  0.960
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.882   5.122   5.530   5.532   6.003   6.757

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8071  0.8272  0.8374  0.8430  0.8550  0.9494

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --out barnswallow_orgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' barnswallow_orgSNP.filtered.recode.vcf > barnswallow_orgSNP.filtered.renamed.vcf # polish the header line to be compatible with bcftools

## Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o orgpan_delly.filtered.vcf.gz orgpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats orgpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,2000)
summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     780.0     804.1    1200.0 5160.0
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02191  0.04000  0.24000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.002957  0.005915  0.009612  0.021915  0.014418  0.280591

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results 
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --out orgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' orgpan_delly.filtered2.recode.vcf > orgpan_delly.filtered2.renamed.vcf # polish the header line to be compatible with bcftools

## Filtering SV vcf files from vg 
cd ${ORGPAN}/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools
module load samtools
module load nf-core
module load bbmap
module load picard
module load bwa
module load boost

## Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  sed 's/refp#0#//g' ${i}_orgSV.vcf > ${i}_orgSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_orgSV2.vcf -Oz -o ${i}_orgSV.sorted.vcf.gz
  bcftools index ${i}_orgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_orgSV2.vcf
done

## Combining separately called SV vcf files
bcftools merge -m all -Oz -o barnswallow_orgSV.merged.vcf.gz --threads ${N} *_orgSV.sorted.vcf.gz 

## Filtering with population-level parameters 
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSV.merged.vcf.gz > vcf-stats.txt

## In R, plotting and summarizing vcf stats
module load r

R
library(tidyverse)
library(ggplot2)
var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0,500)
summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean   3rd Qu.   Max.
# -21.65     12.71     26.29     42.38    46.74 31367.90
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.160    0.875    1.787    2.312 3435.120
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.080    0.267    0.520    0.875    1.250    1.909    2.882    4.750    3435.120

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0000  0.0800   0.2787  0.5200  1.0000
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4626   1.4524   1.6978   1.6670   1.8723   2.2785

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.2621  0.2685  0.2724  0.2787  0.2824  0.3828

quit()

## Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --out barnswallow_orgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Annotating SVs in vcf file
bcftools sort barnswallow_orgSV.filtered.recode.vcf -Oz -o barnswallow_orgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input barnswallow_orgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > barnswallow_orgSV.filtered.popped.vcf # remove large (>100 kb) alleles in MC raphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
vcfwave -I 1000 -t ${N} barnswallow_orgSV.filtered.popped.vcf > barnswallow_orgSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

truvari anno svinfo -m 50 -o barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python3 ${APP}/SVextractor_fromTruvariAnno.py 50 barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)