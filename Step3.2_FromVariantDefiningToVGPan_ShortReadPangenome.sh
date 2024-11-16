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

# Constructing a VG pangenome
module load biocontainers
module load bwa
module load picard
module load boost
module load bcftools
module load htslib

cd ${BASE}
mkdir -p ${VGPAN}/augmented/called
cd ${VGPAN}/augmented

## Running sarek with simulated samples to augment the linear pangenome
### Preparing an input sheet
ls -lrt /scratch/negishi/jeon96/swallow/sim/raw/ | tr -s ' ' | cut -d ' ' -f 9 | grep "primary" | cut -c 1-24 | uniq | sort > primary; \
ls -lrt /scratch/negishi/jeon96/swallow/sim/raw/ | tr -s ' ' | cut -d ' ' -f 9 | grep "alternate" | cut -c 1-26 | uniq | sort > alternate; \
cat primary alternate > sample; \
printf '1\n%.0s' {1..40} > lane ;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*primary_R1.fq.gz | sort > primary_fastq_1;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*alternate_R1.fq.gz | sort > alternate_fastq_1;  \
cat primary_fastq_1 alternate_fastq_1 > fastq_1; \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*primary_R2.fq.gz | sort > primary_fastq_2;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*alternate_R2.fq.gz | sort > alternate_fastq_2;  \
cat primary_fastq_2 alternate_fastq_2 > fastq_2; \
echo "patient,sample,lane,fastq_1,fastq_2" > header; \
paste sample sample lane fastq_1 fastq_2 | tr "\t" "," > content; \
cat header content > samplesheet.csv; \
rm header content sample lane fastq_1 fastq_2 primary alternate primary_fastq_1 primary_fastq_2 alternate_fastq_1 alternate_fastq_2 # generate samplesheet.csv as an input to sarek

### Generating nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
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
#    "multiqc_title": "short-read-pan_swallow"
#}

## Calling SVs using delly2
### Mapping for delly2
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.fa --OUTPUT ${PREFIX}_panref2_sorted.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.fa

cd ${VGPAN}/augmented
mkdir -p ./mapped/tmp/
cd ./mapped/

for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${LINPAN}/${PREFIX}_panref2_sorted.fa  ${SIM}/raw/sim${k}_${i}_${j}_R1.fq.gz ${SIM}/raw/sim${k}_${i}_${j}_R2.fq.gz > /sim${k}_${i}_${j}.sam
    PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT sim${k}_${i}_${j}.sam --OUTPUT sim${k}_${i}_${j}_sorted.bam --SORT_ORDER coordinate
    PicardCommandLine MarkDuplicates --INPUT sim${k}_${i}_${j}_sorted.bam --OUTPUT ./sim${k}_${i}_${j}_sorted.marked.bam --METRICS_FILE sim${k}_${i}_${j}_metrics.txt
    PicardCommandLine BuildBamIndex --INPUT ./sim${k}_${i}_${j}_sorted.marked.bam
    done 
  done
done

### =Running delly2
cd ${VGPAN}/augmented/called/  
mkdir -p ./delly/
cd ./delly/

export OMP_NUM_THREADS=40 # equal or less than the number of samples
samples=""
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    samples+="${VGPAN}/augmented/mapped/sim${k}_${i}_${j}_sorted.marked.bam "
    done
  done
done

${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > sim_delly.vcf 
bcftools view -Ob -o sim_delly.bcf sim_delly.vcf
bcftools index sim_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o sim_delly.filtered.bcf sim_delly.bcf
bcftools view -Ov -o sim_delly.filtered.vcf sim_delly.filtered.bcf 

## Filtering variants (SNPs)
cd ../

### Compressing and indexing each vcf file first
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    bcftools sort ./variant_calling/haplotypecaller/sim${k}_${i}_${j}/sim${k}_${i}_${j}.haplotypecaller.vcf.gz -Oz -o sim${k}_${i}_${j}_SNP.sorted.vcf.gz
    bcftools index sim${k}_${i}_${j}_SNP.sorted.vcf.gz --threads ${N}
    bcftools view sim${k}_${i}_${j}_SNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
    done
  done
done

### Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o sim_SNP.merged.vcf.gz --threads ${N} *_SNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf sim_SNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf sim_SNP.merged.vcf.gz --missing-site
vcftools --gzvcf sim_SNP.merged.vcf.gz --depth
vcftools --gzvcf sim_SNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf sim_SNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sim_SNP.merged.vcf.gz > vcf-stats.txt

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
#  Min.    1st Qu.   Median    Mean     3rd Qu.  Max.
# 30.0     301.0     426.1     449.6    570.1 81304.1
   
var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()   
a + theme_light() + xlim(0,100)        
summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    5.583    7.375    7.971    9.250 6713.750
quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    4.250    5.250    6.000    6.750    7.375    8.000    8.750    9.875    11.750    6713.750

var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.7000  0.9000   0.7685  0.9000  0.9750
                       
ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 5.297   5.318   5.985   7.657   9.841   11.966

ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.7248  0.7286  0.7526  0.7685  0.7984  0.8468

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                     
vcftools --gzvcf sim_SNP.merged.vcf.gz --out sim_SNP.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Filtering variants (SVs)
cd ./delly/
bcftools convert --thread ${N} -Oz -o sim_delly.filtered.vcf.gz sim_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf sim_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf sim_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf sim_delly.filtered.vcf.gz --depth
vcftools --gzvcf sim_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf sim_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sim_delly.filtered.vcf.gz > vcf-stats.txt

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
#  Min.  1st Qu.  Median   Mean   3rd Qu. Max.
# 300     647     1020     1052    1200 10000
   
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
# 0.0000  0.0250  0.1000   0.1026  0.2000  0.2500
                       
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
#   Min.   1st Qu.  Median     Mean    3rd Qu.   Max.
# 0.01619  0.03147  0.05195  0.10262  0.13422  0.32350

quit()

### Filtering vcf file based on the determined conservative cut-off values from above results 
vcftools --gzvcf sim_delly.filtered.vcf.gz --out ../sim_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Constructing a VG pangenome
cd ../

### Preprocessing input files
bcftools convert --thread ${N} -Ov -o sim_delly.filtered2.vcf sim_delly.filtered2.recode.vcf.gz
while read -r line; do
    pair=($line)
    echo "s/${pair[0]}/${pair[1]}/" >> delly_replace.sed
done < delly_sample.list # make a sed script from the list of "search-replace" pairs of "wrong-correct" sample name pairs
sed -f delly_replace.sed sim_delly.filtered2.vcf > sim_delly.filtered2.replaced.vcf # fix sample names

while read -r line; do
    pair=($line)
    echo "s/${pair[0]}/${pair[1]}/" >> snp_replace.sed
done < snp_sample.list # make a sed script from the list of "search-replace" pairs of "wrong-correct" sample name pairs
sed -f snp_replace.sed sim_SNP.filtered.recode.vcf > sim_SNP.filtered.replaced.vcf # fix sample names

bcftools query -l sim_delly.filtered2.replaced.vcf > samples.txt # make a sample list
bcftools view -S samples.txt sim_delly.filtered2.replaced.vcf > sim_delly.filtered2.ordered.vcf # order samples; check the sample header! bcftools view sim_delly.filtered2.ordered.vcf | grep -v "##" | head
cat sim_delly.filtered2.ordered.vcf | grep -v "SVTYPE=BND" | grep -v "SVTYPE=DUP" > sim_delly.filtered3.ordered.vcf # remove unsupported SV types in vg construct (i.e., DUP, BND)
bcftools view -S samples.txt sim_SNP.filtered.replaced.vcf > sim_SNP.filtered.ordered.vcf # order samples; check the sample header! bcftools view sim_SNP.filtered.ordered.vcf | grep -v "##" | head

bcftools convert --thread ${N} -Oz -o sim_SNP.filtered.vcf.gz sim_SNP.filtered.ordered.vcf
bcftools convert --thread ${N} -Oz -o sim_delly.filtered3.vcf.gz sim_delly.filtered3.ordered.vcf
tabix -p vcf sim_SNP.filtered.vcf.gz
tabix -p vcf sim_delly.filtered3.vcf.gz
bcftools concat --allow-overlaps -Oz -o sim_SNP_SV.vcf.gz --threads ${N} sim_SNP.filtered.vcf.gz sim_delly.filtered3.vcf.gz # combine SNP and SV vcf files

bcftools norm -a --atom-overlaps . -m -any --threads ${N} -Oz -o sim_SNP_SV.normed.vcf.gz sim_SNP_SV.vcf.gz # normalize the vcf file      
tabix -p vcf sim_SNP_SV.normed.vcf.gz

### Constructing a graph pangenome using vg and the linear pangenome
${APP}/vg construct -a -f -S -r ${LINPAN}/${PREFIX}_panref2_sorted.fa -v sim_SNP_SV.normed.vcf.gz -m 1000000 > ../${PREFIX}2-pg.vg #| ${APP}/vg convert - | ${APP}/vg mod -x 32 -  # too big inversions exist for Protobuf to stream -> moving node chopping outside of vg construct

${APP}/vg stats -lz ${PREFIX}2-pg.vg # check pangenome basic statistics
#nodes   61690
#edges   88747
#length  1123193383


# Augmenting the graph pangenome (preprocessing steps before alignmening were modified from Sacomandi et al. (2023)'s scripts)
cd ${VGPAN}/

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./augmented/
mkdir -p ./tmp/

${APP}/vg mod -t ${N} -O ${PREFIX}2-pg.vg > ${PREFIX}2_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}2_mod.xg ${PREFIX}2_mod.vg

## Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_mod.vg > ${PREFIX}2_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}2_mod_chopped.xg ${PREFIX}2_mod_chopped.vg

## Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_mod_chopped.vg > ${PREFIX}2_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_mod_chopped_pruned.gcsa ${PREFIX}2_mod_chopped_pruned.vg

## Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}2_mod_chopped.vg -x ${PREFIX}2_mod_chopped_new_L.xg #-L: preserve alt paths in the xg

## Aligning the individuals that were used to build the pangenome; do augmentation steps when reads that were used to build the pangenome are different from reads that will be used to call variants like this case.  
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for k in {1..4}; do 
    cat ${SIM}/raw/sim${k}_${i}_primary_R1.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R1.fq.gz >  ${SIM}/raw/sim${k}_${i}_R1.fq.gz 
    cat ${SIM}/raw/sim${k}_${i}_primary_R2.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R2.fq.gz >  ${SIM}/raw/sim${k}_${i}_R2.fq.gz
    ${APP}/vg map -t ${N} -f ${SIM}/raw/sim${k}_${i}_R1.fq.gz -f ${SIM}/raw/sim${k}_${i}_R2.fq.gz -x ${PREFIX}2_mod_chopped.xg -g ${PREFIX}2_mod_chopped_pruned.gcsa > sim${k}_${i}_aln.gam

## Filtering secondary and ambiguous read mappings out of the GAM of the above step
    ${APP}/vg filter -t ${N} sim${k}_${i}_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${PREFIX}2_mod_chopped.xg > sim${k}_${i}_aln.filtered.gam #-r : minimum score to keep primary alignment; -f: normalize score based on length; -u: use substitution count instead of score; -m: filter reads that don't begin with at least N matches on each end; -q: filter alignments with mapping quality < N; -D: clip back the ends of reads that are ambiguously aligned, up to N bases
  done
done

cat *filtered.gam > combined_filtered.gam 

## Augmenting the graph with all variation from the GAM of the above step
${APP}/vg convert -t ${N} ${PREFIX}2_mod_chopped_new_L.xg -p > ${PREFIX}2.pg
${APP}/vg augment -t ${N} ${PREFIX}2.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_aug.gam > ${PREFIX}2_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 

## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_aug.pg > ${PREFIX}2_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}2_aug_chopped.xg ${PREFIX}2_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_aug_chopped.pg > ${PREFIX}2_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_aug_chopped_pruned.gcsa ${PREFIX}2_aug_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}2_aug_chopped.pg -x ${PREFIX}2_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}2_aug_chopped_new_L.xg -p > ${PREFIX}2_aug_new.pg