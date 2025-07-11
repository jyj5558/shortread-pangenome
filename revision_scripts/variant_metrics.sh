#!/bin/bash
#SBATCH --job-name=variant_metrics
#SBATCH -A highmem
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=jeon96@purdue.edu
#SBATCH --mail-type=END,FAIL

###-----Define frequently used variables-----###
N=64 # number of cores
PREFIX=short-read
PREFIX2=short-swallow-out
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
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic
GENOME_SIZE=1.1g
POPGEN=/scratch/negishi/jeon96/swallow/popgen
BENCHMARK=/scratch/negishi/jeon96/swallow/benchmarking

module purge
module load biocontainers
module load bcftools
module load vcftools
module load samtools

# Minigraph-Cactus pangenome
cd ${MCPAN}/${PREFIX}-pg/aligned_re/

## Mapping metrics
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samtools flagstat -@ ${N} ${i}.${PREFIX}-pg_surject.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
  samtools view -@ ${N} ${i}.${PREFIX}-pg_surject.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
  cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
  cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 96.3364
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 48.9911

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${APP}/vg stats --threads ${N} --alignments ${i}_aug_aln.gam ../${PREFIX}_aug_chopped.xg > ${i}_mapping.vgstats
  total_alignment=`cat ${i}_mapping.vgstats | grep "Total alignments:" | cut -d ":" -f 2`
  total_aligned=`cat ${i}_mapping.vgstats | grep "Total aligned:" | cut -d ":" -f 2`
  printf %.3f "$((10**3 * 100 * ${total_aligned}/${total_alignment}))e-3" >> vgmapping.rate
  echo "" >> vgmapping.rate
  cat ${i}_mapping.vgstats | grep "Mapping quality:" | cut -d ":" -f 2 | cut -d "," -f 1 | cut -d " " -f 3 >> vgmap.q
done

awk '{sum += $1} END {print sum/NR}' vgmapping.rate > vgmapping_rate.average # 97.2451
awk '{sum += $1} END {print sum/NR}' vgmap.q > vgmapq.average # 50.5876

## Number of variants
cd ${MCPAN}/${PREFIX}-pg/called_re/sv/vg/

### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  sed 's/bHirRus1_LinPan#0#//g' ${i}_mcSV.vcf > ${i}_mcSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_mcSV2.vcf -Oz -o ${i}_mcSV.sorted.vcf.gz
  bcftools index ${i}_mcSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_mcSV2.vcf
done

### Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSV.merged.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --out ${PREFIX}_mcSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Annotating SVs in vcf file
bcftools sort ${PREFIX}_mcSV.filtered.recode.vcf -Oz -o ${PREFIX}_mcSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}_mcSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}_mcSV.filtered.popped.vcf # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree.

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

### Filtering SNP vcf files
cd ${MCPAN}/${PREFIX}-pg/called_re/snp

### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_mcSNP.sorted.vcf.gz
  bcftools index ${i}_mcSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining aseparately called SNP vcf files
bcftools merge -m all -Ov -o ${PREFIX}_mcSNP.merged.vcf --threads ${N} *_mcSNP.sorted.vcf.gz 
sed 's/bHirRus1_LinPan#0#//g' ${PREFIX}_mcSNP.merged.vcf > ${PREFIX}_mcSNP.renamed.vcf
bcftools convert --threads ${N} -Oz -o ${PREFIX}_mcSNP.merged.vcf.gz ${PREFIX}_mcSNP.merged.vcf

### Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSNP.merged.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --out ${PREFIX}_mcSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' ${PREFIX}_mcSNP.filtered.recode.vcf > ${PREFIX}_mcSNP.filtered.renamed.vcf

### Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o mcpan_delly.filtered.vcf.gz mcpan_delly.filtered.vcf
sed 's/bHirRus1_LinPan#0#//g' mcpan_delly.filtered.vcf > mcpan_delly.filtered.renamed.vcf

### Filtering with population-level parameters
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats mcpan_delly.filtered.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --out mcpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' mcpan_delly.filtered2.recode.vcf > mcpan_delly.filtered2.renamed.vcf
bcftools sort mcpan_delly.filtered2.renamed.vcf -Ov -o mcpan_delly.filtered2.sorted.vcf    

cd ../../

bcftools stats ./snp/${PREFIX}_mcSNP.renamed.vcf > SNP_nofiltered.count # 45449268
bcftools stats ./snp/${PREFIX}_mcSNP.filtered.renamed.vcf > SNP_filtered.count # 196736
bcftools stats ./sv/vg/${PREFIX}_mcSV.merged.vcf.gz > vgSV_nofiltered.count # 146140
bcftools stats ./sv/vg/${PREFIX}_mcSV.decomp.annot.ext.vcf > vgSV_filtered.count # 211
bcftools stats ./sv/delly/mcpan_delly.filtered.renamed.vcf > dellySV_nofiltered.count # 3300
bcftools stats ./sv/delly/mcpan_delly.filtered2.renamed.vcf > dellySV_filtered.count # 2772


# VG pangenome
cd ${VGPAN}/aligned_re/

## Mapping metrics
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samtools flagstat -@ ${N} ${i}.${PREFIX}2-pg_surject.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
  samtools view -@ ${N} ${i}.${PREFIX}2-pg_surject.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
  cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
  cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 97.7728
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 49.6909

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${APP}/vg stats --threads ${N} --alignments ${i}_aug_aln.gam ../augmented_re/${PREFIX}2_aug_chopped.xg > ${i}_mapping.vgstats
  total_alignment=`cat ${i}_mapping.vgstats | grep "Total alignments:" | cut -d ":" -f 2`
  total_aligned=`cat ${i}_mapping.vgstats | grep "Total aligned:" | cut -d ":" -f 2`
  printf %.3f "$((10**3 * 100 * ${total_aligned}/${total_alignment}))e-3" >> vgmapping.rate
  echo "" >> vgmapping.rate
cat ${i}_mapping.vgstats | grep "Mapping quality:" | cut -d ":" -f 2 | cut -d "," -f 1 | cut -d " " -f 3 >> vgmap.q
done

awk '{sum += $1} END {print sum/NR}' vgmapping.rate > vgmapping_rate.average # 97.7912
awk '{sum += $1} END {print sum/NR}' vgmap.q > vgmapq.average # 50.8172

## Number of variants
cd ${VGPAN}/called_re/sv/vg/

### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools sort ${i}_vgSV.vcf -Oz -o ${i}_vgSV.sorted.vcf.gz
  bcftools index ${i}_vgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_vgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}2_vgSV.merged.vcf.gz --threads ${N} *_vgSV.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}2_vgSV.merged.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --out ${PREFIX}2_vgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Annotating SVs in vcf file
bcftools sort ${PREFIX}2_vgSV.filtered.recode.vcf -Oz -o ${PREFIX}2_vgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}2_vgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}2_vgSV.filtered.popped.vcf # remove large (>100 kb) alleles in graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

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

### Filtering SNP vcf files
cd ${VGPAN}/called_re/snp/

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

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --out ${PREFIX}2_vgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Filtering SV vcf files from delly2
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

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --out vgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

cd ../../

bcftools stats ./snp/${PREFIX}2_vgSNP.merged.vcf.gz > SNP_nofiltered.count # 44442867
bcftools stats ./snp/${PREFIX}2_vgSNP.filtered.recode.vcf > SNP_filtered.count # 158787
bcftools stats ./sv/vg/${PREFIX}2_vgSV.merged.vcf.gz > vgSV_nofiltered.count # 114362
bcftools stats ./sv/vg/${PREFIX}2_vgSV.decomp.annot.ext.vcf > vgSV_filtered.count # 173
bcftools stats ./sv/delly/vgpan_delly.filtered.vcf > dellySV_nofiltered.count # 3321
bcftools stats ./sv/delly/vgpan_delly.filtered2.recode.vcf > dellySV_filtered.count # 2778


# Original pangenome
cd ${ORGPAN}/aligned_re

## Mapping metrics
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  samtools flagstat -@ ${N} ${i}.${PREFIX2}-pg_surject.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
  samtools view -@ ${N} ${i}.${PREFIX2}-pg_surject.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
  cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
  cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 90.1204
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 48.2891

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${APP}/vg stats --threads ${N} --alignments ${i}_org_aln.gam ../${PREFIX2}_mod_chopped.xg > ${i}_mapping.vgstats
  total_alignment=`cat ${i}_mapping.vgstats | grep "Total alignments:" | cut -d ":" -f 2`
  total_aligned=`cat ${i}_mapping.vgstats | grep "Total aligned:" | cut -d ":" -f 2`
  printf %.3f "$((10**3 * 100 * ${total_aligned}/${total_alignment}))e-3" >> vgmapping.rate
  echo "" >> vgmapping.rate
cat ${i}_mapping.vgstats | grep "Mapping quality:" | cut -d ":" -f 2 | cut -d "," -f 1 | cut -d " " -f 3 >> vgmap.q
done

awk '{sum += $1} END {print sum/NR}' vgmapping.rate > vgmapping_rate.average # 99.0533
awk '{sum += $1} END {print sum/NR}' vgmap.q > vgmapq.average # 49.9777
 
## Number of variants
cd ${ORGPAN}/called_re/snp/

### Filtering SNP vcf files
### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools sort ${i}.haplotypecaller.vcf.gz -Oz -o ${i}_orgSNP.sorted.vcf.gz
  bcftools index ${i}_orgSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining separately called SNP vcf files
bcftools merge -m all -Ov -o barnswallow_orgSNP.merged.vcf --threads ${N} *_orgSNP.sorted.vcf.gz 
sed 's/refp#0#//g' barnswallow_orgSNP.merged.vcf > barnswallow_orgSNP.renamed.vcf
bcftools convert --threads ${N} -Oz -o barnswallow_orgSNP.merged.vcf.gz barnswallow_orgSNP.merged.vcf

### Filtering with population-level parameters
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSNP.merged.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --out barnswallow_orgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' barnswallow_orgSNP.filtered.recode.vcf > barnswallow_orgSNP.filtered.renamed.vcf # polish the header line to be compatible with bcftools

### Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o orgpan_delly.filtered.vcf.gz orgpan_delly.filtered.vcf
sed 's/refp#0#//g' orgpan_delly.filtered.vcf > orgpan_delly.filtered.renamed.vcf

### Filtering with population-level parameters
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats orgpan_delly.filtered.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results 
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --out orgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' orgpan_delly.filtered2.recode.vcf > orgpan_delly.filtered2.renamed.vcf # polish the header line to be compatible with bcftools

### Filtering SV vcf files from vg 
cd ../vg/

### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  sed 's/refp#0#//g' ${i}_orgSV.vcf > ${i}_orgSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_orgSV2.vcf -Oz -o ${i}_orgSV.sorted.vcf.gz
  bcftools index ${i}_orgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_orgSV2.vcf
done

### Combining separately called SV vcf files
bcftools merge -m all -Oz -o barnswallow_orgSV.merged.vcf.gz --threads ${N} *_orgSV.sorted.vcf.gz 

### Filtering with population-level parameters 
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSV.merged.vcf.gz > vcf-stats.txt

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                      
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --out barnswallow_orgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Annotating SVs in vcf file
bcftools sort barnswallow_orgSV.filtered.recode.vcf -Oz -o barnswallow_orgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input barnswallow_orgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > barnswallow_orgSV.filtered.popped.vcf # remove large (>100 kb) alleles in MC raphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
vcfwave -I 1000 -t ${N} barnswallow_orgSV.filtered.popped.vcf > barnswallow_orgSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

truvari anno svinfo -m 50 -o barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python3 ${APP}/SVextractor_fromTruvariAnno.py 50 barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)

cd ../../

bcftools stats ./snp/barnswallow_orgSNP.merged.vcf.gz > SNP_nofiltered.count # 44686859
bcftools stats ./snp/barnswallow_orgSNP.filtered.renamed.vcf > SNP_filtered.count # 197589
bcftools stats ./sv/vg/barnswallow_orgSV.merged.vcf.gz > vgSV_nofiltered.count # 172398
bcftools stats ./sv/vg/barnswallow_orgSV.decomp.annot.ext.vcf > vgSV_filtered.count # 517
bcftools stats ./sv/delly/orgpan_delly.filtered.renamed.vcf > dellySV_nofiltered.count # 3497
bcftools stats ./sv/delly/orgpan_delly.filtered2.renamed.vcf > dellySV_filtered.count # 2899



# Variant Benchmarking (based on pre-filtering vcf files)
module purge
module load biocontainers
module load bcftools
module load picard 
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bbmap
module load samtools

cd ${BASE}
mkdir -p ./benchmarking_re/
cd ./benchmarking_re/

## Linear representative genome vs. Linear pangenome
mkdir -p ./ncbi_linpan/
cd ./ncbi_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
bcftools reheader --samples ${BASE}/original/original/SRR_Acc_List.txt ${REF}/called/ncbi_SNP.sorted.filtered.recode.vcf > ncbi_SNP_final_postfilter.vcf # polish sample names in the header
bcftools reheader --samples ${BASE}/original/original/SRR_Acc_List.txt ${REF}/called/ncbi_SNP.sorted.vcf.gz > ncbi_SNP_final_prefilter.vcf # polish sample names in the header
cp ${REF}/called/ncbi_delly.filtered2.recode.vcf ./ncbi_delly_final_postfilter.vcf
cp ${REF}/called/ncbi_delly.filtered.vcf ./ncbi_delly_final_prefilter.vcf

bcftools reheader --samples ${BASE}/original/original/SRR_Acc_List.txt ${LINPAN}/called/linpan_SNP.sorted.filtered.recode.vcf > linpan_SNP_final_postfilter.vcf # polish sample names in the header
bcftools reheader --samples ${BASE}/original/original/SRR_Acc_List.txt ${LINPAN}/called/linpan_SNP.sorted.vcf.gz > linpan_SNP_final_prefilter.vcf # polish sample names in the header
cp ${LINPAN}/called/linpan_delly.filtered2.recode.vcf ./linpan_delly_final_postfilter.vcf
cp ${LINPAN}/called/linpan_delly.filtered.vcf ./linpan_delly_final_prefilter.vcf

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools view -a -I -s ${i} --thread ${N} ncbi_SNP_final_prefilter.vcf > ncbi_SNP_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} ncbi_delly_final_prefilter.vcf > ncbi_delly_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} linpan_SNP_final_prefilter.vcf > linpan_SNP_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} linpan_delly_final_prefilter.vcf > linpan_delly_final_prefilter_${i}.vcf
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools view -a -I -s ${i} --thread ${N} ncbi_SNP_final_postfilter.vcf > ncbi_SNP_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} ncbi_delly_final_postfilter.vcf > ncbi_delly_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} linpan_SNP_final_postfilter.vcf > linpan_SNP_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} linpan_delly_final_postfilter.vcf > linpan_delly_final_postfilter_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o ncbi_SNP_final_prefilter_norm_${i}.vcf.gz ncbi_SNP_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o ncbi_delly_final_prefilter_norm_${i}.vcf.gz ncbi_delly_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o linpan_SNP_final_prefilter_norm_${i}.vcf.gz linpan_SNP_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o linpan_delly_final_prefilter_norm_${i}.vcf.gz linpan_delly_final_prefilter_${i}.vcf
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o ncbi_SNP_final_postfilter_norm_${i}.vcf.gz ncbi_SNP_final_postfilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o ncbi_delly_final_postfilter_norm_${i}.vcf.gz ncbi_delly_final_postfilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o linpan_SNP_final_postfilter_norm_${i}.vcf.gz linpan_SNP_final_postfilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o linpan_delly_final_postfilter_norm_${i}.vcf.gz linpan_delly_final_postfilter_${i}.vcf
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i ncbi_SNP_final_prefilter_norm_${i}.vcf.gz -o ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i linpan_SNP_final_prefilter_norm_${i}.vcf.gz -o linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i ncbi_SNP_final_postfilter_norm_${i}.vcf.gz -o ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i linpan_SNP_final_postfilter_norm_${i}.vcf.gz -o linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

${DEPOT}/apps/rtg-tools-3.12.1/rtg format -o linpan.sdf ${LINPAN}/${PREFIX}_panref2_sorted.fa # format the reference genome the variants are called against to RTG's Sequence Data File format 

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools index -t --threads ${N} ncbi_delly_final_prefilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} linpan_delly_final_prefilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} ncbi_delly_final_postfilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} linpan_delly_final_postfilter_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b ncbi_delly_final_prefilter_norm_${i}.vcf.gz -c linpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b ncbi_delly_final_postfilter_norm_${i}.vcf.gz -c linpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


## Linear representative genome vs. Minigraph-Cactus pangenome
mkdir -p ./ncbi_mcpan/
cd ./ncbi_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

cp ${MCPAN}/${PREFIX}-pg/called_re/snp/${PREFIX}_mcSNP.renamed.vcf ./mcpan_SNP_final_prefilter.vcf 
cp ${MCPAN}/${PREFIX}-pg/called_re/snp/${PREFIX}_mcSNP.filtered.renamed.vcf ./mcpan_SNP_final_postfilter.vcf 
cp ${MCPAN}/${PREFIX}-pg/called_re/sv/delly/mcpan_delly.filtered.renamed.vcf ./mcpan_delly_final_prefilter.vcf
cp ${MCPAN}/${PREFIX}-pg/called_re/sv/delly/mcpan_delly.filtered2.sorted.vcf ./mcpan_delly_final_postfilter.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
picard UpdateVcfSequenceDictionary I=./mcpan_SNP_final_prefilter.vcf O=./mcpan_SNP_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./mcpan_delly_final_prefilter.vcf O=./mcpan_delly_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_SNP_final_prefilter_updt.vcf O=./mcpan_SNP_final_prefilter_sorted.vcf
picard SortVcf I=./mcpan_delly_final_prefilter_updt.vcf O=./mcpan_delly_final_prefilter_sorted.vcf

picard UpdateVcfSequenceDictionary I=./mcpan_SNP_final_postfilter.vcf O=./mcpan_SNP_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./mcpan_delly_final_postfilter.vcf O=./mcpan_delly_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_SNP_final_postfilter_updt.vcf O=./mcpan_SNP_final_postfilter_sorted.vcf
picard SortVcf I=./mcpan_delly_final_postfilter_updt.vcf O=./mcpan_delly_final_postfilter_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools view -a -I -s ${i} --thread ${N} mcpan_SNP_final_prefilter_sorted.vcf > mcpan_SNP_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} mcpan_delly_final_prefilter_sorted.vcf > mcpan_delly_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} mcpan_SNP_final_postfilter_sorted.vcf > mcpan_SNP_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} mcpan_delly_final_postfilter_sorted.vcf > mcpan_delly_final_postfilter_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_SNP_final_prefilter_norm_${i}.vcf.gz mcpan_SNP_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_delly_final_prefilter_norm_${i}.vcf.gz mcpan_delly_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_SNP_final_postfilter_norm_${i}.vcf.gz mcpan_SNP_final_postfilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_delly_final_postfilter_norm_${i}.vcf.gz mcpan_delly_final_postfilter_${i}.vcf
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i mcpan_SNP_final_prefilter_norm_${i}.vcf.gz -o mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i mcpan_SNP_final_postfilter_norm_${i}.vcf.gz -o mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} mcpan_delly_final_prefilter_norm_${i}.vcf.gz
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} mcpan_delly_final_postfilter_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b ncbi_delly_final_prefilter_norm_${i}.vcf.gz -c mcpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b ncbi_delly_final_postfilter_norm_${i}.vcf.gz -c mcpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


## Linear representative genome vs. VG pangenome
mkdir -p ./ncbi_vgpan/
cd ./ncbi_vgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/
mkdir -p ./sv_comp/vg/

cp ${VGPAN}/called_re/snp/${PREFIX}2_vgSNP.merged.vcf.gz ./vgpan_SNP_final_prefilter.vcf.gz 
cp ${VGPAN}/called_re/snp/${PREFIX}2_vgSNP.filtered.recode.vcf ./vgpan_SNP_final_postfilter.vcf 
cp ${VGPAN}/called_re/sv/delly/vgpan_delly.filtered.vcf.gz ./vgpan_delly_final_prefilter.vcf.gz
cp ${VGPAN}/called_re/sv/delly/vgpan_delly.filtered2.recode.vcf ./vgpan_delly_final_postfilter.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
picard UpdateVcfSequenceDictionary I=./vgpan_SNP_final_prefilter.vcf.gz O=./vgpan_SNP_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./vgpan_delly_final_prefilter.vcf.gz O=./vgpan_delly_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_SNP_final_prefilter_updt.vcf O=./vgpan_SNP_final_prefilter_sorted.vcf
picard SortVcf I=./vgpan_delly_final_prefilter_updt.vcf O=./vgpan_delly_final_prefilter_sorted.vcf

picard UpdateVcfSequenceDictionary I=./vgpan_SNP_final_postfilter.vcf O=./vgpan_SNP_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./vgpan_delly_final_postfilter.vcf O=./vgpan_delly_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_SNP_final_postfilter_updt.vcf O=./vgpan_SNP_final_postfilter_sorted.vcf
picard SortVcf I=./vgpan_delly_final_postfilter_updt.vcf O=./vgpan_delly_final_postfilter_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools view -a -I -s ${i} --thread ${N} vgpan_SNP_final_prefilter_sorted.vcf > vgpan_SNP_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} vgpan_delly_final_prefilter_sorted.vcf > vgpan_delly_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} vgpan_SNP_final_postfilter_sorted.vcf > vgpan_SNP_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} vgpan_delly_final_postfilter_sorted.vcf > vgpan_delly_final_postfilter_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_SNP_final_prefilter_norm_${i}.vcf.gz vgpan_SNP_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_delly_final_prefilter_norm_${i}.vcf.gz vgpan_delly_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_SNP_final_postfilter_norm_${i}.vcf.gz vgpan_SNP_final_postfilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_delly_final_postfilter_norm_${i}.vcf.gz vgpan_delly_final_postfilter_${i}.vcf 
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i vgpan_SNP_final_prefilter_norm_${i}.vcf.gz -o vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i vgpan_SNP_final_postfilter_norm_${i}.vcf.gz -o vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} vgpan_delly_final_prefilter_norm_${i}.vcf.gz
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} vgpan_delly_final_postfilter_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b ncbi_delly_final_prefilter_norm_${i}.vcf.gz -c vgpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b ncbi_delly_final_postfilter_norm_${i}.vcf.gz -c vgpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


## Original pangenome vs. Minigraph-Cactus pangenome
mkdir -p ./orgpan_mcpan/
cd ./orgpan_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

cp ${MCPAN}/${PREFIX}-pg/called_re/sv/vg/${PREFIX}_mcSV.merged.vcf.gz ./mcpan_vg_final_prefilter.vcf.gz
cp ${MCPAN}/${PREFIX}-pg/called_re/sv/vg/${PREFIX}_mcSV.decomp.annot.ext.vcf ./mcpan_vg_final_postfilter.vcf
cp ${ORGPAN}/called_re/snp/barnswallow_orgSNP.renamed.vcf ./orgpan_SNP_final_prefilter.vcf 
cp ${ORGPAN}/called_re/snp/barnswallow_orgSNP.filtered.renamed.vcf ./orgpan_SNP_final_postfilter.vcf 
cp ${ORGPAN}/called_re/sv/delly/orgpan_delly.filtered.renamed.vcf ./orgpan_delly_final_prefilter.vcf
cp ${ORGPAN}/called_re/sv/delly/orgpan_delly.filtered2.renamed.vcf ./orgpan_delly_final_postfilter.vcf
cp ${ORGPAN}/called_re/sv/vg/barnswallow_orgSV.merged.vcf.gz  ./orgpan_vg_final_prefilter.vcf.gz
cp ${ORGPAN}/called_re/sv/vg/barnswallow_orgSV.decomp.annot.ext.vcf  ./orgpan_vg_final_postfilter.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 

python3 # run the code below in python
def extract_contigs_from_vcf(vcf_file):
    contigs_dict = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##contig='):
                parts = line.strip().split('<')[1].split('>')[0]
                contig_id = None
                contig_length = None
                for part in parts.split(','):
                    if part.startswith('ID='):
                        contig_id = part.split('=')[1]
                    elif part.startswith('length='):
                        contig_length = int(part.split('=')[1])
                if contig_id and contig_length:
                    contigs_dict[contig_id] = contig_length           
    return contigs_dict

vcf_file="orgpan_vg_final_postfilter.vcf"
contigs = extract_contigs_from_vcf(vcf_file)
print("Contigs:", contigs)

with open("orgSV_chroms.txt", 'w') as f:
    for key, value in contigs.items():
        f.write(str(key) + " " + str(value) + '\n') 

exit()

sort -k2,2 -n orgSV_chroms.txt > sorted_orgSV_chroms.txt
cut -f1,2 ${REF}/${GENOME}.fna.fai > ref_chroms.txt
sort -k2,2 -n ref_chroms.txt > sorted_ref_chroms.txt
join -1 2 -2 2 sorted_orgSV_chroms.txt sorted_ref_chroms.txt > combined_chroms.txt # check if there are duplicates and if so manually curate them (in this case, duplicates at lengths of "18466", "25102", and "65965"; refer to the fai file)
cut -d " " -f2,3 combined_chroms.txt > chr_name_conv.txt

#sed 's/|/\//g' orgpan_vg_final.vcf > orgpan_vg_final_sbst.vcf # vg vcf files include "|" after truvari
#sed 's/|/\//g' mcpan_vg_final.vcf > mcpan_vg_final_sbst.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt mcpan_vg_final_prefilter.vcf.gz -Ov -o mcpan_vg_final_prefilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_vg_final_prefilter.vcf.gz -Ov -o orgpan_vg_final_prefilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_delly_final_prefilter.vcf -Ov -o orgpan_delly_final_prefilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_SNP_final_prefilter.vcf -Ov -o orgpan_SNP_final_prefilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt mcpan_vg_final_postfilter.vcf -Ov -o mcpan_vg_final_postfilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_vg_final_postfilter.vcf -Ov -o orgpan_vg_final_postfilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_delly_final_postfilter.vcf -Ov -o orgpan_delly_final_postfilter_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_SNP_final_postfilter.vcf -Ov -o orgpan_SNP_final_postfilter_annot.vcf

picard UpdateVcfSequenceDictionary I=./mcpan_vg_final_prefilter_annot.vcf O=./mcpan_vg_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_SNP_final_prefilter_annot.vcf O=./orgpan_SNP_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_delly_final_prefilter_annot.vcf O=./orgpan_delly_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_vg_final_prefilter_annot.vcf O=./orgpan_vg_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_vg_final_prefilter_updt.vcf O=./mcpan_vg_final_prefilter_sorted.vcf
picard SortVcf I=./orgpan_SNP_final_prefilter_updt.vcf O=./orgpan_SNP_final_prefilter_sorted.vcf
picard SortVcf I=./orgpan_delly_final_prefilter_updt.vcf O=./orgpan_delly_final_prefilter_sorted.vcf
picard -Xmx500g SortVcf I=./orgpan_vg_final_prefilter_updt.vcf O=./orgpan_vg_final_prefilter_sorted.vcf # used highmem

picard UpdateVcfSequenceDictionary I=./mcpan_vg_final_postfilter_annot.vcf O=./mcpan_vg_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_SNP_final_postfilter_annot.vcf O=./orgpan_SNP_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_delly_final_postfilter_annot.vcf O=./orgpan_delly_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_vg_final_postfilter_annot.vcf O=./orgpan_vg_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_vg_final_postfilter_updt.vcf O=./mcpan_vg_final_postfilter_sorted.vcf
picard SortVcf I=./orgpan_SNP_final_postfilter_updt.vcf O=./orgpan_SNP_final_postfilter_sorted.vcf
picard SortVcf I=./orgpan_delly_final_postfilter_updt.vcf O=./orgpan_delly_final_postfilter_sorted.vcf
picard SortVcf I=./orgpan_vg_final_postfilter_updt.vcf O=./orgpan_vg_final_postfilter_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' mcpan_vg_final_prefilter_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > mcpan_vg_final_prefilter_${i}.vcf # sed step to bypass bcftools error
  bcftools view -a -I -s ${i} --thread ${N} orgpan_SNP_final_prefilter_sorted.vcf > orgpan_SNP_final_prefilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} orgpan_delly_final_prefilter_sorted.vcf > orgpan_delly_final_prefilter_${i}.vcf
  sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' orgpan_vg_final_prefilter_sorted.vcf |bcftools view -a -I -s ${i} --thread ${N} > orgpan_vg_final_prefilter_${i}.vcf
  
  sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' mcpan_vg_final_postfilter_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > mcpan_vg_final_postfilter_${i}.vcf # sed step to bypass bcftools error
  bcftools view -a -I -s ${i} --thread ${N} orgpan_SNP_final_postfilter_sorted.vcf > orgpan_SNP_final_postfilter_${i}.vcf
  bcftools view -a -I -s ${i} --thread ${N} orgpan_delly_final_postfilter_sorted.vcf > orgpan_delly_final_postfilter_${i}.vcf
  sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' orgpan_vg_final_postfilter_sorted.vcf |bcftools view -a -I -s ${i} --thread ${N} > orgpan_vg_final_postfilter_${i}.vcf  
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_SNP_final_prefilter_norm_${i}.vcf.gz orgpan_SNP_final_prefilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_delly_final_prefilter_norm_${i}.vcf.gz orgpan_delly_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_vg_final_prefilter_norm_${i}.vcf.gz mcpan_vg_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_vg_final_prefilter_norm_${i}.vcf.gz orgpan_vg_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_SNP_final_postfilter_norm_${i}.vcf.gz orgpan_SNP_final_postfilter_${i}.vcf
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_delly_final_postfilter_norm_${i}.vcf.gz orgpan_delly_final_postfilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o mcpan_vg_final_postfilter_norm_${i}.vcf.gz mcpan_vg_final_postfilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o orgpan_vg_final_postfilter_norm_${i}.vcf.gz orgpan_vg_final_postfilter_${i}.vcf 
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_mcpan/mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i orgpan_SNP_final_prefilter_norm_${i}.vcf.gz -o orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
  cp ../ncbi_mcpan/mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i orgpan_SNP_final_postfilter_norm_${i}.vcf.gz -o orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_mcpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_mcpan/mcpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} orgpan_delly_final_prefilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} mcpan_vg_final_prefilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} orgpan_vg_final_prefilter_norm_${i}.vcf.gz
  cp ../ncbi_mcpan/mcpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} orgpan_delly_final_postfilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} mcpan_vg_final_postfilter_norm_${i}.vcf.gz
  bcftools index -t --threads ${N} orgpan_vg_final_postfilter_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b orgpan_delly_final_prefilter_norm_${i}.vcf.gz -c mcpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_vg_final_prefilter_norm_${i}.vcf.gz -c mcpan_vg_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_delly_final_postfilter_norm_${i}.vcf.gz -c mcpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_vg_final_postfilter_norm_${i}.vcf.gz -c mcpan_vg_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000
done

cd ../


## Original pangenome vs. VG pangenome
mkdir -p ./orgpan_vgpan/
cd ./orgpan_vgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

cp ${VGPAN}/called_re/sv/vg/${PREFIX}2_vgSV.merged.vcf.gz ./vgpan_vg_final_prefilter.vcf.gz
cp ${VGPAN}/called_re/sv/vg/${PREFIX}2_vgSV.decomp.annot.ext.vcf ./vgpan_vg_final_postfilter.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
zcat vgpan_vg_final_prefilter.vcf.gz | awk '$5 != "." || /^#/' > vgpan_vg_cleaned.vcf # removed 76 variants with ALT = "." that hinders picard steps
picard UpdateVcfSequenceDictionary I=./vgpan_vg_cleaned.vcf O=./vgpan_vg_final_prefilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard -Xmx100g SortVcf I=./vgpan_vg_final_prefilter_updt.vcf O=./vgpan_vg_final_prefilter_sorted.vcf
picard UpdateVcfSequenceDictionary I=./vgpan_vg_final_postfilter.vcf O=./vgpan_vg_final_postfilter_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_vg_final_postfilter_updt.vcf O=./vgpan_vg_final_postfilter_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' vgpan_vg_final_prefilter_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > vgpan_vg_final_prefilter_${i}.vcf # sed step to bypass bcftools error
  sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' vgpan_vg_final_postfilter_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > vgpan_vg_final_postfilter_${i}.vcf # sed step to bypass bcftools error
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_vg_final_prefilter_norm_${i}.vcf.gz vgpan_vg_final_prefilter_${i}.vcf 
  bcftools norm -m -any --thread ${N} -Oz -o vgpan_vg_final_postfilter_norm_${i}.vcf.gz vgpan_vg_final_postfilter_${i}.vcf
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_mcpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../orgpan_mcpan/orgpan_vg_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_vg_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} vgpan_vg_final_prefilter_norm_${i}.vcf.gz
  cp ../orgpan_mcpan/orgpan_vg_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_vg_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  bcftools index -t --threads ${N} vgpan_vg_final_postfilter_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b orgpan_delly_final_prefilter_norm_${i}.vcf.gz -c vgpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_vg_final_prefilter_norm_${i}.vcf.gz -c vgpan_vg_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_delly_final_postfilter_norm_${i}.vcf.gz -c vgpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_vg_final_postfilter_norm_${i}.vcf.gz -c vgpan_vg_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


## Linear representative genome vs. Original pangenome
mkdir -p ./ncbi_orgpan/
cd ./ncbi_orgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/ncbi_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b ncbi_delly_final_prefilter_norm_${i}.vcf.gz -c orgpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b ncbi_delly_final_postfilter_norm_${i}.vcf.gz -c orgpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000
done

cd ../


## Original pangenome vs. Linear pangenome
mkdir -p ./orgpan_linpan/
cd ./orgpan_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/linpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/linpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/orgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b orgpan_delly_final_prefilter_norm_${i}.vcf.gz -c linpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b orgpan_delly_final_postfilter_norm_${i}.vcf.gz -c linpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000
done

cd ../


## VG pangenome vs. Linear pangenome
mkdir -p ./vgpan_linpan/
cd ./vgpan_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_linpan/linpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_linpan/linpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_linpan/linpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b vgpan_delly_final_prefilter_norm_${i}.vcf.gz -c linpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b vgpan_delly_final_postfilter_norm_${i}.vcf.gz -c linpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


## VG pangenome vs. Minigraph-Cactus pangenome
mkdir -p ./vgpan_mcpan/
cd ./vgpan_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../ncbi_mcpan/mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_mcpan/mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_prefilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_prefilter -t linpan.sdf
  ${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_postfilter_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i}_postfilter -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  cp ../orgpan_mcpan/mcpan_vg_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/mcpan_vg_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_mcpan/mcpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_vgpan/vgpan_vg_final_prefilter_norm_${i}.vcf.gz ./
  cp ../orgpan_vgpan/vgpan_vg_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_prefilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_mcpan/mcpan_vg_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_mcpan/mcpan_vg_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_mcpan/mcpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_mcpan/mcpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../orgpan_vgpan/vgpan_vg_final_postfilter_norm_${i}.vcf.gz ./
  cp ../orgpan_vgpan/vgpan_vg_final_postfilter_norm_${i}.vcf.gz.tbi ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz ./
  cp ../ncbi_vgpan/vgpan_delly_final_postfilter_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/original/SRR_Acc_List.txt`; do
  truvari bench -b vgpan_delly_final_prefilter_norm_${i}.vcf.gz -c mcpan_delly_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b vgpan_vg_final_prefilter_norm_${i}.vcf.gz -c mcpan_vg_final_prefilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_prefilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b vgpan_delly_final_postfilter_norm_${i}.vcf.gz -c mcpan_delly_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/delly/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
  truvari bench -b vgpan_vg_final_postfilter_norm_${i}.vcf.gz -c mcpan_vg_final_postfilter_norm_${i}.vcf.gz -o ./sv_comp/vg/${i}_postfilter --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../