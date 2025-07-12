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

# Constructing a VG pangenome
module purge
module load biocontainers
module load bwa
module load picard
module load boost
module load bcftools
module load htslib

## Calling SVs using delly2
### Mapping for delly2
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.fa --OUTPUT ${PREFIX}_panref2_sorted.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.fa

cd ${VGPAN}/augmented

### Running delly2
cd ${VGPAN}/augmented/called/  
mkdir -p ./delly/
cd ./delly/

export OMP_NUM_THREADS=${num_samples} # equal or less than the number of samples
samples=""
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
    samples+="${VGPAN}/augmented/mapped/${i}_sorted.marked.bam "
done

${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > sample_delly.vcf 
bcftools view -Ob -o sample_delly.bcf sample_delly.vcf
bcftools index sample_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o sample_delly.filtered.bcf sample_delly.bcf
bcftools view -Ov -o sample_delly.filtered.vcf sample_delly.filtered.bcf 

## Filtering variants (SNPs)
cd ../

### Compressing and indexing each vcf file first
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_SNP.sorted.vcf.gz
  bcftools index ${i}_SNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_SNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o sample_SNP.merged.vcf.gz --threads ${N} *_SNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf sample_SNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf sample_SNP.merged.vcf.gz --missing-site
vcftools --gzvcf sample_SNP.merged.vcf.gz --depth
vcftools --gzvcf sample_SNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf sample_SNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sample_SNP.merged.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r
Rscript variant_check.R

### Filtering vcf file based on the determined conservative cut-off values from above results                                                                     
vcftools --gzvcf sample_SNP.merged.vcf.gz --out sample_SNP.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### Filtering variants (SVs)
cd ./delly/
bcftools convert --thread ${N} -Oz -o sample_delly.filtered.vcf.gz sample_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf sample_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf sample_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf sample_delly.filtered.vcf.gz --depth
vcftools --gzvcf sample_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf sample_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sample_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r
Rscript variant_check.R


### Filtering vcf file based on the determined conservative cut-off values from above results 
vcftools --gzvcf sample_delly.filtered.vcf.gz --out ../sample_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Constructing a VG pangenome
cd ../

### Preprocessing input files *** modify
bcftools convert --thread ${N} -Ov -o sample_delly.filtered2.vcf sample_delly.filtered2.recode.vcf.gz
#while read -r line; do
#    pair=($line)
#    echo "s/${pair[0]}/${pair[1]}/" >> delly_replace.sed
#done < delly_sample.list # make a sed script from the list of "search-replace" pairs of "wrong-correct" sample name pairs (delly_sample.list)
#sed -f delly_replace.sed sample_delly.filtered2.vcf > sample_delly.filtered2.replaced.vcf # fix sample names

#while read -r line; do
#    pair=($line)
#    echo "s/${pair[0]}/${pair[1]}/" >> snp_replace.sed
#done < snp_sample.list # make a sed script from the list of "search-replace" pairs of "wrong-correct" sample name pairs (snp_sample.list)
#sed -f snp_replace.sed sample_SNP.filtered.recode.vcf > sample_SNP.filtered.replaced.vcf # fix sample names

bcftools query -l sample_delly.filtered2.vcf > samples.txt # make a sample list
bcftools view -S samples.txt sample_delly.filtered2.recode.vcf > sample_delly.filtered2.ordered.vcf # order samples; check the sample header! "bcftools view sample_delly.filtered2.ordered.vcf | grep -v "##" | head"
cat sample_delly.filtered2.ordered.vcf | grep -v "SVTYPE=BND" | grep -v "SVTYPE=DUP" > sample_delly.filtered3.ordered.vcf # remove unsupported SV types in vg construct (i.e., DUP, BND)
bcftools view -S samples.txt sample_SNP.filtered.recode.vcf > sample_SNP.filtered.ordered.vcf # order samples; check the sample header! "bcftools view sample_SNP.filtered.ordered.vcf | grep -v "##" | head"

bcftools convert --thread ${N} -Oz -o sample_SNP.filtered.vcf.gz sample_SNP.filtered.ordered.vcf
bcftools convert --thread ${N} -Oz -o sample_delly.filtered3.vcf.gz sample_delly.filtered3.ordered.vcf
tabix -p vcf sample_SNP.filtered.vcf.gz
tabix -p vcf sample_delly.filtered3.vcf.gz
bcftools concat --allow-overlaps -Oz -o sample_SNP_SV.vcf.gz --threads ${N} sample_SNP.filtered.vcf.gz sample_delly.filtered3.vcf.gz # combine SNP and SV vcf files

bcftools norm -a --atom-overlaps . -m -any --threads ${N} -Oz -o sample_SNP_SV.normed.vcf.gz sample_SNP_SV.vcf.gz # normalize the vcf file      
tabix -p vcf sample_SNP_SV.normed.vcf.gz
