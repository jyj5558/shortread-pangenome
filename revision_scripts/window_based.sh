#!/bin/bash
#SBATCH --job-name=variant_metrics
#SBATCH -A fnrdewoody
#SBATCH -t 14-00:00:00
#SBATCH -N 1
#SBATCH -n 64
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

module load biocontainers
module load bcftools
module load vcftools
module load r

cd /scratch/negishi/jeon96/swallow/
mkdir -p popgen_re/window_based

cd popgen_re/window_based
cp ${BASE}/benchmarking/ncbi_linpan/ncbi_SNP_final.vcf ./
cp ${BASE}/benchmarking/ncbi_linpan/linpan_SNP_final.vcf ./
cp ${BASE}/benchmarking/ncbi_vgpan/vgpan_SNP_final_sorted.vcf ./
cp ${BASE}/benchmarking/ncbi_mcpan/mcpan_SNP_final_sorted.vcf ./
cp ${BASE}/benchmarking/orgpan_mcpan/orgpan_SNP_final_sorted.vcf ./

bcftools convert ncbi_SNP_final.vcf -Oz -o ncbi_SNP_final.vcf.gz
bcftools index ncbi_SNP_final.vcf.gz
bcftools view ncbi_SNP_final.vcf.gz --regions NC_053450.1 -Oz -o ncbi_SNP_NC_053450.vcf.gz
vcftools --gzvcf ncbi_SNP_NC_053450.vcf.gz --TajimaD 50000 --out ncbi_SNP_NC_053450
vcftools --gzvcf ncbi_SNP_NC_053450.vcf.gz --window-pi 50000 --window-pi-step 10000 --out ncbi_SNP_NC_053450

bcftools convert linpan_SNP_final.vcf -Oz -o linpan_SNP_final.vcf.gz
bcftools index linpan_SNP_final.vcf.gz
bcftools view linpan_SNP_final.vcf.gz --regions NC_053450.1 -Oz -o linpan_SNP_NC_053450.vcf.gz
vcftools --gzvcf linpan_SNP_NC_053450.vcf.gz --TajimaD 50000 --out linpan_SNP_NC_053450
vcftools --gzvcf linpan_SNP_NC_053450.vcf.gz --window-pi 50000 --window-pi-step 10000 --out linpan_SNP_NC_053450

bcftools convert vgpan_SNP_final_sorted.vcf -Oz -o vgpan_SNP_final_sorted.vcf.gz
bcftools index vgpan_SNP_final_sorted.vcf.gz
bcftools view vgpan_SNP_final_sorted.vcf.gz --regions NC_053450.1 -Oz -o vgpan_SNP_NC_053450.vcf.gz
vcftools --gzvcf vgpan_SNP_NC_053450.vcf.gz --TajimaD 50000 --out vgpan_SNP_NC_053450
vcftools --gzvcf vgpan_SNP_NC_053450.vcf.gz --window-pi 50000 --window-pi-step 10000 --out vgpan_SNP_NC_053450

bcftools convert mcpan_SNP_final_sorted.vcf -Oz -o mcpan_SNP_final_sorted.vcf.gz
bcftools index mcpan_SNP_final_sorted.vcf.gz
bcftools view mcpan_SNP_final_sorted.vcf.gz --regions NC_053450.1 -Oz -o mcpan_SNP_NC_053450.vcf.gz
vcftools --gzvcf mcpan_SNP_NC_053450.vcf.gz --TajimaD 50000 --out mcpan_SNP_NC_053450
vcftools --gzvcf mcpan_SNP_NC_053450.vcf.gz --window-pi 50000 --window-pi-step 10000 --out mcpan_SNP_NC_053450

bcftools convert orgpan_SNP_final_sorted.vcf -Oz -o orgpan_SNP_final_sorted.vcf.gz
bcftools index orgpan_SNP_final_sorted.vcf.gz
bcftools view orgpan_SNP_final_sorted.vcf.gz --regions NC_053450.1 -Oz -o orgpan_SNP_NC_053450.vcf.gz
vcftools --gzvcf orgpan_SNP_NC_053450.vcf.gz --TajimaD 50000 --out orgpan_SNP_NC_053450
vcftools --gzvcf orgpan_SNP_NC_053450.vcf.gz --window-pi 50000 --window-pi-step 10000 --out orgpan_SNP_NC_053450

R
options(scipen = 999)
library(ggplot2)
library(scales)

windowed_tajima_ncbi <- read.table("ncbi_SNP_NC_053450.Tajima.D", sep="\t", header=TRUE)
str(windowed_tajima_ncbi)
summary(windowed_tajima_ncbi) # Min.   :-1.21993; Max.   : 3.63543
quantile(windowed_tajima_ncbi$TajimaD, probs = c(.95, .99, .999), na.rm=TRUE)
#      95%       99%     99.9%
#0.4054974 0.8709698 2.8285792

svg("ncbi_NC053450_tajima.svg", width = 8, height = 6)
ggplot(windowed_tajima_ncbi, aes(x = BIN_START, y = TajimaD)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Tajima's D") +
  ylim(-1.3, 4.2) +
  theme_classic(base_size = 14)
dev.off()

windowed_tajima_linpan <- read.table("linpan_SNP_NC_053450.Tajima.D", sep="\t", header=TRUE)
str(windowed_tajima_linpan)
summary(windowed_tajima_linpan) # Min.   :-1.21993; Max.   : 3.63543
quantile(windowed_tajima_linpan$TajimaD, probs = c(.95, .99, .999), na.rm=TRUE)
#      95%       99%     99.9%
#0.3986813 0.8500352 2.5158037

svg("linpan_NC053450_tajima.svg", width = 8, height = 6)
ggplot(windowed_tajima_linpan, aes(x = BIN_START, y = TajimaD)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Tajima's D") +
  ylim(-1.3, 4.2) +
  theme_classic(base_size = 14)
dev.off()

windowed_tajima_vgpan <- read.table("vgpan_SNP_NC_053450.Tajima.D", sep="\t", header=TRUE)
str(windowed_tajima_vgpan)
summary(windowed_tajima_vgpan) #  Min.   :-1.11851; Max.   : 3.50728
quantile(windowed_tajima_vgpan$TajimaD, probs = c(.95, .99, .999), na.rm=TRUE)
#      95%       99%     99.9%
#0.448896 1.106265 2.649118

svg("vgpan_NC053450_tajima.svg", width = 8, height = 6)
ggplot(windowed_tajima_vgpan, aes(x = BIN_START, y = TajimaD)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Tajima's D") +
  ylim(-1.3, 4.2) +
  theme_classic(base_size = 14)
dev.off()

windowed_tajima_mcpan <- read.table("mcpan_SNP_NC_053450.Tajima.D", sep="\t", header=TRUE)
str(windowed_tajima_mcpan)
summary(windowed_tajima_mcpan) # Min.   :-1.14547; Max.   : 4.19489
quantile(windowed_tajima_mcpan$TajimaD, probs = c(.95, .99, .999), na.rm=TRUE)
#      95%       99%     99.9%
#0.483467 1.033918 3.329551

svg("mcpan_NC053450_tajima.svg", width = 8, height = 6)
ggplot(windowed_tajima_mcpan, aes(x = BIN_START, y = TajimaD)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Tajima's D") +
  ylim(-1.3, 4.2) +
  theme_classic(base_size = 14)
dev.off()

windowed_tajima_orgpan <- read.table("orgpan_SNP_NC_053450.Tajima.D", sep="\t", header=TRUE)
str(windowed_tajima_orgpan)
summary(windowed_tajima_orgpan) # Min.   :-1.1185; Max.   : 3.3629
quantile(windowed_tajima_orgpan$TajimaD, probs = c(.95, .99, .999), na.rm=TRUE)
#      95%       99%     99.9%
#0.4224632 1.0085906 1.8351232

svg("orgpan_NC053450_tajima.svg", width = 8, height = 6)
ggplot(windowed_tajima_orgpan, aes(x = BIN_START, y = TajimaD)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Tajima's D") +
  ylim(-1.3, 4.2) +
  theme_classic(base_size = 14)
dev.off()


windowed_pi_ncbi <- read.table("ncbi_SNP_NC_053450.windowed.pi", sep="\t", header=TRUE)
str(windowed_pi_ncbi)
summary(windowed_pi_ncbi) # Min.   :0.000001241; Max.   :0.000214747
quantile(windowed_pi_ncbi$PI, probs = c(.95, .99, .999), na.rm=TRUE)
#         95%          99%        99.9%
#6.078755e-05 8.742359e-05 1.509328e-04

svg("ncbi_NC053450_pi.svg", width = 8, height = 6)
ggplot(windowed_pi_ncbi, aes(x = BIN_START, y = PI)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0.0001, linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Nucleotide diversity") +
  ylim(0, 0.00030) +
  theme_classic(base_size = 14)
dev.off()

windowed_pi_linpan <- read.table("linpan_SNP_NC_053450.windowed.pi", sep="\t", header=TRUE)
str(windowed_pi_linpan)
summary(windowed_pi_linpan) # Min.   :0.000001241; Max.   :0.000214747 
quantile(windowed_pi_linpan$PI, probs = c(.95, .99, .999), na.rm=TRUE)
#         95%          99%        99.9%
#6.070738e-05 8.665703e-05 1.509330e-04

svg("linpan_NC053450_pi.svg", width = 8, height = 6)
ggplot(windowed_pi_linpan, aes(x = BIN_START, y = PI)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0.0001, linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Nucleotide diversity") +
  ylim(0, 0.00030) +
  theme_classic(base_size = 14)
dev.off()

windowed_pi_vgpan <- read.table("vgpan_SNP_NC_053450.windowed.pi", sep="\t", header=TRUE)
str(windowed_pi_vgpan)
summary(windowed_pi_vgpan) # Min.   :0.000001241; Max.   :0.000207564 
quantile(windowed_pi_vgpan$PI, probs = c(.95, .99, .999), na.rm=TRUE)
#         95%          99%        99.9%
#4.293738e-05 6.713535e-05 1.418350e-04

svg("vgpan_NC053450_pi.svg", width = 8, height = 6)
ggplot(windowed_pi_vgpan, aes(x = BIN_START, y = PI)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0.0001, linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Nucleotide diversity") +
  ylim(0, 0.00030) +
  theme_classic(base_size = 14)
dev.off()

windowed_pi_mcpan <- read.table("mcpan_SNP_NC_053450.windowed.pi", sep="\t", header=TRUE)
str(windowed_pi_mcpan)
summary(windowed_pi_mcpan) #Min.   :0.000001241; Max.   :0.000276511 
quantile(windowed_pi_mcpan$PI, probs = c(.95, .99, .999), na.rm=TRUE)
#         95%          99%        99.9%
#4.068762e-05 6.185374e-05 1.649521e-04

svg("mcpan_NC053450_pi.svg", width = 8, height = 6)
ggplot(windowed_pi_mcpan, aes(x = BIN_START, y = PI)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0.0001, linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  scale_y_continuous(labels = function(x) sprintf("%.5f", x)) +
  labs(x = "Chromosome position (bp)", y = "Nucleotide diversity") +
  ylim(0, 0.00030) +
  theme_classic(base_size = 14)
dev.off()

windowed_pi_orgpan <- read.table("orgpan_SNP_NC_053450.windowed.pi", sep="\t", header=TRUE)
str(windowed_pi_orgpan)
summary(windowed_pi_orgpan) # Min.   :0.000001241; Max.   :0.000161184
quantile(windowed_pi_orgpan$PI, probs = c(.95, .99, .999), na.rm=TRUE)
#         95%          99%        99.9%
#4.215630e-05 6.042943e-05 1.126380e-04

svg("orgpan_NC053450_pi.svg", width = 8, height = 6)
ggplot(windowed_pi_orgpan, aes(x = BIN_START, y = PI)) + 
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0.0001, linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -2, linetype="dashed", color = "red") +
  labs(x = "Chromosome position (bp)", y = "Nucleotide diversity") +
  ylim(0, 0.00030) +
  theme_classic(base_size = 14)
dev.off()
