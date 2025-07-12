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

# Augmenting the graph pangenome (preprocessing steps before alignmening were modified from Sacomandi et al. (2023)'s scripts)
cd ${VGPAN}/

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./augmented/
mkdir -p ./tmp/

## Aligning the individuals that were used to build the pangenome; do augmentation steps when reads that were used to build the pangenome are different from reads that will be used to call variants like this case.  
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  cat ${SRA}/Nremoved/${i}*_1.fq* >  ${SRA}/Nremoved/${i}_R1.fq.gz 
  cat ${SRA}/Nremoved/${i}*_2.fq* >  ${SRA}/Nremoved/${i}_R2.fq.gz
  ${APP}/vg map -t ${N} -f ${SRA}/Nremoved/${i}_R1.fq.gz -f ${SRA}/Nremoved/${i}_R2.fq.gz -x ${PREFIX}-vg_mod_chopped.xg -g ${PREFIX}-vg_mod_chopped_pruned.gcsa > ${i}_aln.gam

## Filtering secondary and ambiguous read mappings out of the GAM of the above step
  ${APP}/vg filter -t ${N} ${i}_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${PREFIX}-vg_mod_chopped.xg >${i}_aln.filtered.gam #-r : minimum score to keep primary alignment; -f: normalize score based on length; -u: use substitution count instead of score; -m: filter reads that don't begin with at least N matches on each end; -q: filter alignments with mapping quality < N; -D: clip back the ends of reads that are ambiguously aligned, up to N bases
done

cat *filtered.gam > combined_filtered.gam 

## Augmenting the graph with all variation from the GAM of the above step
${APP}/vg convert -t ${N} ${PREFIX}-vg_mod_chopped_new_L.xg -p > ${PREFIX}-vg.pg
${APP}/vg augment -t ${N} ${PREFIX}-vg.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}-vg_aug.gam > ${PREFIX}-vg_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 

## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}-vg_aug.pg > ${PREFIX}-vg_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}-vg_aug_chopped.xg ${PREFIX}-vg_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}-vg_aug_chopped.pg > ${PREFIX}-vg_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}-vg_aug_chopped_pruned.gcsa ${PREFIX}-vg_aug_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}-vg_aug_chopped.pg -x ${PREFIX}-vg_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}-vg_aug_chopped_new_L.xg -p > ${PREFIX}-vg_aug_new.pg
