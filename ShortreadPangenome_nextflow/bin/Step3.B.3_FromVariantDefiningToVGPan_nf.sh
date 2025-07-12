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

cd ${VGPAN}/augmented/called/  

### Constructing a graph pangenome using vg and the linear pangenome
${APP}/vg construct -a -f -S -r ${LINPAN}/${PREFIX}_panref2_sorted.fa -v sample_SNP_SV.normed.vcf.gz -m 1000000 > ../${PREFIX}-vg-pg.vg #| ${APP}/vg convert - | ${APP}/vg mod -x 32 -  # too big inversions exist for Protobuf to stream -> moving node chopping outside of vg construct

${APP}/vg stats -lz ${PREFIX}-vg-pg.vg # check pangenome basic statistics


# Augmenting the graph pangenome (preprocessing steps before alignmening were modified from Sacomandi et al. (2023)'s scripts)
cd ${VGPAN}/

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./augmented/
mkdir -p ./tmp/

${APP}/vg mod -t ${N} -O ${PREFIX}-vg-pg.vg > ${PREFIX}-vg_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}-vg_mod.xg ${PREFIX}-vg_mod.vg

## Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}-vg_mod.vg > ${PREFIX}-vg_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}-vg_mod_chopped.xg ${PREFIX}-vg_mod_chopped.vg

## Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}-vg_mod_chopped.vg > ${PREFIX}-vg_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}-vg_mod_chopped_pruned.gcsa ${PREFIX}-vg_mod_chopped_pruned.vg

## Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}-vg_mod_chopped.vg -x ${PREFIX}-vg_mod_chopped_new_L.xg #-L: preserve alt paths in the xg
