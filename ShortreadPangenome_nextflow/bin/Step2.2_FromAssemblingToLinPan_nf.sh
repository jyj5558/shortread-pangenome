#!/bin/bash

SRA=$1
REF=$2
LINPAN=$3
CONTIG=$4
GENOME=$5
N=$6
APP=$7
MASURCA=$8
GENOME_SIZE=$9
PREFIX=${10}
MASURCA_CONFIG_MAPPED=${11}
MASURCA_CONFIG_UNMAPPED=${12}

# Please note that this script was based on Hu et al. (2020) in Legume Genomics: Methods and Protocols (and its scripts) AND Yao et al. (2015) in Genome Biology

cd ${SRA}/
cd ./mapped

# Assembling mapped and unmapped reads independently, redundantly (MaSuRCA and MegaHit) modifying the approach of Yao et al. (2015)
## MaSuRCA assembly using mapped reads
module purge
module load biocontainers
module load megahit

mkdir -p ./masurca_mapped
mkdir -p ./masurca_unmapped

cd ./masurca_mapped
cp ~/masurca_config_example ./${MASURCA_CONFIG_MAPPED}

FRAG_MEAN=$(awk 'NR==8 {print $6}' ../ihist.picard.txt)
FRAG_SD=$(awk 'NR==8 {print $7}' ../ihist.picard.txt)

# confirm values are non-empty
if [[ -z "$FRAG_MEAN" || -z "$FRAG_SD" ]]; then
  echo "ERROR: Could not extract FRAG_MEAN and FRAG_SD from ihist.picard.txt"
  exit 1
fi


FWD_READS=${SRA}/mapped/mapped_R1.fastq
REV_READS=${SRA}/mapped/mapped_R2.fastq
SGL_READS=${SRA}/mapped/mapped_unpaired.fastq

sed -i "s|PE= pe |PE= pe ${FRAG_MEAN} ${FRAG_SD} ${FWD_READS} ${REV_READS}|" ${MASURCA_CONFIG_MAPPED}
sed -i "s|PE= se |PE= se ${FRAG_MEAN} ${FRAG_SD} ${SGL_READS}|" ${MASURCA_CONFIG_MAPPED}
##

echo "Masurca assembly of mapped reads started"
${MASURCA}/masurca ${MASURCA_CONFIG_MAPPED} # configure the configuration file beforehand following updated MaSuRCA manual
bash assemble.sh
echo "Masurca assembly of mapped reads finished"
cd ../
#

# MaSuRCA assembly using unmapped reads
cd ./masurca_unmapped
cp ~/masurca_config_example ./${MASURCA_CONFIG_UNMAPPED}

# extract fragment statistics
FRAG_MEAN=$(awk 'NR==8 {print $6}' ../ihist.picard.txt)
FRAG_SD=$(awk 'NR==8 {print $7}' ../ihist.picard.txt)

# confirm values are non-empty
if [[ -z "$FRAG_MEAN" || -z "$FRAG_SD" ]]; then
  echo "ERROR: Could not extract FRAG_MEAN and FRAG_SD from ihist.picard.txt"
  exit 1
fi

# set read paths
FWD_READS=${SRA}/mapped/unmapped_R1.fastq
REV_READS=${SRA}/mapped/unmapped_R2.fastq
SGL_READS=${SRA}/mapped/unmapped_unpaired.fastq

sed -i "s|PE= pe |PE= pe ${FRAG_MEAN} ${FRAG_SD} ${FWD_READS} ${REV_READS}|" ${MASURCA_CONFIG_UNMAPPED}
sed -i "s|PE= se |PE= se ${FRAG_MEAN} ${FRAG_SD} ${SGL_READS}|" ${MASURCA_CONFIG_UNMAPPED}

echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca ${MASURCA_CONFIG_UNMAPPED}
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../
