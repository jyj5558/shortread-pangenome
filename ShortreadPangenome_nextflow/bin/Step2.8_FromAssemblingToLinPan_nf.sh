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
MASURCA_CONFIG_UNMAPPED2=${11}

# Please note that this script was based on Hu et al. (2020) in Legume Genomics: Methods and Protocols (and its scripts) AND Yao et al. (2015) in Genome Biology
# Meta-assemblyng unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu

cd ${SRA}/
cd mapped2/
mkdir -p cd ${SRA}/mapped2/canu_unmapped2/

## Canu assmebly using unmapped contigs
cd ${SRA}/mapped2/canu_unmapped2/
echo "Canu assembly of unmapped reads started"
canu -p final_unmapped2 -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped2/masurca_unmapped2/unmapped2_masurca_200.fa ${SRA}/mapped2/megahit_unmapped2/unmapped2_megahit_200.fa 
cat final_unmapped2.contigs.fasta final_unmapped2.unassembled.fasta > final_unmapped2.fa
echo "Canu assembly of unmapped reads finishted"
cd ../

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd canu_${g}/
  sortbyname.sh in=final_${g}.fa out=final_${g}_sorted.fasta length descending overwrite=T
  reformat.sh minlength=500 in=final_${g}_sorted.fasta out=final_${g}_500.fa overwrite=T
  cd ../
done

## Duplicate contig removal
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd ./canu_${g}/
  dedupe.sh in=final_${g}_500.fa out=final_${g}_dedup.fa overwrite=T
  cd ../
done

## Renaming contig headers to avoid redundant names from different fasta files
sed 's/>tig/>unmap2_tig/g' ./canu_unmapped2/final_unmapped2_dedup.fa > ./canu_unmapped2/final_unmapped2_renamed.fa


# Confirming unmapped supercontigs are unmappable and removing mappable ones
mkdir -p ${CONTIG}

module purge
module load biocontainers
module load minimap2
module load bamtools
module load samtools

cd ../ 
cat ./mapped/canu_unmapped/final_unmapped_renamed.fa ./mapped2/canu_unmapped2/final_unmapped2_renamed.fa > ${CONTIG}/final_unmapped_all.fa # gather unmapped contigs
minimap2 -ax asm5 ${REF}/${GENOME}.fna ${CONTIG}/final_unmapped_all.fa > ${CONTIG}/unmaptig_mapped.sam
samtools view -@ ${N} -f 4 ${CONTIG}/unmaptig_mapped.sam > ${CONTIG}/unmaptig.sam # confirmed unmapped contigs
samtools sort -n -@ ${N} -o ${CONTIG}/unmaptig.bam ${CONTIG}/unmaptig.sam 
bamtools convert -in ${CONTIG}/unmaptig.bam -format fasta -out ${CONTIG}/unmaptig.fa


# Confirming mapped contigs are mappable and removing unmappable ones (technical artifacts - kimeric and alternative contigs, etc.)
minimap2 -ax asm5 ${REF}/${GENOME}.fna ./mapped/canu_mapped/final_mapped_renamed.fa > ${CONTIG}/maptig_mapped.sam
samtools sort -@ ${N} -o ${CONTIG}/maptig_mapped.bam ${CONTIG}/maptig_mapped.sam 
samtools flagstat -@ ${N} ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.mappingrate
samtools depth -@ ${N} -a ${CONTIG}/maptig_mapped.bam | awk '{c++;s+=$3}END{print s/c}' > ${CONTIG}/maptig.depth
samtools depth -@ ${N} -a ${CONTIG}/maptig_mapped.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${CONTIG}/maptig.breadth
samtools coverage ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.stat

samtools view -@ ${N} -b -F 4 ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.bam # confirmed mapped contigs
bamtools convert -in ${CONTIG}/maptig.bam -format fasta -out ${CONTIG}/maptig.fa


# Pooling all the supercontigs (mapped or unmapped) to generate a sample-derived contig collection
cd ${CONTIG}/

module purge
module load biocontainers
module load bbmap

## Check duplicates
dedupe.sh in=unmaptig.fa out=unmaptig_dedup.fa overwrite=T # 0 reads were identified
dedupe.sh in=maptig.fa out=maptig_dedup.fa overwrite=T

cat maptig_dedup.fa unmaptig_dedup.fa > sample_contigs2.fa

## Duplicate contig removal
dedupe.sh in=sample_contigs2.fa out=sample_contigs2_dedup.fa overwrite=T


# Updating the original linear representative genome by appending supercontigs from unmapped reads -> Get the linear pangenome by the "iterative map-then-assemble" approach (Hu et al. 2020, Yao et al. 2015)
module purge
module load biocontainers
module load bbmap
module load samtools

## Filtering out small contigs (<10kb)
reformat.sh minlength=10000 in=unmaptig_dedup.fa out=unmaptig_10kb.fa overwrite=T

cat ${REF}/${GENOME}.fna ./unmaptig_10kb.fa > ${PREFIX}_panref2.fa # a.k.a. linear pangenome
 
## Duplicate contig removal
dedupe.sh in=${PREFIX}_panref2.fa out=${PREFIX}_panref2_dedup.fa overwrite=T # 0% duplication

## Sorting the updated linear representative genome (= linear pangenome)
sortbyname.sh in=${PREFIX}_panref2_dedup.fa out=${PREFIX}_panref2_sorted.fa length descending overwrite=T

## Indexing the updated linear representative genome (= linear pangenome)
samtools faidx ${PREFIX}_panref2_sorted.fa  # check contig length distribution to decide the length cutoff later (>10kb can include all the reference contigs + alpha) 

## Assembly quality check
module load biocontainers
module load quast
module load busco
export AUGUSTUS_CONFIG_PATH=$RCAC_SCRATCH/augustus/config

quast.py --threads ${N} -o Quast_panref2 --eukaryote ./${PREFIX}_panref2_sorted.fa
#### this gave an error for -l busco -f -m genome -c ${N} -o Busco_panref2 -l --auto-lineage-euk -i ./${PREFIX}_panref2_sorted.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
busco -f -m genome -c ${N} -o Busco_panref2 -l aves_odb10 -i ./${PREFIX}_panref2_sorted.fa --augustus --augustus_species chicken --limit 5 --datasets_version odb10
