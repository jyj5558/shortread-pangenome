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
module load biocontainers
module load biopython 
module load boost
module load megahit

cd ${SRA}/
cd mapped2/

python3 ${APP}/splitUP.py unmapped2_merged.fastq

mv unmapped2_merged.fastq_R1.fastq unmapped2_R1.fastq
mv unmapped2_merged.fastq_R2.fastq unmapped2_R2.fastq
mv unmapped2_merged.fastq_unpaired.fastq unmapped2_unpaired.fastq

# Assemblyng unmapped reads independently 
## MaSuRCA assembly using unmapped reads
mkdir -p ./masurca_unmapped2

cd ./masurca_unmapped2
cp ~/masurca_config_example ./${MASURCA_CONFIG_UNMAPPED2}
FRAG_MEAN=`awk 'NR==8 {print $6}' ../ihist2.picard.txt`
FRAG_SD=`awk 'NR==8 {print $7}' ../ihist2.picard.txt`
FWD_READS=${SRA}/mapped2/unmapped2_R1.fastq
REV_READS=${SRA}/mapped2/unmapped2_R2.fastq
SGL_READS=${SRA}/mapped2/unmapped2_unpaired.fastq

sed -i "s#PE= pe #PE= pe ${FRAG_MEAN} ${FRAG_SD} ${FWD_READS} ${REV_READS}#" ${MASURCA_CONFIG_UNMAPPED2}
sed -i "s#PE= se #PE= se ${FRAG_MEAN} ${FRAG_SD} ${SGL_READS}#" ${MASURCA_CONFIG_UNMAPPED2}

echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca ${MASURCA_CONFIG_UNMAPPED2} # configure the configuration file beforehand following the updated MaSuRCA manual
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../

## MegaHit assembly using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped2_R1.fastq -2 unmapped2_R2.fastq -r unmapped2_unpaired.fastq -o ./megahit_unmapped2 --out-prefix unmapped2_megahit 
echo "Megahit assembly of unmapped reads finished"

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd masurca_${g}/
  cd CA/
  cat primary.genome.scf.fasta alternative.genome.scf.fasta > ../${g}_masurca.fasta # a limitation in that we cannot distinguish primary and alternate contigs (of an individual or among individuals) since we pooled all the reads among samples 
  cd ../
  sortbyname.sh in=${g}_masurca.fasta out=${g}_masurca_sorted.fasta length descending overwrite=T
  reformat.sh minlength=200 in=${g}_masurca_sorted.fasta out=${g}_masurca_200.fa
  cd ../
  
  cd megahit_${g}/
  cp ${g}_megahit.contigs.fa ${g}_megahit.fasta
  sortbyname.sh in=${g}_megahit.fasta out=${g}_megahit_sorted.fasta length descending overwrite=T
  reformat.sh minlength=200 in=${g}_megahit_sorted.fasta out=${g}_megahit_200.fa
  cd ../
done


# Meta-assemblyng unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu
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
  reformat.sh minlength=500 in=final_${g}_sorted.fasta  out=final_${g}_500.fa
  cd ../
done

## Duplicate contig removal
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd ./canu_${g}/
  dedupe.sh in=final_${g}_500.fa out=final_${g}_dedup.fa
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
dedupe.sh in=unmaptig.fa out=unmaptig_dedup.fa # 0 reads were identified
dedupe.sh in=maptig.fa out=maptig_dedup.fa

cat maptig_dedup.fa unmaptig_dedup.fa > sample_contigs2.fa

## Duplicate contig removal
dedupe.sh in=sample_contigs2.fa out=sample_contigs2_dedup.fa


# Updating the original linear representative genome by appending supercontigs from unmapped reads -> Get the linear pangenome by the "iterative map-then-assemble" approach (Hu et al. 2020, Yao et al. 2015)
module purge
module load biocontainers
module load bbmap
module load samtools

## Filtering out small contigs (<10kb)
reformat.sh minlength=10000 in=unmaptig_dedup.fa out=unmaptig_10kb.fa

cat ${REF}/${GENOME}.fna ./unmaptig_10kb.fa > ${PREFIX}_panref2.fa # a.k.a. linear pangenome
 
## Duplicate contig removal
dedupe.sh in=${PREFIX}_panref2.fa out=${PREFIX}_panref2_dedup.fa # 0% duplication

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
