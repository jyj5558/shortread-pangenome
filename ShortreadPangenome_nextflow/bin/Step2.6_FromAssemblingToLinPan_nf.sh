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

cd ${SRA}/
cd ./mapped

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in mapped unmapped; do
  cd canu_${g}/
  sortbyname.sh in=final_${g}.fa out=final_${g}_sorted.fasta length descending overwrite=T
  reformat.sh minlength=500 in=final_${g}_sorted.fasta out=final_${g}_500.fa
  cd ../
done

## Duplicate contig removal
module purge
module load biocontainers
module load bbmap

for g in mapped unmapped; do
  cd ./canu_${g}/
  dedupe.sh in=final_${g}_500.fa out=final_${g}_dedup.fa
  cd ../
done

## Renaming contig headers to avoid redundant names from different fasta files
sed 's/>tig/>map_tig/g' ./canu_mapped/final_mapped_dedup.fa > ./canu_mapped/final_mapped_renamed.fa
sed 's/>tig/>unmap_tig/g' ./canu_unmapped/final_unmapped_dedup.fa > ./canu_unmapped/final_unmapped_renamed.fa

## Assembly quality check
#module load quast
#module load busco
#export AUGUSTUS_CONFIG_PATH=$RCAC_SCRATCH/augustus/config

### Quast assessment
#echo "Canu mapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_canu_mapped --eukaryote ./canu_mapped/final_mapped_renamed.fa
#echo "Canu mapped assembly Quast finished"

#echo "Canu unmapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_canu_unmapped --eukaryote ./canu_unmapped/final_unmapped_renamed.fa
#echo "Canu unmapped assembly Quast finished"

### Busco assessment # change options accordingly for each species 
#echo "Canu mapped assembly Busco started"
# removed update data flag
#busco -f -m genome -c ${N} -o Busco_canu_mapped -l aves_odb10 -i ./canu_mapped/final_mapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
#busco -f -m genome -c ${N} -o Busco_canu_mapped -l aves_odb10 -i ./canu_mapped/final_mapped_renamed.fa --augustus --augustus_species chicken --limit 5 --datasets_version odb10
#echo "Canu mapped assembly Busco finished"

#echo "Canu unmapped assembly Busco started"
# removed update data flag
#busco -f -m genome -c ${N} -o Busco_canu_unmapped -l aves_odb10 -i ./canu_unmapped/final_unmapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
#busco -f -m genome -c ${N} -o Busco_canu_unmapped -l aves_odb10 -i ./canu_unmapped/final_unmapped_renamed.fa --augustus --augustus_species chicken --limit 5 --datasets_version odb10
#echo "Canu unmapped assembly Busco finished"


# Updating the original linear representative genome by appending all the supercontigs from unmapped reads
cat ${REF}/${GENOME}.fna ./canu_unmapped/final_unmapped_renamed.fa > ${LINPAN}/${PREFIX}_panref1.fa


# Second round of updating the linear representative genome (using unmapped reads only)
# Read mapping using Bowtie2
module purge
module load biocontainers
module load bowtie2
module load samtools
module load bbmap
module load r/4.3

bowtie2-build ${LINPAN}/${PREFIX}_panref1.fa ${LINPAN}/${PREFIX}_panref1 # now use the updated representative genome as the reference

cd ../
mkdir -p mapped2/
cd mapped2/

bowtie2 -p ${N} -I 0 -X 1000 -x ${LINPAN}/${PREFIX}_panref1 -1 ../Nremoved/merged_Nremoved_1.fq -2 ../Nremoved/merged_Nremoved_2.fq --end-to-end --sensitive -S mapped2.sam 2> bowtie2.log

samtools sort -@ ${N} -o mapping2_sorted.bam mapped2.sam

## Average insert size of the mapped paired-end reads for the MaSuRCA assembly
#reformat.sh in=mapped2.sam ihist=ihist2.txt
${APP}/jdk-21.0.1/bin/java -Duser.country=US -Duser.language=en -jar ${APP}/picard.jar CollectInsertSizeMetrics -I mapping2_sorted.bam -O ihist2.picard.txt -H ihist2.picard.pdf -M 0.05 #--DEVIATIONS # picard is better than reformat.sh

### Mapping rate, depth, and breadth
echo "mapping rate is" > ./mapping2_rate.txt
samtools flagstat -@ ${N} ./mapping2_sorted.bam >> ./mapping2_rate.txt 
echo "depth is" > ./mapping2_depth.txt
samtools depth -@ ${N} -a ./mapping2_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> ./mapping2_depth.txt 
echo "breadth is" > ./mapping2_breadth.txt
samtools depth -@ ${N} -a ./mapping2_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./mapping2_breadth.txt


# Extracting unmapped reads only
module purge
module load biocontainers
module load samtools
module load bamtools
module load biopython
module load megahit

samtools view -@ ${N} -b -f 68 mapped2.sam > unmapped2_f.bam # unmapped forward reads
samtools view -@ ${N} -b -f 132 mapped2.sam > unmapped2_r.bam # unmapped reverse reads
samtools merge -@ ${N} unmapped2_merged.bam unmapped2_f.bam unmapped2_r.bam
samtools sort -@ ${N} -n -o unmapped2_sorted.bam unmapped2_merged.bam
bamtools convert -in unmapped2_sorted.bam -out unmapped2_merged.fastq -format fastq

mv unmapped2_merged.fastq_R1.fastq unmapped2_R1.fastq
mv unmapped2_merged.fastq_R2.fastq unmapped2_R2.fastq
mv unmapped2_merged.fastq_unpaired.fastq unmapped2_unpaired.fastq
