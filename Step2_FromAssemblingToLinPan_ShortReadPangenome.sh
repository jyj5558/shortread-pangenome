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

# Please note that this script was based on Hu et al. (2020) in Legume Genomics: Methods and Protocols (and its scripts) AND Yao et al. (2015) in Genome Biology

cd ${SRA}/

# Read mapping using Bowtie2 
module purge
module load biocontainers
module load bowtie2
module load bbmap
module load samtools

bowtie2-build ${REF}/${GENOME}.fna ${REF}/${GENOME}

mkdir -p ./mapped
cd ./mapped

bowtie2 -p ${N} -I 0 -X 1000 -x ${REF}/${GENOME} -1 ../Nremoved/merged_Nremoved_1.fq -2 ../Nremoved/merged_Nremoved_2.fq --end-to-end --sensitive -S mapped.sam 2> bowtie.log

## Summary stats on bam files
samtools sort -@ ${N} -o mapping_sorted.bam mapped.sam

### Average insert size of the mapped paired-end reads for the MaSuRCA assembly
#reformat.sh in=mapping_sorted.bam ihist=ihist.txt
${APP}/jdk-21.0.1/bin/java -Duser.country=US -Duser.language=en -jar ${APP}/picard.jar CollectInsertSizeMetrics -I mapping_sorted.bam -O ihist.picard.txt -H ihist.picard.pdf -M 0.05 #--DEVIATIONS # picard is better than reformat.sh

### Mapping rate, depth, and breadth
echo "mapping rate is" > ./mapping_rate.txt
samtools flagstat -@ ${N} ./mapping_sorted.bam >> ./mapping_rate.txt 
echo "depth is" > ./mapping_depth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> ./mapping_depth.txt 
echo "breadth is" > ./mapping_breadth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./mapping_breadth.txt


# Extracting mapped and unmapped reads (on the linear representative genome) separately
module purge
module load biocontainers
module load samtools
module load bamtools
module load biopython

samtools view -@ ${N} -b -F 4 -f 64 mapped.sam > mapped_f.bam # mapped forward reads
samtools view -@ ${N} -b -F 4 -f 128 mapped.sam > mapped_r.bam # mapped reverse reads
samtools merge -@ ${N} mapped_merged.bam mapped_f.bam mapped_r.bam
samtools sort -@ ${N} -n -o mapped_sorted.bam mapped_merged.bam
bamtools convert -in mapped_sorted.bam -out mapped_merged.fastq -format fastq
python ~/splitUP.py mapped_merged.fastq # splitUP.py script from Hu et al. (2020)
# If the mapped_merged.fastq file is too large, it will cause an out-of-memory issue. In that case, split the file by 4n lines first and then concatenate them later like below.

# cp mapped_merged.fastq mapped_merged_cp.fastq # copy the file for redundancy
# split -l 100000000 --numeric-suffixes=1 --additional-suffix=.fastq --verbose mapped_merged.fastq mapped_merged_

# for i in {01..52}; do
# echo "printing tail ${i}"
# tail -4 mapped_merged_${i}.fastq # manually check if all the files are split at the end of a fastq-formatted-sequence line
# done

# for i in {01..52}; do
# python ~/splitUP.py mapped_merged_${i}.fastq
# done

# For the unpaired, concatenate them and rerun splitUP.py again to avoid artefacts due to split command.

# for i in {01..52}; do
# cat mapped_merged_${i}.fastq_unpaired.fastq >> mapped_merged_53.fastq # Do not forget later that this '53' file is a concatenated file of all '1'-'52' unpaired files.
# done 

# python ~/splitUP.py mapped_merged_53.fastq
# mv mapped_merged_53.fastq_unpaired.fastq mapped_merged.fastq_unpaired.fastq

# for i in {01..53}; do
# cat mapped_merged_${i}.fastq_R1.fastq >> mapped_merged.fastq_R1.fastq
# cat mapped_merged_${i}.fastq_R2.fastq >> mapped_merged.fastq_R2.fastq
# done

samtools view -@ ${N} -b -f 68 mapped.sam > unmapped_f.bam # unmapped forward reads
samtools view -@ ${N} -b -f 132 mapped.sam > unmapped_r.bam # unmapped reverse reads
samtools merge -@ ${N} unmapped_merged.bam unmapped_f.bam unmapped_r.bam
samtools sort -@ ${N} -n -o unmapped_sorted.bam unmapped_merged.bam
bamtools convert -in unmapped_sorted.bam -out unmapped_merged.fastq -format fastq
python ~/splitUP.py unmapped_merged.fastq

mv mapped_merged.fastq_R1.fastq mapped_R1.fastq
mv mapped_merged.fastq_R2.fastq mapped_R2.fastq
mv mapped_merged.fastq_unpaired.fastq mapped_unpaired.fastq
mv unmapped_merged.fastq_R1.fastq unmapped_R1.fastq
mv unmapped_merged.fastq_R2.fastq unmapped_R2.fastq
mv unmapped_merged.fastq_unpaired.fastq unmapped_unpaired.fastq


# Assembling mapped and unmapped reads independently, redundantly (MaSuRCA and MegaHit) modifying the approach of Yao et al. (2015)
## MaSuRCA assembly using mapped reads
module purge
module load biocontainers
module load megahit

MASURCA=/home/jeon96/app/MaSuRCA-4.1.0/bin

mkdir -p ./masurca_mapped
mkdir -p ./masurca_unmapped

cd ./masurca_mapped
echo "Masurca assembly of mapped reads started"
${MASURCA}/masurca masurca_config_mapped.txt # configure the configuration file beforehand following updated MaSuRCA manual
bash assemble.sh
echo "Masurca assembly of mapped reads finished"
cd ../

## MaSuRCA assembly using unmapped reads
cd ./masurca_unmapped
echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca masurca_config_unmapped.txt # configure the configuration file beforehand following updated MaSuRCA manual
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../

## MegaHit assembly using mapped reads
echo "Megahit assembly of mapped reads started"
megahit -t ${N} -1 mapped_R1.fastq -2 mapped_R2.fastq -r mapped_unpaired.fastq -o ./megahit_mapped --out-prefix mapped_megahit 
echo "Megahit assembly of mapped reads finished"

## MegaHit assembly using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped_R1.fastq -2 unmapped_R2.fastq -r unmapped_unpaired.fastq -o ./megahit_unmapped --out-prefix unmapped_megahit 
echo "Megahit assembly of unmapped reads finished"

## Small contig filtering
module purge
module load biocontainers
module load bbmap

for g in mapped unmapped; do
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

## Quast assessment - skip
#module load biocontainers
#module load quast

#echo "Masurca mapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_masurca_mapped --eukaryote ./masurca_mapped/mapped_masurca_200.fa
#echo "Masurca mapped assembly Quast finished"
 
#echo "Masurca unmapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_masurca_unmapped --eukaryote ./masurca_unmapped/unmapped_masurca_200.fa
#echo "Masurca mapped assembly Quast finished"

#echo "Megahit mapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_megahit_mapped --eukaryote ./megahit_mapped/mapped_megahit_200.fa
#echo "Megahit mapped assembly Quast finished"
 
#echo "Megahit unmapped assembly Quast started" 
#quast.py --threads ${N} -o Quast_megahit_unmapped --eukaryote ./megahit_unmapped/unmapped_megahit_200.fa
#echo "Megahit unmapped assembly Quast finished"


# Meta-assembling mapped and unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu

mkdir -p ./canu_mapped/
mkdir -p ./canu_unmapped/

## Canu assmebly using mapped contigs
cd ./canu_mapped
echo "Canu assembly of mapped reads started"
canu -p final_mapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_mapped/mapped_masurca_200.fa ${SRA}/mapped/megahit_mapped/mapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_mapped.contigs.fasta final_mapped.unassembled.fasta > final_mapped.fa
echo "Canu assembly of mapped reads finishted"
cd ../

## Canu assmebly using unmapped contigs
cd ./canu_unmapped
echo "Canu assembly of unmapped reads started"
canu -p final_unmapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_unmapped/unmapped_masurca_200.fa ${SRA}/mapped/megahit_unmapped/unmapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_unmapped.contigs.fasta final_unmapped.unassembled.fasta > final_unmapped.fa
echo "Canu assembly of unmapped reads finishted"
cd ../

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
module load quast
module load busco
export AUGUSTUS_CONFIG_PATH=$RCAC_SCRATCH/augustus/config

### Quast assessment
echo "Canu mapped assembly Quast started" 
quast.py --threads ${N} -o Quast_canu_mapped --eukaryote ./canu_mapped/final_mapped_renamed.fa
echo "Canu mapped assembly Quast finished"

echo "Canu unmapped assembly Quast started" 
quast.py --threads ${N} -o Quast_canu_unmapped --eukaryote ./canu_unmapped/final_unmapped_renamed.fa
echo "Canu unmapped assembly Quast finished"

### Busco assessment # change options accordingly for each species 
echo "Canu mapped assembly Busco started"
busco -f -m genome -c ${N} -o Busco_canu_mapped -l aves_odb10 -i ./canu_mapped/final_mapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
echo "Canu mapped assembly Busco finished"

echo "Canu unmapped assembly Busco started"
busco -f -m genome -c ${N} -o Busco_canu_unmapped -l aves_odb10 -i ./canu_unmapped/final_unmapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
echo "Canu unmapped assembly Busco finished"


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
python ~/splitUP.py unmapped2_merged.fastq

mv unmapped2_merged.fastq_R1.fastq unmapped2_R1.fastq
mv unmapped2_merged.fastq_R2.fastq unmapped2_R2.fastq
mv unmapped2_merged.fastq_unpaired.fastq unmapped2_unpaired.fastq


# Assemblyng unmapped reads independently 
## MaSuRCA assembly using unmapped reads
mkdir -p ./masurca_unmapped2

cd ./masurca_unmapped2
echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca masurca_config_unmapped2.txt # configure the configuration file beforehand following the updated MaSuRCA manual
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
mkdir -p ./canu_unmapped2/

## Canu assmebly using unmapped contigs
cd ./canu_unmapped2/
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
busco -f -m genome -c ${N} -o Busco_panref2 -l --auto-lineage-euk -i ./${PREFIX}_panref2_sorted.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10