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

cd ${SRA}/

# 0. Downloading published sequencing data from the GenBank
## 0.1. short reads
mkdir -p ./raw/
cd ./raw/

module load biocontainers
module load sra-tools
module load cutadapt
module load hifiasm
module load purge_dups

cat ../../SRR_Acc_List.txt | while read g # need the "SRR_Acc_List.txt" file beforehand that includes accession numbers of SRA file to download
do
mkdir ${g}
cd ${g}
prefetch ${g} --max-size 500GB -O ./
fasterq-dump ${g} -e ${N} --progress
find . -name '*.fastq' -exec mv {} ../ \;
cd ../
rm -r ${g}
done

cd ../../

## 0.2. HiFi reads
mkdir -p ./hifi/raw/
cd ./hifi/raw/

cat ../../hifi_SRR_Acc_List.txt | while read g # need the "hifi_SRR_Acc_List.txt" file beforehand that includes accession numbers of SRA file to download
do
mkdir ${g}
cd ${g}
prefetch ${g} --max-size 500GB -O ./
fasterq-dump ${g} -e ${N} --progress
find . -name '*.fastq' -exec mv {} ../ \;
cd ../
rm -r ${g}
done

cd ../

### 0.2.1. decontamination by adapters on HiFi raw reads.
mkdir -p ./cleaned/
cd ./cleaned/

cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588214_trimmed.ccs.fq SRR22588214.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588215_trimmed.ccs.fq SRR22588215.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588216_trimmed.ccs.fq SRR22588216.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588217_trimmed.ccs.fq SRR22588217.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588218_trimmed.ccs.fq SRR22588218.ccs.fq -j 64 --revcomp -e 0.01

cd ../

### 0.2.2. Assembly with Hifiasm.
mkdir -p ./assembled/
cd ./assembled/

hifiasm --primary -o SRR22588214_trimmed.ccs.asm -t 64 SRR22588214_trimmed.ccs.fq
hifiasm --primary -o SRR22588215_trimmed.ccs.asm -t 64 SRR22588215_trimmed.ccs.fq
hifiasm --primary -o SRR22588216_trimmed.ccs.asm -t 64 SRR22588216_trimmed.ccs.fq
hifiasm --primary -o SRR22588217_trimmed.ccs.asm -t 64 SRR22588217_trimmed.ccs.fq
hifiasm --primary -o SRR22588218_trimmed.ccs.asm -t 64 SRR22588218_trimmed.ccs.fq

### 0.2.3. Separate primary and alternate assemblies.

# primary
awk '/^S/{print ">"$2;print $3}' SRR22588214_trimmed.ccs.asm.p_ctg.gfa > SRR22588214_primary.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588215_trimmed.ccs.asm.p_ctg.gfa > SRR22588215_primary.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588216_trimmed.ccs.asm.p_ctg.gfa > SRR22588216_primary.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588217_trimmed.ccs.asm.p_ctg.gfa > SRR22588217_primary.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588218_trimmed.ccs.asm.p_ctg.gfa > SRR22588218_primary.fasta 
# alternate
awk '/^S/{print ">"$2;print $3}' SRR22588214_trimmed.ccs.asm.a_ctg.gfa > SRR22588214_alternate.fasta 
awk '/^S/{print ">"$2;print $3}' SRR22588215_trimmed.ccs.asm.a_ctg.gfa > SRR22588215_alternate.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588216_trimmed.ccs.asm.a_ctg.gfa > SRR22588216_alternate.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588217_trimmed.ccs.asm.a_ctg.gfa > SRR22588217_alternate.fasta
awk '/^S/{print ">"$2;print $3}' SRR22588218_trimmed.ccs.asm.a_ctg.gfa > SRR22588218_alternate.fasta

cd ../raw/

### 0.2.4. run merqury to generate the Meryl database with the _submit_build.sh script included in the package modified according to your cluster. The script can be found here https://github.com/marbl/merqury/blob/master/_submit_build.sh .

bash _submit_build_mod.sh 21 input_SRR22588214.fofn SRR22588214 
bash _submit_build_mod.sh 21 input_SRR22588215.fofn SRR22588215 
bash _submit_build_mod.sh 21 input_SRR22588216.fofn SRR22588216 
bash _submit_build_mod.sh 21 input_SRR22588217.fofn SRR22588217 
bash _submit_build_mod.sh 21 input_SRR22588218.fofn SRR22588218 

cd ../

### 0.2.5. run Genomescope2.0 online (http://qb.cshl.edu/genomescope/genomescope2.0/) with 31 k-mers, uploading the .histo file outputted by the previous script.

### 0.2.6. run purge dups with custom cutoffs for HiFi reads.The scripts can be found here https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups
cd ./assembled/
mkdir -p ./purged/
cd ./purged/

#### 0.2.6.1. add the -xasm20 option for HiFi reads to minimap2.sh and minimap2_self.sh 

#### 0.2.6.2. calculate custom cutoffs starting from the kcov computed by Genomescope2.0.

#SRR22588214
#kcov by Genomescope 2.0 = 10.1 
#value for -m option = kcov*1.5 = 15.15
#value for -u option = (value for -m option)*3 = 45.45

#SRR22588215
#kcov by Genomescope 2.0 = 9.81
#value for -m option = kcov*1.5 = 14.715
#value for -u option = (value for -m option)*3 = 44.145

#SRR22588216
#kcov by Genomescope 2.0 = 12.3
#value for -m option = kcov*1.5 = 18.45
#value for -u option = (value for -m option)*3 = 55.35

#SRR22588217
#kcov by Genomescope 2.0 = 7.47
#value for -m option = kcov*1.5 = 11.205
#value for -u option = (value for -m option)*3 = 33.615

#SRR22588218
#kcov by Genomescope 2.0 = 16.5
#value for -m option = kcov*1.5 = 24.75
#value for -u option = (value for -m option)*3 = 74.25

#### 0.2.6.3. modify the script purge_dups.sh with the custom cutoffs; add the modified minimap2.sh, minimap2_self.sh and purge_dups.sh scripts into _submit_purge_dups.sh and run it
# in the purge_dups.sh script, replace to this when running SRR22588214: calcuts -m 15.15 -u 45.45 PB.stat > cutoffs 2>calcults.log
mkdir -p ./SRR22588214/
mkdir -p ./SRR22588214/primary
mkdir -p ./SRR22588214/alternate
cp ../../raw/input_SRR22588214.fofn ./input.fofn
cp ../SRR22588214_primary.fasta SRR22588214_primary.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588214_primary.fasta input.fofn  # Repeat also for alternate assemblies. Modify the script to accommodate a custome-named fofn file.
mv purged.fa SRR22588214_primary_purged.fasta 
mv SRR22588214.* SRR22588214_* ./SRR22588214/primary
mv *.log PB* cutoffs dups.bed hap.fa  ./SRR22588214/primary

cp ../SRR22588214_alternate.fasta SRR22588214_alternate.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588214_alternate.fasta input.fofn  
mv purged.fa SRR22588214_alternate_purged.fasta
mv SRR22588214.* SRR22588214_* ./SRR22588214/alternate
mv *.log PB* cutoffs dups.bed hap.fa  ./SRR22588214/alternate
mv input.fofn ./SRR22588214/

# in the purge_dups.sh script, replace to this when running SRR22588215: calcuts -m 14.715 -u 44.145 PB.stat > cutoffs 2>calcults.log 
mkdir -p ./SRR22588215/
mkdir -p ./SRR22588215/primary
mkdir -p ./SRR22588215/alternate
cp ../../raw/input_SRR22588215.fofn ./input.fofn  
cp ../SRR22588215_primary.fasta SRR22588215_primary.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588215_primary.fasta input.fofn  # Repeat also for alternate assemblies.
mv purged.fa SRR22588215_primary_purged.fasta 
mv SRR22588215.* SRR22588215_* ./SRR22588215/primary
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588215/primary

cp ../SRR22588215_alternate.fasta SRR22588215_alternate.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588215_alternate.fasta input.fofn  
mv purged.fa SRR22588215_alternate_purged.fasta 
mv SRR22588215.* SRR22588215_* ./SRR22588215/alternate
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588215/alternate
mv input.fofn ./SRR22588215/

# in the purge_dups.sh script, replace to this when running SRR22588216: calcuts -m 18.45 -u 55.35 PB.stat > cutoffs 2>calcults.log   
mkdir -p ./SRR22588216/
mkdir -p ./SRR22588216/primary
mkdir -p ./SRR22588216/alternate
cp ../../raw/input_SRR22588216.fofn ./input.fofn  
cp ../SRR22588216_primary.fasta SRR22588216_primary.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588216_primary.fasta input.fofn  # Repeat also for alternate assemblies.
mv purged.fa SRR22588216_primary_purged.fasta 
mv SRR22588216.* SRR22588216_* ./SRR22588216/primary
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588216/primary

cp ../SRR22588216_alternate.fasta SRR22588216_alternate.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588216_alternate.fasta input.fofn  
mv purged.fa SRR22588216_alternate_purged.fasta 
mv SRR22588216.* SRR22588216_* ./SRR22588216/alternate
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588216/alternate
mv input.fofn ./SRR22588216/

# in the purge_dups.sh script, replace to this when running SRR22588217: calcuts -m 11.205 -u 33.615 PB.stat > cutoffs 2>calcults.log  
mkdir -p ./SRR22588217/
mkdir -p ./SRR22588217/primary
mkdir -p ./SRR22588217/alternate
cp ../../raw/input_SRR22588217.fofn ./input.fofn  
cp ../SRR22588217_primary.fasta SRR22588217_primary.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588217_primary.fasta input.fofn  # Repeat also for alternate assemblies.
mv purged.fa SRR22588217_primary_purged.fasta 
mv SRR22588217.* SRR22588217_* ./SRR22588217/primary
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588217/primary

cp ../SRR22588217_alternate.fasta SRR22588217_alternate.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588217_alternate.fasta input.fofn  
mv purged.fa SRR22588217_alternate_purged.fasta 
mv SRR22588217.* SRR22588217_* ./SRR22588217/alternate
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588217/alternate
mv input.fofn ./SRR22588217/

# in the purge_dups.sh script, replace to this when running SRR22588218: calcuts -m 24.75 -u 74.25 PB.stat > cutoffs 2>calcults.log  
mkdir -p ./SRR22588218/
mkdir -p ./SRR22588218/primary
mkdir -p ./SRR22588218/alternate
cp ../../raw/input_SRR22588218.fofn ./input.fofn  
cp ../SRR22588218_primary.fasta SRR22588218_primary.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588218_primary.fasta input.fofn  # Repeat also for alternate assemblies.
mv purged.fa SRR22588218_primary_purged.fasta 
mv SRR22588218.* SRR22588218_* ./SRR22588218/primary
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588218/primary

cp ../SRR22588218_alternate.fasta SRR22588218_alternate.fasta
bash ${APP}/vgp/_submit_purge_dups.sh SRR22588218_alternate.fasta input.fofn  
mv purged.fa SRR22588218_alternate_purged.fasta 
mv SRR22588218.* SRR22588218_* ./SRR22588218/alternate
mv *.log PB* cutoffs dups.bed hap.fa ./SRR22588218/alternate
mv input.fofn ./SRR22588218/

#0.2.7. run asm_stats.sh to obtain statistics before and after purging. 
${APP}/vgp/stats/asm_stats.sh SRR22588214_primary.fasta 1079274132 c #the genome size is the mean predicted size from Genomescope2.0
${APP}/vgp/stats/asm_stats.sh SRR22588214_primary_purged.fasta 1079274132 c

cd ${BASE}

## 0.3. Simulating short read data from a published pangenome - when working with simulated reads only
#module load biocontainers
#module load bbmap

mkdir -p ./sim/raw/temp/
cd ./sim/raw/

## 0.3.1. Indexing the pangenome haplotypes
#PAN_ID=barn_swallow_pangenome

#vg mod -t ${N} -X 256 ${PAN}/${PAN_ID}.vg > ${PAN}/${PAN_ID}_chopped.vg
#vg index -x ${PAN}/${PAN_ID}.xg -b temp/ --threads ${N} -p ${PAN}/${PAN_ID}_chopped.vg #-g ${PAN}/${PAN_ID}.gcsa 

## 0.3.2. Simulating short reads
NGSNGS=/depot/fnrdewoody/apps/NGSNGS
PURGED=/scratch/negishi/jeon96/swallow/original/hifi/assembled/purged

for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    ${NGSNGS}/ngsngs -i ${PURGED}/${i}/${j}/${i}_${j}_purged.fasta -c 5 -ld Norm,500,50 -seq paired-end -f fq.gz -o sim${k}_${i}_${j} -qs 40 -cl 150 -s ${k} -t ${N} -t2 ${N}
    done 
  done
done # 20 samples in total with 5x coverage for both primary and alternate

#for FQ in `cat ${BASE}/sra/SRR_Acc_List.txt`; do
#seed=$(echo ${FQ} | cut -d "R" -f 3)
#seed1=`expr ${seed} \* 10 + 1`
#seed2=`expr ${seed} \* 10 + 2`
#vg sim -t ${N} -x ${PAN}/${PAN_ID}.xg -N -n 25000000 -l 150 -p 500 -v 50 -a -r -s ${seed1} > ${FQ}_sim1.gam -F ${BASE}/sra/raw/${FQ}_1.fastq # -n 25000000 for ~7x coverage data
#vg sim -t ${N} -x ${PAN}/${PAN_ID}.xg -N -n 25000000 -l 150 -p 500 -v 50 -a -r -s ${seed2} > ${FQ}_sim2.gam -F ${BASE}/sra/raw/${FQ}_2.fastq # -n 25000000 for ~7x coverage data
#vg view --threads ${N} -a ${FQ}_sim1.gam -X -i > ${FQ}_sim1.fastq
#vg view --threads ${N} -a ${FQ}_sim2.gam -X -i > ${FQ}_sim2.fastq

### 0.3.3. Converting to paired fastq files
#reformat.sh in=${FQ}_sim1.fastq out1=${FQ}_sim1_1.fastq.gz out2=${FQ}_sim1_2.fastq.gz
#reformat.sh in=${FQ}_sim2.fastq out1=${FQ}_sim2_1.fastq.gz out2=${FQ}_sim2_2.fastq.gz

done


# 1. Data Quality Control - may be skipped if using simulated reads of high quality without adapter sequences
cd ${SRA}

module purge
module load biocontainers
module load kraken2
module load biopython
module load cutadapt
module load fastqc
module load trim-galore

## 1.1. Initial quality check using Fastqc
mkdir -p ./qc/
cd ./qc/

for g in `ls -lt ../raw/ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  fastqc ../raw/${g}R1_001.fastq --extract --quiet
  fastqc ../raw/${g}R2_001.fastq --extract --quiet
  
  cat ${g}1_fastqc/summary.txt ${g}2_fastqc/summary.txt > ${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}fastqc_summary.txt)
  echo "raw"
  echo "$FAIL"
  rm -r ${g}?_fastqc* 
done
cd ../

## 1.2. Adapter & low quality reads removal using Trim_galore
mkdir -p ./cleaned/
cd ./cleaned/

for g in `ls -lt ../raw/ | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ./ --paired ../raw/${g}R1_001.fastq ../raw/${g}R2_001.fastq
  cat ${g}1.fastq_trimming_report.txt ${g}2.fastq_trimming_report.txt > ${g}fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
done
cd ../

## 1.3. Contamination removal using Kraken 
mkdir -p ./filtered/
cd ./filtered/

DB=/scratch/bell/dewoody/LEPC/JY-test/db/contam_lib # make a "contam_lib" beforehand with potential contaminant sequences
KRAKEN2_NUM_THREADS=64
export KRAKEN2_DB_PATH="/scratch/bell/dewoody/LEPC/JY-test/db/contam_lib/"

mkdir -p ./filtered/
cd ./filtered/

for g in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
  do
  kraken2 --db ${DB} --quick --threads ${KRAKEN2_NUM_THREADS} --output ${g}cseq --paired --unclassified-out ${g}filtered#.fq ${SRA}/cleaned/${g}R1_001_val_1.fq.gz ${SRA}/cleaned/${g}R2_001_val_2.fq.gz --report ${i}_kraken_report.txt --use-names
done
cd ../

## 1.4. N removal using Remove_reads_with_Ns.py script
mkdir -p ./Nremoved/
cd ./Nremoved/

for i in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`
  do
  python ~/Remove_reads_with_Ns.py ../filtered/${i}_*_filtered_1.fq ../filtered/${i}_*_filtered_2.fq
done

for i in `ls -lt ../cleaned | grep "fq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`
  do
  cat ./${i}*_1.fq* >> merged_Nremoved_1.fq
  cat ./${i}*_2.fq* >> merged_Nremoved_2.fq
done
cd ../


# 2. Read mapping using Bowtie2
module purge
module load biocontainers
module load bowtie2
module load bbmap
module load samtools

bowtie2-build ${REF}/${GENOME}.fna ${REF}/${GENOME}

mkdir -p ./mapped
cd ./mapped

bowtie2 -p ${N} -I 0 -X 1000 -x ${REF}/${GENOME} -1 ../Nremoved/merged_Nremoved_1.fq -2 ../Nremoved/merged_Nremoved_2.fq --end-to-end --sensitive -S mapped.sam 2> bowtie.log

## 2.1. Summary stats on bam files
samtools sort -@ ${N} -o mapping_sorted.bam mapped.sam

### 2.1.1. Average insert size of the mapped paired-end reads for the MaSuRCA assembly
#reformat.sh in=mapping_sorted.bam ihist=ihist.txt
${APP}/jdk-21.0.1/bin/java -Duser.country=US -Duser.language=en -jar ${APP}/picard.jar CollectInsertSizeMetrics -I mapping_sorted.bam -O ihist.picard.txt -H ihist.picard.pdf -M 0.05 #--DEVIATIONS # picard is better

### 2.1.2. Mapping rate, depth, and breadth
echo "mapping rate is" > ./mapping_rate.txt
samtools flagstat -@ ${N} ./mapping_sorted.bam >> ./mapping_rate.txt 
echo "depth is" > ./mapping_depth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> ./mapping_depth.txt 
echo "breadth is" > ./mapping_breadth.txt
samtools depth -@ ${N} -a ./mapping_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./mapping_breadth.txt


# 3. Extraction of mapped and unmapped reads (on the linear reference genome) separately
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
python ~/splitUP.py mapped_merged.fastq
# If the mapped_merged.fastq file is too large, it will cause an out-of-memory issue. In that case, split the file by 4n lines first and then concatenate them later.
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


# 4. Independent assemblies of mapped and unmapped reads
## 4.1. MaSuRCA mapping using mapped reads
module purge
module load biocontainers
module load megahit

MASURCA=/home/jeon96/app/MaSuRCA-4.1.0/bin

mkdir -p ./masurca_mapped
mkdir -p ./masurca_unmapped

cd ./masurca_mapped
echo "Masurca assembly of mapped reads started"
${MASURCA}/masurca masurca_config_mapped.txt # configure the configuration file beforehand
bash assemble.sh
echo "Masurca assembly of mapped reads finished"
cd ../

## 4.2. MaSuRCA mapping using unmapped reads
cd ./masurca_unmapped
echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca masurca_config_unmapped.txt # configure the configuration file beforehand
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../

## 4.3. MegaHit mapping using mapped reads
echo "Megahit assembly of mapped reads started"
megahit -t ${N} -1 mapped_R1.fastq -2 mapped_R2.fastq -r mapped_unpaired.fastq -o ./megahit_mapped --out-prefix mapped_megahit 
echo "Megahit assembly of mapped reads finished"

## 4.4. MegaHit mapping using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped_R1.fastq -2 unmapped_R2.fastq -r unmapped_unpaired.fastq -o ./megahit_unmapped --out-prefix unmapped_megahit 
echo "Megahit assembly of unmapped reads finished"

## 4.5. Small contig filtering
module purge
module load biocontainers
module load bbmap
#module load bioawk

for g in mapped unmapped; do
  cd masurca_${g}/
  cd CA/
  cat primary.genome.scf.fasta alternative.genome.scf.fasta > ../${g}_masurca.fasta # a limitation in that we cannot distinguish primary and alternate contigs (of an individual or among individuals) since we pooled all the reads among samples
  cd ../
  sortbyname.sh in=${g}_masurca.fasta out=${g}_masurca_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 200) { print ">"$name; print $seq }}' ${g}_masurca_sorted.fasta  > ${g}_masurca_200.fa
  reformat.sh minlength=200 in=${g}_masurca_sorted.fasta out=${g}_masurca_200.fa
  cd ../
  
  cd megahit_${g}/
  cp ${g}_megahit.contigs.fa ${g}_megahit.fasta
  sortbyname.sh in=${g}_megahit.fasta out=${g}_megahit_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 200) { print ">"$name; print $seq }}' ${g}_megahit_sorted.fasta  > ${g}_megahit_200.fa
  reformat.sh minlength=200 in=${g}_megahit_sorted.fasta out=${g}_megahit_200.fa  
  cd ../
done

## 4.5. Quast assessment - skip
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


# 5. Meta-assemblies of mapped and unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu

mkdir -p ./canu_mapped/
mkdir -p ./canu_unmapped/

## 5.1. Canu assmebly using mapped contigs

cd ./canu_mapped
echo "Canu assembly of mapped reads started"
canu -p final_mapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_mapped/mapped_masurca_200.fa ${SRA}/mapped/megahit_mapped/mapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_mapped.contigs.fasta final_mapped.unassembled.fasta > final_mapped.fa
echo "Canu assembly of mapped reads finishted"
cd ../

## 5.2. Canu assmebly using unmapped contigs
cd ./canu_unmapped
echo "Canu assembly of unmapped reads started"
canu -p final_unmapped -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -stopOnLowCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped/masurca_unmapped/unmapped_masurca_200.fa ${SRA}/mapped/megahit_unmapped/unmapped_megahit_200.fa # use "pacbio-hifi" option since its sequencing error rate is the most similar to shor-read data among long-read data types but with a higher overlap error rate allowed
cat final_unmapped.contigs.fasta final_unmapped.unassembled.fasta > final_unmapped.fa
echo "Canu assembly of unmapped reads finishted"
cd ../

## 5.3. Small contig filtering
module purge
module load biocontainers
module load bbmap
#module load bioawk

for g in mapped unmapped; do
  cd canu_${g}/
  sortbyname.sh in=final_${g}.fa out=final_${g}_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' final_${g}_sorted.fasta  > final_${g}_500.fa
  reformat.sh minlength=500 in=final_${g}_sorted.fasta out=final_${g}_500.fa
  cd ../
done

## 5.4. Duplicate contig removal
module purge
module load biocontainers
module load bbmap

for g in mapped unmapped; do
  cd ./canu_${g}/
  dedupe.sh in=final_${g}_500.fa out=final_${g}_dedup.fa
  cd ../
done

## 5.5. Renaming contig headers to avoid redundant names from different fasta files
sed 's/>tig/>map_tig/g' ./canu_mapped/final_mapped_dedup.fa > ./canu_mapped/final_mapped_renamed.fa
sed 's/>tig/>unmap_tig/g' ./canu_unmapped/final_unmapped_dedup.fa > ./canu_unmapped/final_unmapped_renamed.fa

## 5.6. Assembly quality check
module load quast
module load busco
export AUGUSTUS_CONFIG_PATH=$RCAC_SCRATCH/augustus/config

### 5.6.1. Quast assessment
echo "Canu mapped assembly Quast started" 
quast.py --threads ${N} -o Quast_canu_mapped --eukaryote ./canu_mapped/final_mapped_renamed.fa
echo "Canu mapped assembly Quast finished"

echo "Canu unmapped assembly Quast started" 
quast.py --threads ${N} -o Quast_canu_unmapped --eukaryote ./canu_unmapped/final_unmapped_renamed.fa
echo "Canu unmapped assembly Quast finished"

### 5.6.2. Busco assessment # change options accordingly for each species 
echo "Canu mapped assembly Busco started"
busco -f -m genome -c ${N} -o Busco_canu_mapped -l aves_odb10 -i ./canu_mapped/final_mapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
echo "Canu mapped assembly Busco finished"

echo "Canu unmapped assembly Busco started"
busco -f -m genome -c ${N} -o Busco_canu_unmapped -l aves_odb10 -i ./canu_unmapped/final_unmapped_renamed.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10
echo "Canu unmapped assembly Busco finished"


# Combine contig sets to generate sample-derived supercontig collections - delete
#cat ./canu_mapped/final_mapped_renamed.fa ./canu_unmapped/final_unmapped_renamed.fa > sample_contigs.fa


# 6. Updating the original reference genome by appending supercontigs from unmapped reads
cat ${REF}/${GENOME}.fna ./canu_unmapped/final_unmapped_renamed.fa > ${LINPAN}/${PREFIX}_panref1.fa


# Second round of updating reference genome (using unmapped reads only)
# 7. Read mapping using Bowtie2
module purge
module load biocontainers
module load bowtie2
module load samtools
module load bbmap
module load r/4.3

bowtie2-build ${LINPAN}/${PREFIX}_panref1.fa ${LINPAN}/${PREFIX}_panref1

cd ../
mkdir -p mapped2/
cd mapped2/

bowtie2 -p ${N} -I 0 -X 1000 -x ${LINPAN}/${PREFIX}_panref1 -1 ../Nremoved/merged_Nremoved_1.fq -2 ../Nremoved/merged_Nremoved_2.fq --end-to-end --sensitive -S mapped2.sam 2> bowtie2.log

samtools sort -@ ${N} -o mapping2_sorted.bam mapped2.sam

## 7.1. Average insert size of the mapped paired-end reads for the MaSuRCA assembly
#reformat.sh in=mapped2.sam ihist=ihist2.txt
${APP}/jdk-21.0.1/bin/java -Duser.country=US -Duser.language=en -jar ${APP}/picard.jar CollectInsertSizeMetrics -I mapping2_sorted.bam -O ihist2.picard.txt -H ihist2.picard.pdf -M 0.05 #--DEVIATIONS # picard is better

### 7.1.2. Mapping rate, depth, and breadth
echo "mapping rate is" > ./mapping2_rate.txt
samtools flagstat -@ ${N} ./mapping2_sorted.bam >> ./mapping2_rate.txt 
echo "depth is" > ./mapping2_depth.txt
samtools depth -@ ${N} -a ./mapping2_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> ./mapping2_depth.txt 
echo "breadth is" > ./mapping2_breadth.txt
samtools depth -@ ${N} -a ./mapping2_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./mapping2_breadth.txt


# 8. Extraction of unmapped reads
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


# 9. Independent assemblies of unmapped reads
## 9.1. MaSuRCA mapping using unmapped reads
mkdir -p ./masurca_unmapped2

cd ./masurca_unmapped2
echo "Masurca assembly of unmapped reads started"
${MASURCA}/masurca masurca_config_unmapped2.txt # configure the configuration file beforehand
bash assemble.sh
echo "Masurca assembly of unmapped reads finished"
cd ../

## 9.2. MegaHit mapping using unmapped reads
echo "Megahit assembly of unmapped reads started"
megahit -t ${N} -1 unmapped2_R1.fastq -2 unmapped2_R2.fastq -r unmapped2_unpaired.fastq -o ./megahit_unmapped2 --out-prefix unmapped2_megahit 
echo "Megahit assembly of unmapped reads finished"

## 9.3. Small contig filtering
module purge
module load biocontainers
module load bbmap
#module load bioawk

for g in unmapped2; do
  cd masurca_${g}/
  cd CA/
  cat primary.genome.scf.fasta alternative.genome.scf.fasta > ../${g}_masurca.fasta # a limitation in that we cannot distinguish primary and alternate contigs (of an individual or among individuals) since we pooled all the reads among samples 
  cd ../
  sortbyname.sh in=${g}_masurca.fasta out=${g}_masurca_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 200) { print ">"$name; print $seq }}' ${g}_masurca_sorted.fasta  > ${g}_masurca_200.fa
  reformat.sh minlength=200 in=${g}_masurca_sorted.fasta out=${g}_masurca_200.fa
  cd ../
  
  cd megahit_${g}/
  cp ${g}_megahit.contigs.fa ${g}_megahit.fasta
  sortbyname.sh in=${g}_megahit.fasta out=${g}_megahit_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 200) { print ">"$name; print $seq }}' ${g}_megahit_sorted.fasta  > ${g}_megahit_200.fa
  reformat.sh minlength=200 in=${g}_megahit_sorted.fasta out=${g}_megahit_200.fa
  cd ../
done


# 10. Meta-assemblies of unmapped reads independently into supercontigs
module purge
module load biocontainers
module load canu
mkdir -p ./canu_unmapped2/

## 10.1. Canu assmebly using unmapped contigs
cd ./canu_unmapped2/
echo "Canu assembly of unmapped reads started"
canu -p final_unmapped2 -d ./ genomeSize=${GENOME_SIZE} -trimmed -minReadLength=200 -minOverlapLength=100 -minInputCoverage=0 -correctedErrorRate=0.105 -assemble -pacbio-hifi ${SRA}/mapped2/masurca_unmapped2/unmapped2_masurca_200.fa ${SRA}/mapped2/megahit_unmapped2/unmapped2_megahit_200.fa 
cat final_unmapped2.contigs.fasta final_unmapped2.unassembled.fasta > final_unmapped2.fa
echo "Canu assembly of unmapped reads finishted"
cd ../

## 10.2. Small contig filtering
module purge
module load biocontainers
module load bbmap
#module load bioawk

for g in unmapped2; do
  cd canu_${g}/
  sortbyname.sh in=final_${g}.fa out=final_${g}_sorted.fasta length descending overwrite=T
  #bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' final_${g}_sorted.fasta  > final_${g}_500.fa
  reformat.sh minlength=500 in=final_${g}_sorted.fasta  out=final_${g}_500.fa
  cd ../
done

## 10.3. Duplicate contig removal
module purge
module load biocontainers
module load bbmap

for g in unmapped2; do
  cd ./canu_${g}/
  dedupe.sh in=final_${g}_500.fa out=final_${g}_dedup.fa
  cd ../
done

## 10.4. Renaming contig headers to avoid redundant names from different fasta files
sed 's/>tig/>unmap2_tig/g' ./canu_unmapped2/final_unmapped2_dedup.fa > ./canu_unmapped2/final_unmapped2_renamed.fa


# 11. Confirming unmapped supercontigs are unmappable and removing mappable ones
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


# 12. Confirming mapped contigs are mappable and removing unmappable ones (technical artifacts - kimeric and alternative contigs, etc.)
minimap2 -ax asm5 ${REF}/${GENOME}.fna ./mapped/canu_mapped/final_mapped_renamed.fa > ${CONTIG}/maptig_mapped.sam
samtools sort -@ ${N} -o ${CONTIG}/maptig_mapped.bam ${CONTIG}/maptig_mapped.sam 
samtools flagstat -@ ${N} ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.mappingrate
samtools depth -@ ${N} -a ${CONTIG}/maptig_mapped.bam | awk '{c++;s+=$3}END{print s/c}' > ${CONTIG}/maptig.depth
samtools depth -@ ${N} -a ${CONTIG}/maptig_mapped.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${CONTIG}/maptig.breadth
samtools coverage ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.stat

samtools view -@ ${N} -b -F 4 ${CONTIG}/maptig_mapped.bam > ${CONTIG}/maptig.bam # confirmed mapped contigs
bamtools convert -in ${CONTIG}/maptig.bam -format fasta -out ${CONTIG}/maptig.fa


# 13. Pooling all the supercontigs (mapped or unmapped) to generate sample-derived contig collections
cd ${CONTIG}/

module purge
module load biocontainers
module load bbmap

## 13.1. Check duplicates
dedupe.sh in=unmaptig.fa out=unmaptig_dedup.fa # 0 reads were identified
dedupe.sh in=maptig.fa out=maptig_dedup.fa

#cat ./mapped/sample_contigs.fa ./canu_unmapped2/final_unmapped2_renamed.fa > sample_contigs2.fa - delete
#cat ./mapped/canu_mapped/final_mapped_renamed.fa ./unmaptig.fa > sample_contigs2.fa - delete
cat maptig_dedup.fa unmaptig_dedup.fa > sample_contigs2.fa

## 13.2. Duplicate contig removal
dedupe.sh in=sample_contigs2.fa out=sample_contigs2_dedup.fa


# 14. Updating the original reference genome by appending supercontigs from unmapped reads -> Get the linear pangenome by the "iterative map-then-assemble" approach
module purge
module load biocontainers
module load bbmap
module load samtools
#module load bioawk 

## 14.1. Filtering out small contigs (<10kb)
#bioawk -c fastx '{ if(length($seq) >= 10000) { print ">"$name; print $seq }}' unmaptig_dedup.fa > unmaptig_10kb.fa
reformat.sh minlength=10000 in=unmaptig_dedup.fa out=unmaptig_10kb.fa

cat ${REF}/${GENOME}.fna ./unmaptig_10kb.fa > ${PREFIX}_panref2.fa
 
## 14.2. Duplicate contig removal
dedupe.sh in=${PREFIX}_panref2.fa out=${PREFIX}_panref2_dedup.fa # 0% duplication

## 14.3. Sorting the updated reference
sortbyname.sh in=${PREFIX}_panref2_dedup.fa out=${PREFIX}_panref2_sorted.fa length descending overwrite=T

## 14.4. Indexing the updated reference
samtools faidx ${PREFIX}_panref2_sorted.fa  # check contig length distribution to decide the length cutoff below (>10kb can all the reference contigs + alpha) 

## 14.5. Assembly quality check
module load biocontainers
module load quast
module load busco
export AUGUSTUS_CONFIG_PATH=$RCAC_SCRATCH/augustus/config

quast.py --threads ${N} -o Quast_panref2 --eukaryote ./${PREFIX}_panref2_sorted.fa
busco -f -m genome -c ${N} -o Busco_panref2 -l --auto-lineage-euk -i ./${PREFIX}_panref2_sorted.fa --augustus --augustus_species chicken --limit 5 --update-data --datasets_version odb10


# 15. Retrieving each sample's contigs from the sample-derived supercontig pool
module purge
module load biocontainers
module load bowtie2
module load bioawk
module load samtools
module load minimap2

mkdir -p ./retrieved/
cd ./retrieved/

## 15.1. Mapping each sample's reads onto sample_contigs2.fa
bowtie2-build ../sample_contigs2_dedup.fa ../sample_contigs2

touch all_mappingrate.txt
touch all_depth.txt
touch all_coverage.txt
touch all_covstat.txt

for i in `ls -lt ${SRA}/raw | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do
  bowtie2 -p ${N} -I 0 -X 1000 -x ../sample_contigs2 -1 ${SRA}/Nremoved/${i}_*1.fq_Ns_removed -2 ${SRA}/Nremoved/${i}_*2.fq_Ns_removed --end-to-end --sensitive -S ${i}_mapped.sam 2> ${i}_bowtie.log
  samtools sort -@ ${N} -o ${i}_mapped_sorted.bam ${i}_mapped.sam
  samtools flagstat -@ ${N} ${i}_mapped_sorted.bam | grep "mapped (" | grep -v "primary" | cut -d "(" -f 2 | cut -d "%" -f 1 >> all_mappingrate.txt # avg.rate = ? vs. ? against the original reference
  samtools coverage ${i}_mapped_sorted.bam | sort -k3 -k6 -k7 -k9 -r -g > ${i}_cov_sorted.txt # sort by 6th (breadth) then 7th (depth) 
  cut -f1,3,6,7,9 ${i}_cov_sorted.txt >> all_covstat.txt # keep only length, breadth, depth, then mapq
  samtools depth -@ ${N} -a ${i}_mapped_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> all_depth.txt # avg.depth = ? vs. ? against the original reference
  samtools depth -@ ${N} -a ${i}_mapped_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> all_breadth.txt # avg.coverage = ? vs. ? against the original reference
done
sort -k2 -k3 -k4 -k5 -r -g all_covstat.txt > all_covstat_sorted.txt # sort by length, breadth, depth, then mapq

### Example when using simulated reads is as follows:
#for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
#  for j in primary alternate; do
#    for k in {1..4}; do
#     bowtie2 -p ${N} -I 0 -X 1000 -x ../sample_contigs2 -1 sim${k}_${i}_${j}_R1.fq.gz -2 sim${k}_${i}_${j}_R2.fq.gz --end-to-end --sensitive -S sim${k}_${i}_${j}_mapped.sam 2> sim${k}_${i}_${j}_bowtie.log
#     samtools sort -@ ${N} -o sim${k}_${i}_${j}_mapped_sorted.bam sim${k}_${i}_${j}_mapped.sam
#     samtools flagstat -@ ${N} sim${k}_${i}_${j}_mapped_sorted.bam | grep "mapped (" | grep -v "primary" | cut -d "(" -f 2 | cut -d "%" -f 1 >> all_mappingrate.txt # avg.rate = ? vs. ? against the original reference
#     samtools coverage sim${k}_${i}_${j}_mapped_sorted.bam | sort -k3 -k6 -k7 -k9 -r -g > sim${k}_${i}_${j}_cov_sorted.txt # sort by 6th (breadth) then 7th (depth) 
#     cut -f1,3,6,7,9 sim${k}_${i}_${j}_cov_sorted.txt >> all_covstat.txt # keep only length, breadth, depth, then mapq
#     samtools depth -@ ${N} -a sim${k}_${i}_${j}_mapped_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> all_depth.txt # avg.depth = ? vs. ? against the original reference
#     samtools depth -@ ${N} -a sim${k}_${i}_${j}_mapped_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> all_breadth.txt # avg.coverage = ? vs. ? against the original reference
#    done
#  done
#done
#sort -k2 -k3 -k4 -k5 -r -g all_covstat.txt > all_covstat_sorted.txt 

# In R, as total,
# quantile(breadth, probs=c(0.25,0.5,0.75,0.9,0.95)) 
# = 25%: 64.7833  50%: 79.6269  75%: 86.0457  90%: 90.6313  95%: 93.8471
# summary(breadth)
# = Min.: 0.000  Median: 79.63  Mean: 69.95  Max.: 100.00
# quantile(depth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 1.595200  50%: 2.429270  75%: 2.981620  90%: 3.720740  95%: 4.989538
# summary(depth)
# = Min.: 0.000  Median: 2.429  Mean: 3.859  Max.: 7688.670
# quantile(mapq, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 2.03  50%: 3.62  75%: 5.72  90%: 11.60  95%: 19.20
# summary(mapq)
# = Min.: 0.000  Median: 3.620  Mean: 5.386  Max.: 42.000

# among contigs >100kb 
# quantile(breadth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 82.04980  50%: 84.51305  75%: 86.65110  90%: 88.37187  95%: 89.31033
# summary(breadth)
# = Min.: 55.19  Median: 84.51  Mean: 83.98  Max.: 95.76
# quantile(depth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 2.559728  50%: 2.754490  75%: 2.953157  90%: 3.150981  95%: 3.277524
# summary(depth)
# = Min.: 1.067  Median: 2.754  Mean: 2.749  Max.: 10.306
# quantile(mapq, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 3.72  50%: 4.21  75%: 4.75  90%: 5.27  95%: 5.65
# summary(mapq) = Min.: 1.740  Median: 4.210  Mean: 4.288  Max.: 27.400

# among contigs >10kb  & <=100kb                
# quantile(breadth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 79.91628  50%: 83.73460  75%: 86.88900  90%: 89.37240  95%: 90.75180
# summary(breadth)
# = Min.: 30.30  Median: 83.73  Mean: 82.93  Max.: 100.00
# quantile(depth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 2.39775  50%: 2.69171  75%: 2.99618  90%: 3.31164  95%: 3.53992
# summary(depth)
# = Min.: 0.5069  Median: 2.6917  Mean: 2.7340  Max.: 439.3360
# quantile(mapq, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 2.71  50%: 3.62  75%: 4.62  90%: 5.76  95%: 6.64
# summary(mapq)
# = Min.: 0.906  Median: 3.620  Mean: 3.829  Max.: 41.600

# among contigs <=10kb                
# quantile(breadth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 54.7256  50%: 75.6683  75%: 85.2350  90%: 91.5292  95%: 95.2286
# summary(breadth)
# = Min.: 0.00  Median: 75.67  Mean: 65.40  Max.: 100.00
# quantile(depth, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 1.15265  50%: 2.19497  75%: 2.96989  90%: 4.04638  95%: 6.15076
# summary(depth)
# = Min.: 0.000  Median: 2.195  Mean: 4.253  Max.: 7688.670
# quantile(mapq, probs=c(0.25,0.5,0.75,0.9,0.95))
# = 25%: 1.65  50%: 3.58  75%: 6.44  90%: 14.60  95%: 22.70
# summary(mapq)
# = Min.: 0.000  Median: 3.580  Mean: 5.925  Max.: 42.000    

## 15.2. Saving supercontigs of coverage >80, depth >3, mapq >5 for contigs >100kb; coverage >85, depth >3, mapq >5 for contigs >10kb; coverage >90, depth >3, mapq >5 for contigs <=10kb; the size of the set of supercontig pool is a double of the original reference genome (so depth should be the half of the original 5x and mapq should be around 3 (50% accuracy) if assuming random distribution) 
for i in `ls -lt ${SRA}/raw | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do
  #total_lines_long=$(awk '{if ($3 > 100000) {print $0}}' ${i}_cov_sorted.txt | sort -n | wc -l)
  #ninety_index_long=$(echo "($total_lines_long * 90) / 100" | bc)
  #breadth_long=$(awk '{if ($3 > 100000) {print $0}}' ${i}_cov_sorted.txt | sort -n -k6 | awk -v idx="$ninety_index_long" 'NR == idx {print $6}')
  #depth_long=$(awk '{if ($3 > 100000) {print $0}}' ${i}_cov_sorted.txt | sort -n -k7 | awk -v idx="$ninety_index_long" 'NR == idx {print $7}')
  #mapq_long=$(awk '{if ($3 > 100000) {print $0}}' ${i}_cov_sorted.txt | sort -n -k9 | awk -v idx="$ninety_index_long" 'NR == idx {print $9}')
  awk '{if ($3 > 100000 && $6 > 80 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt > ${i}_contigs.stat
  awk '{if ($3 > 10000 && $3 <= 100000 && $6 > 85 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt >> ${i}_contigs.stat
  awk '{if ($3 <= 10000 && $6 > 90 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt >> ${i}_contigs.stat
  grep -v "#" ${i}_contigs.stat | awk 'BEGIN {FS="\t"}; {print $1}' > ${i}_contigs.region
  cat ${i}_contigs.region >> all_contigs.region
done
sort all_contigs.region | uniq > all_contigs_uniq.region

### Example when using simulated reads is as follows:
#for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
#  for j in primary alternate; do
#    for k in {1..4}; do
#     awk '{if ($3 > 100000 && $6 > 80 && $7 > 3 && $9 > 5) {print $0}}' sim${k}_${i}_${j}_cov_sorted.txt > sim${k}_${i}_${j}_contigs.stat
#     awk '{if ($3 > 10000 && $3 <= 100000 && $6 > 85 && $7 > 3 && $9 > 5) {print $0}}' sim${k}_${i}_${j}_cov_sorted.txt >> sim${k}_${i}_${j}_contigs.stat
#     awk '{if ($3 <= 10000 && $6 > 90 && $7 > 3 && $9 > 5) {print $0}}' sim${k}_${i}_${j}_cov_sorted.txt >> sim${k}_${i}_${j}_contigs.stat
#     grep -v "#" sim${k}_${i}_${j}_contigs.stat | awk 'BEGIN {FS="\t"}; {print $1}' > sim${k}_${i}_${j}_contigs.region
#     cat sim${k}_${i}_${j}_contigs.region >> all_contigs.region
#    done
#  done
#done
#sort all_contigs.region | uniq > all_contigs_uniq.region

samtools faidx ../sample_contigs2_dedup.fa > ../sample_contigs2_dedup.fa.fai
for i in `ls -lt ../raw | grep "fastq.gz" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do
  samtools faidx -r ${i}_contigs.region ../sample_contigs2_dedup.fa > ${i}_contigs.fa
done

## 15.3. Test-mapping the all_contigs onto the ${PREFIX}_panref2 in order to see mapping statistics (can help for forecasting a pangenome structure)
samtools faidx -r all_contigs_uniq.region ../sample_contigs2_dedup.fa > all_contigs.fa

minimap2 -ax asm5 ${LINPAN}/${PREFIX}_panref2_sorted.fa all_contigs.fa > sample_to_panref.sam

samtools sort -@ ${N} -o sample_to_panref.bam sample_to_panref.sam
samtools flagstat -@ ${N} sample_to_panref.bam > sample_to_panref_stat.txt
samtools depth -@ ${N} -a sample_to_panref.bam | awk '{c++;s+=$3}END{print s/c}' > sample_to_panref_depth.txt 
samtools depth -@ ${N} -a sample_to_panref.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > sample_to_panref_breadth.txt
samtools coverage sample_to_panref.bam > sample_to_panref_coverage.txt


# 16. Constructing a Minigraph-Cactus pangenome
cd ${BASE}

export MYBUCKET=${BASE}/mcpan
keep_contigs=$(awk '{print $1}' ${REF}/${GENOME}.fna.fai)

mkdir -p ./mcpan/
cd ./mcpan/

## 16.1. Creating a seqFile containing haplotype file information 
#touch ${PREFIX}Pangenome.txt

#echo -e "# Reference individual:" >> ${PREFIX}Pangenome.txt 
#echo -e "${RefInd}\t${LINPAN}/${PREFIX}_panref2_sorted.fa" >> ${PREFIX}Pangenome.txt

#echo -e "\n# Haploid sample:" >> ${PREFIX}Pangenome.txt
## add a alternate haplotype of the reference individual if there is (consider adding ".1" and ".2" for different haplotypes of the same reference diploid individual. Refer to the Minigraph-Cactus manual.)
#echo -e "bHirRus1_alt\t${REF}/GCA_015227815.3_bHirRus1.alt.v3_genomic.fna" >> ${PREFIX}Pangenome.txt

#for g in `ls -lt ../raw | grep "fastq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq`
#do
#  i=`echo ${g} | cut -d "_" -f 1`
#  echo -e "${i}\t${CONTIG}/retrieved/${i}_contigs.fa" >> ${PREFIX}Pangenome.txt
#done

### Example when using simulated reads is as follows:
#for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
# for j in primary alternate; do
#     for k in {1..4}; do
#     echo -e "sim${k}_${i}_${j}\t${CONTIG}/retrieved/sim${k}_${i}_${j}_contigs.fa" >> ${PREFIX}Pangenome.txt
#   done
# done
#done

# manually confirm the Pangenome.txt file! Check the manuaml at: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md

#module purge
#module load biocontainers
#module load cactus/2.6.5
#alias python='/opt/globus/bin/python3.9'
module load anaconda
conda activate cactus_env # pre-install cactus in conda environment "cactus_env"
source ~/app/venv-cactus-v2.7.0/bin/activate
export PATH=$PATH:/home/jeon96/app/cactus-bin-v2.7.0/bin/

## 16.2. Preprocessing seqFile
mkdir -p Temp/
cp ${PREFIX}Pangenome.txt ${PREFIX}Pangenome.mc.seqfile
cactus-preprocess ./jobstore ${PREFIX}Pangenome.txt ${PREFIX}Pangenome.mc.seqfile --pangenome --workDir ./Temp/ --binariesMode local
echo "preprocessing done"

## 16.3 Making a SV graph with Minigraph in GFA format 
cactus-minigraph ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}.gfa --reference ${RefInd} --maxCores ${N} --workDir Temp/ --binariesMode local
echo "minigraph done"

## 16.4. Mapping each input assembly back to the graph (assembly-to-graph alignments) using Minigraph
cactus-graphmap ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}.gfa ${PREFIX}.paf --outputGAFDir ${MYBUCKET}/${PREFIX}-mc-gaf --outputFasta ${PREFIX}.gfa.fa --reference ${RefInd} --nodeStorage 1000 --delFilter 10000000 --workDir Temp/ --maxMemory 1000G --mapCores 16 --maxCores ${N} --maxNodes 25 --binariesMode local
echo "graphmap done"

## 16.5. Splitting by chromosomes (or scaffolds)
cactus-graphmap-split ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}.gfa ${PREFIX}.paf --outDir ${MYBUCKET}/contigs --otherContig contigOther --refContigs $(for i in $keep_contigs; do echo ${i}; done) --reference ${RefInd} --nodeStorage 1000 --workDir Temp/ --maxMemory 1000G --maxCores ${N} --maxNodes 5 --logFile ${PREFIX}.split.log --binariesMode local
echo "split done"

## 16.6. Creating a Cactus base alignment and a "raw" pangenome graph
cactus-align --batch ./jobstore ./contigs/chromfile.txt ${MYBUCKET}/align --consCores ${N} --maxCores ${N} --maxMemory 1000G --workDir Temp/ --logFile ${PREFIX}.align.log --maxNodes 20 --nodeStorage 1000 --pangenome --maxLen 10000 --reference ${RefInd} --outVG --noLinkImports --noMoveExports --cleanWorkDir=onSuccess --realTimeLogging --binariesMode local
echo "align done"

## 16.7. Creating and indexing the final pangenome graph and produce a VCF ("--filter 4" filters sequences covered by less than 10% of these 40 haplotypes)
cactus-graphmap-join ./jobstore --vg ./align/*.vg --hal ./align/*.hal --outDir ./${PREFIX}-pg --outName ${PREFIX}-pg --reference ${RefInd} --vcf --gbz clip full --gfa --chrom-vg --giraffe clip --filter 4 --clip 10000 --nodeStorage 1000 --maxNodes 20 --workDir Temp/ --logFile ${PREFIX}.join.log --logDebug --indexCores $((N-1)) --maxCores ${N} --maxMemory 1500G --cleanWorkDir=onSuccess --disableCaching --chrom-og --viz --xg --binariesMode local --draw
echo "join" done

${APP}/vg stats -lz ${PREFIX}-pg.gbz # check pangenome basic statistics
#nodes   33092983
#edges   44756088
#length  1151889971

# compare to the length of RefInd in the graph: 1122798926 from the Quast output above

${APP}/vg stats -lz ${PREFIX}-pg.full.gbz # check how much additional sequence is added without clipping
#nodes   33165853
#edges   44837006
#length  1178700242

${APP}/vg stats -lz ${PREFIX}-pg.vg
#nodes	32946821
#edges	44609926
#length	1151889971

# 17. Augmenting the graph pangenome
cd ${MCPAN}/

## 17.1. Converting file type
gunzip ${PREFIX}-pg.gfa.gz 
${APP}/vg convert -g ${PREFIX}-pg.gfa -v > ${PREFIX}-pg.vg

## 17.2. Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./${PREFIX}-pg/
mkdir -p ./tmp/

${APP}/vg mod -t ${N} -O ../${PREFIX}-pg.vg > ${PREFIX}_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}_mod.xg ${PREFIX}_mod.vg

## 17.3. Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_mod.vg > ${PREFIX}_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}_mod_chopped.xg ${PREFIX}_mod_chopped.vg

## 17.4. Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_mod_chopped.vg > ${PREFIX}_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_mod_chopped_pruned.gcsa ${PREFIX}_mod_chopped_pruned.vg

## 17.5. Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}_mod_chopped.vg -x ${PREFIX}_mod_chopped_new_L.xg #-L: preserve alt paths in the xg

## 17.6 Aligning the individuals that were used to build the pangenome; do steps 17.6-17.10 when reads that were used to build the pangenome are different from reads that will be used to call variants  
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for k in {1..4}; do 
    cat ${SIM}/raw/sim${k}_${i}_primary_R1.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R1.fq.gz >  ${SIM}/raw/sim${k}_${i}_R1.fq.gz
    cat ${SIM}/raw/sim${k}_${i}_primary_R2.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R2.fq.gz >  ${SIM}/raw/sim${k}_${i}_R2.fq.gz
    ${APP}/vg map -t ${N} -f ${SIM}/raw/sim${k}_${i}_R1.fq.gz -f ${SIM}/raw/sim${k}_${i}_R2.fq.gz -x ${PREFIX}_mod_chopped.xg -g ${PREFIX}_mod_chopped_pruned.gcsa > sim${k}_${i}_aln.gam

## 17.7 Filtering secondary and ambiguous read mappings out of the GAM of the above step
    ${APP}/vg filter -t ${N} sim${k}_${i}_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${PREFIX}_mod_chopped.xg > sim${k}_${i}_aln.filtered.gam #-r : minimum score to keep primary alignment; -f: normalize score based on length; -u: use substitution count instead of score; -m: filter reads that don't begin with at least N matches on each end; -q: filter alignments with mapping quality < N; -D: clip back the ends of reads that are ambiguously aligned, up to N bases
  done
done

cat *filtered.gam > combined_filtered.gam 

## 17.8. Augmenting the graph with all variation from the GAM of the above step
${APP}/vg convert -t ${N} ${PREFIX}_mod_chopped_new_L.xg -p > ${PREFIX}.pg
${APP}/vg augment -t ${N} ${PREFIX}.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_aug.gam > ${PREFIX}_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 

## 17.9. Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_aug.pg > ${PREFIX}_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}_aug_chopped.xg ${PREFIX}_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_aug_chopped.pg > ${PREFIX}_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_aug_chopped_pruned.gcsa ${PREFIX}_aug_chopped_pruned.pg

## 17.10. Index the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}_aug_chopped.pg -x ${PREFIX}_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}_aug_chopped_new_L.xg -p > ${PREFIX}_aug_new.pg


# 18. Mapping test data reads - when not working with simulated reads, use the reference files based on ${PREFIX}-pg.vg, not ${PREFIX}_aug.pg
cd ${MCPAN}/${PREFIX}-pg/
mkdir -p ./aligned/tmp
cd ./aligned/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## 18.1. Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../${PREFIX}_aug_chopped.xg -g ../${PREFIX}_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## 18.2. Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now)
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam

## 18.3. Augmenting the graph with all variation from the GAM
${APP}/vg augment -t ${N} ../${PREFIX}_aug_new.pg combined_aug_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}_augSV.gam > ${PREFIX}_augSV.pg 

# 18.4. Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_augSV.pg > ${PREFIX}_augSV_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}_augSV_chopped.xg ${PREFIX}_augSV_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_augSV_chopped.pg > ${PREFIX}_augSV_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_augSV_chopped_pruned.gcsa ${PREFIX}_augSV_chopped_pruned.pg

## 18.5. Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}_augSV_chopped.pg -x ${PREFIX}_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}_augSV_chopped_new_L.xg -p > ${PREFIX}_augSV_new.pg


# 19. Variant calling (for SVs) for each individual using vg
cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/sv/vg/

## 19.1. Computing the snarls
${APP}/vg snarls -t ${N} ${PREFIX}_augSV_new.pg > ${PREFIX}_augSV.snarls

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## 19.2. Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}_augSV_chopped.xg -g ./${PREFIX}_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam

## 19.3. Filtering secondary and ambiguous read mappings out of the gam
  ${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam

## 19.4. Computing the support
  ${APP}/vg pack -t ${N} -x ${PREFIX}_augSV_new.pg -g ${i}_augSV_aln.filtered.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
done

## 19.5. Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg call -t ${N} ${PREFIX}_augSV_new.pg -r ${PREFIX}_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 --progress > ../called/sv/vg/${i}_mcSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called
done


# 20. Filtering SV vcf files from vg
cd ${MCPAN}/${PREFIX}-pg/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools
module load samtools
module load nf-core
module load bbmap
module load picard
module load bwa
module load boost

## 20.1. Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  sed 's/bHirRus1_LinPan#0#//g' ${i}_mcSV.vcf > ${i}_mcSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_mcSV2.vcf -Oz -o ${i}_mcSV.sorted.vcf.gz
  bcftools index ${i}_mcSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_mcSV2.vcf
done

## 20.2. Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}_mcSV.merged.vcf.gz --threads ${N} *_mcSV.sorted.vcf.gz 

## 20.2. Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSV.merged.vcf.gz > vcf-stats.txt

## 20.3. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.       1st Qu.    Median     Mean     3rd Qu.   Max.
# -26.684     9.542     21.751     35.119    43.435 23775.500
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.050    0.583    1.570    2.000 6173.960
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.000    0.120    0.320    0.583    1.000    1.545    2.500    4.320    6173.960

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0000  0.0800   0.3031  0.6000  1.0000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4012   1.3013   1.4679   1.4516   1.6741   1.9699

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.2819  0.2897  0.2963  0.3031  0.3061  0.4271

#quit()

## 20.4. Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf ${PREFIX}_mcSV.merged.vcf.gz --out ${PREFIX}_mcSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## 20.5. Annotating SVs in vcf file
bcftools sort ${PREFIX}_mcSV.filtered.recode.vcf -Oz -o ${PREFIX}_mcSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}_mcSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}_mcSV.filtered.popped.vcf # remove large (>100 kb) alleles in MC graphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
PREFIX=short-read # needs to define again inside Apptainer
vcfwave -I 1000 -t ${N} ${PREFIX}_mcSV.filtered.popped.vcf > ${PREFIX}_mcSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
truvari anno svinfo -m 50 -o ${PREFIX}_mcSV.decomposed.annotated.vcf ${PREFIX}_mcSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python ${APP}/SVextractor_fromTruvariAnno.py 50 ${PREFIX}_mcSV.decomposed.annotated.vcf ${PREFIX}_mcSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)


# 21. Variant calling (for SNPs using Sarek and for SVs using delly2) based on surjected bamfiles for each individual 
cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/snp

module purge
module load biocontainers
module load picard
module load bwa
module load boost

## 21.1. Creating a list of reference paths
${APP}/vg paths -x ${MCPAN}/${PREFIX}-pg/${PREFIX}-pg.full.gbz -S ${RefInd} -L > ../${RefInd}.${PREFIX}-pg_paths.txt # use .full graph just to make path lists (used updated version of vg for this - v.1.53)

## 21.2. Projecting each sample's reads to the RefInd
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../${PREFIX}_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam
  ${APP}/vg surject -x ../${PREFIX}_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -F ../${RefInd}.${PREFIX}-pg_paths.txt -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}-pg_surject.bam
done

## 21.3 Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### 21.3.1. Indexing bam files
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}-pg_surject.sorted.bam ${i}.${PREFIX}-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}-pg_surject.sorted.bam   
done

### 21.3.2. Preprocessing input files for sarek
#cd $CLUSTER_SCRATCH
#mkdir -p nf-core
#mkdir -p /scratch/negishi/jeon96/singularity/cache
#cd nf-core  
#nf-core download sarek # download sarek pipeline

cd ${MCPAN}/${PREFIX}-pg/aligned/

cat ${BASE}/original/SRR_Acc_List.txt > sample; \
#printf '1\n%.0s' {1..25} > lane ;  \
ls -1 ${MCPAN}/${PREFIX}-pg/aligned/*surject.sorted.bam > bam; \
ls -1 ${MCPAN}/${PREFIX}-pg/aligned/*surject.sorted.bam.bai > bai; \
echo "patient,sample,lane,bam,bai" > header; \
paste sample sample lane bam bai | tr "\t" "," > content; \
cat header content > mapped.csv; \
rm header content sample lane bam bai # generate mapped.csv as an input to sarek

sed 's/>/>bHirRus1_LinPan#0#/g' ${LINPAN}/${PREFIX}_panref2_sorted.fa > ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa # accommodate "reference contig" names of the surjected bam files
samtools faidx ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa
cat ${PREFIX}_panref2_sorted.renamed.fa.fai | cut -f1 | grep "unmap_tig" | sort > unmap_tig.sorted
cat ${PREFIX}_panref2_sorted.renamed.fa.fai | cut -f1 | grep -v "unmap_tig" | sort > ref_tig.sorted
cat unmap_tig.sorted ref_tig.sorted > bamcontig.list
sortbyname.sh in=${PREFIX}_panref2_sorted.renamed.fa out=${PREFIX}_panref2_resorted.renamed.fa list=bamcontig.list overwrite=T # match the order of contigs of the linpan reference with the order in the bam files for sarek
samtools faidx ${PREFIX}_panref2_resorted.renamed.fa 

# generate nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#    "input": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/aligned\/mapped.csv",
#    "step": "markduplicates",
#    "outdir": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/called\/snp\/",
#    "tools": "haplotypecaller",
#    "skip_tools": "baserecalibrator",
#    "trim_fastq": true,
#    "genome": "custom",
#    "fasta": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa",
#    "fasta_fai": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa.fai",
#    "igenomes_ignore": true,
#    "email": "jeon96@purdue.edu",
#    "multiqc_title": "short-read-pan_swallow"
#}

### 21.3.3. Running sarek pipeline
nextflow run nf-core/sarek -r 3.3.2 -profile singularity -work-dir ${MCPAN}/${PREFIX}-pg/aligned/ -resume -params-file ${MCPAN}/${PREFIX}-pg/aligned/nf-params.json 
#--step markduplicates --input ${MCPAN}/${PREFIX}-pg/aligned/mapped.csv --outdir ${MCPAN}/${PREFIX}-pg/called/snp --tools haplotypecaller --skip_tools baserecalibrator  

## 21.4. Calling SVs using delly2
### 21.4.1. Preprocessing the reference genome and bam files
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.renamed.fa --OUTPUT ${PREFIX}_panref2_sorted.renamed.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.renamed.fa

cd ${MCPAN}/${PREFIX}-pg/aligned/
mkdir -p ../called/sv/delly

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}-pg_surject.sorted.marked.bam
done

### 21.4.2. Running delly2
cd ../called/sv/delly

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${MCPAN}/${PREFIX}-pg/aligned/${i}.${PREFIX}-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa ${samples} > mcpan_delly.vcf 
bcftools view -Ob -o mcpan_delly.bcf mcpan_delly.vcf
bcftools index mcpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o mcpan_delly.filtered.bcf mcpan_delly.bcf
bcftools view -Ov -o mcpan_delly.filtered.vcf mcpan_delly.filtered.bcf


# 22. Filtering vcf files
cd ${MCPAN}/${PREFIX}-pg/called/snp

module purge
module load biocontainers
module load bcftools
module load vcftools
module load picard

## 22.1. Filtering SNP vcf files
### 22.1.1. Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_mcSNP.sorted.vcf.gz
  bcftools index ${i}_mcSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_mcSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### 22.1.2. Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o ${PREFIX}_mcSNP.merged.vcf.gz --threads ${N} *_mcSNP.sorted.vcf.gz 

### 22.1.3. Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}_mcSNP.merged.vcf.gz > vcf-stats.txt

### 22.1.4. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     73.64     121.80     157.47    194.96 148296.00
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.418    6.500 2741.100
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.143    3.000    4.000    4.278    5.000    5.500    6.000    7.000    8.000    2741.100

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8434  0.9600  0.9600
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.962   5.308   5.681   5.702   6.238   6.938

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8071  0.8277  0.8378  0.8434  0.8550  0.9499

#quit()

### 22.1.4. Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf ${PREFIX}_mcSNP.merged.vcf.gz --out ${PREFIX}_mcSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' ${PREFIX}_mcSNP.filtered.recode.vcf > ${PREFIX}_mcSNP.filtered.renamed.vcf

## 22.2. Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o mcpan_delly.filtered.vcf.gz mcpan_delly.filtered.vcf

### 22.2.3. Filtering with population-level parameters
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf mcpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats mcpan_delly.filtered.vcf.gz > vcf-stats.txt

### 22.2.4. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2000)
#summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     840.0     868.4    1200.0 10000.0
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02161  0.04000  0.24000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.004124  0.006186  0.009072  0.021608  0.012371  0.277938

#quit()

vcftools --gzvcf mcpan_delly.filtered.vcf.gz --out mcpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/bHirRus1_LinPan#0#//g' mcpan_delly.filtered2.recode.vcf > mcpan_delly.filtered2.renamed.vcf
bcftools sort mcpan_delly.filtered2.renamed.vcf -Ov -o mcpan_delly.filtered2.sorted.vcf    


# 23. Constructing a VG pangenome
module load biocontainers
module load bwa
module load picard
module load boost
module load bcftools
module load htslib

cd ${BASE}
mkdir -p ${VGPAN}/augmented/called
cd ${VGPAN}/augmented

## 23.1. Running sarek with simulated samples to augment the linear pangenome
### 23.1.1. Preparing an input sheet
ls -lrt /scratch/negishi/jeon96/swallow/sim/raw/ | tr -s ' ' | cut -d ' ' -f 9 | grep "primary" | cut -c 1-24 | uniq | sort > primary; \
ls -lrt /scratch/negishi/jeon96/swallow/sim/raw/ | tr -s ' ' | cut -d ' ' -f 9 | grep "alternate" | cut -c 1-26 | uniq | sort > alternate; \
cat primary alternate > sample; \
printf '1\n%.0s' {1..40} > lane ;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*primary_R1.fq.gz | sort > primary_fastq_1;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*alternate_R1.fq.gz | sort > alternate_fastq_1;  \
cat primary_fastq_1 alternate_fastq_1 > fastq_1; \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*primary_R2.fq.gz | sort > primary_fastq_2;  \
ls -1 /scratch/negishi/jeon96/swallow/sim/raw/*alternate_R2.fq.gz | sort > alternate_fastq_2;  \
cat primary_fastq_2 alternate_fastq_2 > fastq_2; \
echo "patient,sample,lane,fastq_1,fastq_2" > header; \
paste sample sample lane fastq_1 fastq_2 | tr "\t" "," > content; \
cat header content > samplesheet.csv; \
rm header content sample lane fastq_1 fastq_2 primary alternate primary_fastq_1 primary_fastq_2 alternate_fastq_1 alternate_fastq_2 # generate samplesheet.csv as an input to sarek

### 23.1.2. generating nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#    "input": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/aligned\/mapped.csv",
#    "step": "markduplicates",
#    "outdir": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/called\/snp\/",
#    "tools": "haplotypecaller",
#    "skip_tools": "baserecalibrator",
#    "trim_fastq": true,
#    "genome": "custom",
#    "fasta": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa",
#    "fasta_fai": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa.fai",
#    "igenomes_ignore": true,
#    "email": "jeon96@purdue.edu",
#    "multiqc_title": "short-read-pan_swallow"
#}

## 23.2. Calling SVs using delly2
### 23.2.1. Mapping for delly2
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.fa --OUTPUT ${PREFIX}_panref2_sorted.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.fa

cd ${VGPAN}/augmented
mkdir -p ./mapped/tmp/
cd ./mapped/

for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${LINPAN}/${PREFIX}_panref2_sorted.fa  ${SIM}/raw/sim${k}_${i}_${j}_R1.fq.gz ${SIM}/raw/sim${k}_${i}_${j}_R2.fq.gz > /sim${k}_${i}_${j}.sam
    PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT sim${k}_${i}_${j}.sam --OUTPUT sim${k}_${i}_${j}_sorted.bam --SORT_ORDER coordinate
    PicardCommandLine MarkDuplicates --INPUT sim${k}_${i}_${j}_sorted.bam --OUTPUT ./sim${k}_${i}_${j}_sorted.marked.bam --METRICS_FILE sim${k}_${i}_${j}_metrics.txt
    PicardCommandLine BuildBamIndex --INPUT ./sim${k}_${i}_${j}_sorted.marked.bam
    done 
  done
done

### 23.2.2. Running delly2
cd ${VGPAN}/augmented/called/  
mkdir -p ./delly/
cd ./delly/

export OMP_NUM_THREADS=40 # equal or less than the number of samples
samples=""
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    samples+="${VGPAN}/augmented/mapped/sim${k}_${i}_${j}_sorted.marked.bam "
    done
  done
done

${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > sim_delly.vcf 
bcftools view -Ob -o sim_delly.bcf sim_delly.vcf
bcftools index sim_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o sim_delly.filtered.bcf sim_delly.bcf
bcftools view -Ov -o sim_delly.filtered.vcf sim_delly.filtered.bcf 

## 23.3. Filtering variants (SNPs)
cd ../

### 23.3.1. Compressing and indexing each vcf file first
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    bcftools sort ./variant_calling/haplotypecaller/sim${k}_${i}_${j}/sim${k}_${i}_${j}.haplotypecaller.vcf.gz -Oz -o sim${k}_${i}_${j}_SNP.sorted.vcf.gz
    bcftools index sim${k}_${i}_${j}_SNP.sorted.vcf.gz --threads ${N}
    bcftools view sim${k}_${i}_${j}_SNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
    done
  done
done

### 23.3.2. Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o sim_SNP.merged.vcf.gz --threads ${N} *_SNP.sorted.vcf.gz 

### 23.3.3. Filtering with population-level parameters
vcftools --gzvcf sim_SNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf sim_SNP.merged.vcf.gz --missing-site
vcftools --gzvcf sim_SNP.merged.vcf.gz --depth
vcftools --gzvcf sim_SNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf sim_SNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sim_SNP.merged.vcf.gz > vcf-stats.txt

### 23.3.4. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.    1st Qu.   Median    Mean     3rd Qu.  Max.
# 30.0     301.0     426.1     449.6    570.1 81304.1
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    5.583    7.375    7.971    9.250 6713.750
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    4.250    5.250    6.000    6.750    7.375    8.000    8.750    9.875    11.750    6713.750

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.7000  0.9000   0.7685  0.9000  0.9750
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 5.297   5.318   5.985   7.657   9.841   11.966

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.7248  0.7286  0.7526  0.7685  0.7984  0.8468

#quit()

### 23.3.5. Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf sim_SNP.merged.vcf.gz --out sim_SNP.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

### 23.4. Filtering variants (SVs)
cd ./delly/
bcftools convert --thread ${N} -Oz -o sim_delly.filtered.vcf.gz sim_delly.filtered.vcf

### 23.4.1. Filtering with population-level parameters
vcftools --gzvcf sim_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf sim_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf sim_delly.filtered.vcf.gz --depth
vcftools --gzvcf sim_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf sim_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats sim_delly.filtered.vcf.gz > vcf-stats.txt

### 23.4.2. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2500)
#summary(var_qual$qual)
#  Min.  1st Qu.  Median   Mean   3rd Qu. Max.
# 300     647     1020     1052    1200 10000
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu. Median    Mean   3rd Qu.  Max.
# 0.0000  0.0250  0.1000   0.1026  0.2000  0.2500
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#   Min.   1st Qu.  Median     Mean    3rd Qu.   Max.
# 0.01619  0.03147  0.05195  0.10262  0.13422  0.32350

#quit()

### 23.4.3. Filtering vcf file based on the determined cut-off values 
vcftools --gzvcf sim_delly.filtered.vcf.gz --out ../sim_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## 23.2. Constructing a VG pangenome
cd ../

### 23.2.1. Preprocessing input files
bcftools convert --thread ${N} -Ov -o sim_delly.filtered2.vcf sim_delly.filtered2.recode.vcf.gz
while read -r line; do
    pair=($line)
    echo "s/${pair[0]}/${pair[1]}/" >> delly_replace.sed
done < delly_sample.list # make a sed script from the list of "search-replace" pairs
sed -f delly_replace.sed sim_delly.filtered2.vcf > sim_delly.filtered2.replaced.vcf # fix sample names

while read -r line; do
    pair=($line)
    echo "s/${pair[0]}/${pair[1]}/" >> snp_replace.sed
done < snp_sample.list # make a sed script from the list of "search-replace" pairs
sed -f snp_replace.sed sim_SNP.filtered.recode.vcf > sim_SNP.filtered.replaced.vcf # fix sample names

bcftools query -l sim_delly.filtered2.replaced.vcf > samples.txt # make a sample list
bcftools view -S samples.txt sim_delly.filtered2.replaced.vcf > sim_delly.filtered2.ordered.vcf # order samples; check the sample header! bcftools view sim_delly.filtered2.ordered.vcf | grep -v "##" | head
cat sim_delly.filtered2.ordered.vcf | grep -v "SVTYPE=BND" | grep -v "SVTYPE=DUP" > sim_delly.filtered3.ordered.vcf # remove unsupported SV types in vg construct (i.e., DUP, BND)
bcftools view -S samples.txt sim_SNP.filtered.replaced.vcf > sim_SNP.filtered.ordered.vcf # order samples; check the sample header! bcftools view sim_SNP.filtered.ordered.vcf | grep -v "##" | head

bcftools convert --thread ${N} -Oz -o sim_SNP.filtered.vcf.gz sim_SNP.filtered.ordered.vcf
bcftools convert --thread ${N} -Oz -o sim_delly.filtered3.vcf.gz sim_delly.filtered3.ordered.vcf
tabix -p vcf sim_SNP.filtered.vcf.gz
tabix -p vcf sim_delly.filtered3.vcf.gz
bcftools concat --allow-overlaps -Oz -o sim_SNP_SV.vcf.gz --threads ${N} sim_SNP.filtered.vcf.gz sim_delly.filtered3.vcf.gz # combine SNP and SV vcf files

bcftools norm -a --atom-overlaps . -m -any --threads ${N} -Oz -o sim_SNP_SV.normed.vcf.gz sim_SNP_SV.vcf.gz # normalize the vcf file      
tabix -p vcf sim_SNP_SV.normed.vcf.gz

### 23.2.2. Constructing a graph pangenome using vg and the linear pangenome
${APP}/vg construct -a -f -S -r ${LINPAN}/${PREFIX}_panref2_sorted.fa -v sim_SNP_SV.normed.vcf.gz -m 1000000 > ../${PREFIX}2-pg.vg #| ${APP}/vg convert - | ${APP}/vg mod -x 32 -  # too big inversions exist for Protobuf to stream -> moving node chopping outside of vg construct

${APP}/vg stats -lz ${PREFIX}2-pg.vg # check pangenome basic statistics
#nodes   61690
#edges   88747
#length  1123193383


# 24. Augmenting the graph pangenome
cd ${VGPAN}/

## 24.1. Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./augmented/
mkdir -p ./tmp/

${APP}/vg mod -t ${N} -O ${PREFIX}2-pg.vg > ${PREFIX}2_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}2_mod.xg ${PREFIX}2_mod.vg

## 24.2. Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_mod.vg > ${PREFIX}2_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}2_mod_chopped.xg ${PREFIX}2_mod_chopped.vg

## 24.3. Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_mod_chopped.vg > ${PREFIX}2_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_mod_chopped_pruned.gcsa ${PREFIX}2_mod_chopped_pruned.vg

## 24.4. Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}2_mod_chopped.vg -x ${PREFIX}2_mod_chopped_new_L.xg #-L: preserve alt paths in the xg

## 24.5 Aligning the individuals that were used to build the pangenome; do steps 17.6-17.10 when reads that were used to build the pangenome are different from reads that will be used to call variants  
for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for k in {1..4}; do 
    #cat ${SIM}/raw/sim${k}_${i}_primary_R1.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R1.fq.gz >  ${SIM}/raw/sim${k}_${i}_R1.fq.gz
    #cat ${SIM}/raw/sim${k}_${i}_primary_R2.fq.gz ${SIM}/raw/sim${k}_${i}_alternate_R2.fq.gz >  ${SIM}/raw/sim${k}_${i}_R2.fq.gz
    ${APP}/vg map -t ${N} -f ${SIM}/raw/sim${k}_${i}_R1.fq.gz -f ${SIM}/raw/sim${k}_${i}_R2.fq.gz -x ${PREFIX}2_mod_chopped.xg -g ${PREFIX}2_mod_chopped_pruned.gcsa > sim${k}_${i}_aln.gam

## 24.6 Filtering secondary and ambiguous read mappings out of the GAM of the above step
    ${APP}/vg filter -t ${N} sim${k}_${i}_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${PREFIX}2_mod_chopped.xg > sim${k}_${i}_aln.filtered.gam #-r : minimum score to keep primary alignment; -f: normalize score based on length; -u: use substitution count instead of score; -m: filter reads that don't begin with at least N matches on each end; -q: filter alignments with mapping quality < N; -D: clip back the ends of reads that are ambiguously aligned, up to N bases
  done
done

cat *filtered.gam > combined_filtered.gam 

## 24.7. Augmenting the graph with all variation from the GAM of the above step
${APP}/vg convert -t ${N} ${PREFIX}2_mod_chopped_new_L.xg -p > ${PREFIX}2.pg
${APP}/vg augment -t ${N} ${PREFIX}2.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_aug.gam > ${PREFIX}2_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 

## 24.8. Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_aug.pg > ${PREFIX}2_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}2_aug_chopped.xg ${PREFIX}2_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_aug_chopped.pg > ${PREFIX}2_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_aug_chopped_pruned.gcsa ${PREFIX}2_aug_chopped_pruned.pg

## 24.9. Index the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}2_aug_chopped.pg -x ${PREFIX}2_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}2_aug_chopped_new_L.xg -p > ${PREFIX}2_aug_new.pg


# 25. Mapping test data reads - when not working with simulated reads, use the reference files based on ${PREFIX}2-pg.vg, not ${PREFIX}2_aug.pg
cd ${VGPAN}/
mkdir -p ./aligned/tmp
cd ./aligned/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## 25.1. Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ../augmented/${PREFIX}2_aug_chopped.xg -g ../augmented/${PREFIX}2_aug_chopped_pruned.gcsa > ${i}_aug_aln.gam

## 25.2. Filtering secondary and ambiguous read mappings out of the gam (for SV detection for now)
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_aln.filtered.gam

done

cat *_aug_aln.filtered.gam > combined_aug_aln.filtered.gam

## 25.3. Augmenting the graph with all variation from the GAM
${APP}/vg augment -t ${N} ../augmented/${PREFIX}2_aug_new.pg combined_aug_aln.filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}2_augSV.gam > ${PREFIX}2_augSV.pg 

# 25.4. Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}2_augSV.pg > ${PREFIX}2_augSV_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}2_augSV_chopped.xg ${PREFIX}2_augSV_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}2_augSV_chopped.pg > ${PREFIX}2_augSV_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}2_augSV_chopped_pruned.gcsa ${PREFIX}2_augSV_chopped_pruned.pg

## 25.5. Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}2_augSV_chopped.pg -x ${PREFIX}2_augSV_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}2_augSV_chopped_new_L.xg -p > ${PREFIX}2_augSV_new.pg


# 26. Variant calling (for SVs) for each individual using vg
cd ${VGPAN}/aligned/
mkdir -p ../called/sv/vg/

## 26.1. Computing the snarls
${APP}/vg snarls -t ${N} ${PREFIX}2_augSV_new.pg > ${PREFIX}2_augSV.snarls

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do

## 26.2. Aligning the individuals
  ${APP}/vg map -t ${N} -f ${CLEANED_SRA}/${i}_1_val_1.fq -f ${CLEANED_SRA}/${i}_2_val_2.fq -x ./${PREFIX}2_augSV_chopped.xg -g ./${PREFIX}2_augSV_chopped_pruned.gcsa > ${i}_augSV_aln.gam

## 26.3. Filtering secondary and ambiguous read mappings out of the gam
  ${APP}/vg filter -t ${N} ${i}_augSV_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ./${PREFIX}2_augSV_chopped.xg > ${i}_augSV_aln.filtered.gam

## 26.4. Computing the support
  ${APP}/vg pack -t ${N} -x ${PREFIX}2_augSV_new.pg -g ${i}_augSV_aln.filtered.gam -Q 5 -o ${i}_augSV.pack #-Q 5: ignore mapping and base qualitiy < 5
done

## 26.5. Calling variants; run this step using highmem queue (otherwise it can't be finished)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg call -t ${N} ${PREFIX}2_augSV_new.pg -r ${PREFIX}2_augSV.snarls -k ${i}_augSV.pack -s ${i} -a -A -c 50 -C 100000 > ../called/sv/vg/${i}_vgSV.vcf #-a: calling every snarl using the same coordinates and including reference calls; -c: minimum length of variants to be called; -C: maximum length of variants to be called
done


# 27. Filtering SV vcf files from vg
cd ${VGPAN}/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools
module load samtools
module load nf-core
module load bbmap
module load picard
module load bwa
module load boost

## 27.1. Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ${i}_vgSV.vcf -Oz -o ${i}_vgSV.sorted.vcf.gz
  bcftools index ${i}_vgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_vgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

## 27.2. Combining separately called SV vcf files
bcftools merge -m all -Oz -o ${PREFIX}2_vgSV.merged.vcf.gz --threads ${N} *_vgSV.sorted.vcf.gz 

## 27.2. Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}2_vgSV.merged.vcf.gz > vcf-stats.txt

## 27.3. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.       1st Qu.    Median     Mean     3rd Qu.   Max.
# -19.510     9.542     17.835     31.086    36.919 6701.480
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.000    0.600    1.662    2.040 2739.560
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.000    0.040    0.250    0.600    1.000    1.680    2.800    4.877    2739.560

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0400  0.2800   0.4412  0.9600  1.0000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4893   1.4052   1.6702   1.6089   1.8430   2.2308

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.3727  0.4326  0.4374  0.4412  0.4463  0.5403

#quit()

## 27.4. Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf ${PREFIX}2_vgSV.merged.vcf.gz --out ${PREFIX}2_vgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## 27.5. Annotating SVs in vcf file
bcftools sort ${PREFIX}2_vgSV.filtered.recode.vcf -Oz -o ${PREFIX}2_vgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input ${PREFIX}2_vgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > ${PREFIX}2_vgSV.filtered.popped.vcf # remove large (>100 kb) alleles in MC raphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
PREFIX=short-read # needs to define again inside Apptainer
vcfwave -I 1000 -t ${N} ${PREFIX}2_vgSV.filtered.popped.vcf > ${PREFIX}2_vgSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
truvari anno svinfo -m 50 -o ${PREFIX}2_vgSV.decomposed.annotated.vcf ${PREFIX}2_vgSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python ${APP}/SVextractor_fromTruvariAnno.py 50 ${PREFIX}2_vgSV.decomposed.annotated.vcf ${PREFIX}2_vgSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)


# 28. Variant calling (for SNPs using Sarek and for SVs using delly2) based on surjected bamfiles for each individual 
cd ${VGPAN}/aligned/
mkdir -p ../called/snp
#${APP}/vg gbwt -p -E -x ../augmented/${PREFIX}2-pg.vg --gbz-format -g ../augmented/${PREFIX}2-pg.gbz # make a gbz file 

module load biocontainers
module load nf-core

## 28.1. Projecting each sample's reads to the reference paths
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  ${APP}/vg filter -t ${N} ${i}_aug_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -i -x ../augmented/${PREFIX}2_aug_chopped.xg > ${i}_aug_interleaved_aln.filtered.gam
  ${APP}/vg surject -x ../augmented/${PREFIX}2_aug_chopped.xg ${i}_aug_interleaved_aln.filtered.gam --threads ${N} --prune-low-cplx --interleaved -b -N ${i} -R "ID:1 LB:lib1 SM:${i} PL:illumina PU:unit1" > ${i}.${PREFIX}2-pg_surject.bam
done

## 28.2 Calling SNPs (and small SVs < 50 bp) using nf-core Sarek pipeline
### 28.2.1. Indexing bam files
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samtools sort -@ ${N} -o ${i}.${PREFIX}2-pg_surject.sorted.bam ${i}.${PREFIX}2-pg_surject.bam 
  samtools index -@ ${N} -b ${i}.${PREFIX}2-pg_surject.sorted.bam   
done

### 28.2.2. Preprocessing input files for sarek
cat ${BASE}/original/SRR_Acc_List.txt > sample; \
#printf '1\n%.0s' {1..25} > lane ;  \
ls -1 ${VGPAN}/aligned/*surject.sorted.bam > bam; \
ls -1 ${VGPAN}/aligned/*surject.sorted.bam.bai > bai; \
echo "patient,sample,bam,bai" > header; \
paste sample sample bam bai | tr "\t" "," > content; \
cat header content > mapped.csv; \
rm header content sample bam bai # generate mapped.csv as an input to sarek

# generate nf-params.json file as input to sarek at: https://oldsite.nf-co.re/launch
#e.g., 
#{
#    "input": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/aligned\/mapped.csv",
#    "step": "markduplicates",
#    "outdir": "\/scratch\/negishi\/jeon96\/swallow\/mcpan\/short-read-pg\/called\/snp\/",
#    "tools": "haplotypecaller",
#    "skip_tools": "baserecalibrator",
#    "trim_fastq": true,
#    "genome": "custom",
#    "fasta": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa",
#    "fasta_fai": "\/scratch\/negishi\/jeon96\/swallow\/linpan\/${PREFIX}_panref2_resorted.renamed.fa.fai",
#    "igenomes_ignore": true,
#    "email": "jeon96@purdue.edu",
#    "multiqc_title": "short-read-pan_swallow"
#}

### 28.2.3. Running sarek pipeline
nextflow run nf-core/sarek -r 3.3.2 -profile singularity -work-dir ${VGPAN}/aligned/ -resume -params-file ${VGPAN}/aligned/nf-params.json 
#--step markduplicates --input ${VGPAN}/aligned/mapped.csv --outdir ${VGPAN}/called/snp --tools haplotypecaller --skip_tools baserecalibrator  

## 28.3. Calling SVs using delly2
### 28.3.1. Preprocessing the reference genome and bam files
cd ${VGPAN}

cd ${VGPAN}/aligned/
mkdir -p ../called/sv/delly

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  PicardCommandLine MarkDuplicates --INPUT ${i}.${PREFIX}2-pg_surject.sorted.bam --OUTPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam --METRICS_FILE ${i}.${PREFIX}2-pg_surject_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}.${PREFIX}2-pg_surject.sorted.marked.bam
done

### 28.3.2. Running delly2
cd ../called/sv/delly

module load boost

export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${VGPAN}/aligned/${i}.${PREFIX}2-pg_surject.sorted.marked.bam "
done
   
${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > vgpan_delly.vcf 
bcftools view -Ob -o vgpan_delly.bcf vgpan_delly.vcf
bcftools index vgpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o vgpan_delly.filtered.bcf vgpan_delly.bcf
bcftools view -Ov -o vgpan_delly.filtered.vcf vgpan_delly.filtered.bcf


# 29. Filtering vcf files
cd ${VGPAN}/called/snp/

module purge
module load biocontainers
module load bcftools
module load vcftools

## 29.1. Filtering SNP vcf files
### 29.1.1. Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/haplotypecaller/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_vgSNP.sorted.vcf.gz
  bcftools index ${i}_vgSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_vgSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### 29.1.2. Combining aseparately called SNP vcf files
bcftools merge -m all -Oz -o ${PREFIX}2_vgSNP.merged.vcf.gz --threads ${N} *_vgSNP.sorted.vcf.gz 

### 29.1.3. Filtering with population-level parameters
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --missing-site
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --depth
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ${PREFIX}2_vgSNP.merged.vcf.gz > vcf-stats.txt

### 29.1.4. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     75.64     121.84     159.29    194.96 138011.00
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.452    6.500 2599.000
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.250    3.000    4.000    4.333    5.000    5.500    6.000    7.000    8.000    2599.000

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8418  0.9600  0.9600
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.996   5.357   5.743   5.765   6.306   7.005

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8059  0.8259  0.8364  0.8418  0.8537  0.9488

#quit()

### 29.1.4. Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf ${PREFIX}2_vgSNP.merged.vcf.gz --out ${PREFIX}2_vgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## 29.2. Filtering SV vcf files from delly2
cd ../sv/delly/
bcftools convert --thread ${N} -Oz -o vgpan_delly.filtered.vcf.gz vgpan_delly.filtered.vcf

### 29.2.3. Filtering with population-level parameters
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf vgpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats vgpan_delly.filtered.vcf.gz > vcf-stats.txt

### 29.2.4. In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2000)
#summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     780.0     809.2    1200.0 10000.0
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02266  0.04000  0.24000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.004234  0.006543  0.009623  0.022664  0.015397  0.283680

#quit()

vcftools --gzvcf vgpan_delly.filtered.vcf.gz --out vgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all


# Calling and filtering variants from NCBI linear representative genome
# black's codes

### Compressing and indexing each vcf file first (By Andrew Black)
for i in `ls -1 *vcf.gz`; do   bcftools index $i -c -t ; done
bcftools merge *vcf.gz > ncbi_raw_snp.vcf

### Changing name and sorting  (By JongYoon from here)
cd ${REF}/called/

bcftools sort -Oz -o ncbi_SNP.sorted.vcf.gz ncbi_raw_snp.vcf

## Calling SVs using delly2
module load biocontainers
module load picard
module load bwa
module load boost

### Preprocessing the linear pangenome
cd ${REF}/
PicardCommandLine CreateSequenceDictionary --REFERENCE ${GENOME}.fna --OUTPUT ${GENOME}.dict
bwa index -a bwtsw ${GENOME}.fna

### Mapping for delly2
mkdir -p ./mapped/tmp/
cd ./mapped/
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
    bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${REF}/${GENOME}.fna  ${CLEANED_SRA}/${i}_1_val_1.fq ${CLEANED_SRA}/${i}_2_val_2.fq > ${i}.sam
    PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT ${i}.sam --OUTPUT ${i}_sorted.bam --SORT_ORDER coordinate
    PicardCommandLine MarkDuplicates --INPUT ${i}_sorted.bam --OUTPUT ${i}_sorted.marked.bam --METRICS_FILE ${i}_metrics.txt
    PicardCommandLine BuildBamIndex --INPUT ${i}_sorted.marked.bam
done

### Running delly2
cd ../
mkdir -p ./called/
cd ./called/
export OMP_NUM_THREADS=25 #equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
    samples+="${REF}/mapped/${i}_sorted.marked.bam "
done

${DEPOT}/apps/delly/src/delly call -g ${REF}/${GENOME}.fna ${samples} > ncbi_delly.vcf 
bcftools view -Ob -o ncbi_delly.bcf ncbi_delly.vcf
bcftools index ncbi_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o ncbi_delly.filtered.bcf ncbi_delly.bcf
bcftools view -Ov -o ncbi_delly.filtered.vcf ncbi_delly.filtered.bcf

## Filtering variants (SNPs)
module purge
module load biocontainers
module load bcftools
module load vcftools

### Filtering with population-level parameters
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --missing-indv 
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --missing-site
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --depth
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --site-mean-depth
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ncbi_SNP.sorted.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     78.28     129.64      166.88    210.02 209853.00
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    4.000    5.368    6.082    7.000 5288.500
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.500    3.286    4.000    5.000    5.368   6.000    7.000    7.875    9.000    5288.500

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8425  0.9600  0.9600
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 3.180   5.885   6.360   6.398   7.079   7.834

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8096  0.8273  0.8369  0.8425  0.8530  0.9452

#quit()

### Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf ncbi_SNP.sorted.vcf.gz --out ncbi_SNP.sorted.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Filtering variants (SVs)
bcftools convert --thread ${N} -Oz -o ncbi_delly.filtered.vcf.gz ncbi_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --depth
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats ncbi_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2500)
#summary(var_qual$qual)
#  Min.    1st Qu.  Median      Mean   3rd Qu. Max.
# 300.0     360.0     540.0     740.2    1015.5 10000.0
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu. Median    Mean   3rd Qu.  Max.
# 0.0000  0.0000  0.0000   0.0164  0.0400  0.2400
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#    Min.    1st Qu.   Median     Mean     3rd Qu.    Max.
# 0.001764  0.003919  0.006663  0.016390  0.009994  0.234764

#quit()

### Filtering vcf file based on the determined cut-off values 
vcftools --gzvcf ncbi_delly.filtered.vcf.gz --out ncbi_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all


# Calling and filtering variants from linear pangenome
# black's codes

### Compressing and indexing each vcf file first (By Andrew Black)
for i in `ls -1 *vcf.gz`; do   bcftools index $i -c -t ; done
bcftools merge *vcf.gz > linpan_raw_snp.vcf

### Changing names and sorting (By JongYoon from here)
cd ${LINPAN}/called/

bcftools sort -Oz -o linpan_SNP.sorted.vcf.gz linpan_raw_snp.vcf

## Calling SVs using delly2
module load biocontainers
module load picard
module load bwa
module load boost

### Preprocessing the linear pangenome
cd ${LINPAN}
PicardCommandLine CreateSequenceDictionary --REFERENCE ${PREFIX}_panref2_sorted.renamed.fa --OUTPUT ${PREFIX}_panref2_sorted.renamed.dict
bwa index -a bwtsw ${PREFIX}_panref2_sorted.renamed.fa

### Mapping for delly2
mkdir -p ./mapped/tmp/
cd ./mapped/
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bwa mem -t ${N} -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ${LINPAN}/${PREFIX}_panref2_sorted.fa  ${CLEANED_SRA}/${i}_1_val_1.fq ${CLEANED_SRA}/${i}_2_val_2.fq > ${i}.sam
  PicardCommandLine SortSam --TMP_DIR ./tmp/ --INPUT ${i}.sam --OUTPUT ${i}_sorted.bam --SORT_ORDER coordinate
  PicardCommandLine MarkDuplicates --INPUT ${i}_sorted.bam --OUTPUT ${i}_sorted.marked.bam --METRICS_FILE ${i}_metrics.txt
  PicardCommandLine BuildBamIndex --INPUT ${i}_sorted.marked.bam
done

### Running delly2
cd ../
mkdir -p ./called/
cd ./called/
export OMP_NUM_THREADS=25 # equal or less than the number of samples
samples=""
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samples+="${LINPAN}/mapped/${i}_sorted.marked.bam "
done

${DEPOT}/apps/delly/src/delly call -g ${LINPAN}/${PREFIX}_panref2_sorted.fa ${samples} > linpan_delly.vcf 
bcftools view -Ob -o linpan_delly.bcf linpan_delly.vcf
bcftools index linpan_delly.bcf
${DEPOT}/apps/delly/src/delly filter -f germline -o linpan_delly.filtered.bcf linpan_delly.bcf
bcftools view -Ov -o linpan_delly.filtered.vcf linpan_delly.filtered.bcf

## Filtering variants (SNPs)
module purge
module load biocontainers
module load bcftools
module load vcftools

### Filtering with population-level parameters
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --missing-indv 
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --missing-site
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --depth
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --site-mean-depth
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats linpan_SNP.sorted.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     78.28     129.64      166.79    209.64 205949.00
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    4.000    5.333    6.076    7.000 5272.000
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.500    3.250    4.000    5.000    5.333   6.000    7.000    7.857    9.000    5272.000

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.8400  0.9200   0.8427  0.9600  0.9600
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 3.178   5.878   6.361   6.394   7.065   7.833

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8098  0.8274  0.8370  0.8427  0.8531  0.9453

#quit()

### Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf linpan_SNP.sorted.vcf.gz --out linpan_SNP.sorted.filtered --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Filtering variants (SVs)
bcftools convert --thread ${N}-Oz -o linpan_delly.filtered.vcf.gz linpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf linpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf linpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf linpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf linpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf linpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats linpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2500)
#summary(var_qual$qual)
# Min.   1st Qu. Median   Mean 3rd Qu. Max.
# 300     360     540     739    1009 10000
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.    1st Qu.  Median     Mean     3rd Qu.  Max.
# 0.00000  0.00000  0.00000   0.01652  0.04000  0.24000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#    Min.    1st Qu.   Median     Mean     3rd Qu.    Max.
# 0.002148  0.003905  0.006833  0.016525  0.009957  0.236236

#quit()

### Filtering vcf file based on the determined cut-off values 
vcftools --gzvcf linpan_delly.filtered.vcf.gz --out linpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all


# Calling and filtering variants from original pangenome
# Natalie's codes

${APP}/vg stats -lz ${PREFIX}-pg.gbz # check pangenome basic statistics
#nodes   91856637
#edges   125994300
#length   1240496512

# compare to the length of RefInd in the graph: 1122798926 from the Quast output above

${APP}/vg stats -lz ${PREFIX}-pg.full.gbz # check how much additional sequence is added without clipping
#nodes   92163595
#edges   126377766
#length   1360787612

${APP}/vg stats -lz ${PREFIX}-pg.vg
#nodes   91732409
#edges   125870072
#length   1240496512

# Filtering variants from original pangenome
# Natalie's codes

# Filtering vcf files
cd ${ORGPAN}/called/snp/

module purge
module load biocontainers
module load bcftools
module load vcftools

## Filtering SNP vcf files
### Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  bcftools sort ./variant_calling/${i}/${i}.haplotypecaller.vcf.gz -Oz -o ${i}_orgSNP.sorted.vcf.gz
  bcftools index ${i}_orgSNP.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSNP.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
done

### Combining separately called SNP vcf files
bcftools merge -m all -Oz -o barnswallow_orgSNP.merged.vcf.gz --threads ${N} *_orgSNP.sorted.vcf.gz 

### Filtering with population-level parameters
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSNP.merged.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean     3rd Qu.   Max.
# 30.00     76.84     121.84     157.44    194.96 169562.00
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1.000    3.500    5.000    5.383    6.500 3814.730
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 1.000    2.333    3.000    4.000    4.333    5.000    5.500    6.000    7.000    8.000    3814.730

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min. 1st Qu.  Median  Mean 3rd Qu.  Max.
# 0.000  0.840  0.920   0.843  0.960  0.960
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max.
# 2.882   5.122   5.530   5.532   6.003   6.757

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.8071  0.8272  0.8374  0.8430  0.8550  0.9494

#quit()

### Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf barnswallow_orgSNP.merged.vcf.gz --out barnswallow_orgSNP.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' barnswallow_orgSNP.filtered.recode.vcf > barnswallow_orgSNP.filtered.renamed.vcf # polish the header line to be compatible with bcftools

## Filtering SV vcf files from delly2
cd ../delly/
bcftools convert --thread ${N} -Oz -o orgpan_delly.filtered.vcf.gz orgpan_delly.filtered.vcf

### Filtering with population-level parameters
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-indv 
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --missing-site
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-mean-depth
vcftools --gzvcf orgpan_delly.filtered.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats orgpan_delly.filtered.vcf.gz > vcf-stats.txt

### In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,2000)
#summary(var_qual$qual)
#  Min.     1st Qu.   Median    Mean     3rd Qu.  Max.
# 300.0     420.0     780.0     804.1    1200.0 5160.0
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
# Min.  1st Qu.  Median    Mean 3rd Qu.  Max.
# NA      NA       NA       NA    NA     NA

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.     1st Qu.  Median    Mean     3rd Qu.   Max.
# 0.00000  0.00000  0.00000   0.02191  0.04000  0.24000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
# Min. 1st Qu. Median  Mean  3rd Qu.  Max.
# NA    NA      NA      NA    NA      NA

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#   Min.     1st Qu.   Median     Mean     3rd Qu.   Max.
# 0.002957  0.005915  0.009612  0.021915  0.014418  0.280591

#quit()

vcftools --gzvcf orgpan_delly.filtered.vcf.gz --out orgpan_delly.filtered2 --minQ 30 --maf 0.05 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all
sed 's/refp#0#//g' orgpan_delly.filtered2.recode.vcf > orgpan_delly.filtered2.renamed.vcf # polish the header line to be compatible with bcftools

## Filtering SV vcf files from vg 
cd cd ${ORGPAN}/called/sv/vg/

module purge
module load biocontainers
module load bcftools
module load vcftools
module load samtools
module load nf-core
module load bbmap
module load picard
module load bwa
module load boost

## Compressing and indexing each vcf file first
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  sed 's/refp#0#//g' ${i}_orgSV.vcf > ${i}_orgSV2.vcf # polish the header line to be compatible with bcftools
  bcftools sort ${i}_orgSV2.vcf -Oz -o ${i}_orgSV.sorted.vcf.gz
  bcftools index ${i}_orgSV.sorted.vcf.gz --threads ${N}
  bcftools view ${i}_orgSV.sorted.vcf.gz --threads ${N} | grep -v "##" | wc -l
  rm ${i}_orgSV2.vcf
done

## Combining separately called SV vcf files (By Natalie Allen until here)
bcftools merge -m all -Oz -o barnswallow_orgSV.merged.vcf.gz --threads ${N} *_orgSV.sorted.vcf.gz 

## Filtering with population-level parameters (By JongYoon from here)
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-indv 
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --missing-site
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-mean-depth
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --site-quality
cat out.lmiss | awk '!/CHR/' | awk '$6 > 0.1' | cut -f1,2 >> badloci
bcftools stats barnswallow_orgSV.merged.vcf.gz > vcf-stats.txt

## In R, plotting and summarizing vcf stats
module load r

#R
#library(tidyverse)
#library(ggplot2)
#var_qual <- read_delim("./out.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
#a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#a + theme_light() + xlim(0,500)
#summary(var_qual$qual)
#  Min.     1st Qu.    Median     Mean   3rd Qu.   Max.
# -21.65     12.71     26.29     42.38    46.74 31367.90
   
#var_depth <- read_delim("./out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
#a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()   
#a + theme_light() + xlim(0,100)        
#summary(var_depth$mean_depth)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.160    0.875    1.787    2.312 3435.120
#quantile(var_depth$mean_depth, probs = seq(0, 1, 1/10))
#   0%      10%      20%      30%      40%      50%      60%      70%      80%      90%      100%  
# 0.000    0.000    0.080    0.267    0.520    0.875    1.250    1.909    2.882    4.750    3435.120

#var_miss <- read_delim("./out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(var_miss$fmiss)
#  Min.   1st Qu.  Median   Mean  3rd Qu.   Max.
# 0.0000  0.0000  0.0800   0.2787  0.5200  1.0000
                       
#ind_depth <- read_delim("./out.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
#a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_depth$depth)
#  Min.    1st Qu.  Median    Mean    3rd Qu.    Max.
# 0.4626   1.4524   1.6978   1.6670   1.8723   2.2785

#ind_miss  <- read_delim("./out.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light()
#summary(ind_miss$fmiss)
#  Min.   1st Qu.  Median  Mean   3rd Qu.   Max.
# 0.2621  0.2685  0.2724  0.2787  0.2824  0.3828

#quit()

## Filtering vcf file based on the determined cut-off values                                                                      
vcftools --gzvcf barnswallow_orgSV.merged.vcf.gz --out barnswallow_orgSV.filtered --remove-filtered lowad --minQ 30 --maf 0.05 --min-meanDP 5 --max-meanDP 15 --minDP 3 --max-missing 0.8 --exclude-positions badloci --recode --recode-INFO-all

## Annotating SVs in vcf file
bcftools sort barnswallow_orgSV.filtered.recode.vcf -Oz -o barnswallow_orgSV.filtered.recode.vcf.gz
${APP}/vcfbub --input barnswallow_orgSV.filtered.recode.vcf.gz --max-allele-length 100000 --max-level 0 > barnswallow_orgSV.filtered.popped.vcf # remove large (>100 kb) alleles in MC raphs. This removes all non-top-level bubbles from the VCF unless they were nested inside a top-level bubble with a reference length exceeding 100 kb; that is, top-level bubbles longer than that are replaced by their child nodes in the snarl tree

singularity run ${APP}/vcflib-1.09.simg
N=128
vcfwave -I 1000 -t ${N} barnswallow_orgSV.filtered.popped.vcf > barnswallow_orgSV.decomposed.vcf # decompose complex alleles into primitive ones
exit

module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
truvari anno svinfo -m 50 -o barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomposed.vcf # annotate SV info with the minimum allele size of 50
python ${APP}/SVextractor_fromTruvariAnno.py 50 barnswallow_orgSV.decomposed.annotated.vcf barnswallow_orgSV.decomp.annot.ext.vcf # extract alleles that have truvari annotation as SVs (minimum size of 50)


# Comparing assembly sizes, mapping metrics, number of variants
module purge  
module load biocontainers
module load bioawk
module load samtools
module load bcftools

## NCBI linear representative genome
cd ${REF}/

### Assembly size
bioawk -c fastx '{print $name length($seq)}' ${GENOME}.fna | awk '{sum += $1} END {print sum}' > ncbi_asm.length # 1106045383 bp

### Mapping metrics
cd ./mapped/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
samtools flagstat -@ ${N} ${i}.sam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
samtools view -@ ${N} ${i}.sam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 98.8812
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 50.3805

### Number of variants
cd ../called/

bcftools stats ncbi_SNP.sorted.filtered.recode.vcf > SNP.count # 202707
bcftools stats ncbi_delly.filtered2.recode.vcf > SV.count # 1272

## Linear pangenome
cd ${LINPAN}/

### Assembly size
bioawk -c fastx '{print $name length($seq)}' ${PREFIX}_panref2_sorted.fa | awk '{sum += $1} END {print sum}' > linpan_asm.length # 1122830712 bp

### Mapping metrics
cd ./mapped/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
samtools flagstat -@ ${N} ${i}.sam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
samtools view -@ ${N} ${i}.sam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 99.1236
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 50.5236

### Number of variants
cd ../called/

bcftools stats linpan_SNP.sorted.filtered.recode.vcf > SNP.count # 202969
bcftools stats linpan_delly.filtered2.recode.vcf > SV.count # 1292

## Minigraph-Cactus pangenome
cd ${MCPAN}/${PREFIX}-pg/aligned/

### Assembly size
${APP}/vg stats -p ${N} -l ${PREFIX}-pg.vg > mcpan_asm.length # 1151889971 bp

### Mapping metrics
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
samtools flagstat -@ ${N} ${i}.${PREFIX}-pg_surject.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
samtools view -@ ${N} ${i}.${PREFIX}-pg_surject.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 99.68
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 58.9984

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${APP}/vg stats --threads ${N} --alignments ${i}_augSV_aln.gam ${PREFIX}_augSV_chopped.xg > ${i}_mapping.vgstats
total_alignment=`cat ${i}_mapping.vgstats | grep "Total alignments:" | cut -d ":" -f 2`
total_aligned=`cat ${i}_mapping.vgstats | grep "Total aligned:" | cut -d ":" -f 2`
printf %.3f "$((10**3 * 100 * ${total_aligned}/${total_alignment}))e-3" >> vgmapping.rate
echo "" >> vgmapping.rate
cat ${i}_mapping.vgstats | grep "Mapping quality:" | cut -d ":" -f 2 | cut -d "," -f 1 | cut -d " " -f 3 >> vgmap.q
done

awk '{sum += $1} END {print sum/NR}' vgmapping.rate > vgmapping_rate.average 
awk '{sum += $1} END {print sum/NR}' vgmap.q > vgmapq.average

### Number of variants
cd ../called/

bcftools stats ./snp/${PREFIX}_mcSNP.filtered.recode.vcf > SNP.count # 134744
bcftools stats ./sv/vg/${PREFIX}_mcSV.decomp.annot.ext.vcf > vgSV.count # 13
bcftools stats ./sv/delly/mcpan_delly.filtered2.recode.vcf > dellySV.count # 2028

## VG pangenome
cd ${VGPAN}/augmented/

### Assembly size
${APP}/vg stats -p ${N} -l ${PREFIX}2-pg.vg > vgpan_asm.length # 1123193383 bp

### Mapping metrics
cd ${VGPAN}/aligned/

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
samtools flagstat -@ ${N} ${i}.${PREFIX}2-pg_surject.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
samtools view -@ ${N} ${i}.${PREFIX}2-pg_surject.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 100!
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 59.1715

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${APP}/vg stats --threads ${N} --alignments ${i}_augSV_aln.gam ${PREFIX}2_augSV_chopped.xg > ${i}_mapping.vgstats
total_alignment=`cat ${i}_mapping.vgstats | grep "Total alignments:" | cut -d ":" -f 2`
total_aligned=`cat ${i}_mapping.vgstats | grep "Total aligned:" | cut -d ":" -f 2`
printf %.3f "$((10**3 * 100 * ${total_aligned}/${total_alignment}))e-3" >> vgmapping.rate
echo "" >> vgmapping.rate
cat ${i}_mapping.vgstats | grep "Mapping quality:" | cut -d ":" -f 2 | cut -d "," -f 1 | cut -d " " -f 3 >> vgmap.q
done

awk '{sum += $1} END {print sum/NR}' vgmapping.rate > vgmapping_rate.average # 93.79 
awk '{sum += $1} END {print sum/NR}' vgmap.q > vgmapq.average # 50.47 

### Number of variants
cd ../called/

bcftools stats ./snp/${PREFIX}2_vgSNP.filtered.recode.vcf > SNP.count # 142958
bcftools stats ./sv/vg/${PREFIX}2_vgSV.decomp.annot.ext.vcf > vgSV.count # 6
bcftools stats ./sv/delly/vgpan_delly.filtered2.recode.vcf > dellySV.count # 2149

## Original pangenome
cd ${ORGPAN}/

### Assembly size
${APP}/vg stats -p ${N} -l barnswallow-pg.vg > orgpan_asm.length # 1240496512 bp

### Mapping metrics
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
samtools flagstat -@ ${N} ${i}.barnswallow-pg_surjected.bam | awk -F "[(|%]" 'NR == 7 {print "'$i' Mean mapping rate = ", $2}' > ${i}_mapping.metrics
samtools view -@ ${N} ${i}.barnswallow-pg_surjected.bam | awk '{sum+=$5} END {print "'$i' Mean MAPQ = ",sum/NR}' >> ${i}_mapping.metrics
cat ${i}_mapping.metrics | grep "mapping rate" | cut -d "=" -f 2 >> mapping.rate
cat ${i}_mapping.metrics | grep "MAPQ" | cut -d "=" -f 2 >> map.q
done

awk '{sum += $1} END {print sum/NR}' mapping.rate > mapping_rate.average # 97.54
awk '{sum += $1} END {print sum/NR}' map.q > mapq.average # 57.71

### Number of variants
cd ./called/

bcftools stats ./snp/barnswallow_orgSNP.filtered.recode.vcf > SNP.count # 133785
bcftools stats ./sv/vg/barnswallow_orgSV.decomp.annot.ext.vcf > vgSV.count # 15
bcftools stats ./sv/delly/orgpan_delly.filtered2.recode.vcf > dellySV.count # 2247


# Variant benchmarking (per HPRC paper)
module purge
module load biocontainers
module load bcftools

cd ${BASE}
mkdir -p ./benchmarking/
cd ./benchmarking/

## Linear representative genome vs. Linear pangenome
mkdir -p ./ncbi_linpan/
cd ./ncbi_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
bcftools reheader --samples ${BASE}/original/SRR_Acc_List.txt ${REF}/called/ncbi_SNP.sorted.filtered.recode.vcf > ncbi_SNP_final.vcf # polish sample names in the header
cp ${REF}/called/ncbi_delly.filtered2.recode.vcf ./ncbi_delly_final.vcf

bcftools reheader --samples ${BASE}/original/SRR_Acc_List.txt ${LINPAN}/called/linpan_SNP.sorted.filtered.recode.vcf > linpan_SNP_final.vcf # polish sample names in the header
cp ${LINPAN}/called/linpan_delly.filtered2.recode.vcf ./linpan_delly_final.vcf

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools view -a -I -s ${i} --thread ${N} ncbi_SNP_final.vcf > ncbi_SNP_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} ncbi_delly_final.vcf > ncbi_delly_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} linpan_SNP_final.vcf > linpan_SNP_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} linpan_delly_final.vcf > linpan_delly_final_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools norm -m -any --thread ${N} -Oz -o ncbi_SNP_final_norm_${i}.vcf.gz ncbi_SNP_final_${i}.vcf > 
bcftools norm -m -any --thread ${N} -Oz -o ncbi_delly_final_norm_${i}.vcf.gz ncbi_delly_final_${i}.vcf
bcftools norm -m -any --thread ${N} -Oz -o linpan_SNP_final_norm_${i}.vcf.gz linpan_SNP_final_${i}.vcf
bcftools norm -m -any --thread ${N} -Oz -o linpan_delly_final_norm_${i}.vcf.gz linpan_delly_final_${i}.vcf
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i ncbi_SNP_final_norm_${i}.vcf.gz -o ncbi_SNP_final_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i linpan_SNP_final_norm_${i}.vcf.gz -o linpan_SNP_final_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

${DEPOT}/apps/rtg-tools-3.12.1/rtg format -o linpan.sdf ${LINPAN}/${PREFIX}_panref2_sorted.fa # format the reference genome the variants are called against to RTG's Sequence Data File format 

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools index -t --threads ${N} ncbi_delly_final_norm_${i}.vcf.gz
bcftools index -t --threads ${N} linpan_delly_final_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b ncbi_delly_final_norm_${i}.vcf.gz -c linpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Linear representative genome vs. Minigraph-Cactus pangenome
mkdir -p ./ncbi_mcpan/
cd ./ncbi_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

cp ${MCPAN}/${PREFIX}-pg/called/snp/${PREFIX}_mcSNP.filtered.renamed.vcf ./mcpan_SNP_final.vcf 
cp ${MCPAN}/${PREFIX}-pg/called/sv/delly/mcpan_delly.filtered2.sorted.vcf ./mcpan_delly_final.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
module load picard 

picard UpdateVcfSequenceDictionary I=./mcpan_SNP_final.vcf O=./mcpan_SNP_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./mcpan_delly_final.vcf O=./mcpan_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_SNP_final_updt.vcf O=./mcpan_SNP_final_sorted.vcf
picard SortVcf I=./mcpan_delly_final_updt.vcf O=./mcpan_delly_final_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools view -a -I -s ${i} --thread ${N} mcpan_SNP_final_sorted.vcf > mcpan_SNP_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} mcpan_delly_final_sorted.vcf > mcpan_delly_final_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools norm -m -any --thread ${N} -Oz -o mcpan_SNP_final_norm_${i}.vcf.gz mcpan_SNP_final_${i}.vcf
bcftools norm -m -any --thread ${N} -Oz -o mcpan_delly_final_norm_${i}.vcf.gz mcpan_delly_final_${i}.vcf
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i mcpan_SNP_final_norm_${i}.vcf.gz -o mcpan_SNP_final_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz.tbi ./
bcftools index -t --threads ${N} mcpan_delly_final_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b ncbi_delly_final_norm_${i}.vcf.gz -c mcpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Linear representative genome vs. VG pangenome
mkdir -p ./ncbi_vgpan/
cd ./ncbi_vgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/
mkdir -p ./sv_comp/vg/

cp ${VGPAN}/called/snp/${PREFIX}2_vgSNP.filtered.recode.vcf ./vgpan_SNP_final.vcf 
cp ${VGPAN}/called/sv/delly/vgpan_delly.filtered2.recode.vcf ./vgpan_delly_final.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
module load biocontainers
module load picard 
 
picard UpdateVcfSequenceDictionary I=./vgpan_SNP_final.vcf O=./vgpan_SNP_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./vgpan_delly_final.vcf O=./vgpan_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_SNP_final_updt.vcf O=./vgpan_SNP_final_sorted.vcf
picard SortVcf I=./vgpan_delly_final_updt.vcf O=./vgpan_delly_final_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools view -a -I -s ${i} --thread ${N} vgpan_SNP_final_sorted.vcf > vgpan_SNP_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} vgpan_delly_final_sorted.vcf > vgpan_delly_final_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools norm -m -any --thread ${N} -Oz -o vgpan_SNP_final_norm_${i}.vcf.gz vgpan_SNP_final_${i}.vcf
bcftools norm -m -any --thread ${N} -Oz -o vgpan_delly_final_norm_${i}.vcf.gz vgpan_delly_final_${i}.vcf 
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i vgpan_SNP_final_norm_${i}.vcf.gz -o vgpan_SNP_final_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz.tbi ./
bcftools index -t --threads ${N} vgpan_delly_final_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b ncbi_delly_final_norm_${i}.vcf.gz -c vgpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Original pangenome vs. Minigraph-Cactus pangenome
mkdir -p ./orgpan_mcpan/
cd ./orgpan_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

cp ${MCPAN}/${PREFIX}-pg/called/sv/vg/${PREFIX}_mcSV.decomp.annot.ext.vcf ./mcpan_vg_final.vcf
cp ${ORGPAN}/called/snp/barnswallow_orgSNP.filtered.renamed.vcf ./orgpan_SNP_final.vcf 
cp ${ORGPAN}/called/sv/delly/orgpan_delly.filtered2.renamed.vcf ./orgpan_delly_final.vcf
cp ${ORGPAN}/called/sv/delly/barnswallow_orgSV.decomp.annot.ext.vcf  ./orgpan_vg_final.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
module load biocontainers
module load bbmap
module load samtools
module load picard 

#python
#def extract_contigs_from_vcf(vcf_file):
#    contigs_dict = {}
#    with open(vcf_file, 'r') as f:
#        for line in f:
#            if line.startswith('##contig='):
#                parts = line.strip().split('<')[1].split('>')[0]
#                contig_id = None
#                contig_length = None
#                for part in parts.split(','):
#                    if part.startswith('ID='):
#                        contig_id = part.split('=')[1]
#                    elif part.startswith('length='):
#                        contig_length = int(part.split('=')[1])
#                if contig_id and contig_length:
#                    contigs_dict[contig_id] = contig_length           
#    return contigs_dict

#vcf_file="orgpan_vg_final.vcf"
#contigs = extract_contigs_from_vcf(vcf_file)
#print("Contigs:", contigs)

#with open("orgSV_chroms.txt", 'w') as f:
#    for key, value in contigs.items():
#        f.write(str(key) + " " + str(value) + '\n') 

sort -k2,2 -n orgSV_chroms.txt > sorted_orgSV_chroms.txt
cut -f1,2 ${REF}/${GENOME}.fna.fai > ref_chroms.txt
sort -k2,2 -n ref_chroms.txt > sorted_ref_chroms.txt
join -1 2 -2 2 sorted_orgSV_chroms.txt sorted_ref_chroms.txt > combined_chroms.txt # check if there are duplicates and if so manually curate them
cut -d " " -f2,3 combined_chroms.txt > chr_name_conv.txt

#sed 's/|/\//g' orgpan_vg_final.vcf > orgpan_vg_final_sbst.vcf # vg vcf files include "|" after truvari
#sed 's/|/\//g' mcpan_vg_final.vcf > mcpan_vg_final_sbst.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt mcpan_vg_final.vcf -Ov -o mcpan_vg_final_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_vg_final.vcf -Ov -o orgpan_vg_final_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_delly_final.vcf -Ov -o orgpan_delly_final_annot.vcf
bcftools annotate --threads ${N} --rename-chrs chr_name_conv.txt orgpan_SNP_final.vcf -Ov -o orgpan_SNP_final_annot.vcf

picard UpdateVcfSequenceDictionary I=./mcpan_vg_final_annot.vcf O=./mcpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_SNP_final_annot.vcf O=./orgpan_SNP_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_delly_final_annot.vcf O=./orgpan_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard UpdateVcfSequenceDictionary I=./orgpan_vg_final_annot.vcf O=./orgpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_vg_final_updt.vcf O=./mcpan_vg_final_sorted.vcf
picard SortVcf I=./orgpan_SNP_final_updt.vcf O=./orgpan_SNP_final_sorted.vcf
picard SortVcf I=./orgpan_delly_final_updt.vcf O=./orgpan_delly_final_sorted.vcf
picard SortVcf I=./orgpan_vg_final_updt.vcf O=./orgpan_vg_final_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' mcpan_vg_final_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > mcpan_vg_final_${i}.vcf # sed step to bypass bcftools error
bcftools view -a -I -s ${i} --thread ${N} orgpan_SNP_final_sorted.vcf > orgpan_SNP_final_${i}.vcf
bcftools view -a -I -s ${i} --thread ${N} orgpan_delly_final_sorted.vcf > orgpan_delly_final_${i}.vcf
sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' orgpan_vg_final_sorted.vcf |bcftools view -a -I -s ${i} --thread ${N} > orgpan_vg_final_${i}.vcf
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools norm -m -any --thread ${N} -Oz -o orgpan_SNP_final_norm_${i}.vcf.gz orgpan_SNP_final_${i}.vcf
bcftools norm -m -any --thread ${N} -Oz -o orgpan_delly_final_norm_${i}.vcf.gz orgpan_delly_final_${i}.vcf 
bcftools norm -m -any --thread ${N} -Oz -o mcpan_vg_final_norm_${i}.vcf.gz mcpan_vg_final_${i}.vcf 
bcftools norm -m -any --thread ${N} -Oz -o orgpan_vg_final_norm_${i}.vcf.gz orgpan_vg_final_${i}.vcf 
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_mcpan/mcpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_mcpan/mcpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfdecompose --break-mnps --break-indels -i orgpan_SNP_final_norm_${i}.vcf.gz -o orgpan_SNP_final_norm_decomp_${i}.vcf.gz # the multi-nucleotide polymorphisms and complex indels were further decomposed into SNPs and simple indels
done

cp -r ../ncbi_mcpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_mcpan/mcpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_mcpan/mcpan_delly_final_norm_${i}.vcf.gz.tbi ./
bcftools index -t --threads ${N} orgpan_delly_final_norm_${i}.vcf.gz
bcftools index -t --threads ${N} mcpan_vg_final_norm_${i}.vcf.gz
bcftools index -t --threads ${N} orgpan_vg_final_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b orgpan_delly_final_norm_${i}.vcf.gz -c mcpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
truvari bench -b orgpan_vg_final_norm_${i}.vcf.gz -c mcpan_vg_final_norm_${i}.vcf.gz -o ./sv_comp/vg/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Original pangenome vs. VG pangenome
mkdir -p ./orgpan_vgpan/
cd ./orgpan_vgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

cp ${VGPAN}/called/sv/vg/${PREFIX}2_vgSV.decomp.annot.ext.vcf ./vgpan_vg_final.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 
module load biocontainers
module load picard 
 
picard UpdateVcfSequenceDictionary I=./vgpan_vg_final.vcf O=./vgpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_vg_final_updt.vcf O=./vgpan_vg_final_sorted.vcf

### Converting multi-sample vcf files to per-sample vcf files for ease of variant comparison
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
sed 's/##INFO=<ID=AT,Number=R/##INFO=<ID=AT,Number=1/g' vgpan_vg_final_sorted.vcf | bcftools view -a -I -s ${i} --thread ${N} > vgpan_vg_final_${i}.vcf # sed step to bypass bcftools error
done

### Spliting multi-allelic sites into bi-allelic records
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
bcftools norm -m -any --thread ${N} -Oz -o vgpan_vg_final_norm_${i}.vcf.gz vgpan_vg_final_${i}.vcf 
done

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_mcpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_norm_decomp_${i}.vcf.gz -c vgpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../orgpan_mcpan/orgpan_vg_final_norm_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_vg_final_norm_${i}.vcf.gz.tbi ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz.tbi ./
bcftools index -t --threads ${N} vgpan_vg_final_norm_${i}.vcf.gz
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b orgpan_delly_final_norm_${i}.vcf.gz -c vgpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
truvari bench -b orgpan_vg_final_norm_${i}.vcf.gz -c vgpan_vg_final_norm_${i}.vcf.gz -o ./sv_comp/vg/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Linear representative genome vs. Original pangenome
mkdir -p ./ncbi_orgpan/
cd ./ncbi_orgpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b ncbi_SNP_final_norm_decomp_${i}.vcf.gz -c orgpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_linpan/ncbi_delly_final_norm_${i}.vcf.gz.tbi ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b ncbi_delly_final_norm_${i}.vcf.gz -c orgpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## Original pangenome vs. Linear pangenome
mkdir -p ./orgpan_linpan/
cd ./orgpan_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/linpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_linpan/linpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b orgpan_SNP_final_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/linpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_linpan/linpan_delly_final_norm_${i}.vcf.gz.tbi ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz ./
cp ../orgpan_mcpan/orgpan_delly_final_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b orgpan_delly_final_norm_${i}.vcf.gz -c linpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## VG pangenome vs. Linear pangenome
mkdir -p ./vgpan_linpan/
cd ./vgpan_linpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/linpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_linpan/linpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_norm_decomp_${i}.vcf.gz -c linpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_linpan/linpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_linpan/linpan_delly_final_norm_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b vgpan_delly_final_norm_${i}.vcf.gz -c linpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../

## VG pangenome vs. Minigraph-Cactus pangenome
mkdir -p ./vgpan_mcpan/
cd ./vgpan_mcpan/
mkdir -p ./snp_comp/
mkdir -p ./sv_comp/vg/
mkdir -p ./sv_comp/delly/

### Comparing with the truth SNP call set (for each individual)
for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../ncbi_mcpan/mcpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_mcpan/mcpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_SNP_final_norm_decomp_${i}.vcf.gz.tbi ./
done

cp -r ../ncbi_linpan/linpan.sdf ./

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
${DEPOT}/apps/rtg-tools-3.12.1/rtg vcfeval -T ${N} --output-mode annotate --all-records --ref-overlap --no-roc -b vgpan_SNP_final_norm_decomp_${i}.vcf.gz -c mcpan_SNP_final_norm_decomp_${i}.vcf.gz -o ./snp_comp/${i} -t linpan.sdf
done

### Comparing with the truth SV call set (for each individual)
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bcftools

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
cp ../orgpan_mcpan/mcpan_vg_final_norm_${i}.vcf.gz ./
cp ../orgpan_mcpan/mcpan_vg_final_norm_${i}.vcf.gz.tbi ./
cp ../ncbi_mcpan/mcpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_mcpan/mcpan_delly_final_norm_${i}.vcf.gz.tbi ./
cp ../orgpan_vgpan/vgpan_vg_final_norm_${i}.vcf.gz ./
cp ../orgpan_vgpan/vgpan_vg_final_norm_${i}.vcf.gz.tbi ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz ./
cp ../ncbi_vgpan/vgpan_delly_final_norm_${i}.vcf.gz.tbi ./
done

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
truvari bench -b vgpan_delly_final_norm_${i}.vcf.gz -c mcpan_delly_final_norm_${i}.vcf.gz -o ./sv_comp/delly/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
truvari bench -b vgpan_vg_final_norm_${i}.vcf.gz -c mcpan_vg_final_norm_${i}.vcf.gz -o ./sv_comp/vg/${i} --pick multi -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 100000 
done

cd ../


# Conducting population genetic analyses
module purge
module load biocontainers
module load angsd
#module load pcangsd
module load vcftools
module load bcftools
module load bedtools
module load htslib
module load r
module load anaconda
module load picard
module load ngsld
module load samtools

cd ${BASE}
mkdir -p ./popgen/
cd ./popgen/

## Linear representative genome
mkdir -p ./ncbi/
cd ./ncbi/
cp ${BASE}/benchmarking/ncbi_linpan/ncbi_SNP_final.vcf ./
cp ${BASE}/benchmarking/ncbi_linpan/ncbi_delly_final.vcf ./

### Formating SV vcf file for ANGSD
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ ncbi_delly_final.vcf > ncbi_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # but see if GL field exists. If not, do this before picard steps

grep -v ^\# ncbi_delly_final_simpl.vcf | wc -l # 1272
bcftools view --max-alleles 2 -Ov -o ncbi_delly_final_bi.vcf --threads ${N} ncbi_delly_final_simpl.vcf # keep bi-allelic loci only
grep -v ^\# ncbi_delly_final_bi.vcf | wc -l #1272

bcftools query -f '%CHROM %POS %REF %ALT\n' ncbi_delly_final_bi.vcf > ncbi_delly_final_bi.variants # keep info about SV & extract it
bcftools query -f '%CHROM\t%POS\n' ncbi_delly_final_bi.vcf > ncbi_delly_final_bi.chrpos
head -n 2 ncbi_delly_final_bi.variants

Rscript ${HOME}/extract_SV_info_graph.r ncbi_delly_final_bi.variants
head -n 2 ncbi_delly_final_bi.variants.bed

python3 ${HOME}/fasta_extract_flanking_regions_jyj.py ${REF}/${GENOME}.fna ncbi_delly_final_bi.chrpos 1 ncbi_delly_final_bi.variants.ref # extract ref allele at position for dummy vcf
grep ^"#" ncbi_delly_final_bi.vcf | grep -v ^\#\#"contig=<ID=" > ncbi_delly_final_bi.variants.header # keep header
grep -v ^"#" ncbi_delly_final_bi.vcf > ncbi_delly_final_bi.withoutheader # keep without header
Rscript ${HOME}/make_dummy_vcf_snp_jyj.r ncbi_delly_final_bi.withoutheader ncbi_delly_final_bi.variants.ref # edit in REF
cat ncbi_delly_final_bi.variants.header ncbi_delly_final_bi.withoutheader.withdummyREFALT > ncbi_delly_final_bi.forangsd.vcf # format vcf for ANGSD

### Generating GL file or beagle file for each chromosome (modified from Merot 2023)
ls ${REF}/mapped/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${REF}/${GENOME}.fna.fai -o ncbi_delly_final_bi_rehead.forangsd.vcf --threads ${N} ncbi_delly_final_bi.forangsd.vcf
bgzip ncbi_delly_final_bi_rehead.forangsd.vcf 
tabix ncbi_delly_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl ncbi_delly_final_bi_rehead.forangsd.vcf.gz -nind 20 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

### Nucleotide diversity
realSFS -P ${N} out_snp.saf.idx -maxiter 100 -fold 1 > out_snp.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_snp.saf.idx -sfs out_snp.sfs -P ${N} -outname out_snp # calculate thetas for each site

thetaStat do_stat out_snp.thetas.idx -win 50000 -step 50000  -outnames theta_snp.thetasWindow.gz # sliding window estimate
awk '{print $4,$5,$14}' theta_snp.thetasWindow.gz.pestPG > Thetas
sumSites=$(awk 'BEGIN{s=0;}{s=s+$3;}END{print s;}' Thetas) # get number of sites
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas) # get mean Watterson's Theta
meanW=$(awk "BEGIN {print $sumW/$sumSites}") # 0.00954444
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas) # get mean Tajima's Theta
meanN=$(awk "BEGIN {print $sumN/$sumSites}") # 0.0054623

realSFS -P ${N} out_sv.saf.idx -maxiter 100 -fold 1 > out_sv.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_sv.saf.idx -sfs out_sv.sfs -P ${N} -outname out_sv # calculate thetas for each site

### SFS
#R
#norm <- function(x) x/sum(x) # function to normalize 
#sfs <- (scan("out_snp.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("ncbi_snp_sfs_plot.svg", width=5, height=5)
#barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

#sfs <- (scan("out_sv.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("ncbi_sv_sfs_plot.svg", width=5, height=5)
#barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

### ROH
angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_ncbi_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_ncbi_PLraw.txt ROH_ncbi_PLresult.txt ${REF}/${GENOME}.fna.fai # parse raw ROH file

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca # generate SNP beagle

#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand ncbi_SNP_final.vcf -Oz -o ncbi_SNP_final_ld.vcf.gz # prune linkage disequilibrium
#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand ncbi_delly_final.vcf -Oz -o ncbi_delly_final_ld.vcf.gz 

picard UpdateVcfSequenceDictionary I=./ncbi_delly_final_simpl.vcf O=./ncbi_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./ncbi_delly_final_updt.vcf O=./ncbi_delly_final_sorted.vcf
bcftools +tag2tag ncbi_delly_final_sorted.vcf -- -r --pl-to-gl > ncbi_delly_final_sortedGL.vcf

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | grep -v "unmap" | while read i # generate SV beagle from here
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf ncbi_delly_final_sortedGL.vcf --BEAGLE-GL --chr ${i} --out ncbi_delly_${i}
done 

head -n 1 ncbi_delly_NC_053450.1.BEAGLE.GL > ncbi_delly.beagleheader # get the header

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | grep -v "unmap" | while read i
do 
echo "append beagle from ${i}"
tail -n +2 ncbi_delly_${i}.BEAGLE.GL >> ncbi_delly.beagle
done

wc -l ncbi_delly.beagle # should be equal to the nb of variants minus non-biallelic variants in the vcf
cat ncbi_delly_final.vcf | grep -v "#" | wc -l
cat ncbi_delly_final.vcf | grep -v "#" | cut -f5 | grep "," | wc -l
rm ncbi_delly*log
rm ncbi_delly*BEAGLE.GL

cat ncbi_delly.beagle | cut -f 1 | sed 's/:/\t/g'| gzip > sv_pca.pos.gz # prepare a pos file 
sed -i 's/:/_/g' ncbi_delly.beagle # replace the : with _ in the position information
sed -i 's/\[/#/g' ncbi_delly.beagle # replace the [ with # to avoid file reading issue 
sed -i 's/#/_/g' ncbi_delly.beagle # replace the # with _ to avoid file reading issue 

Rscript ${HOME}/normalize_beagle.r ncbi_delly.beagle ncbi_delly.norm.beagle  # normalize the genotype likelihoods as otherwise too small numbers are counted as zero and cannot be divided

cat ncbi_delly.beagleheader ncbi_delly.norm.beagle > ncbi_delly.ready.beagle #add the  proper header again and zip
head ncbi_delly.ready.beagle
gzip ncbi_delly.ready.beagle

zcat snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file 
ngsLD --geno snp_pca_subsampled.geno.beagle.gz --pos snp_pca_subsampled.pos.gz --probs --n_ind 25 --n_sites 922896 --max_kb_dist 1000 --n_threads ${N} --out out_snp.ngsld # run "zcat snp_pca_subsampled.pos.gz | wc -l" to get n_sites
cat out_snp.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5547 156034184 # check by "cat out_snp.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_snp.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat out_snp.ngsld | grep "NC_053450.1" > out_snp_1stchr.ngsld
#cat out_snp.ngsld | grep "NC_053480.1" > out_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
cat out_snp_1stchr.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5000 1005000
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_snp_1stchr.ngsld > out_snp_1stchr_header.ngsld
echo "${PWD}/out_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, LD decays quickly

zcat ncbi_delly.ready.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz # prepare a geno file 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 1272 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites # it seems there is a duplicated position - let's remove one

zcat sv_pca.pos.gz | awk '{if ($2 in seen) {print} else {seen[$2]=1}}' # check duplicates -> "NW_026690477.1  7341" is duplicated; 

zcat ncbi_delly.ready.beagle.gz | awk '!seen[$1]++' > ncbi_delly.uniq.beagle
cat ncbi_delly.uniq.beagle | tail -n +2 | cut -f 1 | cut -f 1,2 -d "_" > sv_pca.pos1
cat ncbi_delly.uniq.beagle | tail -n +2 | cut -f 1 | cut -f 3 -d "_" > sv_pca.pos2
paste -d '\t' sv_pca.pos1 sv_pca.pos2 | gzip > sv_pca.pos.gz
head ncbi_delly.uniq.beagle
gzip ncbi_delly.uniq.beagle
zcat ncbi_delly.uniq.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 1271 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites
cat out_sv.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 152582 155811490 # check by "cat out_sv.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_sv.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat decay.header out_sv.ngsld > out_sv_header.ngsld
echo "${PWD}/out_sv_header.ngsld" > decay_sv.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_sv.ld --max_kb_dist 1000 --out fit_SV_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, pruning seems not needed

${APP}/prune_graph/target/release/prune_graph --in out_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_snp_unlinked.ngsld
sed 's/:/_/g' out_snp_unlinked.ngsld > out_snp_unlinked

${APP}/prune_graph/target/release/prune_graph --in out_sv.ngsld --weight-field "column_7" --weight-filter "column_3 <= 25000 && column_7 >= 0.5" --out out_sv_unlinked.ngsld
sed 's/:/_/g' out_sv_unlinked.ngsld > out_sv_unlinked

zcat snp_pca.beagle.gz | head -n 1 > ncbi_SNP_ld.beagleheader
zcat snp_pca.beagle.gz | grep -Fwf out_snp_unlinked > ncbi_SNP_ld.beagle # keep only unlinked loci
cat ncbi_SNP_ld.beagleheader ncbi_SNP_ld.beagle | gzip > ncbi_SNP_ld.beagle.gz
zcat ncbi_SNP_ld.beagle.gz | head
cat out_snp_unlinked | wc -l # check the number of loci if matching
zcat ncbi_SNP_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

zcat ncbi_delly.uniq.beagle.gz | head -n 1 > ncbi_delly_ld.beagleheader
zcat ncbi_delly.uniq.beagle.gz | grep -Fwf out_sv_unlinked > ncbi_delly_ld.beagle # keep only unlinked loci
cat ncbi_delly_ld.beagleheader ncbi_delly_ld.beagle | gzip > ncbi_delly_ld.beagle.gz
zcat ncbi_delly_ld.beagle.gz | head
cat out_sv_unlinked | wc -l # check the number of loci if matching
zcat ncbi_delly_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

sed 's/:/\t/g' out_snp_unlinked.ngsld > out_snp_unlinked.tab
vcftools --positions out_snp_unlinked.tab --vcf ncbi_SNP_final.vcf --recode --recode-INFO-all --out ncbi_SNP_ld # generate ld-pruned vcf file for downstrean usage

sed 's/:/\t/g' out_sv_unlinked.ngsld > out_sv_unlinked.tab
vcftools --positions out_sv_unlinked.tab --vcf ncbi_delly_final_sortedGL.vcf --recode --recode-INFO-all --out ncbi_delly_ld # generate ld-pruned vcf file for downstrean usage
awk '!seen[$1,$2]++' ncbi_delly_ld.recode.vcf > ncbi_delly_ld.uniq.vcf

pcangsd -b ncbi_SNP_ld.beagle.gz -o ncbi_snp_pca -t ${N}
pcangsd -b ncbi_delly_ld.beagle.gz -o ncbi_sv_pca -t ${N}

R
library(ggplot2)
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("ncbi_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.3%, PC2: 4.3%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("ncbi_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.3%)", y = "Principal Component 2 (4.3%)") + theme_classic(base_size=14)
dev.off()

C <- as.matrix(read.table("ncbi_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.1%, PC2: 6.0%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("ncbi_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.1%)", y = "Principal Component 2 (6.0%)") + theme_classic(base_size=14)
dev.off()

### Selection
pcangsd -b ncbi_SNP_ld.beagle.gz  -t ${N} -o ncbi_snp_sel --selection --sites_save #selection
pcangsd -b ncbi_delly_ld.beagle.gz  -t ${N} -o ncbi_sv_sel --selection --sites_save #selection

zcat ncbi_SNP_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f1,2 -d "_" > snp_chr
zcat ncbi_SNP_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f3 -d "_" > snp_pos
zcat ncbi_SNP_ld.beagle.gz | tail -n +2 | cut -f1 > snp_var

zcat ncbi_delly_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f1,2 -d "_" > sv_chr
zcat ncbi_delly_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f3 -d "_" > sv_pos
zcat ncbi_delly_ld.beagle.gz | tail -n +2 | cut -f1 > sv_var

R
library(qqman)
library(qvalue)
D <- as.matrix(read.table("ncbi_snp_sel.selection"))
chr <- as.matrix(read.table("snp_chr"))
pos <- as.matrix(read.table("snp_pos"))
var <- as.matrix(read.table("snp_var"))
sites <- as.matrix(read.table("ncbi_snp_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:538) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
png("ncbi_snp_manh.png", res=300, units="in", width=25, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

D <- as.matrix(read.table("ncbi_sv_sel.selection"))
chr <- as.matrix(read.table("sv_chr"))
pos <- as.matrix(read.table("sv_pos"))
var <- as.matrix(read.table("sv_var"))
sites <- as.matrix(read.table("ncbi_sv_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:126) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
svg("ncbi_sv_manh.svg", width=5, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

# Fst outlier
vcftools --vcf ncbi_SNP_ld.recode.vcf --out ncbi_SNP_ld.filtered --max-missing 0.88 --max-alleles 2 --recode --recode-INFO-all # filter out non-biallelic loci and loci with missing data proportion > 3/25*100
bcftools annotate --set-id '%CHROM:%POS' -Ov -o ncbi_SNP_ld.filtered.setID.vcf --threads ${N} ncbi_SNP_ld.filtered.recode.vcf  
bcftools annotate --set-id '%CHROM:%POS' -Ov -o ncbi_delly_ld.uniq.setID.vcf --threads ${N} ncbi_delly_ld.uniq.vcf

R
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)
obj.vcfR <- read.vcfR("ncbi_SNP_ld.filtered.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "snp_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "snp_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK # need to redefine functions as suggested in: https://github.com/whitlock/OutFLANK/issues/19
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # 4 outliers (NC_053456.1:11518650, NC_053456.1:26878073, NC_053457.1:2496665, NC_053459.1:14512222)

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "snp_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "snp_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)

obj.vcfR <- read.vcfR("ncbi_delly_ld.uniq.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "sv_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "sv_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # 1 outlier (NC_053456.1:23954480)

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "sv_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)
quit()

cd ../

## Linear pangenome
mkdir -p ./linpan/
cd ./linpan/
cp ${BASE}/benchmarking/ncbi_linpan/linpan_SNP_final.vcf ./
cp ${BASE}/benchmarking/ncbi_linpan/linpan_delly_final.vcf ./

### Formating SV vcf file for ANGSD
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ linpan_delly_final.vcf > linpan_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # but see if GL field exists. If not, do this before picard steps

grep -v ^\# linpan_delly_final_simpl.vcf | wc -l # 1292
bcftools view --max-alleles 2 -Ov -o linpan_delly_final_bi.vcf --threads ${N} linpan_delly_final_simpl.vcf # keep bi-allelic loci only
grep -v ^\# linpan_delly_final_bi.vcf | wc -l #1292

bcftools query -f '%CHROM %POS %REF %ALT\n' linpan_delly_final_bi.vcf > linpan_delly_final_bi.variants # keep info about SV & extract it
bcftools query -f '%CHROM\t%POS\n' linpan_delly_final_bi.vcf > linpan_delly_final_bi.chrpos
head -n 2 linpan_delly_final_bi.variants

Rscript ${HOME}/extract_SV_info_graph.r linpan_delly_final_bi.variants
head -n 2 linpan_delly_final_bi.variants.bed

python3 ${HOME}/fasta_extract_flanking_regions_jyj.py ${LINPAN}/${PREFIX}_panref2_sorted.fa linpan_delly_final_bi.chrpos 1 linpan_delly_final_bi.variants.ref # extract ref allele at position for dummy vcf
grep ^"#" linpan_delly_final_bi.vcf | grep -v ^\#\#"contig=<ID=" > linpan_delly_final_bi.variants.header # keep header
grep -v ^"#" linpan_delly_final_bi.vcf > linpan_delly_final_bi.withoutheader # keep without header
Rscript ${HOME}/make_dummy_vcf_snp_jyj.r linpan_delly_final_bi.withoutheader linpan_delly_final_bi.variants.ref # edit in REF
cat linpan_delly_final_bi.variants.header linpan_delly_final_bi.withoutheader.withdummyREFALT > linpan_delly_final_bi.forangsd.vcf # format vcf for ANGSD

### Generating GL file or beagle file for each chromosome (modified from Merot 2023)
ls /${LINPAN}/mapped/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o linpan_delly_final_bi_rehead.forangsd.vcf --threads ${N} linpan_delly_final_bi.forangsd.vcf
bgzip linpan_delly_final_bi_rehead.forangsd.vcf 
tabix linpan_delly_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl linpan_delly_final_bi_rehead.forangsd.vcf.gz -nind 20 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

### Nucleotide diversity
realSFS -P ${N} out_snp.saf.idx -maxiter 100 -fold 1 > out_snp.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_snp.saf.idx -sfs out_snp.sfs -P ${N} -outname out_snp # calculate thetas for each site

thetaStat do_stat out_snp.thetas.idx -win 50000 -step 50000  -outnames theta_snp.thetasWindow.gz # sliding window estimate
awk '{print $4,$5,$14}' theta_snp.thetasWindow.gz.pestPG > Thetas
sumSites=$(awk 'BEGIN{s=0;}{s=s+$3;}END{print s;}' Thetas) # get number of sites
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas) # get mean Watterson's Theta
meanW=$(awk "BEGIN {print $sumW/$sumSites}") # 0.00954436
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas) # get mean Tajima's Theta
meanN=$(awk "BEGIN {print $sumN/$sumSites}") # 0.00546215

realSFS -P ${N} out_sv.saf.idx -maxiter 100 -fold 1 > out_sv.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_sv.saf.idx -sfs out_sv.sfs -P ${N} -outname out_sv # calculate thetas for each site

### SFS
#R
#norm <- function(x) x/sum(x) # function to normalize 
#sfs <- (scan("out_snp.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("linpan_snp_sfs_plot.svg", width=5, height=5)
#barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

#sfs <- (scan("out_sv.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("linpan_sv_sfs_plot.svg", width=5, height=5)
#barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_linpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_linpan_PLraw.txt ROH_linpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca # generate SNP beagle

#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand linpan_SNP_final_sorted_GL.vcf.gz -Oz -o linpan_SNP_final_sorted_GL_ld.vcf.gz # prune linkage disequilibrium
#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand linpan_delly_final.vcf -Oz -o linpan_delly_final_ld.vcf.gz 

picard UpdateVcfSequenceDictionary I=./linpan_delly_final_simpl.vcf O=./linpan_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./linpan_delly_final_updt.vcf O=./linpan_delly_final_sorted.vcf
bcftools +tag2tag linpan_delly_final_sorted.vcf -- -r --pl-to-gl > linpan_delly_final_sortedGL.vcf

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i # generate SV beagle from here
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf linpan_delly_final_sorted.vcf --BEAGLE-GL --chr ${i} --out linpan_delly_${i}
done 

head -n 1 linpan_delly_NC_053450.1.BEAGLE.GL > linpan_delly.beagleheader # get the header

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i
do 
echo "append beagle from ${i}"
tail -n +2 linpan_delly_${i}.BEAGLE.GL >> linpan_delly.beagle
done

wc -l linpan_delly.beagle # should be equal to the nb of variants minus non-biallelic variants in the vcf
cat linpan_delly_final.vcf | grep -v "#" | wc -l
cat linpan_delly_final.vcf | grep -v "#" | cut -f5 | grep "," | wc -l
rm linpan_delly*log
rm linpan_delly*BEAGLE.GL

cat linpan_delly.beagle | cut -f 1 | sed 's/:/\t/g'| gzip > sv_pca.pos.gz # prepare a pos file 
sed -i 's/:/_/g' linpan_delly.beagle # replace the : with _ in the position information
sed -i 's/\[/#/g' linpan_delly.beagle # replace the [ with # to avoid file reading issue 
sed -i 's/#/_/g' linpan_delly.beagle # replace the # with _ to avoid file reading issue 

Rscript ${HOME}/normalize_beagle.r linpan_delly.beagle linpan_delly.norm.beagle  # normalize the genotype likelihoods as otherwise too small numbers are counted as zero and cannot be divided

cat linpan_delly.beagleheader linpan_delly.norm.beagle > linpan_delly.ready.beagle #add the  proper header again and zip
head linpan_delly.ready.beagle
gzip linpan_delly.ready.beagle

zcat snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file 
ngsLD --geno snp_pca_subsampled.geno.beagle.gz --pos snp_pca_subsampled.pos.gz --probs --n_ind 25 --n_sites 922896 --max_kb_dist 1000 --n_threads ${N} --out out_snp.ngsld # run "zcat snp_pca_subsampled.pos.gz | wc -l" to get n_sites
cat out_snp.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5547 156034184 # check by "cat out_snp.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_snp.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat out_snp.ngsld | grep "NC_053450.1" > out_snp_1stchr.ngsld
#cat out_snp.ngsld | grep "NC_053480.1" > out_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
cat out_snp_1stchr.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5000 1005000
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_snp_1stchr.ngsld > out_snp_1stchr_header.ngsld
echo "${PWD}/out_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, LD decays quickly

zcat linpan_delly.ready.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz # prepare a geno file 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 1292 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites # it seems there are duplicated positions - let's remove one

zcat sv_pca.pos.gz | awk '{if ($2 in seen) {print} else {seen[$2]=1}}' # check duplicates -> "NW_026690477.1  7341" and "unmap_tig00006430  19409" is duplicated; 

zcat linpan_delly.ready.beagle.gz | awk '!seen[$1]++' > linpan_delly.uniq.beagle
cat linpan_delly.uniq.beagle | tail -n +2 | cut -f 1 | cut -f 1,2 -d "_" > sv_pca.pos1
cat linpan_delly.uniq.beagle | tail -n +2 | cut -f 1 | cut -f 3 -d "_" > sv_pca.pos2
paste -d '\t' sv_pca.pos1 sv_pca.pos2 | gzip > sv_pca.pos.gz
head linpan_delly.uniq.beagle
gzip linpan_delly.uniq.beagle
zcat linpan_delly.uniq.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 1291 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites
cat out_sv.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 152582 155811490 # check by "cat out_sv.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_sv.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat decay.header out_sv.ngsld > out_sv_header.ngsld
echo "${PWD}/out_sv_header.ngsld" > decay_sv.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_sv.ld --max_kb_dist 1000 --out fit_SV_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, pruning seems not needed

${APP}/prune_graph/target/release/prune_graph --in out_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_snp_unlinked.ngsld
sed 's/:/_/g' out_snp_unlinked.ngsld > out_snp_unlinked

${APP}/prune_graph/target/release/prune_graph --in out_sv.ngsld --weight-field "column_7" --weight-filter "column_3 <= 25000 && column_7 >= 0.5" --out out_sv_unlinked.ngsld
sed 's/:/_/g' out_sv_unlinked.ngsld > out_sv_unlinked

zcat snp_pca.beagle.gz | head -n 1 > linpan_SNP_ld.beagleheader
zcat snp_pca.beagle.gz | grep -Fwf out_snp_unlinked > linpan_SNP_ld.beagle # keep only unlinked loci
cat linpan_SNP_ld.beagleheader linpan_SNP_ld.beagle | gzip > linpan_SNP_ld.beagle.gz
zcat linpan_SNP_ld.beagle.gz | head
cat out_snp_unlinked | wc -l # check the number of loci if matching
zcat linpan_SNP_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

zcat linpan_delly.uniq.beagle.gz | head -n 1 > linpan_delly_ld.beagleheader
zcat linpan_delly.uniq.beagle.gz | grep -Fwf out_sv_unlinked > linpan_delly_ld.beagle # keep only unlinked loci
cat linpan_delly_ld.beagleheader linpan_delly_ld.beagle | gzip > linpan_delly_ld.beagle.gz
zcat linpan_delly_ld.beagle.gz | head
cat out_sv_unlinked | wc -l # check the number of loci if matching
zcat linpan_delly_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

sed 's/:/\t/g' out_snp_unlinked.ngsld > out_snp_unlinked.tab
vcftools --positions out_snp_unlinked.tab --vcf linpan_SNP_final.vcf --recode --recode-INFO-all --out linpan_SNP_ld # generate ld-pruned vcf file for downstrean usage

sed 's/:/\t/g' out_sv_unlinked.ngsld > out_sv_unlinked.tab
vcftools --positions out_sv_unlinked.tab --vcf linpan_delly_final_sortedGL.vcf --recode --recode-INFO-all --out linpan_delly_ld # generate ld-pruned vcf file for downstrean usage
awk '!seen[$1,$2]++' linpan_delly_ld.recode.vcf > linpan_delly_ld.uniq.vcf

pcangsd -b linpan_SNP_ld.beagle.gz -o linpan_snp_pca -t ${N}
pcangsd -b linpan_delly_ld.beagle.gz -o linpan_sv_pca -t ${N}

R
library(ggplot2)
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("linpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.2%, PC2: 4.3%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("linpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.2%)", y = "Principal Component 2 (4.3%)") + theme_classic(base_size=14)
dev.off()

C <- as.matrix(read.table("linpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.0%, PC2: 5.9%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("linpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.0%)", y = "Principal Component 2 (5.9%)") + theme_classic(base_size=14)
dev.off()
quit()

### Selection
pcangsd -b linpan_SNP_ld.beagle.gz  -t ${N} -o linpan_snp_sel --selection --sites_save #selection
pcangsd -b linpan_delly_ld.beagle.gz  -t ${N} -o linpan_sv_sel --selection --sites_save #selection

zcat linpan_SNP_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f1,2 -d "_" > snp_chr
zcat linpan_SNP_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f3 -d "_" > snp_pos
zcat linpan_SNP_ld.beagle.gz | tail -n +2 | cut -f1 > snp_var

zcat linpan_delly_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f1,2 -d "_" > sv_chr
zcat linpan_delly_ld.beagle.gz | tail -n +2 | cut -f1 | cut -f3 -d "_" > sv_pos
zcat linpan_delly_ld.beagle.gz | tail -n +2 | cut -f1 > sv_var

R
library(qqman)
library(qvalue)
D <- as.matrix(read.table("linpan_snp_sel.selection"))
chr <- as.matrix(read.table("snp_chr"))
pos <- as.matrix(read.table("snp_pos"))
var <- as.matrix(read.table("snp_var"))
sites <- as.matrix(read.table("linpan_snp_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:843) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
png("linpan_snp_manh.png", res=300, units="in", width=25, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

D <- as.matrix(read.table("linpan_sv_sel.selection"))
chr <- as.matrix(read.table("sv_chr"))
pos <- as.matrix(read.table("sv_pos"))
var <- as.matrix(read.table("sv_var"))
sites <- as.matrix(read.table("linpan_sv_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:142) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
svg("linpan_sv_manh.svg", width=5, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()
quit()

# Fst outlier
vcftools --vcf linpan_SNP_ld.recode.vcf --out linpan_SNP_ld.filtered --max-missing 0.88 --max-alleles 2 --recode --recode-INFO-all # filter out non-biallelic loci and loci with missing data proportion > 3/25*100
bcftools annotate --set-id '%CHROM:%POS' -Ov -o linpan_SNP_ld.filtered.setID.vcf --threads ${N} linpan_SNP_ld.filtered.recode.vcf  
bcftools annotate --set-id '%CHROM:%POS' -Ov -o linpan_delly_ld.uniq.setID.vcf --threads ${N} linpan_delly_ld.uniq.vcf

R
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)
obj.vcfR <- read.vcfR("linpan_SNP_ld.filtered.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "snp_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "snp_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK # need to redefine functions as suggested in: https://github.com/whitlock/OutFLANK/issues/19
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outliers

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "snp_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "snp_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)

obj.vcfR <- read.vcfR("linpan_delly_ld.uniq.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "sv_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "sv_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # 1 outlier (NC_053456.1:23954480)

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "sv_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)
quit()

cd ../

## Minigraph-Cactus pangenome
mkdir -p ./mcpan/
cd ./mcpan/
cp ${BASE}/benchmarking/ncbi_mcpan/mcpan_SNP_final_sorted.vcf ./
cp ${BASE}/benchmarking/ncbi_mcpan/mcpan_delly_final_sorted.vcf ./
cp ${MCPAN}/${PREFIX}-pg/called/sv/vg/${PREFIX}_mcSV.filtered.popped.vcf ./mcpan_vg_final.vcf

### Formating a merged SV vcf file for ANGSD
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ mcpan_vg_final_sorted.vcf > mcpan_vg_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # but see if GL field exists. If not, do this before picard steps
picard UpdateVcfSequenceDictionary I=./mcpan_vg_final_simpl.vcf O=./mcpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./mcpan_vg_final_updt.vcf O=./mcpan_vg_final_sorted.vcf
bcftools +tag2tag mcpan_vg_final_sorted.vcf -- -r --pl-to-gl > mcpan_vg_final_sortedGL.vcf

bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ mcpan_delly_final_sorted.vcf > mcpan_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # but see if GL field exists. If not, do this before picard steps
bcftools +tag2tag mcpan_delly_final_simpl.vcf -- -r --pl-to-gl > mcpan_delly_final_sortedGL.vcf

bcftools query -f '%CHROM\t%POS0\t%END\n' mcpan_delly_final_sortedGL.vcf > mcpan_delly.bed # convert vcf region to bed
bcftools query -f '%CHROM\t%POS0\t%END\n' mcpan_vg_final_sortedGL.vcf > mcpan_vg.bed # convert vcf region to bed

cut -f1,2 ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai > linpan.genome # create genome file for bedtools
bedtools slop -i mcpan_delly.bed -g linpan.genome -b 100 > mcpan_delly_padded.bed # pad 100 bp 
bedtools slop -i mcpan_vg.bed -g linpan.genome -b 100 > mcpan_vg_padded.bed # pad 100 bp
bedtools intersect -a mcpan_delly_padded.bed -b mcpan_vg_padded.bed > sv_overlapped.bed # nothing. If there is, overlapped region should not be duplicated in the below step
bedtools intersect -header -a mcpan_delly_final_sortedGL.vcf -b sv_overlapped.bed -v > mcpan_delly_final_excluded.vcf # just keep this nonnecessary code for reproducibility

bcftools view -Oz -o mcpan_delly_final_excluded.vcf.gz --threads ${N} mcpan_delly_final_excluded.vcf
bcftools view -Oz -o mcpan_vg_final_sortedGL.vcf.gz --threads ${N} mcpan_vg_final_sortedGL.vcf
bcftools index mcpan_delly_final_excluded.vcf.gz
bcftools index mcpan_vg_final_sortedGL.vcf.gz
bcftools concat --allow-overlaps -Ov -o mcpan_SV_final.vcf --threads ${N} mcpan_delly_final_excluded.vcf.gz mcpan_vg_final_sortedGL.vcf.gz # merge two SV vcf files

grep -v ^\# mcpan_SV_final.vcf | wc -l # 2051
bcftools view --max-alleles 2 -Ov -o mcpan_SV_final_bi.vcf --threads ${N} mcpan_SV_final.vcf # keep bi-allelic loci only
grep -v ^\# mcpan_SV_final_bi.vcf | wc -l #2049

bcftools query -f '%CHROM %POS %REF %ALT\n' mcpan_SV_final_bi.vcf > mcpan_SV_final_bi.variants # keep info about SV & extract it
bcftools query -f '%CHROM\t%POS\n' mcpan_SV_final_bi.vcf > mcpan_SV_final_bi.chrpos
head -n 2 mcpan_SV_final_bi.variants

Rscript ${HOME}/extract_SV_info_graph.r mcpan_SV_final_bi.variants
head -n 2 mcpan_SV_final_bi.variants.bed

python3 ${HOME}/fasta_extract_flanking_regions_jyj.py ${LINPAN}/${PREFIX}_panref2_sorted.fa mcpan_SV_final_bi.chrpos 1 mcpan_SV_final_bi.variants.ref # extract ref allele at position for dummy vcf
grep ^"#" mcpan_SV_final_bi.vcf | grep -v ^\#\#"contig=<ID=" > mcpan_SV_final_bi.variants.header # keep header
grep -v ^"#" mcpan_SV_final_bi.vcf > mcpan_SV_final_bi.withoutheader # keep without header
Rscript ${HOME}/make_dummy_vcf_snp_jyj.r mcpan_SV_final_bi.withoutheader mcpan_SV_final_bi.variants.ref # edit in REF
cat mcpan_SV_final_bi.variants.header mcpan_SV_final_bi.withoutheader.withdummyREFALT > mcpan_SV_final_bi.forangsd.vcf # format vcf for ANGSD

### generate GL or beagle file for each chromosome (modified from Merot 2023)
ls ${MCPAN}/${PREFIX}-pg/aligned/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o mcpan_SV_final_bi_rehead.forangsd.vcf --threads ${N} mcpan_SV_final_bi.forangsd.vcf
bgzip  mcpan_SV_final_bi_rehead.forangsd.vcf
tabix mcpan_SV_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl mcpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 20 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

### Nucleotide diversity
realSFS -P ${N} out_snp.saf.idx -maxiter 100 -fold 1 > out_snp.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_snp.saf.idx -sfs out_snp.sfs -P ${N} -outname out_snp # calculate thetas for each site

thetaStat do_stat out_snp.thetas.idx -win 50000 -step 10000  -outnames theta_snp.thetasWindow.gz # sliding window estimate
awk '{print $4,$5,$14}' theta_snp.thetasWindow.gz.pestPG > Thetas
sumSites=$(awk 'BEGIN{s=0;}{s=s+$3;}END{print s;}' Thetas) # get number of sites
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas) # get mean Watterson's Theta
meanW=$(awk "BEGIN {print $sumW/$sumSites}") # 0.00949753
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas) # get mean Tajima's Theta
meanN=$(awk "BEGIN {print $sumN/$sumSites}") # 0.00540129

realSFS -P ${N} out_sv.saf.idx -maxiter 100 -fold 1 > out_sv.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_sv.saf.idx -sfs out_sv.sfs -P ${N} -outname out_sv # calculate thetas for each site

### SFS
#R
#norm <- function(x) x/sum(x) # function to normalize 
#sfs <- (scan("out_snp.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("mcpan_snp_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

#sfs <- (scan("out_sv.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("mcpan_sv_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_mcpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_mcpan_PLraw.txt ROH_mcpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa.fai # parse raw ROH file

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca

#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand mcpan_SNP_final_sorted_GL.vcf.gz -Oz -o mcpan_SNP_final_sorted_GL_ld.vcf.gz # prune linkage disequilibrium
#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand mcpan_SV_final.vcf -Oz -o mcpan_SV_final_ld.vcf.gz 

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i 
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf mcpan_SV_final.vcf --BEAGLE-GL --chr ${i} --out mcpan_SV_${i}
done 

head -n 1 mcpan_SV_NC_053450.1.BEAGLE.GL > mcpan_SV.beagleheader # get the headerhead 

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i
do 
echo "append beagle from ${i}"
tail -n +2 mcpan_SV_${i}.BEAGLE.GL >> mcpan_SV.beagle
done

cat mcpan_SV_final.vcf | grep -v "#" | wc -l
wc -l mcpan_SV.beagle # should be equal to the nb of variants minus non-biallelic variants in the vcf
cat mcpan_SV_final.vcf | grep -v "#" | cut -f5 | grep "," | wc -l
rm mcpan_SV*log
rm mcpan_SV*BEAGLE.GL

cat mcpan_SV.beagle | cut -f 1 | sed 's/:/\t/g'| gzip > sv_pca.pos.gz # prepare a pos file 
sed -i 's/:/_/g' mcpan_SV.beagle
sed -i 's/\[/#/g' mcpan_SV.beagle # replace the [ with # to avoid file reading issue 
sed -i 's/#/_/g' mcpan_SV.beagle # replace the # with _ to avoid file reading issue 

Rscript ${HOME}/normalize_beagle.r mcpan_SV.beagle mcpan_SV.norm.beagle  # normalize the genotype likelihoods as otherwise too small numbers are counted as zero and cannot be divided

cat mcpan_SV.beagleheader mcpan_SV.norm.beagle > mcpan_SV.ready.beagle #add the  proper header again and zip
head mcpan_SV.ready.beagle
gzip mcpan_SV.ready.beagle

zcat snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file 
ngsLD --geno snp_pca_subsampled.geno.beagle.gz --pos snp_pca_subsampled.pos.gz --probs --n_ind 25 --n_sites 850446 --max_kb_dist 1000 --n_threads ${N} --out out_snp.ngsld # run "zcat snp_pca_subsampled.pos.gz | wc -l" to get n_sites
cat out_snp.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5040 156021255 # check by "cat out_snp.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_snp.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat out_snp.ngsld | grep "NC_053450.1" > out_snp_1stchr.ngsld
#cat out_snp.ngsld | grep "NC_053480.1" > out_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
cat out_snp_1stchr.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5000 1005000
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_snp_1stchr.ngsld > out_snp_1stchr_header.ngsld
sed -i 's/bHirRus1_LinPan#0#//g' out_snp_1stchr_header.ngsld # avoid "scan()" error in R below
echo "${PWD}/out_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, LD decays quickly

zcat mcpan_SV.ready.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz # prepare a geno file 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 2049 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites
cat out_sv.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 279064 155346731 # check by "cat out_sv.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_sv.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat decay.header out_sv.ngsld > out_sv_header.ngsld
echo "${PWD}/out_sv_header.ngsld" > decay_sv.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_sv.ld --max_kb_dist 1000 --out fit_SV_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, pruning seems not needed

${APP}/prune_graph/target/release/prune_graph --in out_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_snp_unlinked.ngsld
sed 's/:/_/g' out_snp_unlinked.ngsld > out_snp_unlinked

${APP}/prune_graph/target/release/prune_graph --in out_sv.ngsld --weight-field "column_7" --weight-filter "column_3 <= 25000 && column_7 >= 0.5" --out out_sv_unlinked.ngsld
sed 's/:/_/g' out_sv_unlinked.ngsld > out_sv_unlinked

zcat snp_pca.beagle.gz | head -n 1 > mcpan_SNP_ld.beagleheader
zcat snp_pca.beagle.gz | grep -Fwf out_snp_unlinked > mcpan_SNP_ld.beagle # keep only unlinked loci
sed -i 's/bHirRus1_LinPan#0#//g' mcpan_SNP_ld.beagle # polish chromosome names
cat mcpan_SNP_ld.beagleheader mcpan_SNP_ld.beagle | gzip > mcpan_SNP_ld.beagle.gz
zcat mcpan_SNP_ld.beagle.gz | head
cat out_snp_unlinked | wc -l # check the number of loci if matching
zcat mcpan_SNP_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

zcat mcpan_SV.ready.beagle.gz | head -n 1 > mcpan_SV_ld.beagleheader
zcat mcpan_SV.ready.beagle.gz | grep -Fwf out_sv_unlinked > mcpan_SV_ld.beagle # keep only unlinked loci
cat mcpan_SV_ld.beagleheader mcpan_SV_ld.beagle | gzip > mcpan_SV_ld.beagle.gz
zcat mcpan_SV_ld.beagle.gz | head
cat out_sv_unlinked | wc -l # check the number of loci if matching
zcat mcpan_SV_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

sed 's/:/\t/g' out_snp_unlinked.ngsld > out_snp_unlinked.tab
sed -i 's/bHirRus1_LinPan#0#//g' out_snp_unlinked.tab # polish chromosome names
vcftools --positions out_snp_unlinked.tab --vcf mcpan_SNP_final_sorted.vcf --recode --recode-INFO-all --out mcpan_SNP_ld # generate ld-pruned vcf file for downstrean usage

sed 's/:/\t/g' out_sv_unlinked.ngsld > out_sv_unlinked.tab
vcftools --positions out_sv_unlinked.tab --vcf mcpan_SV_final.vcf --recode --recode-INFO-all --out mcpan_SV_ld # generate ld-pruned vcf file for downstrean usage

pcangsd -b mcpan_SNP_ld.beagle.gz -o mcpan_snp_pca -t ${N}
pcangsd -b mcpan_SV_ld.beagle.gz -o mcpan_sv_pca -t ${N}

R
library(ggplot2)
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("mcpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.2%, PC2: 4.3%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("mcpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.2%)", y = "Principal Component 2 (4.3%)") + theme_classic(base_size=14)
dev.off()

C <- as.matrix(read.table("mcpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.5%, PC2: 5.2%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("mcpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.5%)", y = "Principal Component 2 (5.2%)") + theme_classic(base_size=14)
dev.off()
quit()

### Selection
pcangsd -b mcpan_SNP_ld.beagle.gz -t ${N} -o mcpan_snp_sel --selection --sites_save # selection
pcangsd -b mcpan_SV_ld.beagle.gz -t ${N} -o mcpan_sv_sel --selection --sites_save # selection
#pcangsd -b mcpan_SV.ready.beagle.gz -o mcpan_sv_sel_nold -t ${N} --selection --sites_save

cut -f1 mcpan_SNP_ld.beagle | cut -f1,2 -d "_" > snp_chr
cut -f1 mcpan_SNP_ld.beagle | cut -f3 -d "_" > snp_pos
cut -f1 mcpan_SNP_ld.beagle > snp_var

cut -f1 mcpan_SV_ld.beagle | cut -f1,2 -d "_" > sv_chr
cut -f1 mcpan_SV_ld.beagle | cut -f3 -d "_" > sv_pos
cut -f1 mcpan_SV_ld.beagle > sv_var

R
library(qqman)
library(qvalue)
D <- as.matrix(read.table("mcpan_snp_sel.selection"))
chr <- as.matrix(read.table("snp_chr"))
pos <- as.matrix(read.table("snp_pos"))
var <- as.matrix(read.table("snp_var"))
sites <- as.matrix(read.table("mcpan_snp_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:749) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
png("mcpan_snp_manh.png", res=300, units="in", width=25, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

D <- as.matrix(read.table("mcpan_sv_sel.selection"))
chr <- as.matrix(read.table("sv_chr"))
pos <- as.matrix(read.table("sv_pos"))
var <- as.matrix(read.table("sv_var"))
sites <- as.matrix(read.table("mcpan_sv_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:95) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
svg("mcpan_sv_manh.svg", width=5, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()
quit()

### Fst outlier
vcftools --vcf mcpan_SNP_ld.recode.vcf --out mcpan_SNP_ld.filtered --max-missing 0.88 --max-alleles 2 --recode --recode-INFO-all # filter out non-biallelic loci and loci with missing data proportion > 3/25*100
bcftools annotate --set-id '%CHROM:%POS' -Ov -o mcpan_SNP_ld.filtered.setID.vcf --threads ${N} mcpan_SNP_ld.filtered.recode.vcf  
bcftools annotate --set-id '%CHROM:%POS' -Ov -o mcpan_SV_ld.setID.vcf --threads ${N} mcpan_SV_ld.recode.vcf

R
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)
obj.vcfR <- read.vcfR("mcpan_SNP_ld.filtered.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "snp_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "snp_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK # need to redefine functions as suggested in: https://github.com/whitlock/OutFLANK/issues/19
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "snp_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "snp_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)


obj.vcfR <- read.vcfR("mcpan_SV_ld.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "sv_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "sv_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "sv_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)
quit()

cd ../

## VG pangenome
mkdir -p ./vgpan/
cd ./vgpan/
cp ${BASE}/benchmarking/ncbi_vgpan/vgpan_delly_final_sorted.vcf ./
cp ${VGPAN}/called/sv/vg/${PREFIX}2_vgSV.filtered.popped.vcf ./vgpan_vg_final.vcf
cp ${BASE}/benchmarking/ncbi_vgpan/vgpan_SNP_final_sorted.vcf ./

### Formating a merged SV vcf file for ANGSD
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ vgpan_vg_final.vcf > vgpan_vg_final_simpl.vcf
picard UpdateVcfSequenceDictionary I=./vgpan_vg_final_simpl.vcf O=./vgpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./vgpan_vg_final_updt.vcf O=./vgpan_vg_final_sorted.vcf
bcftools +tag2tag vgpan_vg_final_sorted.vcf -- -r --pl-to-gl > vgpan_vg_final_sortedGL.vcf

bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ vgpan_delly_final_sorted.vcf > vgpan_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # but see if GL field exists. If not, do this before picard steps
bcftools +tag2tag vgpan_delly_final_simpl.vcf -- -r --pl-to-gl > vgpan_delly_final_sortedGL.vcf

bcftools query -f '%CHROM\t%POS0\t%END\n' vgpan_delly_final_sortedGL.vcf > vgpan_delly.bed # convert vcf region to bed
bcftools query -f '%CHROM\t%POS0\t%END\n' vgpan_vg_final_sortedGL.vcf > vgpan_vg.bed # convert vcf region to bed

cut -f1,2 ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai > linpan.genome # create genome file for bedtools
bedtools slop -i vgpan_delly.bed -g linpan.genome -b 100 > vgpan_delly_padded.bed # pad 100 bp 
bedtools slop -i vgpan_vg.bed -g linpan.genome -b 100 > vgpan_vg_padded.bed # pad 100 bp
bedtools intersect -a vgpan_delly_padded.bed -b vgpan_vg_padded.bed > sv_overlapped.bed # NC_053450.1     61447612        61447898; overlapped region should not be duplicated in the below step
bedtools intersect -header -a vgpan_delly_final_sortedGL.vcf -b sv_overlapped.bed -v > vgpan_delly_final_excluded.vcf # bedtools intersect -header -a vgpan_delly_final_sortedGL.vcf -b sv_overlapped.bed -v > vgpan_delly_final_excluded.vcf

bcftools view -Oz -o vgpan_delly_final_excluded.vcf.gz --threads ${N} vgpan_delly_final_excluded.vcf
bcftools view -Oz -o vgpan_vg_final_sortedGL.vcf.gz --threads ${N} vgpan_vg_final_sortedGL.vcf # bcftools view -Oz -o vgpan_vg_final_sortedGL.vcf.gz --threads ${N} vgpan_g_final_sortedGL.vcf
bcftools index vgpan_delly_final_excluded.vcf.gz
bcftools index vgpan_vg_final_sortedGL.vcf.gz 
bcftools concat --allow-overlaps -Ov -o vgpan_SV_final.vcf --threads ${N} vgpan_delly_final_excluded.vcf.gz vgpan_vg_final_sortedGL.vcf.gz # merge two SV vcf files # bcftools concat --allow-overlaps -Ov -o vgpan_SV_final.vcf --threads ${N} vgpan_delly_final_excluded.vcf.gz vgpan_vg_final_sortedGL.vcf.gz

grep -v ^\# vgpan_SV_final.vcf | wc -l # 2169
bcftools view --max-alleles 2 -Ov -o vgpan_SV_final_bi.vcf --threads ${N} vgpan_SV_final.vcf # keep bi-allelic loci only
grep -v ^\# vgpan_SV_final_bi.vcf | wc -l #2167

bcftools query -f '%CHROM %POS %REF %ALT\n' vgpan_SV_final_bi.vcf > vgpan_SV_final_bi.variants # keep info about SV & extract it
bcftools query -f '%CHROM\t%POS\n' vgpan_SV_final_bi.vcf > vgpan_SV_final_bi.chrpos
head -n 2 vgpan_SV_final_bi.variants

Rscript ${HOME}/extract_SV_info_graph.r vgpan_SV_final_bi.variants
head -n 2 vgpan_SV_final_bi.variants.bed

python3 ${HOME}/fasta_extract_flanking_regions_jyj.py ${LINPAN}/${PREFIX}_panref2_sorted.fa vgpan_SV_final_bi.chrpos 1 vgpan_SV_final_bi.variants.ref # extract ref allele at position for dummy vcf
grep ^"#" vgpan_SV_final_bi.vcf | grep -v ^\#\#"contig=<ID=" > vgpan_SV_final_bi.variants.header # keep header
grep -v ^"#" vgpan_SV_final_bi.vcf > vgpan_SV_final_bi.withoutheader # keep without header
Rscript ${HOME}/make_dummy_vcf_snp_jyj.r vgpan_SV_final_bi.withoutheader vgpan_SV_final_bi.variants.ref # edit in REF
cat vgpan_SV_final_bi.variants.header vgpan_SV_final_bi.withoutheader.withdummyREFALT > vgpan_SV_final_bi.forangsd.vcf # format vcf for ANGSD

### generate GL or beagle file for each chromosome (modified from Merot 2023)
ls ${VGPAN}/aligned/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o vgpan_SV_final_bi_rehead.forangsd.vcf --threads ${N} vgpan_SV_final_bi.forangsd.vcf
bgzip  vgpan_SV_final_bi_rehead.forangsd.vcf
tabix vgpan_SV_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl vgpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 20 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

### Nucleotide diversity
realSFS -P ${N} out_snp.saf.idx -maxiter 100 -fold 1 > out_snp.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_snp.saf.idx -sfs out_snp.sfs -P ${N} -outname out_snp # calculate thetas for each site

thetaStat do_stat out_snp.thetas.idx -win 50000 -step 10000  -outnames theta_snp.thetasWindow.gz # sliding window estimate
awk '{print $4,$5,$14}' theta_snp.thetasWindow.gz.pestPG > Thetas
sumSites=$(awk 'BEGIN{s=0;}{s=s+$3;}END{print s;}' Thetas) # get number of sites
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas) # get mean Watterson's Theta
meanW=$(awk "BEGIN {print $sumW/$sumSites}") # 0.00955816
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas) # get mean Tajima's Theta
meanN=$(awk "BEGIN {print $sumN/$sumSites}") # 0.00546272

realSFS -P ${N} out_sv.saf.idx -maxiter 100 -fold 1 > out_sv.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_sv.saf.idx -sfs out_sv.sfs -P ${N} -outname out_sv # calculate thetas for each site

### SFS
#R
#norm <- function(x) x/sum(x) # function to normalize 
#sfs <- (scan("out_snp.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("vgpan_snp_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

#sfs <- (scan("out_sv.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("vgpan_sv_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_vgpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_vgpan_PLraw.txt ROH_vgpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i 
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf vgpan_SV_final.vcf --BEAGLE-GL --chr ${i} --out vgpan_SV_${i}
done 

head -n 1 vgpan_SV_NC_053450.1.BEAGLE.GL > vgpan_SV.beagleheader # get the headerhead 

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i
do 
echo "append beagle from ${i}"
tail -n +2 vgpan_SV_${i}.BEAGLE.GL >> vgpan_SV.beagle
done

wc -l vgpan_SV.beagle # should be equal to the nb of variants minus non-biallelic variants in the vcf
cat vgpan_SV_final.vcf | grep -v "#" | wc -l
cat vgpan_SV_final.vcf | grep -v "#" | cut -f5 | grep "," | wc -l
rm vgpan_SV*log
rm vgpan_SV*BEAGLE.GL

cat vgpan_SV.beagle | cut -f 1 | sed 's/:/\t/g'| gzip > sv_pca.pos.gz # prepare a pos file 
sed -i 's/:/_/g' vgpan_SV.beagle
sed -i 's/\[/#/g' vgpan_SV.beagle # replace the [ with # to avoid file reading issue 
sed -i 's/#/_/g' vgpan_SV.beagle # replace the # with _ to avoid file reading issue 

Rscript ${HOME}/normalize_beagle.r vgpan_SV.beagle vgpan_SV.norm.beagle  # normalize the genotype likelihoods as otherwise too small numbers are counted as zero and cannot be divided

cat vgpan_SV.beagleheader vgpan_SV.norm.beagle > vgpan_SV.ready.beagle #add the  proper header again and zip
head vgpan_SV.ready.beagle
gzip vgpan_SV.ready.beagle

zcat snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file
ngsLD --geno snp_pca_subsampled.geno.beagle.gz --pos snp_pca_subsampled.pos.gz --probs --n_ind 25 --n_sites 876525 --max_kb_dist 1000 --n_threads ${N} --out out_snp.ngsld # run "zcat snp_pca_subsampled.pos.gz | wc -l" to get n_sites
cat out_snp.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 4950 156034054 # check by "cat out_snp.ngsld | cut -f1 | grep "NC_053450.1" | cut -f2 -d ":" | sort -n | head" and "cat out_snp.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail"; takes too long -> limit regions
cat out_snp.ngsld | grep "NC_053450.1" > out_snp_1stchr.ngsld
cat out_snp.ngsld | grep "NC_053480.1" > out_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
cat out_snp_1stchr.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5000 1005000
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_snp_1stchr.ngsld > out_snp_1stchr_header.ngsld
echo "${PWD}/out_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, LD decays quickly

zcat vgpan_SV.ready.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz # prepare a geno file 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 2167 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites
cat out_sv.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 808685 155682736 # check by "cat out_sv.ngsld | cut -f1 | grep "NC_053450.1" | cut -f2 -d ":" | sort -n | head" and "cat out_sv.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat decay.header out_sv.ngsld > out_sv_header.ngsld
echo "${PWD}/out_sv_header.ngsld" > decay_sv.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_sv.ld --max_kb_dist 1000 --out fit_SV_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, pruning seems not needed

${APP}/prune_graph/target/release/prune_graph --in out_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_snp_unlinked.ngsld
sed 's/:/_/g' out_snp_unlinked.ngsld > out_snp_unlinked

${APP}/prune_graph/target/release/prune_graph --in out_sv.ngsld --weight-field "column_7" --weight-filter "column_3 <= 25000 && column_7 >= 0.5" --out out_sv_unlinked.ngsld
sed 's/:/_/g' out_sv_unlinked.ngsld > out_sv_unlinked

zcat snp_pca.beagle.gz | head -n 1 > vgpan_SNP_ld.beagleheader
zcat snp_pca.beagle.gz | grep -Fwf out_snp_unlinked > vgpan_SNP_ld.beagle # keep only unlinked loci
cat vgpan_SNP_ld.beagleheader vgpan_SNP_ld.beagle | gzip > vgpan_SNP_ld.beagle.gz
zcat vgpan_SNP_ld.beagle.gz | head
cat out_snp_unlinked | wc -l # check the number of loci if matching
zcat vgpan_SNP_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

zcat vgpan_SV.ready.beagle.gz | head -n 1 > vgpan_SV_ld.beagleheader
zcat vgpan_SV.ready.beagle.gz | grep -Fwf out_sv_unlinked > vgpan_SV_ld.beagle # keep only unlinked loci
cat vgpan_SV_ld.beagleheader vgpan_SV_ld.beagle | gzip > vgpan_SV_ld.beagle.gz
zcat vgpan_SV_ld.beagle.gz | head
cat out_sv_unlinked | wc -l # check the number of loci if matching
zcat vgpan_SV_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

sed 's/:/\t/g' out_snp_unlinked.ngsld > out_snp_unlinked.tab
vcftools --positions out_snp_unlinked.tab --vcf vgpan_SNP_final_sorted.vcf --recode --recode-INFO-all --out vgpan_SNP_ld # generate ld-pruned vcf file for downstrean usage

sed 's/:/\t/g' out_sv_unlinked.ngsld > out_sv_unlinked.tab
vcftools --positions out_sv_unlinked.tab --vcf vgpan_SV_final.vcf --recode --recode-INFO-all --out vgpan_SV_ld # generate ld-pruned vcf file for downstrean usage

pcangsd -b vgpan_SNP_ld.beagle.gz -o vgpan_SNP_pca -t ${N}
pcangsd -b vgpan_SV_ld.beagle.gz -o vgpan_SV_pca -t ${N}

R
library(ggplot2)
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("vgpan_SNP_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.2%, PC2: 4.3% 
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("vgpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.2%)", y = "Principal Component 2 (4.3%)") + theme_classic(base_size=14)
dev.off()

C <- as.matrix(read.table("vgpan_SV_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 9.2%, PC2: 5.2%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("vgpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (9.2%)", y = "Principal Component 2 (5.2%)") + theme_classic(base_size=14)
dev.off()
quit()

### Selection
pcangsd -b vgpan_SNP_ld.beagle.gz -t ${N} -o vgpan_snp_sel --selection --sites_save # selection
pcangsd -b vgpan_SV_ld.beagle.gz -t ${N} -o vgpan_sv_sel --selection --sites_save # selection

cut -f1 vgpan_SNP_ld.beagle | cut -f1,2 -d "_" > snp_chr
cut -f1 vgpan_SNP_ld.beagle | cut -f3 -d "_" > snp_pos
cut -f1 vgpan_SNP_ld.beagle > snp_var

cut -f1 vgpan_SV_ld.beagle | cut -f1,2 -d "_" > sv_chr
cut -f1 vgpan_SV_ld.beagle | cut -f3 -d "_" > sv_pos
cut -f1 vgpan_SV_ld.beagle > sv_var

R
library(qqman)
library(qvalue)
D <- as.matrix(read.table("vgpan_snp_sel.selection"))
chr <- as.matrix(read.table("snp_chr"))
pos <- as.matrix(read.table("snp_pos"))
var <- as.matrix(read.table("snp_var"))
sites <- as.matrix(read.table("vgpan_snp_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:755) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
png("vgpan_snp_manh.png", res=300, units="in", width=25, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

D <- as.matrix(read.table("vgpan_sv_sel.selection"))
chr <- as.matrix(read.table("sv_chr"))
pos <- as.matrix(read.table("sv_pos"))
var <- as.matrix(read.table("sv_var"))
sites <- as.matrix(read.table("vgpan_sv_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:79) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
svg("vgpan_sv_manh.svg", width=5, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()
quit()

### Fst outlier
vcftools --vcf vgpan_SNP_ld.recode.vcf --out vgpan_SNP_ld.filtered --max-missing 0.88 --max-alleles 2 --recode --recode-INFO-all # filter out non-biallelic loci and loci with missing data proportion > 3/25*100
bcftools annotate --set-id '%CHROM:%POS' -Ov -o vgpan_SNP_ld.filtered.setID.vcf --threads ${N} vgpan_SNP_ld.filtered.recode.vcf  
bcftools annotate --set-id '%CHROM:%POS' -Ov -o vgpan_SV_ld.setID.vcf --threads ${N} vgpan_SV_ld.recode.vcf

R
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)
obj.vcfR <- read.vcfR("vgpan_SNP_ld.filtered.setID.vcf") 
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "snp_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "snp_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK # need to redefine functions as suggested in: https://github.com/whitlock/OutFLANK/issues/19
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # 2 outliers (NC_053451.1:80943444, NC_053453.1:12223757)

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "snp_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "snp_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)


obj.vcfR <- read.vcfR("vgpan_SV_ld.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "sv_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "sv_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos) 
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp) 
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "sv_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)
quit()

cd ../

## Original pangenome
mkdir -p ./orgpan/
cd ./orgpan/
cp ${BASE}/benchmarking/orgpan_mcpan/orgpan_SNP_final_sorted.vcf ./
cp ${BASE}/benchmarking/orgpan_mcpan/orgpan_delly_final_sorted.vcf ./
cp ${ORGPAN}/called/sv/vg/barnswallow_orgSV.filtered.popped.vcf ./orgpan_vg_final.vcf
bcftools annotate --threads ${N} --rename-chrs ${BASE}/benchmarking/orgpan_mcpan/chr_name_conv.txt orgpan_vg_final.vcf -Ov -o orgpan_vg_final_annot.vcf

### Formating a merged SV vcf file for ANGSD
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ orgpan_vg_final.vcf > orgpan_vg_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format
picard UpdateVcfSequenceDictionary I=./orgpan_vg_final_simpl.vcf O=./orgpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./orgpan_vg_final_updt.vcf O=./orgpan_vg_final_sorted.vcf
bcftools +tag2tag orgpan_vg_final_sorted.vcf -- -r --pl-to-gl > orgpan_vg_final_sortedGL.vcf

bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ orgpan_delly_final_sorted.vcf > orgpan_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format
bcftools +tag2tag orgpan_delly_final_simpl.vcf -- -r --pl-to-gl > orgpan_delly_final_sortedGL.vcf

bcftools query -f '%CHROM\t%POS0\t%END\n' orgpan_delly_final_sortedGL.vcf > orgpan_delly.bed # convert vcf region to bed
bcftools query -f '%CHROM\t%POS0\t%END\n' orgpan_vg_final_sortedGL.vcf > orgpan_vg.bed # convert vcf region to bed

cut -f1,2 ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai > linpan.genome # create genome file for bedtools
bedtools slop -i orgpan_delly.bed -g linpan.genome -b 100 > orgpan_delly_padded.bed # pad 100 bp 
bedtools slop -i orgpan_vg.bed -g linpan.genome -b 100 > orgpan_vg_padded.bed # pad 100 bp
bedtools intersect -a orgpan_delly_padded.bed -b orgpan_vg_padded.bed > sv_overlapped.bed # nothing. If there is, overlapped region should not be duplicated in the below step
bedtools intersect -header -a orgpan_delly_final_sortedGL.vcf -b sv_overlapped.bed -v > orgpan_delly_final_excluded.vcf # just keep this nonnecessary code for reproducibility

bcftools view -Oz -o orgpan_delly_final_excluded.vcf.gz --threads ${N} orgpan_delly_final_excluded.vcf
bcftools view -Oz -o orgpan_vg_final_sortedGL.vcf.gz --threads ${N} orgpan_vg_final_sortedGL.vcf
bcftools index orgpan_delly_final_excluded.vcf.gz
bcftools index orgpan_vg_final_sortedGL.vcf.gz
bcftools concat --allow-overlaps -Ov -o orgpan_SV_final.vcf --threads ${N} orgpan_delly_final_excluded.vcf.gz orgpan_vg_final_sortedGL.vcf.gz # merge two SV vcf files

grep -v ^\# orgpan_SV_final.vcf | wc -l # 2271
bcftools view --max-alleles 2 -Ov -o orgpan_SV_final_bi.vcf --threads ${N} orgpan_SV_final.vcf # keep bi-allelic loci only
grep -v ^\# orgpan_SV_final_bi.vcf | wc -l #2266

bcftools query -f '%CHROM %POS %REF %ALT\n' orgpan_SV_final_bi.vcf > orgpan_SV_final_bi.variants # keep info about SV & extract it
bcftools query -f '%CHROM\t%POS\n' orgpan_SV_final_bi.vcf > orgpan_SV_final_bi.chrpos
head -n 2 orgpan_SV_final_bi.variants

Rscript ${HOME}/extract_SV_info_graph.r orgpan_SV_final_bi.variants
head -n 2 orgpan_SV_final_bi.variants.bed

python3 ${HOME}/fasta_extract_flanking_regions_jyj.py ${LINPAN}/${PREFIX}_panref2_sorted.fa orgpan_SV_final_bi.chrpos 1 orgpan_SV_final_bi.variants.ref # extract ref allele at position for dummy vcf
grep ^"#" orgpan_SV_final_bi.vcf | grep -v ^\#\#"contig=<ID=" > orgpan_SV_final_bi.variants.header # keep header
grep -v ^"#" orgpan_SV_final_bi.vcf > orgpan_SV_final_bi.withoutheader # keep without header
Rscript ${HOME}/make_dummy_vcf_snp_jyj.r orgpan_SV_final_bi.withoutheader orgpan_SV_final_bi.variants.ref # edit in REF
cat orgpan_SV_final_bi.variants.header orgpan_SV_final_bi.withoutheader.withdummyREFALT > orgpan_SV_final_bi.forangsd.vcf # format vcf for ANGSD

### generate GL or beagle file for each chromosome (modified from Merot 2023)
cd ${ORGPAN}/aligned/

while read -r line; do
    pair=($line)
    echo "s/${pair[0]}/${pair[1]}/" >> chr_replace.sed
done < ${BENCHMARK}/orgpan_mcpan/chr_name_conv.txt # make a sed script from the list of "search-replace" pairs
sort -n chr_replace.sed > chr_replace.sorted.sed

for i in `cat ${BASE}/original/SRR_Acc_List.txt`; do
  samtools view -h ${i}.barnswallow-pg_surject.sorted.marked.bam | sed 's/refp#0#//g' > ${i}.barnswallow-pg_surject.sorted.marked.replaced.sam
  samtools view -H ${i}.barnswallow-pg_surject.sorted.marked.replaced.sam > temp_header
  sed -f chr_replace.sorted.sed temp_header > modified_header
  samtools view -bS ${i}.barnswallow-pg_surject.sorted.marked.replaced.sam > ${i}.barnswallow-pg_surject.sorted.marked.replaced.bam
  samtools reheader modified_header ${i}.barnswallow-pg_surject.sorted.marked.replaced.bam > ${i}.barnswallow-pg_surject.sorted.marked.rehead.bam
done

cd ${POPGEN}/orgpan/
ls ${ORGPAN}/aligned/*.marked.rehead.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o orgpan_SV_final_bi_rehead.forangsd.vcf --threads ${N} orgpan_SV_final_bi.forangsd.vcf
bgzip  orgpan_SV_final_bi_rehead.forangsd.vcf
tabix orgpan_SV_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl orgpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 20 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

### Nucleotide diversity
realSFS -P ${N} out_snp.saf.idx -maxiter 100 -fold 1 > out_snp.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_snp.saf.idx -sfs out_snp.sfs -P ${N} -outname out_snp # calculate thetas for each site

thetaStat do_stat out_snp.thetas.idx -win 50000 -step 10000  -outnames theta_snp.thetasWindow.gz # sliding window estimate
awk '{print $4,$5,$14}' theta_snp.thetasWindow.gz.pestPG > Thetas
sumSites=$(awk 'BEGIN{s=0;}{s=s+$3;}END{print s;}' Thetas) # get number of sites
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas) # get mean Watterson's Theta
meanW=$(awk "BEGIN {print $sumW/$sumSites}") # 0.00954522
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas) # get mean Tajima's Theta
meanN=$(awk "BEGIN {print $sumN/$sumSites}") # 0.00544404

realSFS -P ${N} out_sv.saf.idx -maxiter 100 -fold 1 > out_sv.sfs # obtain ML estimate of SFS using the folded realSFS
realSFS saf2theta out_sv.saf.idx -sfs out_sv.sfs -P ${N} -outname out_sv # calculate thetas for each site

### SFS
#R
#norm <- function(x) x/sum(x) # function to normalize 
#sfs <- (scan("out_snp.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("orgpan_snp_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

#sfs <- (scan("out_sv.sfs")[2:26]) # read data
#sfs<-norm(sfs)  # the variable categories of the sfs 
#svg("orgpan_sv_sfs_plot.svg", width=5, height=5)
#barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
#dev.off()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_mcpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_mcpan_PLraw.txt ROH_mcpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca

#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand mcpan_SNP_final_sorted_GL.vcf.gz -Oz -o mcpan_SNP_final_sorted_GL_ld.vcf.gz # prune linkage disequilibrium
#bcftools +prune -w 100kb -m 0.5 -n 1 -N rand mcpan_SV_final.vcf -Oz -o mcpan_SV_final_ld.vcf.gz 

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i 
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf orgpan_SV_final.vcf --BEAGLE-GL --chr ${i} --out orgpan_SV_${i}
done 

head -n 1 orgpan_SV_NC_053450.1.BEAGLE.GL > orgpan_SV.beagleheader # get the headerhead 

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i
do 
echo "append beagle from ${i}"
tail -n +2 orgpan_SV_${i}.BEAGLE.GL >> orgpan_SV.beagle
done

wc -l orgpan_SV.beagle # should be equal to the nb of variants minus non-biallelic variants in the vcf
cat orgpan_SV_final.vcf | grep -v "#" | wc -l
cat orgpan_SV_final.vcf | grep -v "#" | cut -f5 | grep "," | wc -l
rm orgpan_SV*log
rm orgpan_SV*BEAGLE.GL

cat orgpan_SV.beagle | cut -f 1 | sed 's/:/\t/g'| gzip > sv_pca.pos.gz # prepare a pos file 
sed -i 's/:/_/g' orgpan_SV.beagle
sed -i 's/\[/#/g' orgpan_SV.beagle # replace the [ with # to avoid file reading issue 
sed -i 's/#/_/g' orgpan_SV.beagle # replace the # with _ to avoid file reading issue 

Rscript ${HOME}/normalize_beagle.r orgpan_SV.beagle orgpan_SV.norm.beagle  # normalize the genotype likelihoods as otherwise too small numbers are counted as zero and cannot be divided

cat orgpan_SV.beagleheader orgpan_SV.norm.beagle > orgpan_SV.ready.beagle #add the  proper header again and zip
head orgpan_SV.ready.beagle
gzip orgpan_SV.ready.beagle

zcat snp_pca.beagle.gz | tail -n +2 | awk 'NR % 10 == 0' | cut -f 4- | gzip  > snp_pca_subsampled.geno.beagle.gz # prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file
zcat snp_pca.mafs.gz | tail -n +2 | cut -f 1,2 | awk 'NR % 10 == 0' | sed 's/:/_/g'| gzip > snp_pca_subsampled.pos.gz # prepare a pos file by subsampling one SNP in every 10 SNPs in the beagle file 
ngsLD --geno snp_pca_subsampled.geno.beagle.gz --pos snp_pca_subsampled.pos.gz --probs --n_ind 25 --n_sites 868674 --max_kb_dist 1000 --n_threads ${N} --out out_snp.ngsld # run "zcat snp_pca_subsampled.pos.gz | wc -l" to get n_sites
cat out_snp.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5695 156035288 # check by "cat out_snp.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_snp.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat out_snp.ngsld | grep "NC_053450.1" > out_snp_1stchr.ngsld
#cat out_snp.ngsld | grep "NC_053480.1" > out_snp_1mbchr.ngsld # this is the chromosome with its length of 1Mb
cat out_snp_1stchr.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 5000 1005000
echo -e "pos1\tpos2\tdist\tr2_ExpG\tD\tDp\tr2" > decay.header
cat decay.header out_snp_1stchr.ngsld > out_snp_1stchr_header.ngsld
echo "${PWD}/out_snp_1stchr_header.ngsld" > decay_snp.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_snp.ld --max_kb_dist 1000 --out fit_SNP_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, LD decays quickly

zcat orgpan_SV.ready.beagle.gz | tail -n +2 | cut -f 4- | gzip  > sv_pca.geno.beagle.gz # prepare a geno file 
ngsLD --geno sv_pca.geno.beagle.gz --pos sv_pca.pos.gz --probs --n_ind 25 --n_sites 2266 --max_kb_dist 0 --n_threads ${N} --out out_sv.ngsld # run "zcat sv_pca.pos.gz | wc -l" to get n_sites
cat out_sv.ngsld | bash ${HOME}/LD_blocks.sh NC_053450.1 808685 155682736 # check by "cat out_sv.ngsld | grep "NC_053450.1" | cut -f1 | cut -f2 -d ":" | sort -n | head" and "cat out_sv.ngsld | grep "NC_053450.1" | cut -f2 | cut -f2 -d ":" | sort -n | tail" 
cat decay.header out_sv.ngsld > out_sv_header.ngsld
echo "${PWD}/out_sv_header.ngsld" > decay_sv.ld
Rscript --vanilla --slave ${HOME}/fit_LDdecay.R --ld_files decay_sv.ld --max_kb_dist 1000 --out fit_SV_LDdecay.pdf --n_ind 25 --fit_level 100 --plot_x_lim 100 --plot_size 5,6 --plot_axis_scale free
# check the r2 and LD decay pdf above - but in this case, pruning seems not needed

${APP}/prune_graph/target/release/prune_graph --in out_snp.ngsld --weight-field "column_7" --weight-filter "column_3 <= 100000 && column_7 >= 0.5" --out out_snp_unlinked.ngsld
sed 's/:/_/g' out_snp_unlinked.ngsld > out_snp_unlinked

${APP}/prune_graph/target/release/prune_graph --in out_sv.ngsld --weight-field "column_7" --weight-filter "column_3 <= 25000 && column_7 >= 0.5" --out out_sv_unlinked.ngsld
sed 's/:/_/g' out_sv_unlinked.ngsld > out_sv_unlinked

zcat snp_pca.beagle.gz | head -n 1 > orgpan_SNP_ld.beagleheader
zcat snp_pca.beagle.gz | grep -Fwf out_snp_unlinked > orgpan_SNP_ld.beagle # keep only unlinked loci
cat orgpan_SNP_ld.beagleheader orgpan_SNP_ld.beagle | gzip > orgpan_SNP_ld.beagle.gz
zcat orgpan_SNP_ld.beagle.gz | head
cat out_snp_unlinked | wc -l # check the number of loci if matching
zcat orgpan_SNP_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

zcat orgpan_SV.ready.beagle.gz | head -n 1 > orgpan_SV_ld.beagleheader
zcat orgpan_SV.ready.beagle.gz | grep -Fwf out_sv_unlinked > orgpan_SV_ld.beagle # keep only unlinked loci
cat orgpan_SV_ld.beagleheader orgpan_SV_ld.beagle | gzip > orgpan_SV_ld.beagle.gz
zcat orgpan_SV_ld.beagle.gz | head
cat out_sv_unlinked | wc -l # check the number of loci if matching
zcat orgpan_SV_ld.beagle.gz | tail -n +2 | wc -l # check the number of loci if matching

sed 's/:/\t/g' out_snp_unlinked.ngsld > out_snp_unlinked.tab
vcftools --positions out_snp_unlinked.tab --vcf orgpan_SNP_final_sorted.vcf --recode --recode-INFO-all --out orgpan_SNP_ld # generate ld-pruned vcf file for downstrean usage

sed 's/:/\t/g' out_sv_unlinked.ngsld > out_sv_unlinked.tab
vcftools --positions out_sv_unlinked.tab --vcf orgpan_SV_final.vcf --recode --recode-INFO-all --out orgpan_SV_ld # generate ld-pruned vcf file for downstrean usage

pcangsd -b orgpan_SNP_ld.beagle.gz -o orgpan_snp_pca -t ${N}
pcangsd -b orgpan_SV_ld.beagle.gz -o orgpan_sv_pca -t ${N}

R
library(ggplot2)
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("orgpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.2%, PC2: 4.4%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("orgpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.2%)", y = "Principal Component 2 (4.4%)") + theme_classic(base_size=14)
dev.off()

C <- as.matrix(read.table("orgpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.9%, PC2: 5.1%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("orgpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.9%)", y = "Principal Component 2 (5.1%)") + theme_classic(base_size=14)
dev.off()
quit()

### Selection
pcangsd -b orgpan_SNP_ld.beagle.gz -t ${N} -o orgpan_snp_sel --selection --sites_save # selection
pcangsd -b orgpan_SV_ld.beagle.gz -t ${N} -o orgpan_sv_sel --selection --sites_save # selection

cut -f1 orgpan_SNP_ld.beagle | cut -f1,2 -d "_" > snp_chr
cut -f1 orgpan_SNP_ld.beagle | cut -f3 -d "_" > snp_pos
cut -f1 orgpan_SNP_ld.beagle > snp_var

cut -f1 orgpan_SV_ld.beagle | cut -f1,2 -d "_" > sv_chr
cut -f1 orgpan_SV_ld.beagle | cut -f3 -d "_" > sv_pos
cut -f1 orgpan_SV_ld.beagle > sv_var

R
library(qqman)
library(qvalue)
D <- as.matrix(read.table("orgpan_snp_sel.selection"))
chr <- as.matrix(read.table("snp_chr"))
pos <- as.matrix(read.table("snp_pos"))
var <- as.matrix(read.table("snp_var"))
sites <- as.matrix(read.table("orgpan_snp_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:405) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
png("orgpan_snp_manh.png", res=300, units="in", width=25, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()

D <- as.matrix(read.table("orgpan_sv_sel.selection"))
chr <- as.matrix(read.table("sv_chr"))
pos <- as.matrix(read.table("sv_pos"))
var <- as.matrix(read.table("sv_var"))
sites <- as.matrix(read.table("orgpan_sv_sel.sites"))
chr_pos_var_sites <- as.data.frame(cbind(chr,pos,var,sites))
colnames(chr_pos_var_sites) <- c("CHR","POS","VAR","SITES")
used <- chr_pos_var_sites[chr_pos_var_sites$SITES==1,]
p <- pchisq(D, 1, lower.tail=FALSE) # obtain p-values from PC-based selection scan
q <- qvalue(p)$qvalues
df <- as.data.frame(cbind(used,p,q))
colnames(df)[5] <- "p"
colnames(df)[6] <- "q"
chr_names <- c(unique(chr))
chr_numbers <- c(1:69) # according to length(chr_names)
df$chr_num <- chr_numbers[match(df$CHR, chr_names)]
df$POS <- as.numeric(df$POS)
alpha <- 0.1
outliers <- which(q < alpha)
length(outliers) # 0
svg("orgpan_sv_manh.svg", width=5, height=5)
manhattan(df, chr="chr_num", bp="POS", snp="VAR", p="q", col=c("blue4", "orange3"), ylab=expression(-log[10](italic(q))), chrlabs = c(1:length(unique(chr))), suggestiveline = -log10(0.1))
dev.off()
quit()

### Fst outlier
vcftools --vcf orgpan_SNP_ld.recode.vcf --out orgpan_SNP_ld.filtered --max-missing 0.88 --max-alleles 2 --recode --recode-INFO-all # filter out non-biallelic loci and loci with missing data proportion > 3/25*100
bcftools annotate --set-id '%CHROM:%POS' -Ov -o orgpan_SNP_ld.filtered.setID.vcf --threads ${N} orgpan_SNP_ld.filtered.recode.vcf   
bcftools annotate --set-id '%CHROM:%POS' -Ov -o orgpan_SV_ld.setID.vcf --threads ${N} orgpan_SV_ld.recode.vcf

R
library(OutFLANK)
library(vcfR)
library(ggplot2)
library(dplyr)
obj.vcfR <- read.vcfR("orgpan_SNP_ld.filtered.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "snp_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "snp_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK # need to redefine functions as suggested in: https://github.com/whitlock/OutFLANK/issues/19
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # 1 outlier NC_053459.1:14512222

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "snp_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "snp_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)


obj.vcfR <- read.vcfR("orgpan_SV_ld.setID.vcf")
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) 
chr_pos$position<-as.numeric(as.character(chr_pos$position))
write.table(chr_pos, "sv_ld_snp_pos.txt", sep="\t", quote=F, row.names=F)
geno <- extract.gt(obj.vcfR) # extract and format genotype matrix
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno)) # an empty matrix, (9 stands for missing data)
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dim(G)
G[1:10,1:10] # overview of our data and its first 10 rows/10 columns
write.table(G, "sv_ld_geno_matrix.txt", sep="\t", col.names=F, row.names=F) # will be useful later

info_samples <- read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=T)
head(info_samples)
pop_vector <- info_samples$POP
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector) # FST matrix with OutFLANK
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=25, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim$results$pvaluesRightTail)
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1)

P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos)
P1_pos[P1_pos$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers 
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos, "sv_outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)

pop_comp <- info_samples$ID[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia']
geno_comp <- geno[, colnames(geno) %in% pop_comp]
G_comp <- matrix(9, nrow = nrow(geno_comp), ncol = ncol(geno_comp)) # an empty matrix, (9 stands for missing data)
G_comp[geno_comp %in% c("0/0", "0|0")] <- 0
G_comp[geno_comp  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_comp[geno_comp %in% c("1/1", "1|1")] <- 2
table(as.vector(G_comp))
dim(G_comp)
G_comp[1:10,1:10] # overview of our data and its first 10 rows/10 columns
info_comp <- info_samples[info_samples$POP == 'Morocco' | info_samples$POP == 'Russia',]
head(info_comp)
comp_vector <- info_comp$POP
my_fst_comp <- MakeDiploidFSTMat(t(G_comp), locusNames = id_snp, popNames = comp_vector) # FST matrix with OutFLANK
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=17, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
str(out_trim_comp)
OutFLANKResultsPlotter(out_trim_comp, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) # have a look at the results
hist(out_trim_comp$results$pvaluesRightTail)
P1_comp <- pOutlierFinderChiSqNoCorr(my_fst_comp, Fstbar = out_trim_comp$FSTNoCorrbar, dfInferred = out_trim_comp$dfInferred, qthreshold = 0.05, Hmin=0.1) # run OutFLANK on all SNPs of the data, corrected by the trim dataset
head(P1_comp)

P1_pos_comp<-left_join(P1_comp, chr_pos, by=c("LocusName"="id_snp"))
dim(P1_pos_comp)
P1_pos_comp[P1_pos_comp$OutlierFlag==TRUE,] # no outlier

ggplot(P1_pos_comp, aes(x=He, y=FST, colour=OutlierFlag)) + geom_point() + theme_classic() # look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
ggplot(P1_pos_comp, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ geom_point() + theme_classic() + facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x") + labs(x = "position (in MB)") 
write.table(P1_pos_comp, "sv_outflank_fst_outliers_MorRus.txt", sep="\t", row.names=F, quote=F)
quit()

cd ../