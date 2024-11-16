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

# Downloading published sequencing data from the GenBank
## short reads
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

## HiFi reads
mkdir -p ./hifi/raw/
cd ./hifi/raw/

cat ../../hifi_SRR_Acc_List.txt | while read g # need the "hifi_SRR_Acc_List.txt" file beforehand that includes accession numbers of SRA file to download (SRR22588214-SRR22588218)
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

### Decontaminating adapters on HiFi raw reads.
mkdir -p ./cleaned/
cd ./cleaned/

cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588214_trimmed.ccs.fq SRR22588214.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588215_trimmed.ccs.fq SRR22588215.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588216_trimmed.ccs.fq SRR22588216.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588217_trimmed.ccs.fq SRR22588217.ccs.fq -j 64 --revcomp -e 0.01
cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" --discard-trimmed -o SRR22588218_trimmed.ccs.fq SRR22588218.ccs.fq -j 64 --revcomp -e 0.01

cd ../

### Assembly with Hifiasm.
mkdir -p ./assembled/
cd ./assembled/

hifiasm --primary -o SRR22588214_trimmed.ccs.asm -t 64 SRR22588214_trimmed.ccs.fq
hifiasm --primary -o SRR22588215_trimmed.ccs.asm -t 64 SRR22588215_trimmed.ccs.fq
hifiasm --primary -o SRR22588216_trimmed.ccs.asm -t 64 SRR22588216_trimmed.ccs.fq
hifiasm --primary -o SRR22588217_trimmed.ccs.asm -t 64 SRR22588217_trimmed.ccs.fq
hifiasm --primary -o SRR22588218_trimmed.ccs.asm -t 64 SRR22588218_trimmed.ccs.fq

### Separating primary and alternate assemblies.

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

### Running merqury to generate the Meryl database with the "_submit_build.sh" script included in the package modified according to your cluster. The script can be found here https://github.com/marbl/merqury/blob/master/_submit_build.sh .

bash _submit_build_mod.sh 21 input_SRR22588214.fofn SRR22588214 
bash _submit_build_mod.sh 21 input_SRR22588215.fofn SRR22588215 
bash _submit_build_mod.sh 21 input_SRR22588216.fofn SRR22588216 
bash _submit_build_mod.sh 21 input_SRR22588217.fofn SRR22588217 
bash _submit_build_mod.sh 21 input_SRR22588218.fofn SRR22588218 

cd ../

### Running Genomescope2.0 online (http://qb.cshl.edu/genomescope/genomescope2.0/) with 31 k-mers, uploading the .histo file outputted by the previous script.

### Running purge dups with custom cutoffs for HiFi reads. The scripts can be found here https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups
cd ./assembled/
mkdir -p ./purged/
cd ./purged/

#### Adding the -xasm20 option for HiFi reads to "minimap2.sh" and "minimap2_self.sh" that can be found in the above link (see below also)

#### Calculating custom cutoffs starting from the kcov computed by Genomescope2.0.

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

#### Modifying the script purge_dups.sh with the custom cutoffs; add the modified "minimap2.sh", "minimap2_self.sh" and "purge_dups.sh" scripts into "_submit_purge_dups.sh" and run it
# in the "purge_dups.sh" script, replace to the following parameters when running SRR22588214: calcuts -m 15.15 -u 45.45 PB.stat > cutoffs 2>calcults.log
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

# in the "purge_dups.sh" script, replace to the following parameters when running SRR22588215: calcuts -m 14.715 -u 44.145 PB.stat > cutoffs 2>calcults.log 
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

# in the "purge_dups.sh" script, replace to the following parameters when running SRR22588216: calcuts -m 18.45 -u 55.35 PB.stat > cutoffs 2>calcults.log   
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

# in the "purge_dups.sh" script, replace to the following parameters when running SRR22588217: calcuts -m 11.205 -u 33.615 PB.stat > cutoffs 2>calcults.log  
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

# in the "purge_dups.sh" script, replace to the following parameters when running SRR22588218: calcuts -m 24.75 -u 74.25 PB.stat > cutoffs 2>calcults.log  
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

# Running "asm_stats.sh" to obtain statistics before and after purging 
${APP}/vgp/stats/asm_stats.sh SRR22588214_primary.fasta 1079274132 c # the genome size is the mean predicted size from Genomescope2.0
${APP}/vgp/stats/asm_stats.sh SRR22588214_primary_purged.fasta 1079274132 c

cd ${BASE}

## Simulating short read data from the published pangenome - do this when working with simulated reads only
module load biocontainers
module load bbmap

mkdir -p ./sim/raw/temp/
cd ./sim/raw/

## Simulating short reads
NGSNGS=/depot/fnrdewoody/apps/NGSNGS
PURGED=/scratch/negishi/jeon96/swallow/original/hifi/assembled/purged

for i in SRR22588214 SRR22588215 SRR22588216 SRR22588217 SRR22588218; do
  for j in primary alternate; do
    for k in {1..4}; do
    ${NGSNGS}/ngsngs -i ${PURGED}/${i}/${j}/${i}_${j}_purged.fasta -c 5 -ld Norm,500,50 -seq paired-end -f fq.gz -o sim${k}_${i}_${j} -qs 40 -cl 150 -s ${k} -t ${N} -t2 ${N}
    done 
  done
done # 4 simulated samples per real sample with 5x coverage for both primary and alternate (20 simulated samples in total; 4 * 5 * 2 = 40x coverage simulated in total from the samples used to build the published pangenome)

done
