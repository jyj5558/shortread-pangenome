#!/bin/bash
#SBATCH --job-name=variantcall
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

###-----Define frequently used variables-----###
N=128 # number of cores
BASE=/scratch/negishi/allen715/pangenome/forNA2 # e.g., /scratch/negishi/jeon96/swallow
APP=/home/allen715/vgtoolkit_1.53 # e.g., /home/jeon96/app
PREFIX=short-swallow-out # please define PREFIX as whatever you which
ORGPAN=/scratch/negishi/allen715/pangenome/forNA2/out # e.g., /scratch/negishi/jeon96/swallow/mcpan
REF=/scratch/negishi/allen715/pangenome/forNA2/ref # e.g., scratch/negishi/jeon96/swallow/original/ref 
CLEANED_SRA=/scratch/negishi/allen715/pangenome/forNA2/cleaned_sra # e.g., ${BASE}/original/sra/cleaned
DEPOT=/depot/fnrdewoody
RefInd=refp 
REF=refp
OUT=/scratch/negishi/allen715/pangenome/forNA2/out
REFDIR=/scratch/negishi/allen715/pangenome/forNA2/ref/
GENOME=GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed

# Please note that this script is a stitched version of Natalie M. Allen's individual scripts (modified based on Secomandi et al. (2023) in Cell Reports). Visit https://github.com/nataliemallen/Pangenome for details.
# Download reference (primary and alternate) and all SRR files listed in "seqfile.txt" (can be found on Natalie's github)

# rename_ref.sh: renames fasta headings of reference files. primary should be "chr", alt should be "contig" (must change in the script below for each reference file)
cd /scratch/negishi/allen715/pangenome/forNA2/ref/

# Define input file
genome_file="GCF_015227805.2_bHirRus1.pri.v3_genomic.fna" 
# Define output file
genome_output="GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed.fasta" 

# Function to rename chromosomes in a genome file
rename_chromosomes() {
    input_file=$1
    output_file=$2

    # Use awk to rename chromosomes
    awk '/^>/{print ">chr" ++i; next}{print}' < "$input_file" > "$output_file"
}

# Rename chromosomes in the genome
rename_chromosomes "$genome_file" "$genome_output"

echo "Chromosome renaming completed. Check $genome_output"

# Define input file
genome_file="GCA_015227815.3_bHirRus1.alt.v3_genomic.fna" 
# Define output file
genome_output="GCA_015227815.3_bHirRus1.alt.v3_genomic_renamed.fasta" 

# Function to rename chromosomes in a genome file
rename_chromosomes() {
    input_file=$1
    output_file=$2

    # Use awk to rename chromosomes
    awk '/^>/{print ">contig" ++i; next}{print}' < "$input_file" > "$output_file"
}

# Rename chromosomes in the genome
rename_chromosomes "$genome_file" "$genome_output"

echo "Chromosome renaming completed. Check $genome_output"


# index.sh: indexes the genome file
module load biocontainers
module load samtools

cd /scratch/negishi/allen715/pangenome/forNA2/ref/

samtools faidx GCF_015227805.2_bHirRus1.pri.v3_genomic_renamed.fasta

samtools faidx GCA_015227815.3_bHirRus1.alt.v3_genomic_renamed.fasta


# rename_fasta.sh: renames contigs 
# Array of specific fasta file paths
files=("/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588214/alternate/SRR22588214_alternate_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588214/primary/SRR22588214_primary_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588215/alternate/SRR22588215_alternate_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588215/primary/SRR22588215_primary_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588216/alternate/SRR22588216_alternate_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588216/primary/SRR22588216_primary_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588217/alternate/SRR22588217_alternate_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588217/primary/SRR22588217_primary_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588218/alternate/SRR22588218_alternate_purged.fasta" "/scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588218/primary/SRR22588218_primary_purged.fasta")

# Loop through the array of files
for file_path in "${files[@]}"; do
    # Check if the fasta file exists
    if [ -e "$file_path" ]; then
        # Create a temporary file to store the modified content
        temp_file="${file_path}.tmp"

        # Counter for contig renaming
        count=1

        # Read the original fasta file and modify contig names
        while IFS= read -r line; do
            # Check if the line starts with ">"
            if [[ $line == ">"* ]]; then
                # Rename the contig according to the specified convention
                echo ">A1_pri$count" >> "$temp_file"
                ((count++))
            else
                # Keep other lines unchanged
                echo "$line" >> "$temp_file"
            fi
        done < "$file_path"

        # Replace the original file with the modified content
        mv "$temp_file" "$file_path"

        echo "Contigs in $file_path renamed successfully."
    else
        echo "Fasta file $file_path not found."
    fi
done


# mask.sh: masks SRR genome files
module purge
module load biocontainers
module load repeatmasker
module load bedtools
export PATH=$PATH:/home/allen715/winmasker/ncbi_cxx--25_2_0/GCC850-ReleaseMT64/bin/

#214 primary
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588214/primary/

#windowmasker
windowmasker -mk_counts -in SRR22588214_primary_purged.fasta -out stage1_SRR22588214_primary.counts
windowmasker -ustat stage1_SRR22588214_primary.counts -in SRR22588214_primary_purged.fasta -outfmt fasta -out SRR22588214_primary_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588214_primary_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588214_primary_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588214_primary_purged_winmask.fasta -bed output.bed -fo SRR22588214_primary_purged_masked.fasta

#214 alternate
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588214/alternate/

#windowmasker
windowmasker -mk_counts -in SRR22588214_alternate_purged.fasta -out stage1_SRR22588214_alternate.counts
windowmasker -ustat stage1_SRR22588214_alternate.counts -in SRR22588214_alternate_purged.fasta -outfmt fasta -out SRR22588214_alternate_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588214_alternate_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588214_alternate_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588214_alternate_purged_winmask.fasta -bed output.bed -fo SRR22588214_alternate_purged_masked.fasta

########

#215 primary
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588215/primary/

#windowmasker
windowmasker -mk_counts -in SRR22588215_primary_purged.fasta -out stage1_SRR22588215_primary.counts
windowmasker -ustat stage1_SRR22588215_primary.counts -in SRR22588215_primary_purged.fasta -outfmt fasta -out SRR22588215_primary_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588215_primary_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588215_primary_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588215_primary_purged_winmask.fasta -bed output.bed -fo SRR22588215_primary_purged_masked.fasta

#215 alternate
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588215/alternate/

#windowmasker
windowmasker -mk_counts -in SRR22588215_alternate_purged.fasta -out stage1_SRR22588215_alternate.counts
windowmasker -ustat stage1_SRR22588215_alternate.counts -in SRR22588215_alternate_purged.fasta -outfmt fasta -out SRR22588215_alternate_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588215_alternate_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588215_alternate_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588215_alternate_purged_winmask.fasta -bed output.bed -fo SRR22588215_alternate_purged_masked.fasta

########

#216 primary
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588216/primary/

#windowmasker
windowmasker -mk_counts -in SRR22588216_primary_purged.fasta -out stage1_SRR22588216_primary.counts
windowmasker -ustat stage1_SRR22588216_primary.counts -in SRR22588216_primary_purged.fasta -outfmt fasta -out SRR22588216_primary_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588216_primary_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588216_primary_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588216_primary_purged_winmask.fasta -bed output.bed -fo SRR22588216_primary_purged_masked.fasta

#216 alternate
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588216/alternate/

#windowmasker
windowmasker -mk_counts -in SRR22588216_alternate_purged.fasta -out stage1_SRR22588216_alternate.counts
windowmasker -ustat stage1_SRR22588216_alternate.counts -in SRR22588216_alternate_purged.fasta -outfmt fasta -out SRR22588216_alternate_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588216_alternate_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588216_alternate_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588216_alternate_purged_winmask.fasta -bed output.bed -fo SRR22588216_alternate_purged_masked.fasta

########

#217 primary
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588217/primary/

#windowmasker
windowmasker -mk_counts -in SRR22588217_primary_purged.fasta -out stage1_SRR22588217_primary.counts
windowmasker -ustat stage1_SRR22588217_primary.counts -in SRR22588217_primary_purged.fasta -outfmt fasta -out SRR22588217_primary_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588217_primary_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588217_primary_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588217_primary_purged_winmask.fasta -bed output.bed -fo SRR22588217_primary_purged_masked.fasta

#217 alternate
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588217/alternate/

#windowmasker
windowmasker -mk_counts -in SRR22588217_alternate_purged.fasta -out stage1_SRR22588217_alternate.counts
windowmasker -ustat stage1_SRR22588217_alternate.counts -in SRR22588217_alternate_purged.fasta -outfmt fasta -out SRR22588217_alternate_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588217_alternate_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588217_alternate_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588217_alternate_purged_winmask.fasta -bed output.bed -fo SRR22588217_alternate_purged_masked.fasta

########

#218 primary
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588218/primary/

#windowmasker
windowmasker -mk_counts -in SRR22588218_primary_purged.fasta -out stage1_SRR22588218_primary.counts
windowmasker -ustat stage1_SRR22588218_primary.counts -in SRR22588218_primary_purged.fasta -outfmt fasta -out SRR22588218_primary_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588218_primary_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588218_primary_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588218_primary_purged_winmask.fasta -bed output.bed -fo SRR22588218_primary_purged_masked.fasta

#218 alternate
cd /scratch/negishi/allen715/pangenome/forNA2/purged_hifi/SRR22588218/alternate/

#windowmasker
windowmasker -mk_cSRR22588218ounts -in SRR22588218_alternate_purged.fasta -out stage1_SRR22588218_alternate.counts
windowmasker -ustat stage1_SRR22588218_alternate.counts -in SRR22588218_alternate_purged.fasta -outfmt fasta -out SRR22588218_alternate_purged_winmask.fasta

#repeatmasker
RepeatMasker -pa 32 -xsmall -species aves SRR22588218_alternate_purged.fasta -dir .

awk 'BEGIN{OFS="\t"} {if($1 !~ /^#/) print $5, $6-1, $7, $11, $0}' SRR22588218_alternate_purged.fasta.out | tail -n +4 > output.bed

##bedtools
bedtools maskfasta -soft -fi SRR22588218_alternate_purged_winmask.fasta -bed output.bed -fo SRR22588218_alternate_purged_masked.fasta


# pangenome.sh: constructs the original pangenome
module purge
module load anaconda
conda activate cactus_env
source /home/allen715/cactus/cactus-bin-v2.7.0/venv-cactus-v2.7.0/bin/activate

cactus-minigraph ./jobstore seqfile.txt short-swallow-gfa.gfa --reference $REF --binariesMode local
cactus-graphmap ./jobstore seqfile.txt short-swallow-gfa.gfa short-swallow-paf.paf --outputFasta short-swallow-paf.gfa.fa --reference $REF --binariesMode local --delFilter 10000000 # This option is needed to filter out potentially spurious minigraph-based structural variants based on the developer’s tutorial.
keep_contigs=$(awk '{print $1}' ${REFDIR}/${GENOME}.fasta.fai) #${GENOME}.fna.fai file is the reference genome’s index file.
cactus-graphmap-split ./jobstore seqfile.txt short-swallow-gfa.gfa short-swallow-paf.paf --reference $REF --outDir $OUT --binariesMode local --otherContig contigOther --refContigs $(for i in $keep_contigs; do echo ${i}; done) # These two options are needed to distinguish major contigs from the reference genome and the other ones that are not contained in the reference genome.
cactus-align ./jobstore $OUT/chromfile.txt $OUT/chrom-alignments --batch --pangenome --reference $REF --outVG --maxLen 10000 --binariesMode local
cactus-graphmap-join ./jobstore --vg $OUT/chrom-alignments/*.vg --hal $OUT/chrom-alignments/*.hal --outDir $OUT --outName short-swallow-out --reference $REF --binariesMode local --gbz clip full --gfa --filter 1 --clip 10000 --giraffe clip --vcf --chrom-vg --chrom-og --viz --xg --draw # “--gbz clip full” and “--gfa” options are needed for downstream steps. “--filter 1” and “--clip 10000” options are needed to filter out spurious results (a parameter for “—filter” should be the 10% of the number of total haplotypes in the analysis. So as you are using 12 haplotypes, it should be the integer of “1.2”.). The other output options are not necessary but can be useful for visualization and stat comparison.


# Processing the graph pangenome
cd ${ORGPAN}/

## Converting file type
gunzip ${PREFIX}-pg.gfa.gz 
${APP}/vg convert -g ${PREFIX}-pg.gfa -v > ${PREFIX}-pg.vg

## check graph genome stats; run this step on highmem queue 
# ${APP}/vg paths -M -x ${PREFIX}-pg.vg > paths_table.txt
# ${APP}/vg stats -lz ${PREFIX}-pg.vg
# ${APP}/vg paths -E -v ${PREFIX}-pg.vg > paths_lengths.txt
# ${APP}/vg stats -F ${PREFIX}-pg.vg

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./${PREFIX}-pg/
mkdir -p ./tmp/
 
${APP}/vg mod -t ${N} -O ./${PREFIX}-pg.vg > ${PREFIX}_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}_mod.xg ${PREFIX}_mod.vg
 
## Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}_mod.vg > ${PREFIX}_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}_mod_chopped.xg ${PREFIX}_mod_chopped.vg
 
## Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}_mod_chopped.vg > ${PREFIX}_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}_mod_chopped_pruned.gcsa ${PREFIX}_mod_chopped_pruned.vg
 
## Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}_mod_chopped.vg -x ${PREFIX}_mod_chopped_new_L.xg #-L: preserve alt paths in the xg