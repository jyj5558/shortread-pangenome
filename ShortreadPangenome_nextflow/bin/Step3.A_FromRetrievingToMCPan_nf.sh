#!/bin/bash

SRA=$1
REF=$2
LINPAN=$3
CONTIG=$4
MCPAN=$5
N=$6
PREFIX=$7
GENOME=$8
RefInd=$9
AltInd=${10}
AltHap=${11}
APP=${12}
CACTUS=${13}

cd ${CONTIG}

# Retrieving each sample's contigs from the sample-derived supercontig collection
module purge
module load biocontainers
module load bowtie2
module load bioawk
module load samtools
module load minimap2

mkdir -p ./retrieved/
cd ./retrieved/

## Mapping each sample's original reads onto sample_contigs2.fa (i.e., sample-derived supercontig collection)
bowtie2-build ../sample_contigs2_dedup.fa ../sample_contigs2

touch all_mappingrate.txt
touch all_depth.txt
touch all_coverage.txt
touch all_covstat.txt

for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  bowtie2 -p ${N} -I 0 -X 1000 -x ../sample_contigs2 -1 ${SRA}/Nremoved/${i}_*1.fq_Ns_removed -2 ${SRA}/Nremoved/${i}_*2.fq_Ns_removed --end-to-end --sensitive -S ${i}_mapped.sam 2> ${i}_bowtie.log
  samtools sort -@ ${N} -o ${i}_mapped_sorted.bam ${i}_mapped.sam
  samtools flagstat -@ ${N} ${i}_mapped_sorted.bam | grep "mapped (" | grep -v "primary" | cut -d "(" -f 2 | cut -d "%" -f 1 >> all_mappingrate.txt
  samtools coverage ${i}_mapped_sorted.bam | sort -k3 -k6 -k7 -k9 -r -g > ${i}_cov_sorted.txt # sort by 6th (breadth) then 7th (depth) 
  cut -f1,3,6,7,9 ${i}_cov_sorted.txt >> all_covstat.txt # keep only length, breadth, depth, then mapq
  samtools depth -@ ${N} -a ${i}_mapped_sorted.bam | awk '{c++;s+=$3}END{print s/c}' >> all_depth.txt
  samtools depth -@ ${N} -a ${i}_mapped_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> all_breadth.txt
done
sort -k2 -k3 -k4 -k5 -r -g all_covstat.txt > all_covstat_sorted.txt # sort by length, breadth, depth, then mapq


# Optional: Check distribution of breadth, depth, and mapq by length in R
# Below is my example

# R
## as total,
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

## among contigs >100kb 
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

## among contigs >10kb  & <=100kb                
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

## among contigs <=10kb                
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

## Saving supercontigs of coverage >80, depth >3, mapq >5 for contigs >100kb; coverage >85, depth >3, mapq >5 for contigs >10kb; coverage >90, depth >3, mapq >5 for contigs <=10kb; 
# the fasta size of the set of supercontig collection is a double of the original linear representative genome in my test case (so depth should be the half of the original 5x and mapq should be around 3 (50% accuracy) when assuming random distribution) 
# If needed, modify the cutoff values accordingly

for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  awk '{if ($3 > 100000 && $6 > 80 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt > ${i}_contigs.stat
  awk '{if ($3 > 10000 && $3 <= 100000 && $6 > 85 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt >> ${i}_contigs.stat
  awk '{if ($3 <= 10000 && $6 > 90 && $7 > 3 && $9 > 5) {print $0}}' ${i}_cov_sorted.txt >> ${i}_contigs.stat
  grep -v "#" ${i}_contigs.stat | awk 'BEGIN {FS="\t"}; {print $1}' > ${i}_contigs.region
  cat ${i}_contigs.region >> all_contigs.region
done
sort all_contigs.region | uniq > all_contigs_uniq.region

samtools faidx ../sample_contigs2_dedup.fa > ../sample_contigs2_dedup.fa.fai
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  samtools faidx -r ${i}_contigs.region ../sample_contigs2_dedup.fa > ${i}_contigs.fa
sone

## Test-mapping the all_contigs onto the linear pangenome (i.e., ${PREFIX}_panref2) in order to see mapping statistics (can help imagining the pangenome structure)
samtools faidx -r all_contigs_uniq.region ../sample_contigs2_dedup.fa > all_contigs.fa

minimap2 -ax asm5 ${LINPAN}/${PREFIX}_panref2_sorted.fa all_contigs.fa > sample_to_panref.sam

samtools sort -@ ${N} -o sample_to_panref.bam sample_to_panref.sam
samtools flagstat -@ ${N} sample_to_panref.bam > sample_to_panref_stat.txt
samtools depth -@ ${N} -a sample_to_panref.bam | awk '{c++;s+=$3}END{print s/c}' > sample_to_panref_depth.txt 
samtools depth -@ ${N} -a sample_to_panref.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > sample_to_panref_breadth.txt
samtools coverage sample_to_panref.bam > sample_to_panref_coverage.txt


# Constructing a Minigraph-Cactus pangenome
mkdir -p ${MCPAN}
cd ${MCPAN}

export MYBUCKET=${MCPAN}
keep_contigs=$(awk '{print $1}' ${REF}/${GENOME}.fna.fai)

## Creating a seqFile containing haplotype file information 
touch ${PREFIX}Pangenome.txt

echo -e "# Reference individual:" >> ${PREFIX}Pangenome.txt 
echo -e "${RefInd}\t${LINPAN}/${PREFIX}_panref2_sorted.fa" >> ${PREFIX}Pangenome.txt

echo -e "\n# Haploid sample:" >> ${PREFIX}Pangenome.txt
echo -e "${AltInd}}\t${REF}/${AltHap}" >> ${PREFIX}Pangenome.txt ## add an alternate haplotype of the reference individual if there is (consider adding ".1" and ".2" for different haplotypes of the same reference diploid individual. Refer to the Minigraph-Cactus manual.)

for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
 echo -e "${i}\t${CONTIG}/retrieved/${i}_contigs.fa" >> ${PREFIX}Pangenome.txt
done

# manually confirm the Pangenome.txt file! Check the manuaml at: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md

module load anaconda
conda activate cactus_env # pre-install cactus in conda environment "cactus_env"
source ~/app/venv-cactus-v2.7.0/bin/activate
export PATH=$PATH:${CACTUS}

## Preprocessing seqFile
mkdir -p Temp/
cp ${PREFIX}Pangenome.txt ${PREFIX}Pangenome.mc.seqfile
cactus-preprocess ./jobstore ${PREFIX}Pangenome.txt ${PREFIX}Pangenome.mc.seqfile --pangenome --workDir ./Temp/ --binariesMode local
echo "preprocessing done"

## Making a SV graph with Minigraph in GFA format 
cactus-minigraph ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}-mc.gfa --reference ${RefInd} --maxCores ${N} --workDir Temp/ --binariesMode local
echo "minigraph done"

## Mapping each input assembly back to the graph (assembly-to-graph alignments) using Minigraph
cactus-graphmap ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}-mc.gfa ${PREFIX}-mc.paf --outputGAFDir ${MYBUCKET}/${PREFIX}-mc-gaf --outputFasta ${PREFIX}-mc.gfa.fa --reference ${RefInd} --nodeStorage 1000 --delFilter 10000000 --workDir Temp/ --maxMemory 1000G --mapCores 16 --maxCores ${N} --maxNodes 25 --binariesMode local
echo "graphmap done"

## Splitting by chromosomes (or scaffolds)
cactus-graphmap-split ./jobstore ${PREFIX}Pangenome.mc.seqfile ${PREFIX}-mc.gfa ${PREFIX}-mc.paf --outDir ${MYBUCKET}/contigs --otherContig contigOther --refContigs $(for i in $keep_contigs; do echo ${i}; done) --reference ${RefInd} --nodeStorage 1000 --workDir Temp/ --maxMemory 1000G --maxCores ${N} --maxNodes 5 --logFile ${PREFIX}-mc.split.log --binariesMode local
echo "split done"

## Creating a Cactus base alignment and a "raw" pangenome graph
cactus-align --batch ./jobstore ./contigs/chromfile.txt ${MYBUCKET}/align --consCores ${N} --maxCores ${N} --maxMemory 1000G --workDir Temp/ --logFile ${PREFIX}-mc.align.log --maxNodes 20 --nodeStorage 1000 --pangenome --maxLen 10000 --reference ${RefInd} --outVG --noLinkImports --noMoveExports --cleanWorkDir=onSuccess --realTimeLogging --binariesMode local
echo "align done"

## Creating and indexing the final pangenome graph and produce a VCF ("--filter 4" filters sequences covered by less than 10% of these 40 haplotypes)
cactus-graphmap-join ./jobstore --vg ./align/*.vg --hal ./align/*.hal --outDir ./${PREFIX}-mc-pg --outName ${PREFIX}-mc-pg --reference ${RefInd} --vcf --gbz clip full --gfa --chrom-vg --giraffe clip --filter 4 --clip 10000 --nodeStorage 1000 --maxNodes 20 --workDir Temp/ --logFile ${PREFIX}-mc.join.log --logDebug --indexCores $((N-1)) --maxCores ${N} --maxMemory 1500G --cleanWorkDir=onSuccess --disableCaching --chrom-og --viz --xg --binariesMode local --draw
echo "join" done

${APP}/vg stats -lz ${PREFIX}-mc-pg.gbz # check pangenome basic statistics

# compare to the length of RefInd in the graph: 1122798926 from the previous Quast output

${APP}/vg stats -lz ${PREFIX}-mc-pg.full.gbz # check how much additional sequence is added without clipping

${APP}/vg stats -lz ${PREFIX}-mc-pg.vg # after converting file type as below


# Augmenting the output graph pangenome (preprocessing steps before alignmening were modified from Sacomandi et al. (2023)'s scripts)
cd ${MCPAN}/

## Converting file type
gunzip ${PREFIX}-mc-pg.gfa.gz 
${APP}/vg convert -g ${PREFIX}-mc-pg.gfa -v > ${PREFIX}-mc-pg.vg

## Modifying the pangenome to flip nodes' strands to reduce the number of times paths change strands
cd ./${PREFIX}-mc-pg/
mkdir -p ./tmp/

${APP}/vg mod -t ${N} -O ../${PREFIX}-mc-pg.vg > ${PREFIX}-mc_mod.vg
${APP}/vg index -p -t ${N} -x ${PREFIX}-mc_mod.xg ${PREFIX}-mc_mod.vg

## Chopping nodes in the graph so they are not more than 256 bp and index the graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}-mc_mod.vg > ${PREFIX}-mc_mod_chopped.vg
${APP}/vg index -t ${N} -x ${PREFIX}-mc_mod_chopped.xg ${PREFIX}-mc_mod_chopped.vg

## Pruning the graph with kmer size 45 and index the graph
${APP}/vg prune -t ${N} -k 45 ${PREFIX}-mc_mod_chopped.vg > ${PREFIX}-mc_mod_chopped_pruned.vg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}-mc_mod_chopped_pruned.gcsa ${PREFIX}-mc_mod_chopped_pruned.vg

## Indexing the graph with -L option
${APP}/vg index -t ${N} -L -b ./tmp ${PREFIX}-mc_mod_chopped.vg -x ${PREFIX}-mc_mod_chopped_new_L.xg #-L: preserve alt paths in the xg

## Aligning the individuals that were used to build the pangenome; do augmentation steps when reads that were used to build the pangenome are different from reads that will be used to call variants like this case. 
for i in `ls -lt ${SRA}/raw | grep -Ei "fastq|fq" | tr -s ' ' | cut -d " " -f 9 | cut -d "R" -f 1 | uniq | cut -d "_" -f 1`; do # this can only be correctly applied when forward and reverse sequences are represented as "R"1 and "R"2
  cat ${SRA}/Nremoved/${i}*_1.fq* >  ${SRA}/Nremoved/${i}_R1.fq.gz
  cat ${SRA}/Nremoved/${i}*_2.fq* >  ${SRA}/Nremoved/${i}_R2.fq.gz
  ${APP}/vg map -t ${N} -f ${SRA}/Nremoved/${i}_R1.fq.gz -f ${SRA}/Nremoved/${i}_R2.fq.gz -x ${PREFIX}-mc_mod_chopped.xg -g ${PREFIX}-mc_mod_chopped_pruned.gcsa > ${i}_aln.gam

## Filtering secondary and ambiguous read mappings out of the GAM of the above step
  ${APP}/vg filter -t ${N} ${i}_aln.gam -r 0.90 -fu -m 1 -q 15 -D 999 -x ${PREFIX}-mc_mod_chopped.xg > ${i}_aln.filtered.gam #-r : minimum score to keep primary alignment; -f: normalize score based on length; -u: use substitution count instead of score; -m: filter reads that don't begin with at least N matches on each end; -q: filter alignments with mapping quality < N; -D: clip back the ends of reads that are ambiguously aligned, up to N bases
done

cat *filtered.gam > combined_filtered.gam 

## Augmenting the graph with all variation from the GAM of the above step
${APP}/vg convert -t ${N} ${PREFIX}-mc_mod_chopped_new_L.xg -p > ${PREFIX}-mc.pg
${APP}/vg augment -t ${N} ${PREFIX}-mc.pg combined_filtered.gam -s -m 3 -q 5 -Q 5 -A ${PREFIX}-mc_aug.gam > ${PREFIX}-mc_aug.pg #-s: safely ignore alignments to nodes outside the graph; -m 3: minimum coverage of 3; -q & -Q 5: filtering out mappings and bases with quality < 5 

## Indexing the augmented graph
${APP}/vg mod -t ${N} -X 256 ${PREFIX}-mc_aug.pg > ${PREFIX}-mc_aug_chopped.pg
${APP}/vg index -t ${N} -x ${PREFIX}-mc_aug_chopped.xg ${PREFIX}-mc_aug_chopped.pg
${APP}/vg prune -t ${N} -k 45 ${PREFIX}-mc_aug_chopped.pg > ${PREFIX}-mc_aug_chopped_pruned.pg
${APP}/vg index -t ${N} -b ./tmp -p -g ${PREFIX}-mc_aug_chopped_pruned.gcsa ${PREFIX}-mc_aug_chopped_pruned.pg

## Indexing the augmented graph with -L option
${APP}/vg index -t ${N} -L -b /tmp ${PREFIX}-mc_aug_chopped.pg -x ${PREFIX}-mc_aug_chopped_new_L.xg #-L: preserve alt paths in the xg
${APP}/vg convert -t ${N} ${PREFIX}-mc_aug_chopped_new_L.xg -p > ${PREFIX}-mc_aug_new.pg
