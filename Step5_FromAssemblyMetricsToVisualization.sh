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





# SequenceTubeMap Visualization
cd /scratch/negishi/jeon96/swallow/tubemap/
cp ${MCPAN}/${PREFIX}-pg/${PREFIX}-pg.vg
cp ${VGPAN}/augmented/${PREFIX}2-pg.vg ./
cp ${VGPAN}/augmented/called/sim_SNP_SV.normed.vcf.gz ./${PREFIX}2-pg.vcf.gz
cp ${VGPAN}/augmented/called/sim_SNP_SV.normed.vcf.gz.tbi ./${PREFIX}2-pg.vcf.gz.tbi

PATH=${APP}:$PATH # the directory containing the vg executable needs to be added to the environment path

## Generating vg file's index
cd ${APP}/sequenceTubeMap/scripts/ # install SequenceTubeMap beforehand and modify the src/config.json to have your input files' directory as the "dataPath"
./prepare_vg.sh /scratch/negishi/jeon96/swallow/tubemap/${PREFIX}-pg.vg
./prepare_vg.sh /scratch/negishi/jeon96/swallow/tubemap/${PREFIX}2-pg.vg

## Preparing subgraphs in advance
cd /scratch/negishi/jeon96/swallow/tubemap/
${APP}/sequenceTubeMap/scripts/prepare_chunks.sh -x ${PREFIX}-pg.vg.xg -r bHirRus1_LinPan#0#NC_053459.1:17272192-17276215 -d 'CAMK2N2' -o ${PREFIX}_chunk-CAMK2N2 >> ${PREFIX}-CAMK2N2.bed
${APP}/sequenceTubeMap/scripts/prepare_chunks.sh -x ${PREFIX}2-pg.vg.xg -h ${PREFIX}2-pg.vg.gbwt -r NC_053459.1:17272192-17276215 -d 'CAMK2N2' -o ${PREFIX}2_chunk-CAMK2N2 >> ${PREFIX}2-CAMK2N2.bed

## Executing SequenceTubeMap
# open ThinLinc and in its terminal, run "PATH=<your path to the directory containing vg>:$PATH" then "npm run serve" in the 'sequenceTubeMap' folder
# then open a web brower and go to "localhost:3000"
# select "custom" from the Data arrow and choose the bed file wanted

#${APP}/vg index -t ${N} -x ${PREFIX}2-pg.vg.xg ${PREFIX}2-pg.vg
#${APP}/vg chunk -t ${N} -x ${PREFIX}2-pg.vg.xg -p NC_053459.1:17272192-17272292 -c 2 > ${PREFIX}2_CAMK2N2.vg
#${APP}/vg index -t ${N} -x ${PREFIX}2_CAMK2N2.vg.xg ${PREFIX}2_CAMK2N2.vg

#${APP}/vg index -t ${N} -x ${PREFIX}-pg.vg.xg ${PREFIX}-pg.vg
#${APP}/vg chunk -t ${N} -x ${PREFIX}-pg.vg.xg -p bHirRus1_LinPan#0#NC_053459.1:17272192-17272292 -c 2 > ${PREFIX}_CAMK2N2.vg
#${APP}/vg index -t ${N} -x ${PREFIX}_CAMK2N2.vg.xg ${PREFIX}_CAMK2N2.vg

#${APP}/vg convert -t ${N} --gfa-out ${PREFIX}_CAMK2N2.vg > ${PREFIX}_CAMK2N2.gfa
#odgi build -g ${PREFIX}_CAMK2N2.gfa -o ${PREFIX}_CAMK2N2.og
#odgi sort -i ${PREFIX}_CAMK2N2.og --threads ${N} --optimize -P -o ${PREFIX}_CAMK2N2.sorted.og
#odgi viz -i ${PREFIX}_CAMK2N2.sorted.og -o ${PREFIX}_CAMK2N2.png

#${APP}/vg convert -t ${N} --gfa-out ${PREFIX}2_CAMK2N2.vg > ${PREFIX}2_CAMK2N2.gfa 
#odgi build -g ${PREFIX}2_CAMK2N2.gfa -o ${PREFIX}2_CAMK2N2.og
#odgi sort -i ${PREFIX}2_CAMK2N2.og --threads ${N} --optimize -P -Y -o ${PREFIX}2_CAMK2N2.sorted.og
#odgi viz -i ${PREFIX}2_CAMK2N2.sorted.og -o ${PREFIX}2_CAMK2N2.png
