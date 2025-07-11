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


# Conducting population genetic analyses
module purge
module load biocontainers
module load angsd
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



## Linear reference genome
mkdir -p ./ncbi/
cd ./ncbi/
cp ${BASE}/benchmarking/ncbi_linpan/ncbi_SNP_final.vcf ./
cp ${BASE}/benchmarking/ncbi_linpan/ncbi_delly_final.vcf ./

### Formating SV vcf file for ANGSD (modified based on Mérot et al. (2023) in Molecular Ecology)
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

### Generating GL file or beagle file for each chromosome (modified based on Mérot et al. (2023) in Molecular Ecology)
ls ${REF}/mapped/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${REF}/${GENOME}.fna.fai -o ncbi_delly_final_bi_rehead.forangsd.vcf --threads ${N} ncbi_delly_final_bi.forangsd.vcf
bgzip ncbi_delly_final_bi_rehead.forangsd.vcf 
tabix ncbi_delly_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl ncbi_delly_final_bi_rehead.forangsd.vcf.gz -nind 25 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

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
R
norm <- function(x) x/sum(x) # function to normalize 
sfs <- (scan("out_snp.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("ncbi_snp_sfs_plot.svg", width=5, height=5)
barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()

sfs <- (scan("out_sv.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("ncbi_sv_sfs_plot.svg", width=5, height=5)
barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()
quit()

### ROH
angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_ncbi_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_ncbi_PLraw.txt ROH_ncbi_PLresult.txt ${REF}/${GENOME}.fna.fai # parse raw ROH file

### pairwise Fst
cat ./bam.filelist | grep -E "SRR14127742|SRR14131105|SRR14131116|SRR14131127|SRR14131128" > ./bam.mongolia
cat ./bam.filelist | grep -E "SRR14127744|SRR14131090|SRR14131091" > ./bam.china
cat ./bam.filelist | grep -E "SRR14127745|SRR14127746|SRR14127747|SRR14127748|SRR14127749|SRR14127750|SRR14127751|SRR14127752|SRR14127753" > ./bam.russia
cat ./bam.filelist | grep -E "SRR14131088|SRR14131089|SRR14131092|SRR14131094|SRR14131123|SRR14131124|SRR14131125|SRR14131126" > ./bam.morocco

angsd -bam ./bam.mongolia -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 4 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out mongolia_snp
realSFS mongolia_snp.saf.idx -P ${N} -fold 1 > mongolia_snp.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -bam ./bam.china -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 2 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out china_snp
realSFS china_snp.saf.idx -P ${N} -fold 1 > china_snp.sfs

angsd -bam ./bam.russia -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 7 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out russia_snp
realSFS russia_snp.saf.idx -P ${N} -fold 1 > russia_snp.sfs

angsd -bam ./bam.morocco -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -minQ 30 -minInd 6 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out morocco_snp
realSFS morocco_snp.saf.idx -P ${N} -fold 1 > morocco_snp.sfs

realSFS mongolia_snp.saf.idx china_snp.saf.idx -P ${N} > mongolia.china_snp.ml # calculate the 2D SFS
realSFS mongolia_snp.saf.idx russia_snp.saf.idx -P ${N} > mongolia.russia_snp.ml
realSFS mongolia_snp.saf.idx morocco_snp.saf.idx -P ${N} > mongolia.morocco_snp.ml
realSFS china_snp.saf.idx russia_snp.saf.idx -P ${N} > china.russia_snp.ml
realSFS china_snp.saf.idx morocco_snp.saf.idx -P ${N} > china.morocco_snp.ml
realSFS russia_snp.saf.idx morocco_snp.saf.idx -P ${N} > russia.morocco_snp.ml

realSFS fst index mongolia_snp.saf.idx china_snp.saf.idx -sfs mongolia.china_snp.ml -fstout mongolia.china_snp # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_snp.saf.idx russia_snp.saf.idx -sfs mongolia.russia_snp.ml -fstout mongolia.russia_snp
realSFS fst index mongolia_snp.saf.idx morocco_snp.saf.idx -sfs mongolia.morocco_snp.ml -fstout mongolia.morocco_snp
realSFS fst index china_snp.saf.idx russia_snp.saf.idx -sfs china.russia_snp.ml -fstout china.russia_snp
realSFS fst index china_snp.saf.idx morocco_snp.saf.idx -sfs china.morocco_snp.ml -fstout china.morocco_snp
realSFS fst index russia_snp.saf.idx morocco_snp.saf.idx -sfs russia.morocco_snp.ml -fstout russia.morocco_snp

realSFS fst stats mongolia.china_snp.fst.idx # global pairwise estimates # FST.Unweight[nObs:904619153]:0.071057 Fst.Weight:0.084810
realSFS fst stats mongolia.russia_snp.fst.idx # FST.Unweight[nObs:820113826]:0.029463 Fst.Weight:0.052551
realSFS fst stats mongolia.morocco_snp.fst.idx # FST.Unweight[nObs:877515038]:0.042866 Fst.Weight:0.056386
realSFS fst stats china.russia_snp.fst.idx # FST.Unweight[nObs:832310773]:0.028314 Fst.Weight:0.066978
realSFS fst stats china.morocco_snp.fst.idx # FST.Unweight[nObs:897956109]:0.054698 Fst.Weight:0.070026
realSFS fst stats russia.morocco_snp.fst.idx # FST.Unweight[nObs:820046403]:0.029892 Fst.Weight:0.035464

bcftools view -s SRR14127742,SRR14131105,SRR14131116,SRR14131127,SRR14131128 ncbi_delly_final_bi_rehead.forangsd.vcf.gz > ncbi_delly_mongolia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127744,SRR14131090,SRR14131091 ncbi_delly_final_bi_rehead.forangsd.vcf.gz > ncbi_delly_china_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127745,SRR14127746,SRR14127747,SRR14127748,SRR14127749,SRR14127750,SRR14127751,SRR14127752,SRR14127753 ncbi_delly_final_bi_rehead.forangsd.vcf.gz > ncbi_delly_russia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14131088,SRR14131089,SRR14131092,SRR14131094,SRR14131123,SRR14131124,SRR14131125,SRR14131126 ncbi_delly_final_bi_rehead.forangsd.vcf.gz > ncbi_delly_morocco_bi_rehead.forangsd.vcf.gz

angsd -vcf-gl ncbi_delly_mongolia_bi_rehead.forangsd.vcf.gz -nind 5 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -dosaf 1 -doMajorMinor 5 -P ${N} -out mongolia_sv
realSFS mongolia_sv.saf.idx -P ${N} -fold 1 > mongolia_sv.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -vcf-gl ncbi_delly_china_bi_rehead.forangsd.vcf.gz -nind 3 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -dosaf 1 -doMajorMinor 5 -P ${N} -out china_sv
realSFS china_sv.saf.idx -P ${N} -fold 1 > china_sv.sfs 

angsd -vcf-gl ncbi_delly_russia_bi_rehead.forangsd.vcf.gz -nind 9 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -dosaf 1 -doMajorMinor 5 -P ${N} -out russia_sv
realSFS russia_sv.saf.idx -P ${N} -fold 1 > russia_sv.sfs 

angsd -vcf-gl ncbi_delly_morocco_bi_rehead.forangsd.vcf.gz -nind 8 -fai ${REF}/${GENOME}.fna.fai -anc ${REF}/${GENOME}.fna -dosaf 1 -doMajorMinor 5 -P ${N} -out morocco_sv
realSFS morocco_sv.saf.idx -P ${N} -fold 1 > morocco_sv.sfs 

realSFS mongolia_sv.saf.idx china_sv.saf.idx -P ${N} > mongolia.china_sv.ml # calculate the 2D SFS
realSFS mongolia_sv.saf.idx russia_sv.saf.idx -P ${N} > mongolia.russia_sv.ml
realSFS mongolia_sv.saf.idx morocco_sv.saf.idx -P ${N} > mongolia.morocco_sv.ml
realSFS china_sv.saf.idx russia_sv.saf.idx -P ${N} > china.russia_sv.ml
realSFS china_sv.saf.idx morocco_sv.saf.idx -P ${N} > china.morocco_sv.ml
realSFS russia_sv.saf.idx morocco_sv.saf.idx -P ${N} > russia.morocco_sv.ml

realSFS fst index mongolia_sv.saf.idx china_sv.saf.idx -sfs mongolia.china_sv.ml -fstout mongolia.china_sv # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_sv.saf.idx russia_sv.saf.idx -sfs mongolia.russia_sv.ml -fstout mongolia.russia_sv
realSFS fst index mongolia_sv.saf.idx morocco_sv.saf.idx -sfs mongolia.morocco_sv.ml -fstout mongolia.morocco_sv
realSFS fst index china_sv.saf.idx russia_sv.saf.idx -sfs china.russia_sv.ml -fstout china.russia_sv
realSFS fst index china_sv.saf.idx morocco_sv.saf.idx -sfs china.morocco_sv.ml -fstout china.morocco_sv
realSFS fst index russia_sv.saf.idx morocco_sv.saf.idx -sfs russia.morocco_sv.ml -fstout russia.morocco_sv

realSFS fst stats mongolia.china_sv.fst.idx # global pairwise estimates # FST.Unweight[nObs:1268]:0.049627 Fst.Weight:0.046567
realSFS fst stats mongolia.russia_sv.fst.idx # FST.Unweight[nObs:1268]:0.034028 Fst.Weight:0.031241
realSFS fst stats mongolia.morocco_sv.fst.idx # FST.Unweight[nObs:1268]:0.040300 Fst.Weight:0.039127
realSFS fst stats china.russia_sv.fst.idx # FST.Unweight[nObs:1268]:0.041218 Fst.Weight:0.040243
realSFS fst stats china.morocco_sv.fst.idx # FST.Unweight[nObs:1268]:0.037791 Fst.Weight:0.036081
realSFS fst stats russia.morocco_sv.fst.idx # FST.Unweight[nObs:1268]:0.024782 Fst.Weight:0.023147

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${REF}/${GENOME}.fna -ref ${REF}/${GENOME}.fna -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca # generate SNP beagle

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

cat ncbi_delly.beagleheader ncbi_delly.norm.beagle > ncbi_delly.ready.beagle #add the proper header again and zip
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
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("ncbi_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.25%, PC2: 4.35%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("ncbi_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.25%)", y = "Principal Component 2 (4.35%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()

C <- as.matrix(read.table("ncbi_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.13%, PC2: 6.03%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("ncbi_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.13%)", y = "Principal Component 2 (6.03%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()
quit()

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
quit()

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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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

### Formating SV vcf file for ANGSD (modified based on Mérot et al. (2023) in Molecular Ecology)
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

### Generating GL file or beagle file for each chromosome (modified based on Mérot et al. (2023) in Molecular Ecology)
ls /${LINPAN}/mapped/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o linpan_delly_final_bi_rehead.forangsd.vcf --threads ${N} linpan_delly_final_bi.forangsd.vcf
bgzip linpan_delly_final_bi_rehead.forangsd.vcf 
tabix linpan_delly_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl linpan_delly_final_bi_rehead.forangsd.vcf.gz -nind 25 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

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
R
norm <- function(x) x/sum(x) # function to normalize 
sfs <- (scan("out_snp.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("linpan_snp_sfs_plot.svg", width=5, height=5)
barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()

sfs <- (scan("out_sv.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("linpan_sv_sfs_plot.svg", width=5, height=5)
barplot(sfs,xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()
quit()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_linpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_linpan_PLraw.txt ROH_linpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### pairwise Fst
cat ./bam.filelist | grep -E "SRR14127742|SRR14131105|SRR14131116|SRR14131127|SRR14131128" > ./bam.mongolia
cat ./bam.filelist | grep -E "SRR14127744|SRR14131090|SRR14131091" > ./bam.china
cat ./bam.filelist | grep -E "SRR14127745|SRR14127746|SRR14127747|SRR14127748|SRR14127749|SRR14127750|SRR14127751|SRR14127752|SRR14127753" > ./bam.russia
cat ./bam.filelist | grep -E "SRR14131088|SRR14131089|SRR14131092|SRR14131094|SRR14131123|SRR14131124|SRR14131125|SRR14131126" > ./bam.morocco

angsd -bam ./bam.mongolia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 4 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out mongolia_snp
realSFS mongolia_snp.saf.idx -P ${N} -fold 1 > mongolia_snp.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -bam ./bam.china -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 2 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out china_snp
realSFS china_snp.saf.idx -P ${N} -fold 1 > china_snp.sfs

angsd -bam ./bam.russia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 7 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out russia_snp
realSFS russia_snp.saf.idx -P ${N} -fold 1 > russia_snp.sfs

angsd -bam ./bam.morocco -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 6 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out morocco_snp
realSFS morocco_snp.saf.idx -P ${N} -fold 1 > morocco_snp.sfs

realSFS mongolia_snp.saf.idx china_snp.saf.idx -P ${N} > mongolia.china_snp.ml # calculate the 2D SFS
realSFS mongolia_snp.saf.idx russia_snp.saf.idx -P ${N} > mongolia.russia_snp.ml
realSFS mongolia_snp.saf.idx morocco_snp.saf.idx -P ${N} > mongolia.morocco_snp.ml
realSFS china_snp.saf.idx russia_snp.saf.idx -P ${N} > china.russia_snp.ml
realSFS china_snp.saf.idx morocco_snp.saf.idx -P ${N} > china.morocco_snp.ml
realSFS russia_snp.saf.idx morocco_snp.saf.idx -P ${N} > russia.morocco_snp.ml

realSFS fst index mongolia_snp.saf.idx china_snp.saf.idx -sfs mongolia.china_snp.ml -fstout mongolia.china_snp # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_snp.saf.idx russia_snp.saf.idx -sfs mongolia.russia_snp.ml -fstout mongolia.russia_snp
realSFS fst index mongolia_snp.saf.idx morocco_snp.saf.idx -sfs mongolia.morocco_snp.ml -fstout mongolia.morocco_snp
realSFS fst index china_snp.saf.idx russia_snp.saf.idx -sfs china.russia_snp.ml -fstout china.russia_snp
realSFS fst index china_snp.saf.idx morocco_snp.saf.idx -sfs china.morocco_snp.ml -fstout china.morocco_snp
realSFS fst index russia_snp.saf.idx morocco_snp.saf.idx -sfs russia.morocco_snp.ml -fstout russia.morocco_snp

realSFS fst stats mongolia.china_snp.fst.idx # global pairwise estimates # FST.Unweight[nObs:908846629]:0.071064 Fst.Weight:0.084865
realSFS fst stats mongolia.russia_snp.fst.idx # FST.Unweight[nObs:823297647]:0.029476 Fst.Weight:0.052580
realSFS fst stats mongolia.morocco_snp.fst.idx # FST.Unweight[nObs:881360295]:0.042870 Fst.Weight:0.056429
realSFS fst stats china.russia_snp.fst.idx # FST.Unweight[nObs:835642436]:0.028328 Fst.Weight:0.067021
realSFS fst stats china.morocco_snp.fst.idx # FST.Unweight[nObs:902184547]:0.054691 Fst.Weight:0.070080
realSFS fst stats russia.morocco_snp.fst.idx # FST.Unweight[nObs:823247825]:0.029897 Fst.Weight:0.035491

bcftools view -s SRR14127742,SRR14131105,SRR14131116,SRR14131127,SRR14131128 linpan_delly_final_bi_rehead.forangsd.vcf.gz > linpan_delly_mongolia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127744,SRR14131090,SRR14131091 linpan_delly_final_bi_rehead.forangsd.vcf.gz > linpan_delly_china_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127745,SRR14127746,SRR14127747,SRR14127748,SRR14127749,SRR14127750,SRR14127751,SRR14127752,SRR14127753 linpan_delly_final_bi_rehead.forangsd.vcf.gz > linpan_delly_russia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14131088,SRR14131089,SRR14131092,SRR14131094,SRR14131123,SRR14131124,SRR14131125,SRR14131126 linpan_delly_final_bi_rehead.forangsd.vcf.gz > linpan_delly_morocco_bi_rehead.forangsd.vcf.gz

angsd -vcf-gl linpan_delly_mongolia_bi_rehead.forangsd.vcf.gz -nind 5 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out mongolia_sv
realSFS mongolia_sv.saf.idx -P ${N} -fold 1 > mongolia_sv.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -vcf-gl linpan_delly_china_bi_rehead.forangsd.vcf.gz -nind 3 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out china_sv
realSFS china_sv.saf.idx -P ${N} -fold 1 > china_sv.sfs 

angsd -vcf-gl linpan_delly_russia_bi_rehead.forangsd.vcf.gz -nind 9 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out russia_sv
realSFS russia_sv.saf.idx -P ${N} -fold 1 > russia_sv.sfs 

angsd -vcf-gl linpan_delly_morocco_bi_rehead.forangsd.vcf.gz -nind 8 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out morocco_sv
realSFS morocco_sv.saf.idx -P ${N} -fold 1 > morocco_sv.sfs 

realSFS mongolia_sv.saf.idx china_sv.saf.idx -P ${N} > mongolia.china_sv.ml # calculate the 2D SFS
realSFS mongolia_sv.saf.idx russia_sv.saf.idx -P ${N} > mongolia.russia_sv.ml
realSFS mongolia_sv.saf.idx morocco_sv.saf.idx -P ${N} > mongolia.morocco_sv.ml
realSFS china_sv.saf.idx russia_sv.saf.idx -P ${N} > china.russia_sv.ml
realSFS china_sv.saf.idx morocco_sv.saf.idx -P ${N} > china.morocco_sv.ml
realSFS russia_sv.saf.idx morocco_sv.saf.idx -P ${N} > russia.morocco_sv.ml

realSFS fst index mongolia_sv.saf.idx china_sv.saf.idx -sfs mongolia.china_sv.ml -fstout mongolia.china_sv # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_sv.saf.idx russia_sv.saf.idx -sfs mongolia.russia_sv.ml -fstout mongolia.russia_sv
realSFS fst index mongolia_sv.saf.idx morocco_sv.saf.idx -sfs mongolia.morocco_sv.ml -fstout mongolia.morocco_sv
realSFS fst index china_sv.saf.idx russia_sv.saf.idx -sfs china.russia_sv.ml -fstout china.russia_sv
realSFS fst index china_sv.saf.idx morocco_sv.saf.idx -sfs china.morocco_sv.ml -fstout china.morocco_sv
realSFS fst index russia_sv.saf.idx morocco_sv.saf.idx -sfs russia.morocco_sv.ml -fstout russia.morocco_sv

realSFS fst stats mongolia.china_sv.fst.idx # global pairwise estimates # FST.Unweight[nObs:1288]:0.050093 Fst.Weight:0.047354
realSFS fst stats mongolia.russia_sv.fst.idx # FST.Unweight[nObs:1288]:0.034367 Fst.Weight:0.031867
realSFS fst stats mongolia.morocco_sv.fst.idx # FST.Unweight[nObs:1288]:0.040422 Fst.Weight:0.039238
realSFS fst stats china.russia_sv.fst.idx # FST.Unweight[nObs:1288]:0.040697 Fst.Weight:0.039991
realSFS fst stats china.morocco_sv.fst.idx # FST.Unweight[nObs:1288]:0.038274 Fst.Weight:0.036810
realSFS fst stats russia.morocco_sv.fst.idx # FST.Unweight[nObs:1288]:0.024553 Fst.Weight:0.022950

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca # generate SNP beagle

picard UpdateVcfSequenceDictionary I=./linpan_delly_final_simpl.vcf O=./linpan_delly_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./linpan_delly_final_updt.vcf O=./linpan_delly_final_sorted.vcf
bcftools +tag2tag linpan_delly_final_sorted.vcf -- -r --pl-to-gl > linpan_delly_final_sortedGL.vcf

cat ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai | cut -f1 | while read i # generate SV beagle from here
do 
echo "convert vcf to beagle for ${i}"
vcftools --vcf linpan_delly_final_sortedGL.vcf --BEAGLE-GL --chr ${i} --out linpan_delly_${i}
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
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("linpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.24%, PC2: 4.33%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("linpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.24%)", y = "Principal Component 2 (4.33%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()

C <- as.matrix(read.table("linpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.01%, PC2: 5.93%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("linpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.01%)", y = "Principal Component 2 (5.93%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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

### Formating a merged SV vcf file for ANGSD (modified based on Mérot et al. (2023) in Molecular Ecology)
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

### Generating GL or beagle file for each chromosome (modified based on Mérot et al. (2023) in Molecular Ecology)
ls ${MCPAN}/${PREFIX}-pg/aligned/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o mcpan_SV_final_bi_rehead.forangsd.vcf --threads ${N} mcpan_SV_final_bi.forangsd.vcf
bgzip  mcpan_SV_final_bi_rehead.forangsd.vcf
tabix mcpan_SV_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl mcpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 25 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

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
R
norm <- function(x) x/sum(x) # function to normalize 
sfs <- (scan("out_snp.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("mcpan_snp_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()

sfs <- (scan("out_sv.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("mcpan_sv_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()
quit()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_mcpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_mcpan_PLraw.txt ROH_mcpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa.fai # parse raw ROH file

### pairwise Fst
cat ./bam.filelist | grep -E "SRR14127742|SRR14131105|SRR14131116|SRR14131127|SRR14131128" > ./bam.mongolia
cat ./bam.filelist | grep -E "SRR14127744|SRR14131090|SRR14131091" > ./bam.china
cat ./bam.filelist | grep -E "SRR14127745|SRR14127746|SRR14127747|SRR14127748|SRR14127749|SRR14127750|SRR14127751|SRR14127752|SRR14127753" > ./bam.russia
cat ./bam.filelist | grep -E "SRR14131088|SRR14131089|SRR14131092|SRR14131094|SRR14131123|SRR14131124|SRR14131125|SRR14131126" > ./bam.morocco

angsd -bam ./bam.mongolia -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 4 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out mongolia_snp
realSFS mongolia_snp.saf.idx -P ${N} -fold 1 > mongolia_snp.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -bam ./bam.china -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 2 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out china_snp
realSFS china_snp.saf.idx -P ${N} -fold 1 > china_snp.sfs

angsd -bam ./bam.russia -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 7 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out russia_snp
realSFS russia_snp.saf.idx -P ${N} -fold 1 > russia_snp.sfs

angsd -bam ./bam.morocco -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -minQ 30 -minInd 6 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out morocco_snp
realSFS morocco_snp.saf.idx -P ${N} -fold 1 > morocco_snp.sfs

realSFS mongolia_snp.saf.idx china_snp.saf.idx -P ${N} > mongolia.china_snp.ml # calculate the 2D SFS
realSFS mongolia_snp.saf.idx russia_snp.saf.idx -P ${N} > mongolia.russia_snp.ml
realSFS mongolia_snp.saf.idx morocco_snp.saf.idx -P ${N} > mongolia.morocco_snp.ml
realSFS china_snp.saf.idx russia_snp.saf.idx -P ${N} > china.russia_snp.ml
realSFS china_snp.saf.idx morocco_snp.saf.idx -P ${N} > china.morocco_snp.ml
realSFS russia_snp.saf.idx morocco_snp.saf.idx -P ${N} > russia.morocco_snp.ml

realSFS fst index mongolia_snp.saf.idx china_snp.saf.idx -sfs mongolia.china_snp.ml -fstout mongolia.china_snp # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_snp.saf.idx russia_snp.saf.idx -sfs mongolia.russia_snp.ml -fstout mongolia.russia_snp
realSFS fst index mongolia_snp.saf.idx morocco_snp.saf.idx -sfs mongolia.morocco_snp.ml -fstout mongolia.morocco_snp
realSFS fst index china_snp.saf.idx russia_snp.saf.idx -sfs china.russia_snp.ml -fstout china.russia_snp
realSFS fst index china_snp.saf.idx morocco_snp.saf.idx -sfs china.morocco_snp.ml -fstout china.morocco_snp
realSFS fst index russia_snp.saf.idx morocco_snp.saf.idx -sfs russia.morocco_snp.ml -fstout russia.morocco_snp

realSFS fst stats mongolia.china_snp.fst.idx # global pairwise estimates # FST.Unweight[nObs:873944665]:0.071115 Fst.Weight:0.084861
realSFS fst stats mongolia.russia_snp.fst.idx # FST.Unweight[nObs:785684371]:0.029787 Fst.Weight:0.052580
realSFS fst stats mongolia.morocco_snp.fst.idx # FST.Unweight[nObs:845697138]:0.042728 Fst.Weight:0.056361
realSFS fst stats china.russia_snp.fst.idx # FST.Unweight[nObs:798242932]:0.029064 Fst.Weight:0.066969
realSFS fst stats china.morocco_snp.fst.idx # FST.Unweight[nObs:867016092]:0.054718 Fst.Weight:0.069956
realSFS fst stats russia.morocco_snp.fst.idx # FST.Unweight[nObs:785643777]:0.029949 Fst.Weight:0.035467

bcftools view -s SRR14127742,SRR14131105,SRR14131116,SRR14131127,SRR14131128 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_mongolia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127744,SRR14131090,SRR14131091 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_china_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127745,SRR14127746,SRR14127747,SRR14127748,SRR14127749,SRR14127750,SRR14127751,SRR14127752,SRR14127753 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_russia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14131088,SRR14131089,SRR14131092,SRR14131094,SRR14131123,SRR14131124,SRR14131125,SRR14131126 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_morocco_bi_rehead.forangsd.vcf.gz

angsd -vcf-gl mcpan_SV_mongolia_bi_rehead.forangsd.vcf.gz -nind 5 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out mongolia_sv
realSFS mongolia_sv.saf.idx -P ${N} -fold 1 > mongolia_sv.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -vcf-gl mcpan_SV_china_bi_rehead.forangsd.vcf.gz -nind 3 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out china_sv
realSFS china_sv.saf.idx -P ${N} -fold 1 > china_sv.sfs 

angsd -vcf-gl mcpan_SV_russia_bi_rehead.forangsd.vcf.gz -nind 9 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out russia_sv
realSFS russia_sv.saf.idx -P ${N} -fold 1 > russia_sv.sfs 

angsd -vcf-gl mcpan_SV_morocco_bi_rehead.forangsd.vcf.gz -nind 8 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out morocco_sv
realSFS morocco_sv.saf.idx -P ${N} -fold 1 > morocco_sv.sfs 

realSFS mongolia_sv.saf.idx china_sv.saf.idx -P ${N} > mongolia.china_sv.ml # calculate the 2D SFS
realSFS mongolia_sv.saf.idx russia_sv.saf.idx -P ${N} > mongolia.russia_sv.ml
realSFS mongolia_sv.saf.idx morocco_sv.saf.idx -P ${N} > mongolia.morocco_sv.ml
realSFS china_sv.saf.idx russia_sv.saf.idx -P ${N} > china.russia_sv.ml
realSFS china_sv.saf.idx morocco_sv.saf.idx -P ${N} > china.morocco_sv.ml
realSFS russia_sv.saf.idx morocco_sv.saf.idx -P ${N} > russia.morocco_sv.ml

realSFS fst index mongolia_sv.saf.idx china_sv.saf.idx -sfs mongolia.china_sv.ml -fstout mongolia.china_sv # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_sv.saf.idx russia_sv.saf.idx -sfs mongolia.russia_sv.ml -fstout mongolia.russia_sv
realSFS fst index mongolia_sv.saf.idx morocco_sv.saf.idx -sfs mongolia.morocco_sv.ml -fstout mongolia.morocco_sv
realSFS fst index china_sv.saf.idx russia_sv.saf.idx -sfs china.russia_sv.ml -fstout china.russia_sv
realSFS fst index china_sv.saf.idx morocco_sv.saf.idx -sfs china.morocco_sv.ml -fstout china.morocco_sv
realSFS fst index russia_sv.saf.idx morocco_sv.saf.idx -sfs russia.morocco_sv.ml -fstout russia.morocco_sv

realSFS fst stats mongolia.china_sv.fst.idx # global pairwise estimates # FST.Unweight[nObs:2048]:0.067214 Fst.Weight:0.072245
realSFS fst stats mongolia.russia_sv.fst.idx # FST.Unweight[nObs:2048]:0.050682 Fst.Weight:0.053544
realSFS fst stats mongolia.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.052869 Fst.Weight:0.055217
realSFS fst stats china.russia_sv.fst.idx # FST.Unweight[nObs:2048]:0.052437 Fst.Weight:0.062566
realSFS fst stats china.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.056151 Fst.Weight:0.061786
realSFS fst stats russia.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.034044 Fst.Weight:0.035484

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.renamed.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca

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
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("mcpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.22%, PC2: 4.34%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("mcpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.22%)", y = "Principal Component 2 (4.34%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()

C <- as.matrix(read.table("mcpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.49%, PC2: 5.17%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("mcpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.49%)", y = "Principal Component 2 (5.17%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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

### Formating a merged SV vcf file for ANGSD (modified based on Mérot et al. (2023) in Molecular Ecology)
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

### Generating GL or beagle file for each chromosome (modified based on Mérot et al. (2023) in Molecular Ecology)
ls ${VGPAN}/aligned/*.marked.bam > ./bam.filelist
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_snp

bcftools reheader --fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -o vgpan_SV_final_bi_rehead.forangsd.vcf --threads ${N} vgpan_SV_final_bi.forangsd.vcf
bgzip  vgpan_SV_final_bi_rehead.forangsd.vcf
tabix vgpan_SV_final_bi_rehead.forangsd.vcf.gz
angsd -vcf-gl vgpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 25 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

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
R
norm <- function(x) x/sum(x) # function to normalize 
sfs <- (scan("out_snp.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("vgpan_snp_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()

sfs <- (scan("out_sv.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("vgpan_sv_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()
quit()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_vgpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_vgpan_PLraw.txt ROH_vgpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### pairwise Fst
cat ./bam.filelist | grep -E "SRR14127742|SRR14131105|SRR14131116|SRR14131127|SRR14131128" > ./bam.mongolia
cat ./bam.filelist | grep -E "SRR14127744|SRR14131090|SRR14131091" > ./bam.china
cat ./bam.filelist | grep -E "SRR14127745|SRR14127746|SRR14127747|SRR14127748|SRR14127749|SRR14127750|SRR14127751|SRR14127752|SRR14127753" > ./bam.russia
cat ./bam.filelist | grep -E "SRR14131088|SRR14131089|SRR14131092|SRR14131094|SRR14131123|SRR14131124|SRR14131125|SRR14131126" > ./bam.morocco

angsd -bam ./bam.mongolia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 4 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out mongolia_snp
realSFS mongolia_snp.saf.idx -P ${N} -fold 1 > mongolia_snp.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -bam ./bam.china -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 2 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out china_snp
realSFS china_snp.saf.idx -P ${N} -fold 1 > china_snp.sfs

angsd -bam ./bam.russia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 7 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out russia_snp
realSFS russia_snp.saf.idx -P ${N} -fold 1 > russia_snp.sfs

angsd -bam ./bam.morocco -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 6 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out morocco_snp
realSFS morocco_snp.saf.idx -P ${N} -fold 1 > morocco_snp.sfs

realSFS mongolia_snp.saf.idx china_snp.saf.idx -P ${N} > mongolia.china_snp.ml # calculate the 2D SFS
realSFS mongolia_snp.saf.idx russia_snp.saf.idx -P ${N} > mongolia.russia_snp.ml
realSFS mongolia_snp.saf.idx morocco_snp.saf.idx -P ${N} > mongolia.morocco_snp.ml
realSFS china_snp.saf.idx russia_snp.saf.idx -P ${N} > china.russia_snp.ml
realSFS china_snp.saf.idx morocco_snp.saf.idx -P ${N} > china.morocco_snp.ml
realSFS russia_snp.saf.idx morocco_snp.saf.idx -P ${N} > russia.morocco_snp.ml

realSFS fst index mongolia_snp.saf.idx china_snp.saf.idx -sfs mongolia.china_snp.ml -fstout mongolia.china_snp # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_snp.saf.idx russia_snp.saf.idx -sfs mongolia.russia_snp.ml -fstout mongolia.russia_snp
realSFS fst index mongolia_snp.saf.idx morocco_snp.saf.idx -sfs mongolia.morocco_snp.ml -fstout mongolia.morocco_snp
realSFS fst index china_snp.saf.idx russia_snp.saf.idx -sfs china.russia_snp.ml -fstout china.russia_snp
realSFS fst index china_snp.saf.idx morocco_snp.saf.idx -sfs china.morocco_snp.ml -fstout china.morocco_snp
realSFS fst index russia_snp.saf.idx morocco_snp.saf.idx -sfs russia.morocco_snp.ml -fstout russia.morocco_snp

realSFS fst stats mongolia.china_snp.fst.idx # global pairwise estimates # FST.Unweight[nObs:880200852]:0.071122 Fst.Weight:0.085175
realSFS fst stats mongolia.russia_snp.fst.idx # FST.Unweight[nObs:793171804]:0.029713 Fst.Weight:0.052771 
realSFS fst stats mongolia.morocco_snp.fst.idx # FST.Unweight[nObs:852179958]:0.042682 Fst.Weight:0.056572
realSFS fst stats china.russia_snp.fst.idx # FST.Unweight[nObs:805424126]:0.028959 Fst.Weight:0.067220
realSFS fst stats china.morocco_snp.fst.idx # FST.Unweight[nObs:872807902]:0.054683 Fst.Weight:0.070181
realSFS fst stats russia.morocco_snp.fst.idx # FST.Unweight[nObs:792883934]:0.029935 Fst.Weight:0.035623

bcftools view -s SRR14127742,SRR14131105,SRR14131116,SRR14131127,SRR14131128 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_mongolia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127744,SRR14131090,SRR14131091 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_china_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127745,SRR14127746,SRR14127747,SRR14127748,SRR14127749,SRR14127750,SRR14127751,SRR14127752,SRR14127753 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_russia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14131088,SRR14131089,SRR14131092,SRR14131094,SRR14131123,SRR14131124,SRR14131125,SRR14131126 mcpan_SV_final_bi_rehead.forangsd.vcf.gz > mcpan_SV_morocco_bi_rehead.forangsd.vcf.gz

angsd -vcf-gl mcpan_SV_mongolia_bi_rehead.forangsd.vcf.gz -nind 5 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out mongolia_sv
realSFS mongolia_sv.saf.idx -P ${N} -fold 1 > mongolia_sv.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -vcf-gl mcpan_SV_china_bi_rehead.forangsd.vcf.gz -nind 3 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out china_sv
realSFS china_sv.saf.idx -P ${N} -fold 1 > china_sv.sfs 

angsd -vcf-gl mcpan_SV_russia_bi_rehead.forangsd.vcf.gz -nind 9 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out russia_sv
realSFS russia_sv.saf.idx -P ${N} -fold 1 > russia_sv.sfs 

angsd -vcf-gl mcpan_SV_morocco_bi_rehead.forangsd.vcf.gz -nind 8 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out morocco_sv
realSFS morocco_sv.saf.idx -P ${N} -fold 1 > morocco_sv.sfs 

realSFS mongolia_sv.saf.idx china_sv.saf.idx -P ${N} > mongolia.china_sv.ml # calculate the 2D SFS
realSFS mongolia_sv.saf.idx russia_sv.saf.idx -P ${N} > mongolia.russia_sv.ml
realSFS mongolia_sv.saf.idx morocco_sv.saf.idx -P ${N} > mongolia.morocco_sv.ml
realSFS china_sv.saf.idx russia_sv.saf.idx -P ${N} > china.russia_sv.ml
realSFS china_sv.saf.idx morocco_sv.saf.idx -P ${N} > china.morocco_sv.ml
realSFS russia_sv.saf.idx morocco_sv.saf.idx -P ${N} > russia.morocco_sv.ml

realSFS fst index mongolia_sv.saf.idx china_sv.saf.idx -sfs mongolia.china_sv.ml -fstout mongolia.china_sv # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_sv.saf.idx russia_sv.saf.idx -sfs mongolia.russia_sv.ml -fstout mongolia.russia_sv
realSFS fst index mongolia_sv.saf.idx morocco_sv.saf.idx -sfs mongolia.morocco_sv.ml -fstout mongolia.morocco_sv
realSFS fst index china_sv.saf.idx russia_sv.saf.idx -sfs china.russia_sv.ml -fstout china.russia_sv
realSFS fst index china_sv.saf.idx morocco_sv.saf.idx -sfs china.morocco_sv.ml -fstout china.morocco_sv
realSFS fst index russia_sv.saf.idx morocco_sv.saf.idx -sfs russia.morocco_sv.ml -fstout russia.morocco_sv

realSFS fst stats mongolia.china_sv.fst.idx # global pairwise estimates # FST.Unweight[nObs:2048]:0.067214 Fst.Weight:0.072245
realSFS fst stats mongolia.russia_sv.fst.idx # FST.Unweight[nObs:2048]:0.050682 Fst.Weight:0.053544
realSFS fst stats mongolia.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.052869 Fst.Weight:0.055217
realSFS fst stats china.russia_sv.fst.idx # FST.Unweight[nObs:2048]:0.052437 Fst.Weight:0.062566
realSFS fst stats china.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.056151 Fst.Weight:0.061786
realSFS fst stats russia.morocco_sv.fst.idx # FST.Unweight[nObs:2048]:0.034044 Fst.Weight:0.035484

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
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("vgpan_SNP_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.21%, PC2: 4.33% 
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("vgpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.21%)", y = "Principal Component 2 (4.33%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()

C <- as.matrix(read.table("vgpan_SV_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 9.17%, PC2: 5.16%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("vgpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (9.17%)", y = "Principal Component 2 (5.16%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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

### Formating a merged SV vcf file for ANGSD (modified based on Mérot et al. (2023) in Molecular Ecology)
bcftools annotate -x ^FORMAT/GT,FORMAT/GL,FORMAT/GQ orgpan_vg_final_annot.vcf > orgpan_vg_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format
picard UpdateVcfSequenceDictionary I=./orgpan_vg_final_simpl.vcf O=./orgpan_vg_final_updt.vcf SEQUENCE_DICTIONARY=${LINPAN}/${PREFIX}_panref2_sorted.dict
picard SortVcf I=./orgpan_vg_final_updt.vcf O=./orgpan_vg_final_sorted.vcf
bcftools +tag2tag orgpan_vg_final_sorted.vcf -- -r --pl-to-gl > orgpan_vg_final_sortedGL.vcf

bcftools annotate -x ^FORMAT/GT,FORMAT/PL,FORMAT/GQ orgpan_delly_final_sorted.vcf > orgpan_delly_final_simpl.vcf # remove unnecessary FORMAT fields to match the vcf format # "FORMAT/PL", not "FORMAT/GL", used here
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

### Generating GL or beagle file for each chromosome (modified based on Mérot et al. (2023) in Molecular Ecology)
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
angsd -vcf-gl orgpan_SV_final_bi_rehead.forangsd.vcf.gz -nind 25 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -domaf 1 -dosaf 1 -doMajorMinor 5 -P ${N} -out out_sv

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
R
norm <- function(x) x/sum(x) # function to normalize 
sfs <- (scan("out_snp.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("orgpan_snp_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()

sfs <- (scan("out_sv.sfs")[2:26]) # read data
sfs<-norm(sfs)  # the variable categories of the sfs 
svg("orgpan_sv_sfs_plot.svg", width=5, height=5)
barplot(sfs, xlab="Allele frequency",names=1:length(sfs),ylab="Proportions",main="Site Frequency Spectrum plot",col='blue')
dev.off()
quit()

### ROH
angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -snp_pval 1e-6 -dobcf 1 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out out_roh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' out_roh.bcf | bgzip -c > out_roh.freqs.tab.gz
tabix -s1 -b2 -e2 out_roh.freqs.tab.gz
bcftools roh --AF-file out_roh.freqs.tab.gz --output ROH_mcpan_PLraw.txt --threads ${N} out_roh.bcf # ROH estimation

python3 ${HOME}/ROHparser-pg.py ROH_mcpan_PLraw.txt ROH_mcpan_PLresult.txt ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai # parse raw ROH file

### pairwise Fst
cat ./bam.filelist | grep -E "SRR14127742|SRR14131105|SRR14131116|SRR14131127|SRR14131128" > ./bam.mongolia
cat ./bam.filelist | grep -E "SRR14127744|SRR14131090|SRR14131091" > ./bam.china
cat ./bam.filelist | grep -E "SRR14127745|SRR14127746|SRR14127747|SRR14127748|SRR14127749|SRR14127750|SRR14127751|SRR14127752|SRR14127753" > ./bam.russia
cat ./bam.filelist | grep -E "SRR14131088|SRR14131089|SRR14131092|SRR14131094|SRR14131123|SRR14131124|SRR14131125|SRR14131126" > ./bam.morocco

angsd -bam ./bam.mongolia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 4 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out mongolia_snp
realSFS mongolia_snp.saf.idx -P ${N} -fold 1 > mongolia_snp.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -bam ./bam.china -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 2 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out china_snp
realSFS china_snp.saf.idx -P ${N} -fold 1 > china_snp.sfs

angsd -bam ./bam.russia -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 7 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out russia_snp
realSFS russia_snp.saf.idx -P ${N} -fold 1 > russia_snp.sfs

angsd -bam ./bam.morocco -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -minQ 30 -minInd 6 -docounts 1 -domajorminor 5 -dosaf 1 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out morocco_snp
realSFS morocco_snp.saf.idx -P ${N} -fold 1 > morocco_snp.sfs

realSFS mongolia_snp.saf.idx china_snp.saf.idx -P ${N} > mongolia.china_snp.ml # calculate the 2D SFS
realSFS mongolia_snp.saf.idx russia_snp.saf.idx -P ${N} > mongolia.russia_snp.ml
realSFS mongolia_snp.saf.idx morocco_snp.saf.idx -P ${N} > mongolia.morocco_snp.ml
realSFS china_snp.saf.idx russia_snp.saf.idx -P ${N} > china.russia_snp.ml
realSFS china_snp.saf.idx morocco_snp.saf.idx -P ${N} > china.morocco_snp.ml
realSFS russia_snp.saf.idx morocco_snp.saf.idx -P ${N} > russia.morocco_snp.ml

realSFS fst index mongolia_snp.saf.idx china_snp.saf.idx -sfs mongolia.china_snp.ml -fstout mongolia.china_snp # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_snp.saf.idx russia_snp.saf.idx -sfs mongolia.russia_snp.ml -fstout mongolia.russia_snp
realSFS fst index mongolia_snp.saf.idx morocco_snp.saf.idx -sfs mongolia.morocco_snp.ml -fstout mongolia.morocco_snp
realSFS fst index china_snp.saf.idx russia_snp.saf.idx -sfs china.russia_snp.ml -fstout china.russia_snp
realSFS fst index china_snp.saf.idx morocco_snp.saf.idx -sfs china.morocco_snp.ml -fstout china.morocco_snp
realSFS fst index russia_snp.saf.idx morocco_snp.saf.idx -sfs russia.morocco_snp.ml -fstout russia.morocco_snp

realSFS fst stats mongolia.china_snp.fst.idx # global pairwise estimates # FST.Unweight[nObs:873481558]:0.071104 Fst.Weight:0.085555
realSFS fst stats mongolia.russia_snp.fst.idx # FST.Unweight[nObs:787784820]:0.029710 Fst.Weight:0.052998
realSFS fst stats mongolia.morocco_snp.fst.idx # FST.Unweight[nObs:845857125]:0.042680 Fst.Weight:0.056830
realSFS fst stats china.russia_snp.fst.idx # FST.Unweight[nObs:799872132]:0.028943 Fst.Weight:0.067522
realSFS fst stats china.morocco_snp.fst.idx # FST.Unweight[nObs:866019601]:0.054670 Fst.Weight:0.070522
realSFS fst stats russia.morocco_snp.fst.idx # FST.Unweight[nObs:787455919]:0.029935 Fst.Weight:0.035796

bcftools view -s SRR14127742,SRR14131105,SRR14131116,SRR14131127,SRR14131128 orgpan_SV_final_bi_rehead.forangsd.vcf.gz > orgpan_SV_mongolia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127744,SRR14131090,SRR14131091 orgpan_SV_final_bi_rehead.forangsd.vcf.gz > orgpan_SV_china_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14127745,SRR14127746,SRR14127747,SRR14127748,SRR14127749,SRR14127750,SRR14127751,SRR14127752,SRR14127753 orgpan_SV_final_bi_rehead.forangsd.vcf.gz > orgpan_SV_russia_bi_rehead.forangsd.vcf.gz
bcftools view -s SRR14131088,SRR14131089,SRR14131092,SRR14131094,SRR14131123,SRR14131124,SRR14131125,SRR14131126 orgpan_SV_final_bi_rehead.forangsd.vcf.gz > orgpan_SV_morocco_bi_rehead.forangsd.vcf.gz

angsd -vcf-gl orgpan_SV_mongolia_bi_rehead.forangsd.vcf.gz -nind 5 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out mongolia_sv
realSFS mongolia_sv.saf.idx -P ${N} -fold 1 > mongolia_sv.sfs # calculate the 1D SFS from allele freq likelihoods

angsd -vcf-gl orgpan_SV_china_bi_rehead.forangsd.vcf.gz -nind 3 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out china_sv
realSFS china_sv.saf.idx -P ${N} -fold 1 > china_sv.sfs 

angsd -vcf-gl orgpan_SV_russia_bi_rehead.forangsd.vcf.gz -nind 9 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out russia_sv
realSFS russia_sv.saf.idx -P ${N} -fold 1 > russia_sv.sfs 

angsd -vcf-gl orgpan_SV_morocco_bi_rehead.forangsd.vcf.gz -nind 8 -fai ${LINPAN}/${PREFIX}_panref2_sorted.fa.fai -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -dosaf 1 -doMajorMinor 5 -P ${N} -out morocco_sv
realSFS morocco_sv.saf.idx -P ${N} -fold 1 > morocco_sv.sfs 

realSFS mongolia_sv.saf.idx china_sv.saf.idx -P ${N} > mongolia.china_sv.ml # calculate the 2D SFS
realSFS mongolia_sv.saf.idx russia_sv.saf.idx -P ${N} > mongolia.russia_sv.ml
realSFS mongolia_sv.saf.idx morocco_sv.saf.idx -P ${N} > mongolia.morocco_sv.ml
realSFS china_sv.saf.idx russia_sv.saf.idx -P ${N} > china.russia_sv.ml
realSFS china_sv.saf.idx morocco_sv.saf.idx -P ${N} > china.morocco_sv.ml
realSFS russia_sv.saf.idx morocco_sv.saf.idx -P ${N} > russia.morocco_sv.ml

realSFS fst index mongolia_sv.saf.idx china_sv.saf.idx -sfs mongolia.china_sv.ml -fstout mongolia.china_sv # index sample so same sites are analyzed for each pop
realSFS fst index mongolia_sv.saf.idx russia_sv.saf.idx -sfs mongolia.russia_sv.ml -fstout mongolia.russia_sv
realSFS fst index mongolia_sv.saf.idx morocco_sv.saf.idx -sfs mongolia.morocco_sv.ml -fstout mongolia.morocco_sv
realSFS fst index china_sv.saf.idx russia_sv.saf.idx -sfs china.russia_sv.ml -fstout china.russia_sv
realSFS fst index china_sv.saf.idx morocco_sv.saf.idx -sfs china.morocco_sv.ml -fstout china.morocco_sv
realSFS fst index russia_sv.saf.idx morocco_sv.saf.idx -sfs russia.morocco_sv.ml -fstout russia.morocco_sv

realSFS fst stats mongolia.china_sv.fst.idx # global pairwise estimates # FST.Unweight[nObs:2265]:0.074479 Fst.Weight:0.082557
realSFS fst stats mongolia.russia_sv.fst.idx # FST.Unweight[nObs:2265]:0.056421 Fst.Weight:0.061847
realSFS fst stats mongolia.morocco_sv.fst.idx # FST.Unweight[nObs:2265]:0.055592 Fst.Weight:0.060582
realSFS fst stats china.russia_sv.fst.idx # FST.Unweight[nObs:2265]:0.055748 Fst.Weight:0.067572
realSFS fst stats china.morocco_sv.fst.idx # FST.Unweight[nObs:2265]:0.059210 Fst.Weight:0.067290
realSFS fst stats russia.morocco_sv.fst.idx # FST.Unweight[nObs:2265]:0.035650 Fst.Weight:0.037628

### PCAngsd 
conda activate pcangsd-1.2

angsd -bam ./bam.filelist -anc ${LINPAN}/${PREFIX}_panref2_sorted.fa -ref ${LINPAN}/${PREFIX}_panref2_sorted.fa -GL 1 -doglf 2 -minmaf 0.05 -snp_pval 1e-6 -dopost 1 -minQ 30 -minInd 20 -setMinDepth 100 -setMaxDepth 375 -docounts 1 -domajorminor 5 -domaf 1 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -C 50 -P ${N} -out snp_pca

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
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
meta <- as.matrix(read.table("/scratch/negishi/jeon96/swallow/original/swallow_metadata.txt", header=TRUE))
C <- as.matrix(read.table("orgpan_snp_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 7.22%, PC2: 4.36%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("orgpan_snp_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (7.22%)", y = "Principal Component 2 (4.36%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
dev.off()

C <- as.matrix(read.table("orgpan_sv_pca.cov"))
e <- eigen(C)
pca_result <- prcomp(C, scale = TRUE)
summary(pca_result) # PC1: 8.87%, PC2: 5.14%
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)
svg("orgpan_sv_pca.svg", width=6, height=5)
ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) + geom_point(size = 3) + labs(x = "Principal Component 1 (8.87%)", y = "Principal Component 2 (5.14%)") + theme_classic(base_size=14) + scale_color_manual(values = cb_palette)
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim <- OutFLANK(my_fst, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=4, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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
out_trim_comp <- OutFLANK(my_fst_comp, LeftTrimFraction = 0.1, RightTrimFraction = 0.1, NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1) # run outFLANK on pruned SNPs
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