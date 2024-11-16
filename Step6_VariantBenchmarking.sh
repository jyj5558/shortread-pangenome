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


# Variant benchmarking (modified based on Liao et al. (2023) in Nature "A draft human pangenome reference")
module purge
module load biocontainers
module load bcftools
module load picard 
module load anaconda
module load use.own
module load conda-env/mypackages-py3.9.13
module load bbmap
module load samtools

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
cp ${ORGPAN}/called/sv/vg/barnswallow_orgSV.decomp.annot.ext.vcf  ./orgpan_vg_final.vcf

### Updating vcf headers and sorting variants following the order of linear pangenome contigs 

python3 # run the code below in python
def extract_contigs_from_vcf(vcf_file):
    contigs_dict = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##contig='):
                parts = line.strip().split('<')[1].split('>')[0]
                contig_id = None
                contig_length = None
                for part in parts.split(','):
                    if part.startswith('ID='):
                        contig_id = part.split('=')[1]
                    elif part.startswith('length='):
                        contig_length = int(part.split('=')[1])
                if contig_id and contig_length:
                    contigs_dict[contig_id] = contig_length           
    return contigs_dict

vcf_file="orgpan_vg_final.vcf"
contigs = extract_contigs_from_vcf(vcf_file)
print("Contigs:", contigs)

with open("orgSV_chroms.txt", 'w') as f:
    for key, value in contigs.items():
        f.write(str(key) + " " + str(value) + '\n') 

exit()

sort -k2,2 -n orgSV_chroms.txt > sorted_orgSV_chroms.txt
cut -f1,2 ${REF}/${GENOME}.fna.fai > ref_chroms.txt
sort -k2,2 -n ref_chroms.txt > sorted_ref_chroms.txt
join -1 2 -2 2 sorted_orgSV_chroms.txt sorted_ref_chroms.txt > combined_chroms.txt # check if there are duplicates and if so manually curate them (in this case, duplicates at lengths of "18466", "25102", and "65965")
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
