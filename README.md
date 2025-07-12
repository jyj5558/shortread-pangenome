# Short-read Pangenome
The scripts are optimized for paired-end short reads to be run in Purdue's RCAC cluster system. 
Details may be different from your computing system. You need to check your system's requirement and modify accordingly.

Step 0 was used for generating and processing the benchmarking data for the paper (TBA).

Steps 1-3 are used to build each short-read pangenome (2 for linear pangenome, then from there, 3A for MC graph pangenome or 3B for VG graph pangenome).

Steps 4-7 were used for variant calling, assembly comparison, variant benchmarking, and population genetic analyses using an independent empirical test dataset in the same paper.


The revision_script subdirectory includes scripts which were used during the manuscript revision. 
Using these, 'pre-filtering' and 'post-filtering' variants were newly generated (without pangenome augmentation using an independent empirical test dataset) and a couple of simple window-based genetic analyses were added.

## The main directory includes individual scripts for each step with appropriate orders.
To use each individual script, modify directory paths, slurm options, and how to load or call each program in the bash scripts

## The sub Nextflow directory includes modified individual scripts to work in the nextflow pipeline, nextflow scripts, and other necessary files.
To use the automated nextflow pipe to build short-read pangenomes:

(1) Modify 'params.csv' and 'nextflow.config' files accordingly. 

(2) Choose what short-read pangenome you will build (e.g., MC Pangenome) then check the main .nf file (e.g., main_ShortreadPangenome_MCPan.nf) and nf.sh file (e.g., ShortreadPangenome_MCPan_nf.sh). Fill in or modify directory paths (e.g., 'params.csv = ""', 'cd ~') and clusterOptions according to your computing system.

(3) In the 'bin' directory, there are bash scripts which will be called by the nextflow script. Check all the relevant ones for your desired pangenome and modify how to load or call each program in the bash scripts.

(4) run the nf.sh file; e.g., "sbatch ShortreadPangenome_MCPan_nf.sh"

Contact: jeon96@purdue.edu for any inquiries.
