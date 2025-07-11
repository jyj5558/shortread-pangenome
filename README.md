# shortread-pangenome

## The main directory includes individual scripts for each step with appropriate orders.

## The sub Nextflow directory includes modified individual scripts to work in the nextflow pipeline, nextflow scripts, and other necessary files.
To use the automated nextflow pipe to build short-read pangenomes:

(1) Modify 'params.csv' and 'nextflow.config' files accordingly. 

(2) Choose what short-read pangenome you will build (e.g., MC Pangenome) then check the main .nf file (e.g., main_ShortreadPangenome_MCPan.nf). Fill in or modify directory paths (e.g., 'params.csv = ""', 'cd ~') and clusterOptions according to your computing system.

(3) In the 'bin' directory, there are bash scripts which will be called by the nextflow script. Check all the relevant ones for your desired pangenome and modify how to load or call each program in the bash scripts.

Contact: jeon96@purdue.edu for any inquiries.
