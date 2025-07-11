#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH -t 14-00:00:00

# Loads Nextflow
module load nextflow

# Change to a good place to store the nextflow 'work' directory with logged info
cd /scratch/negishi/allen715/shortread_pan/

# Run Nextflow
# nextflow run main_ShortreadPangenome_MCPan.nf -with-report report.html -with-timeline timeline.html -with-trace

nextflow run main_ShortreadPangenome_MCPan.nf -entry step2_only -with-report report.html -with-timeline timeline.html -with-trace


# -resume
# !! PLEASE LOOK AT EVERY SCRIPT THOROUGHLY BEFORE RUNNING THE NEXTFLOW PIPELINE AND MODIFY THEM ACCORDINGLY !!"
# SCRIPTS HAVE BEEN WRITTEN TO BE RUN UNDER PURDUE CLUSTER SYSTEM
