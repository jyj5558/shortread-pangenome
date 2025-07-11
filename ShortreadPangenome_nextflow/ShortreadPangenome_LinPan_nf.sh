#!/bin/bash
#SBATCH -A standby
#SBATCH -t 0-00:10:00

# Loads Nextflow
module load nextflow

# Change to a good place to store the nextflow 'work' directory with logged info
cd /scratch/negishi/jeon96/shortread_pan/

# Run Nextflow
nextflow run main_ShortreadPangenome_LinPan.nf  

# !! PLEASE LOOK AT EVERY SCRIPT THOROUGHLY BEFORE RUNNING THE NEXTFLOW PIPELINE AND MODIFY THEM ACCORDINGLY !!"
# SCRIPTS HAVE BEEN WRITTEN TO BE RUN UNDER PURDUE CLUSTER SYSTEM
