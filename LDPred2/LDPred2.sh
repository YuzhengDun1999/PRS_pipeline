#!/usr/bin/env bash
#$ -N mb_ref
#$ -cwd
#$ -pe local 6 
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jma59@jhu.edu
#$ -t 1-22

# Change name after -N.
# Change number of nodes after -pe local. Not necessary to change.
# Change memory size request after -l. Since we request 6 cores, for each task, we have 6*10=60G.
# Change email address after -M.

# ----------------------- Linux commands -----------------------
module load conda_R
Rscript /dcs04/nilanjan/data/jialan/LDPred2/LDPred2.R ${SGE_TASK_ID} 1  
# Change directory here. ${SGE_TASK_ID} do not need to be changed. Example: Rscript /x/xx/xxx/LDPred2.R ${SGE_TASK_ID}. This specifies chromosome 1-22.
# 1 is the number k in R script, which specifies race.
# races = c('EUR','AFR','ASIAN','AMR'). 
# For example, if your race is Africa, change the number to 2. That is,
# Rscript /dcs04/nilanjan/data/jialan/LDPred2/LDPred2.R ${SGE_TASK_ID} 2