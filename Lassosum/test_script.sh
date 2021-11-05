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
# Submit this bash script to the CLUSTER by "qsub -cwd test_script.sh"

# --bfile: Individual binary files
# --extract: # SNP id to keep (rs2013030)
# --keep: Individual id to keep (familiy ID; individual ID) 
# --score: 52: column number of SNP id (rs2013030)  53: effect allele A1 (A/T/C/G)
# **IMPORTANT CHANGE --score-col-nums: 11: column of the best-fit effect sizes (calculated betas). CHANGE this to index_max from tuning Rscript.
# --out: Output: polygenic risk scores. Rows: invidual. Columns: effect size version.

/dcl01/chatterj/data/jialan/Tool/plink2 \
--threads 6 \
--rm-dup exclude-all \
--bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3 \
--extract /dcs04/nilanjan/data/jialan/1000G/EUR/snp_plink_chr${SGE_TASK_ID}.txt \
--keep /dcs04/nilanjan/data/jialan/Genotype/test_id_14.txt \ 
--score /dcs04/nilanjan/data/jialan/1000G/EUR/plink_beta/beta-use-chr${SGE_TASK_ID}.txt 52 53 header \
ignore-dup-ids \
--score-col-nums 11 \
--out /dcs04/nilanjan/data/jialan/LDPred2/prs_test/prs_test_chr${SGE_TASK_ID}

/dcl01/chatterj/data/jialan/Tool/plink2 \
--threads 6 \
--rm-dup exclude-all \
--bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3 \
--extract /dcs04/nilanjan/data/jialan/1000G/EUR/snp_plink_chr${SGE_TASK_ID}.txt \
--keep /dcs04/nilanjan/data/jialan/Genotype/val_id_14.txt \
--score /dcs04/nilanjan/data/jialan/1000G/EUR/plink_beta/beta-use-chr${SGE_TASK_ID}.txt 52 53 header \
ignore-dup-ids \
--score-col-nums 11 \
--out /dcs04/nilanjan/data/jialan/LDPred2/prs_val/prs_val_chr${SGE_TASK_ID}