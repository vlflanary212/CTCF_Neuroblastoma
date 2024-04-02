#!/bin/bash

##### SLURM #####
#SBATCH --job-name=GSE101295
#SBATCH --partition=express
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3
#SBATCH --mem=128G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########### AUTHORS #########
# Victoria L. Flanary

##### WORKING DIRECTORY #####
wd="/data/project/sen-lab/CTCF_Neuroblastoma"

##### ARRAY #####
sample_list=$wd"/doc/chip-seq/GSE101295/SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$sample_list")
echo "sample:$sample"

##### PACKAGES ######
module load SRA-Toolkit/3.0.0-centos_linux64

##### COMMANDS #####
# download fastq files for specified samples from the SRR_Acc_List.txt
fasterq-dump $sample --outdir $wd/"dat/chip-seq/GSE101295/fastq/" 

##### END #####
echo "done"
