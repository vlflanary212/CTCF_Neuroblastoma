#!/bin/bash

##### SLURM #####
#SBATCH --job-name=peaks
#SBATCH --partition=express
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=1-6
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

##### VARIABLES #####
# directories
dat=$wd/"dat/chip-seq/GSE101295"  #data storage directory
bam_dir=$dat/"bam"  #bam storage directory
peak_dir=$dat/"peaks"  #peaks storage directory

##### PACKAGES #####
module load MACS2/2.1.1.20160309-foss-2016a-Python-2.7.11

##### COMMANDS #####
macs2 callpeak --treatment $bam_dir/"$sample"_final.bam --name "$sample" --outdir $peak_dir --gsize hs --nomodel --shift 0 --extsize 200 --keep-dup all --call-summits --pvalue 0.01 

##### END #####
echo "done"
