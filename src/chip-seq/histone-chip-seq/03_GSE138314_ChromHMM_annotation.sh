#!/bin/bash 

##### SLURM #####
#SBATCH --job-name=chromhmm
#SBATCH --partition=long
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

##### WORKING DIRECTORY #####
wd="/data/project/sen-lab/CTCF_Neuroblastoma"

# directories
bam="/data/project/sen-lab/CTCF_Neuroblastoma/dat/chip-seq/GSE138314/bam"
in="/data/user/home/bajones1/Documents/ChromHMM"
bin="$in/GSE138314/binarized"
out="$in/GSE138314/learnModelOut"

# cellmarkfiletable
list="$in/GSE138314/cellmarkfiletable.txt"

# load module and conda environment with ChromHmm installed
module load Anaconda3
conda activate env

# go to ChromHmm directory
cd $in

# Binarize bam files
java -jar ChromHMM.jar BinarizeBam CHROMSIZES/hg38.txt $bam $list $bin

# Learn model
java -jar ChromHMM.jar LearnModel -noautoopen -p 4 $bin $out 10 hg38


# Finish 
echo "done" 

