#!/bin/bash 
#SBATCH -N 1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=4           
#SBATCH --mem=120Gb                 
#SBATCH -t 6:00:00                   
#SBATCH --partition=short        
#SBATCH -e slurm-%j.err 
#SBATCH --output slurm-%j.out 


# directories
bam="/data/project/sen-lab/CTCF_Neuroblastoma/dat/chip-seq/combined/bam"
in="/data/user/home/bajones1/Documents/ChromHMM"
bin="$in/SK_N_AS/binarized"
out="$in/SK_N_AS/learnModelOut"

# cellmarkfiletable
list="$in/SK_N_AS/cellmarkfiletable.txt"

# load module and conda environment with ChromHmm installed
module load Anaconda3
conda activate env

# go to ChromHmm directory
cd $in

# Binarize bam files
java -jar ChromHMM.jar BinarizeBam CHROMSIZES/hg38.txt $bam $list $bin

# Learn model
java -jar ChromHMM.jar LearnModel -noautoopen -p 8 $bin $out 10 hg38


# Finish 
echo "done" 

