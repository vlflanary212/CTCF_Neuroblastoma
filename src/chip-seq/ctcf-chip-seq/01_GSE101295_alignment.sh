#!/bin/bash

##### SLURM #####
#SBATCH --job-name=align
#SBATCH --partition=short
#SBATCH --time=03:30:00
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
fastq_dir=$dat/"fastq"  #fastq storage directory
bam_dir=$dat/"bam"  #bam storage directory
bigwig_dir=$dat/"bigwig" #bigwig storage directory

# files
ref_genome="/data/project/sen-lab/genome/hg38/bwa/Homo_sapiens_assembly38.fasta"  #reference genome
fastq="$fastq_dir/$sample".fastq  #individual fastq files

##### PACKAGES ######
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.9-foss-2018b
module load deepTools/3.3.1-foss-2018b-Python-3.6.6

##### COMMANDS #####
# map reads and make bam files
bwa mem -t 10 $ref_genome $fastq | samtools view -S -b -h -F 4 -q 20 > $bam_dir/"$sample"_mapped.bam

# sort by coordinates
samtools sort -@ 8 -o $bam_dir/"$sample"_sort.bam $bam_dir/"$sample"_mapped.bam

# remove duplicates
samtools rmdup -s $bam_dir/"$sample"_sort.bam $bam_dir/"$sample"_final.bam

# index 
samtools index $bam_dir/"$sample"_final.bam

# create bigwig file
bamCoverage -b $bam_dir/"$sample"_final.bam -o $bigwig_dir/"$sample".bw --binSize 200 --normalizeUsing RPKM --effectiveGenomeSize 2913022398 

##### END #####
echo "done"
