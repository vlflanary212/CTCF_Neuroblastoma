#!/bin/bash

##### SLURM #####
#SBATCH --job-name=align
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=1-2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########### AUTHORS #########
# Victoria L. Flanary and Brianna A. Jones

##### WORKING DIRECTORY #####
wd="/data/project/sen-lab/CTCF_Neuroblastoma" 

##### ARRAY #####
sample_list=$wd"/doc/atac-seq/GSE108517/SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$sample_list")
echo "sample:$sample"

##### VARIABLES #####
# directories
dat=$wd/"dat/atac-seq/GSE108517"  #data storage directory
fastq_dir=$dat/"fastq"  #fastq storage directory
bam_dir=$dat/"bam"  #bam storage directory
bigwig_dir=$dat/"bigwig" #bigwig storage directory

# reference genome
ref_genome="/data/project/sen-lab/genome/hg38/bwa/Homo_sapiens_assembly38.fasta"

# fastq files (2/sample due to paired-end sequencing)
fastq1="$fastq_dir/$sample"_1.fastq 
fastq2="$fastq_dir/$sample"_2.fastq

##### PACKAGES ######
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.9-foss-2018b
module load deepTools/3.3.1-foss-2018b-Python-3.6.6

##### COMMANDS #####
# map reads and make bam files
bwa mem -M -t 10 $ref_genome $fastq1 $fastq2 | samtools view -S -b -h -F 4 -q 20 > $bam_dir/"$sample"_mapped.bam

# sort by coordinates
samtools sort -@ 8 -o $bam_dir/"$sample"_sort.bam $bam_dir/"$sample"_mapped.bam

# remove duplicates
samtools rmdup -S $bam_dir/"$sample"_sort.bam $bam_dir/"$sample"_rmdup.bam

# index 
samtools index $bam_dir/"$sample"_rmdup.bam

# remove mitochondrial reads
samtools idxstats $bam_dir/"$sample"_rmdup.bam | cut -f 1 | grep -v chrM | xargs samtools view -b $bam_dir/"$sample"_rmdup.bam -> $bam_dir/"$sample"_final.bam

# index again
samtools index $bam_dir/"$sample"_final.bam

# create bigwig file
bamCoverage -b $bam_dir/"$sample"_final.bam -o $bigwig_dir/"$sample".bw --binSize 200 --normalizeUsing None --effectiveGenomeSize 2913022398 

##### END #####
echo "done"
