#!/bin/bash

#SBATCH --partition=Draco
#SBATCH --job-name=mapping-Klebs-phages
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=0GB
#SBATCH --time=10-0
#SBATCH -o slurm-%x-%j.out
#SBATCH --mail-type=END,FAIL,REQUEUE

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# Track the time it takes to run the script
SECONDS=0

# Load dependencies
module load anaconda3
eval "$(conda shell.bash hook)"

# CPUs per node
CPUS=$SLURM_CPUS_ON_NODE


echo "====================================================="
echo "Gathering stats                             : $(date)"
echo "====================================================="
# .fasta input file
input1=$1
input2=$2
ref_gen=$3
conda activate bwa
bwa index -p {basename $input} $input

bwa mem -t $CPUS $ref_gen $1 $2 | samtools sort -@2 -o $(basename $1 .fasta)_sorted.bam -
samtools index $(basename $1 .fasta)_sorted.bam
samtools flagstat $(basename $1 .fasta)_sorted.bam > $(basename $1 .fasta)_mapping_stats.txt

conda activate freebayes
freebayes -f $ref_gen -p 1 $(basename $1 .fasta)_sorted.bam > $(basename $1 .fasta)_FB.vcf

tabix -p vcf $(basename $1 .fasta)_FB.vcf

#igv
#snpeff


echo "===================================================="
echo "              ~Mapping in Progress~                 "
echo "===================================================="

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
echo ""


