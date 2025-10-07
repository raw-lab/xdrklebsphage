#!/bin/bash

#SBATCH --partition=Draco
#SBATCH --job-name=metaome-Klebs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=450GB
#SBATCH --time=15-00
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

# CPUs per node
CPUS=$SLURM_CPUS_ON_NODE
eval "$(conda shell.bash hook)"
module load anaconda3
for i in ./phages/*.fasta
do
	conda activate /projects/raw_lab/envs/binning/metaomestats
	countAssembly.py -i 100 --fasta $i
	echo "MetaomeStats has completed with $i"
done


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
echo ""

