#!/bin/bash

#SBATCH --partition=Draco
#SBATCH --job-name=MC-checkV-metaome-Klebs
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

for i in ./phages/*.fasta
do
	module load anaconda3
	conda activate /projects/raw_lab/envs/MetaCerberus-1.4.0
	metacerberus.py --nanopore --super $i --dir-out ./phages/results/metacerberus/$i-MC --db-path /projects/raw_lab/databases/MetaCerberus/
	conda deactivate
	echo "MetaCerberus has completed $i"
	conda activate /projects/raw_lab/envs/binning/checkv
	checkv end_to_end $i ./phages/results/checkv/$i-checkv -d /projects/raw_lab/envs/binning/checkv/checkv-db-v1.5 -t 36
	echo "CheckV has completed $i"
done

echo "====================================================="
echo "Gathering stats                             : $(date)"
echo "====================================================="

echo "===================================================="
echo "              Running ~MetaCerberus~                "
echo "===================================================="



echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
echo ""

