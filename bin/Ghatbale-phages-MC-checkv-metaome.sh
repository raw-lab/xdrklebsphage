#!/bin/bash

#SBATCH --partition=Draco
#SBATCH --job-name=MC-checkV-metaome-Klebs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250GB
#SBATCH --time=1-08
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

for i in ./Klebsiella-phages-Ghatbale2024/*.fasta
do
	module load anaconda3
	conda activate /projects/raw_lab/envs/MetaCerberus-1.4.0
	metacerberus.py --illumina --super $i --dir-out ./Ghatbale_phages_results/metacerberus/$(basename $i .fasta)-MC --db-path /projects/raw_lab/databases/MetaCerberus/
	conda deactivate
	echo "MetaCerberus has completed $i"
	conda activate /projects/raw_lab/envs/binning/checkv
	checkv end_to_end $i ./Ghatbale_phages_results/checkv/$(basename $i .fasta)-checkv -d /projects/raw_lab/envs/binning/checkv/checkv-db-v1.5 -t 36
	echo "CheckV has completed $i"
	conda deactivate
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

