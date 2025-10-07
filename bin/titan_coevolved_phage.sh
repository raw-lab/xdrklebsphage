#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=titan_Klebs_phages
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400GB
#SBATCH --time=3-0
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
conda activate /projects/raw_lab/envs/titanomics-0.2


#CPUs per node
CPUS=$SLURM_CPUS_ON_NODE

echo "===================================================="
echo "              Running ~Titan~                "
echo "===================================================="

titanomics.py --cpus $CPUS --illumina --super /projects/raw_lab/projects/Klebs/Klebsiella-phage-coevolution/TrnPhage*/* --dir_out /projects/raw_lab/projects/Klebs/coevolved_phage_results/titan/

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
echo ""

