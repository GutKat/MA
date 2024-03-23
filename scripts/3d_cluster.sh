#!/bin/bash
#SBATCH --job-name=cluster_name
#SBATCH --time=5:00:00    # Set maximum wall time to 1 hour
#SBATCH --mem=4G          # Request 4GB of memory
#SBATCH --mail-user=katrin.gutenbrunner@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/scr/aldea/kgutenbrunner/SLURM/SLURM_jobs/cluster_name.%j.out   # Set the name for the standard output file
#SBATCH --error=/scr/aldea/kgutenbrunner/SLURM/SLURM_jobs/cluster_name.%j.err    # Set the name for the standard error file
#SBATCH --cpus-per-task=16

# Load any required modules here (if needed)
# module load python

# Define the paths and filenames
WORK_DIR=/scr/aldea/kgutenbrunner/simRNA/name

# Move to the working directory
cd $WORK_DIR

conda activate ir_env

./status.sh
./cluster.sh
./pp_cluster3.sh
