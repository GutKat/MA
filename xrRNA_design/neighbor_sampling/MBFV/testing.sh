#!/bin/bash
#SBATCH --job-name=neighbor_sampling
#SBATCH --time=24:00:00    # Set maximum wall time to 24 hour
#SBATCH --mem=2G          # Request 2GB of memory
#SBATCH --mail-user=katrin.gutenbrunner@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/scr/aldea/kgutenbrunner/working/SLURM/SLURM_jobs/neighbor_sampling.%j.out   # Set the name for the standard output file
#SBATCH --error=/scr/aldea/kgutenbrunner/working/SLURM/SLURM_jobs/neighbor_sampling.%j.err    # Set the name for the standard error file
#SBATCH --cpus-per-task=16

WORK_DIR=/scr/aldea/kgutenbrunner/working/xrRNA_design/neighbor_sampling/katrin_design
PYTHON_SCRIPT=$WORK_DIR/design.py

# Move to the working directory
cd $WORK_DIR
python $PYTHON_SCRIPT

