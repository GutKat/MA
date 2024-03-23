#!/bin/bash
#SBATCH --job-name=seq_design_mc_500000
#SBATCH --time=5:00:00    # Set maximum wall time to 1 hour
#SBATCH --mem=4G          # Request 4GB of memory
#SBATCH --mail-user=katrin.gutenbrunner@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/scr/aldea/kgutenbrunner/SLURM/SLURM_jobs/seq_design_mc_500000.%j.out   # Set the name for the standard output file
#SBATCH --error=/scr/aldea/kgutenbrunner/SLURM/SLURM_jobs/seq_design_mc_500000.%j.err    # Set the name for the standard error file
#SBATCH --cpus-per-task=16

# Load any required modules here (if needed)
# module load python

# Define the paths and filenames
WORK_DIR=/scr/aldea/kgutenbrunner/xrRNA_design/TBFV_design
PYTHON_SCRIPT=$WORK_DIR/design.py
OUTPUT_FILE=$WORK_DIR/seqs/design_output.fa

# Move to the working directory
cd $WORK_DIR

# Run the Python script
python $PYTHON_SCRIPT > $OUTPUT_FILE
