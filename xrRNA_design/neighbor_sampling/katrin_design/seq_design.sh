#!/bin/bash
#SBATCH --job-name=neighbor_sampling_design
#SBATCH --time=5:00:00    # Set maximum wall time to 1 hour
#SBATCH --mem=4G          # Request 4GB of memory
#SBATCH --mail-user=katrin.gutenbrunner@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/scr/aldea/kgutenbrunner/working/SLURM/SLURM_jobs/neighbor_sampling_design.%j.out   # Set the name for the standard output file
#SBATCH --error=/scr/aldea/kgutenbrunner/working/SLURM/SLURM_jobs/neighbor_sampling_design.%j.err    # Set the name for the standard error file
#SBATCH --cpus-per-task=16

WORK_DIR=/scr/aldea/kgutenbrunner/working/xrRNA_design/neighbor_sampling/katrin_design/
PYTHON_SCRIPT=$WORK_DIR/design.py
OUTPUT_FILE=$WORK_DIR/design_output.fa
cd $WORK_DIR
python $PYTHON_SCRIPT > $OUTPUT_FILE

