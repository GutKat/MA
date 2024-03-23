#!/bin/bash
#SBATCH --job-name=qrna_name
#SBATCH --time=15:00:00    # Set maximum wall time to 5 hour
#SBATCH --mem=8G          # Request 8GB of memory
#SBATCH --mail-user=katrin.gutenbrunner@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/scr/aldea/kgutenbrunner/simRNA/name/qrna_name_%j.out
#SBATCH --error=/scr/aldea/kgutenbrunner/simRNA/name/qrna_name_%j.err
#SBATCH --cpus-per-task=16

SIMNAME="name"
BASEDIR="/scr/aldea/kgutenbrunner/simRNA/name"


cd $BASEDIR
QRNA="${BASEDIR}/QRNAs"

if [ ! -d $QRNA ]; then
        mkdir $QRNA
        echo "QRNAs folder was created"
fi

# Get inputs using ls command
inputs=($(ls ${BASEDIR}/all_thrs5.00A_clust*/all_thrs5.00A_clust*-000001_AA.pdb))

# Loop through each input
for input in "${inputs[@]}"; do

    # Extract cluster name and output_name
    cluster_name=$(echo "$input" | sed 's/.*_\(clust[0-9][0-9]\).*/\1/')
    output_file="${QRNA}/${cluster_name}_qrna.pdb"

    echo "Running QRNA command with input: $input; output: $output_file"

    # Run the QRNA command for 20 minutes
    timeout 20m QRNA -i "$input" -o "$output_file"

    # Check if timeout occurred
    if [ $? -eq 124 ]; then
        echo "Timeout occurred. Continuing to the next input."
    fi

done

echo "All iterations completed."
