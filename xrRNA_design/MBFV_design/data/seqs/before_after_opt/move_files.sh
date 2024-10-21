#!/bin/bash

# Set the source and destination directories
src_dir="/scr/aldea/kgutenbrunner/working/xrRNA_design/MBFV_design/data/seqs/before_after_opt/before_after_opt_2"
dest_dir="/scr/aldea/kgutenbrunner/working/xrRNA_design/MBFV_design/data/seqs/before_after_opt"

# Loop through all files matching the pattern in the source directory
for file in "$src_dir"/designs_before_after_opt_[1-9].csv "$src_dir"/designs_before_after_opt_[1-4][0-9].csv "$src_dir"/designs_before_after_opt_50.csv; do
    # Extract the base filename (without directory)
    filename=$(basename "$file")
    
    # Extract the number part (xz) from the filename
    num=$(echo "$filename" | grep -o -E '[0-9]+')

    # Add 50 to the number
    new_num=$((num + 50))

    # Form the new filename with the new number
    new_filename=$(echo "$filename" | sed "s/[0-9]\+/$new_num/")

    # Move the file to the destination with the new filename
    mv "$file" "$dest_dir/$new_filename"
    echo "Moved $file to $dest_dir/$new_filename"
done

