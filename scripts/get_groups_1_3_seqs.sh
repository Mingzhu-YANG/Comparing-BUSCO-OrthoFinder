#!/bin/bash

# Define directory
SOURCE_DIR="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/5_Orthogroup_Sequences_63451"   
TARGET_DIR="$SOURCE_DIR/OF_1_3_seqs"

# Create target dir (if not exist)
mkdir -p "$TARGET_DIR"

# Loop through all the fa file
for file in "$SOURCE_DIR"/*.fa; do
    # Caculate the number of seqs in each file
    count=$(grep -c "^>" "$file")

    # Check number of seqs, whether they are in between the defined number
    if [ "$count" -ge 1 ] && [ "$count" -le 3 ]; then
        # Move file to target folder
        mv "$file" "$TARGET_DIR/"
        # echo "移动文件: $file"
    fi
done

echo "Done！"
