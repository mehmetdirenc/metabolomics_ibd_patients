#!/bin/bash

# Define the quantification directory
quantification_dir="/mnt/lustre/home/mager/magmu818/results/mahana/ms/2024_results/quantification"

# Loop over each folder in the quantification directory
for folder in "$quantification_dir"/*; do
    # Ensure it's a directory
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder") # Get the folder name

        # Find the .featureXML file in the folder
        feature_file=$(find "$folder/output" -name "*.featureXML")

        if [ -n "$feature_file" ]; then
            sample_name=$(basename "$feature_file" .featureXML) # Extract sample_name

            # Define the output paths
            output_mzTab="/mnt/lustre/home/mager/magmu818/results/mahana/ms/2024_results/compounds/${sample_name}.mzTab"
            output_annotation="/mnt/lustre/home/mager/magmu818/results/mahana/ms/2024_results/compounds/${sample_name}.featureXML"

            # Construct and execute the AccurateMassSearch command
            AccurateMassSearch \
                -in "$feature_file" \
                -out "$output_mzTab" \
                -out_annotation "$output_annotation" \
                -algorithm:ionization_mode negative
        else
            echo "No .featureXML file found in $folder/output"
        fi
    fi
done