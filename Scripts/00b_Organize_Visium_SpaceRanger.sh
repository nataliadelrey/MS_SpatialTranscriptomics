#!/bin/bash
# Organizing data in the Space Ranger directory structure
# Author: Natalia del Rey Díez
# Date: 05/03/2024

################################################################################
#--------------------------- Set working directory ----------------------------#
################################################################################
# Directory with ST data downloaded from GEO
ST_dir="../ST_GSE279181/GSE279181_RAW"

# Create output directory
output_dir="../ST_GSE279181/GSE279181_SpaceRanger"
mkdir -p ../ST_GSE279181/GSE279181_SpaceRanger


################################################################################
#----------------------- Organize downloaded GEO files ------------------------#
################################################################################
for file in "$ST_dir"/*; do
    if [[ -f "$file" ]]; then
       
        filename=$(basename "$file")
        
        # Extract the sample prefix (GSMxxxx_COxx or GSMxxxx_MSxx)
        if [[ $filename =~ ^(GSM[0-9]+_(CO|MS)[0-9]+[A-Z]?)_ ]]; then
            sample="${BASH_REMATCH[1]}"
            sample_dir="$output_dir/$sample/outs"

            mkdir -p "$sample_dir/spatial"
            mkdir -p "$sample_dir/filtered_feature_bc_matrix"

            case "$file" in
            # Copy files to the spatial folder
            *"tissue_hires_image.png"*) cp "$file" "$sample_dir/spatial/tissue_hires_image.png";;
            *"tissue_lowres_image.png"*) cp "$file" "$sample_dir/spatial/tissue_lowres_image.png";;
            *"aligned_fiducials.jpg"*) cp "$file" "$sample_dir/spatial/aligned_fiducials.jpg";;
            *"detected_tissue_image.jpg"*) cp "$file" "$sample_dir/spatial/detected_tissue_image.jpg";;
            *"scalefactors_json.json"*) cp "$file" "$sample_dir/spatial/scalefactors_json.json";;
            *"tissue_positions_list.csv"*) cp "$file" "$sample_dir/spatial/tissue_positions_list.csv";;
            *"spatial_enrichment.csv"*) cp "$file" "$sample_dir/spatial/spatial_enrichment.csv";;

            # Copy files to the filtered_feature_bc_matrix folder
            *"barcodes.tsv"*) cp "$file" "$sample_dir/filtered_feature_bc_matrix/barcodes.tsv";;
            *"features.tsv"*) cp "$file" "$sample_dir/filtered_feature_bc_matrix/features.tsv";;
            *"matrix.mtx"*) cp "$file" "$sample_dir/filtered_feature_bc_matrix/matrix.mtx";;
            esac
        fi
    fi
done

echo "Organization completed in: $output_dir"
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#