#!/bin/bash
# Download data from GEO
# Author: Natalia del Rey Díez
# Date: 04/03/2024

################################################################################
#--------------------------- Download data from GEO ---------------------------#
################################################################################
# Spatial transcriptomics
mkdir -p ../ST_GSE279181
url_ST="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279181/suppl/GSE279181_RAW.tar"
wget -c "$url_ST" -O ../ST_GSE279181/GSE279181_RAW.tar

# snRNA-Seq
mkdir -p ../snRNA_GSE279180
url_snRNA="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279180/suppl/GSE279180_sn_atlas.h5ad"
wget -c "$url_snRNA" -O ../snRNA_GSE279180/GSE279180_sn_atlas.h5ad


################################################################################
#--------------------------- Unzip downloaded files ---------------------------#
################################################################################
# Spatial transcriptomics
mkdir -p ../ST_GSE279181/GSE279181_RAW
tar -xvf ../ST_GSE279181/GSE279181_RAW.tar -C ../ST_GSE279181/GSE279181_RAW
gunzip -d ../ST_GSE279181/GSE279181_RAW/*.gz
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#