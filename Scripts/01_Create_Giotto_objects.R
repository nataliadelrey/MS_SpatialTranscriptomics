# Creating Giotto objects in R
# Author: Natalia del Rey Díez
# Date: 10/03/2024

################################################################################
#------------------------------ Install packages ------------------------------#
################################################################################
# Giotto
if(!"pak" %in% installed.packages()) {
  install.packages("pak")
}

if(!"Giotto" %in% installed.packages()) {
  pak::pkg_install("drieslab/Giotto")
}

## Python environment for Giotto
genv_exists <- Giotto::checkGiottoEnvironment()
if(!genv_exists){
  Giotto::installGiottoEnvironment()
}

# ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")}

# patchwork
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")}


################################################################################
#----------------------------- Loading packages -------------------------------#
################################################################################
library(Giotto)
library(ggplot2)
library(patchwork)


################################################################################
#--------------------------- Set working directory ----------------------------#
################################################################################
GSE279181_dir <- "../ST_GSE279181/GSE279181_SpaceRanger"

# Sample subfolders
samples <- list.dirs(GSE279181_dir, full.names = TRUE, recursive = FALSE)


################################################################################
#--------------------------- Create Giotto objects ----------------------------#
################################################################################
giotto_objects <- list()

for (sample in samples) {
  # Directory with Space Ranger outputs
  visium_dir <- file.path(sample, "outs")
  
  # Create the Giotto object
  giotto_object <- createGiottoVisiumObject(visium_dir = visium_dir, 
                                            expr_data = "filter", 
                                            gene_column_index = 2) # GENE SYMBOL
  
  sample_name <- strsplit(basename(sample), "_")[[1]][2]  
  giotto_objects[[sample_name]] <- giotto_object
}


################################################################################
#---------------------------- Add sample metadata -----------------------------#
################################################################################
# Read file with metadata of all samples
metadata <- read.csv(file = "../sample_metadata.tsv",
                     header = TRUE, sep = "\t")

for (obj_name in names(giotto_objects)) { 
  metadata_filtered <- metadata[metadata$Biomaterial.ID == obj_name, ]
  
  g <- giotto_objects[[obj_name]]
  
  m1 <- getCellMetadata(g, output = "data.table")
  m2 <- data.frame(m1, metadata_filtered)
  mobj = createCellMetaObj(m2)
  
  g <- Giotto::setCellMetadata(gobject = g, x = mobj)
  giotto_objects[[obj_name]] <- g
}


################################################################################
#---------------- Subset on spots that were covered by tissue -----------------#
################################################################################
giotto_objects <- lapply(giotto_objects, function(g) {
  Giotto::subsetGiotto(g,
                       cell_ids = pDataDT(g)[in_tissue == 1]$cell_ID)
})

# Plots
MS497T <- Giotto::spatPlot2D(gobject = giotto_objects$MS497T, 
                             show_image = TRUE, 
                             point_size = 2,
                             point_alpha = 0.4, 
                             title = "MS497T")

MS549T <- Giotto::spatPlot2D(gobject = giotto_objects$MS549T, 
                             show_image = TRUE, 
                             point_size = 2,
                             point_alpha = 0.5,
                             title = "MS549T")

spotsvacios = MS497T + plot_spacer() + MS549T + 
  plot_layout(widths = c(1, 0.4, 1))

ggsave("../Graficos/spotsvacios.png",
       spotsvacios,
       width = 12,
       height = 5,
       dpi = 300)


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
for (name in names(giotto_objects)) {
  Giotto::saveGiotto(gobject = giotto_objects[[name]], 
                     foldername = name, 
                     dir = "../RData/ST_GSE279181/01_GiottoObjects",
                     method = "RDS", 
                     overwrite = TRUE, 
                     export_image = TRUE,
                     image_filetype = "PNG",
                     include_feat_coord = TRUE,
                     verbose = TRUE)}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#