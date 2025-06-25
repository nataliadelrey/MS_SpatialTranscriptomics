# ST Deconvolution
# Author: Natalia del Rey Díez
# Date: 24/04/2024


################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)
library(patchwork)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
# ST
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/06_Leiden/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)

# snRNA-Seq
giotto_SC <- Giotto::loadGiotto(path_to_folder = "../RData/snRNA_GSE279180/02_Processed/giotto_SC",
                                reconnect_giottoImage = TRUE,
                                init_gobject = TRUE,
                                verbose = TRUE)

# Marker Genes
load("../RData/ST_GSE279181/06_Leiden/scran_markers.rda")


################################################################################
#---------------------- Signature Matrix for SpatialDWLS ----------------------#
################################################################################
# Top 50 marker genes per cluster
topgenes_scran <- scran_markers[, head(.SD, 50), by = "cluster"]$feats


DWLS_matrix = Giotto::makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(giotto_SC,
                         values = "normalized",
                         output = "matrix"),
  cell_type = pDataDT(giotto_SC)$celltype,
  sign_gene = topgenes_scran)
  

################################################################################
#------------------------------- Deconvolution --------------------------------#
################################################################################
combo_filt = Giotto::runDWLSDeconv(combo_filt,
                                   expression_values = "default",
                                   cluster_column = "leiden_kNN.050",
                                   n_cell = 5,
                                   sign_matrix = DWLS_matrix)


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/07_Deconvolucion",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
cell_color <- Giotto::getDistinctColors(9)
names(cell_color) <- unique(pDataDT(giotto_SC)$celltype)

# CO37
CO37 = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(combo_filt)[list_ID == "CO37"]$cell_ID)
Giotto::spatDeconvPlot(CO37,  
                       show_image = FALSE,
                       radius = 32,
                       cell_color_code = cell_color,
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "Deconv/spatDeconv",
                                         save_name = "CO37",
                                         base_width = 5,   
                                         base_height = 6,  
                                         units = "in",      
                                         dpi = 300))
Giotto::spatCellPlot(CO37,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.5,
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "CO37",                            
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))

# CO85
CO85 = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(combo_filt)[list_ID == "CO85"]$cell_ID)
Giotto::spatDeconvPlot(CO85,  
                       show_image = FALSE,
                       radius = 32,          
                       cell_color_code = cell_color,
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "Deconv/spatDeconv",
                                         save_name = "CO85",
                                         base_width = 5,   
                                         base_height = 6,  
                                         units = "in",      
                                         dpi = 300))
Giotto::spatCellPlot(CO85,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.5, 
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "CO85",                                    
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))

# MS197U
MS197U = Giotto::subsetGiotto(combo_filt, 
                             cell_ids = pDataDT(combo_filt)[list_ID == "MS197U"]$cell_ID)
Giotto::spatDeconvPlot(MS197U,  
                       show_image = FALSE,
                       radius = 32,    
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       cell_color_code = cell_color,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                        save_folder = "Deconv/spatDeconv",
                                        save_name = "MS197U",
                                        base_width = 5,   
                                        base_height = 6,  
                                        units = "in",      
                                        dpi = 300))
Giotto::spatCellPlot(MS197U,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.5,
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "MS197U",                                  
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))

# MS377T
MS377T = Giotto::subsetGiotto(combo_filt, 
                              cell_ids = pDataDT(combo_filt)[list_ID == "MS377T"]$cell_ID)
Giotto::spatDeconvPlot(MS377T,  
                       show_image = FALSE,
                       radius = 31,          
                       cell_color_code = cell_color,
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "Deconv/spatDeconv",
                                         save_name = "MS377T",
                                         base_width = 5,   
                                         base_height = 6,  
                                         units = "in",      
                                         dpi = 300))
Giotto::spatCellPlot(MS377T,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.5,
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "MS377T",
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))


# MS549H
MS549H = Giotto::subsetGiotto(combo_filt, 
                              cell_ids = pDataDT(combo_filt)[list_ID == "MS549H"]$cell_ID)
Giotto::spatDeconvPlot(MS549H,
                       show_image = FALSE,
                       radius = 65,          
                       cell_color_code = cell_color,
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "Deconv/spatDeconv",
                                         save_name = "MS549H",
                                         base_width = 5,   
                                         base_height = 6,  
                                         units = "in",      
                                         dpi = 300))
Giotto::spatCellPlot(MS549H,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.6,
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "MS549H",
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))

# MS497I
MS497I = Giotto::subsetGiotto(combo_filt, 
                              cell_ids = pDataDT(combo_filt)[list_ID == "MS497I"]$cell_ID)
Giotto::spatDeconvPlot(MS497I,  
                       show_image = FALSE,
                       radius = 65,          
                       cell_color_code = cell_color,
                       axis_text = 12,
                       axis_title = 12,
                       legend_text = 16,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "Deconv/spatDeconv",
                                         save_name = "MS497I",
                                         base_width = 5,   
                                         base_height = 6,  
                                         units = "in",      
                                         dpi = 300))
Giotto::spatCellPlot(MS497I,
                     spat_enr_names = "DWLS",
                     cell_annotation_values = c("OL", "AS", "MG"), 
                     gradient_style = "divergent",
                     gradient_midpoint = 0.5,
                     cell_color_gradient = c("#440154FF", "#21908CFF", "#FDE725FF"),
                     gradient_limits = c(0.0, 1.0),
                     cow_n_col = 1,
                     coord_fix_ratio = 1,
                     show_legend = TRUE,
                     legend_text = 12,
                     axis_text = 10,
                     axis_title = 12,
                     point_size = 1.5,
                     save_plot = TRUE,
                     save_param = list(save_dir = "../Graficos",
                                       save_folder = "Deconv/",
                                       save_name = "MS497I",
                                       base_width = 5,   
                                       base_height = 12,  
                                       units = "in",      
                                       dpi = 300))
# CTRL
CTRL = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(combo_filt)[Lesion.type == "CTRL"]$cell_ID)
list_ids <- unique(CTRL@cell_metadata$cell$rna$list_ID)
metadata = pDataDT(CTRL)

plots_CTRL <- lapply(list_ids, function(SAMPLE) {
  cell_ids <- metadata[list_ID == SAMPLE]$cell_ID
  gobject <- subsetGiotto(combo_filt, cell_ids = cell_ids)
  Giotto::spatDeconvPlot(gobject,  
                         show_image = FALSE,
                         radius = 32,
                         title = SAMPLE)})

CTRL <- wrap_plots(plots_CTRL, ncol = 6)
ggsave("../Graficos/Deconv/CTRL.png",
       plot = CTRL,
       width = 18,
       height = 4,
       dpi = 300)

# MS-CA
CA = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(combo_filt)[Lesion.type == "EM-CA"]$cell_ID)
list_ids <- unique(CA@cell_metadata$cell$rna$list_ID)
metadata = pDataDT(CA)

plots_CA <- lapply(list_ids, function(SAMPLE) {
  cell_ids <- metadata[list_ID == SAMPLE]$cell_ID
  gobject <- subsetGiotto(combo_filt, cell_ids = cell_ids)
  Giotto::spatDeconvPlot(gobject,  
                         show_image = FALSE,
                         radius = 35,
                         title = SAMPLE)})

CA <- wrap_plots(plots_CA, ncol = 6)
ggsave("../Graficos/Deconv/CA.png",
       plot = CA,
       width = 18,
       height = 8,
       dpi = 300)

# MS-CI
CI = Giotto::subsetGiotto(combo_filt, 
                          cell_ids = pDataDT(combo_filt)[Lesion.type == "EM-CI"]$cell_ID)
list_ids <- unique(CI@cell_metadata$cell$rna$list_ID)
metadata = pDataDT(CI)

plots_CI <- lapply(list_ids, function(SAMPLE) {
  cell_ids <- metadata[list_ID == SAMPLE]$cell_ID
  gobject <- subsetGiotto(combo_filt, cell_ids = cell_ids)
  Giotto::spatDeconvPlot(gobject,  
                         show_image = FALSE,
                         radius = 65,
                         title = SAMPLE)})

CI <- wrap_plots(plots_CI, ncol = 4)
ggsave("../Graficos/Deconv/CI.png",
       plot = CI,
       width = 12,
       height = 4,
       dpi = 300)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#