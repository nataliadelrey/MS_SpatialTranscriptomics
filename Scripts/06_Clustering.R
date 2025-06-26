# ST Clustering
# Author: Natalia del Rey Díez
# Date: 14/04/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)
library(patchwork)
library(ggplot2)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/05b_UMAP/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#---------------------------- Harmony Integration -----------------------------#
################################################################################
combo_filt <- Giotto::runGiottoHarmony(combo_filt, 
                                       vars_use = c("list_ID", 
                                                    "Batch.Visium..ST."),
                                       dim_reduction_to_use = "pca",
                                       dim_reduction_name = "pca_default_loess",
                                       dimensions_to_use = 1:10,
                                       name = "harmony",
                                       plot_convergence = TRUE)


################################################################################
#----------------------------- UMAP with Harmony ------------------------------#
################################################################################
combo_filt <- Giotto::runUMAP(combo_filt, 
                              expression_values = "default",
                              reduction = "cells",
                              name = "umap_harmony",
                              dim_reduction_to_use = "harmony",
                              dim_reduction_name = "harmony") 


################################################################################
#------------------------------ Nearest Network -------------------------------#
################################################################################
#kNN
combo_filt <- Giotto::createNearestNetwork(gobject = combo_filt,
                                           type = "kNN",
                                           dim_reduction_to_use = "harmony", 
                                           dim_reduction_name = "harmony",
                                           name = "kNN.harmony",
                                           expression_values = "default",
                                           dimensions_to_use = 1:15, 
                                           k = 15) 

#sNN
combo_filt <- Giotto::createNearestNetwork(gobject = combo_filt,
                                           type = "sNN",
                                           dim_reduction_to_use = "harmony", 
                                           dim_reduction_name = "harmony",
                                           name = "sNN.harmony",
                                           expression_values = "default",
                                           dimensions_to_use = 1:15, 
                                           k = 15) 


################################################################################
#----------------------------- Leiden clustering ------------------------------#
################################################################################
#### kNN Network ####
# Res 0.15
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_kNN.015",
                                      nn_network_to_use = "kNN",
                                      network_name = "kNN.harmony",
                                      resolution = 0.15, 
                                      n_iterations = -1000)
# Res 0.25
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_kNN.025",
                                      nn_network_to_use = "kNN",
                                      network_name = "kNN.harmony",
                                      resolution = 0.25, 
                                      n_iterations = -1000)
# Res 0.50
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_kNN.050",
                                      nn_network_to_use = "kNN",
                                      network_name = "kNN.harmony",
                                      resolution = 0.50, 
                                      n_iterations = -1000)


#### sNN Network ####
# Res 0.15
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_sNN.015",
                                      nn_network_to_use = "sNN",
                                      network_name = "sNN.harmony",
                                      resolution = 0.15, 
                                      n_iterations = -1000)
# Res 0.25
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_sNN.025",
                                      nn_network_to_use = "sNN",
                                      network_name = "sNN.harmony",
                                      resolution = 0.25, 
                                      n_iterations = -1000)
# Res 0.50
combo_filt <- Giotto::doLeidenCluster(gobject = combo_filt, 
                                      name = "leiden_sNN.050",
                                      nn_network_to_use = "sNN",
                                      network_name = "sNN.harmony",
                                      resolution = 0.50, 
                                      n_iterations = -1000)


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/06_Leiden",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
# kNN VS sNN for an example sample (MS377T)
MS377T = Giotto::subsetGiotto(combo_filt, 
                              cell_ids = pDataDT(combo_filt)[list_ID == "MS377T"]$cell_ID)

cell_color <- Giotto::getDistinctColors(24)
names(cell_color) <- as.character(1:length(cell_color))


# leiden_kNN.015
kNN_015 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_kNN.015",
                             cell_color_code = cell_color,
                             point_size = 1.3,        
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "kNN - Resolución 0,15") 
# leiden_kNN.025
kNN_025 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_kNN.025",
                             cell_color_code = cell_color,
                             point_size = 1.3,          
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "kNN - Resolución 0,25") 
# leiden_kNN.050
kNN_050 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_kNN.050",
                             cell_color_code = cell_color,
                             point_size = 1.3,          
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "kNN - Resolución 0,50") 
# leiden_sNN.015
sNN_015 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_sNN.015",
                             cell_color_code = cell_color,
                             point_size = 1.3,            
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "sNN - Resolución 0,15") 
# leiden_sNN.025
sNN_025 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_sNN.025",
                             cell_color_code = cell_color,
                             point_size = 1.3,            
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "sNN - Resolución 0,25") 
# leiden_sNN.050
sNN_050 = Giotto::spatPlot2D(gobject = MS377T,
                             cell_color = "leiden_sNN.050",
                             cell_color_code = cell_color,
                             point_size = 1.3, 
                             point_alpha = 0.85,
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "sNN - Resolución 0,50") 

MS377T_plot <- (sNN_015 | sNN_025 | sNN_050) / 
  plot_spacer() /
  (kNN_015 | kNN_025 | kNN_050) + plot_layout(heights = c(1, 0.1, 1)) 

ggsave("../Graficos/Leiden/MS377T_Leiden.png",
       plot = MS377T_plot,
       width = 14,
       height = 8,
       dpi = 300)


################################################################################
#---------------------------------- plotUMAP ----------------------------------#
################################################################################
# kNN 0,15
Giotto::plotUMAP(combo_filt, 
                 dim_reduction_name = "umap_harmony",
                 cell_color = "leiden_kNN.015",
                 cell_color_code = cell_color,
                 title = "kNN - Resolución 0,15",
                 point_alpha = 0.85,
                 axis_text = 12, axis_title = 15, label_size = 6,
                 show_center_label = TRUE,
                 legend_text = 14,
                 legend_symbol_size = 6,
                 save_plot = TRUE,
                 save_param = list(save_dir = "../Graficos",
                                   save_folder = "Leiden/UMAP/combo_filt",
                                   save_name = "leiden_kNN.015"))
# kNN 0,25
Giotto::plotUMAP(combo_filt, 
                 dim_reduction_name = "umap_harmony",
                 cell_color = "leiden_kNN.025",
                 cell_color_code = cell_color,
                 title = "kNN - Resolución 0,25",
                 point_alpha = 0.85,
                 axis_text = 12, axis_title = 15, label_size = 6,
                 show_center_label = TRUE,
                 legend_text = 14,
                 legend_symbol_size = 6,
                 save_plot = TRUE, 
                 save_param = list(save_dir = "../Graficos",
                                   save_folder = "Leiden/UMAP/combo_filt",
                                   save_name = "leiden_kNN.025"))

# kNN 0,50
Giotto::plotUMAP(combo_filt, 
                 dim_reduction_name = "umap_harmony",
                 cell_color = "leiden_kNN.050",
                 cell_color_code = cell_color,
                 title = "kNN - Resolución 0,50",
                 point_alpha = 0.85,
                 show_center_label = TRUE,
                 axis_text = 12, axis_title = 15, label_size = 6, 
                 legend_text = 14,
                 legend_symbol_size = 6,
                 save_plot = TRUE, 
                 save_param = list(save_dir = "../Graficos",
                                   save_folder = "Leiden/UMAP/combo_filt",
                                   save_name = "leiden_kNN.050"))


################################################################################
#--------------------------------- Comparison ---------------------------------#
################################################################################
kNN050 = Giotto::spatPlot2D(gobject = combo_filt,
                           cell_color = "leiden_kNN.050",
                           cell_color_code = cell_color,
                           group_by = "list_ID", 
                           group_by_subset = c("CO37", "CO40", "MS411", 
                                               "MS197U", "MS497I", "MS549T"),
                           point_size = 1.5,
                           point_alpha = 0.9,
                           cow_n_col = 3,
                           axis_text = 7, 
                           save_plot = TRUE, 
                           save_param = list(save_dir = "../Graficos",
                                             save_folder = "Leiden/",
                                             save_name = "leiden_kNN.050")) 

kNN050_all = Giotto::spatPlot2D(gobject = combo_filt,
                                cell_color = "leiden_kNN.050",
                                cell_color_code = cell_color,
                                group_by = "list_ID", 
                                point_size = 1,
                                point_alpha = 0.9,
                                cow_n_col = 6,
                                axis_text = 7, 
                                save_plot = TRUE, 
                                save_param = list(save_dir = "../Graficos",
                                                  save_folder = "Leiden/",
                                                  save_name = "leiden_kNN.050_all",
                                                  base_width = 18,
                                                  base_height = 10)) 


################################################################################
#----------------------------- Output Directory -------------------------------#
################################################################################
dir_path <- "../RData/ST_GSE279181/06_Leiden"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)}


################################################################################
#-------------------------------- Marker Genes --------------------------------#
################################################################################
scran_markers <- Giotto::findMarkers_one_vs_all(gobject = combo_filt,
                                                method = "scran",
                                                expression_values = "default",
                                                cluster_column = "leiden_kNN.050",
                                                min_feats = 5)

save(scran_markers, file = "../RData/ST_GSE279181/06_Leiden/scran_markers.rda")


topgenes_scran <- scran_markers[, head(.SD, 3), by = "cluster"]$feats

Giotto::plotMetaDataHeatmap(gobject = combo_filt,
                            expression_values = "default",
                            selected_feats = unique(topgenes_scran),
                            metadata_cols = "leiden_kNN.050",
                            show_values = "zscores",
                            x_text_size = 10, y_text_size = 10,
                            custom_cluster_order = c("1", "2", "3", "4", "5", 
                                                     "6", "7", "8", "9", "10",
                                                     "11", "12"),
                            save_plot = TRUE, 
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "Leiden/",
                                              save_name = "scran_markers"))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#