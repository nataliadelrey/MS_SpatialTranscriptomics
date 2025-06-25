# ST Dimension Reduction - UMAP
# Author: Natalia del Rey Díez
# Date: 10/04/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/05a_PCA/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#--------------------- Log-normalization by library size ----------------------#
################################################################################
#### Covariance groups ####
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             dim_reduction_name = "pca_default_groups", 
                             name = "umap_default_groups",
                             feats_to_use = "hvf_groups",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 

#### Loess #### 
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             dim_reduction_name = "pca_default_loess",
                             name = "umap_default_loess",
                             feats_to_use = "hvf_loess",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 
#### Pearson #### 
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             dim_reduction_name = "pca_default_pearson",
                             name = "umap_default_pearson",
                             feats_to_use = "hvf_pearson",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 


################################################################################
#---------------------- Pearson Residuals Normalization -----------------------#
################################################################################
#### Covariance groups ####
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             dim_reduction_name = "pca_pearson_groups",
                             name = "umap_pearson_groups",
                             feats_to_use = "hvf_groups",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 

#### Loess #### 
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             dim_reduction_name = "pca_pearson_loess",
                             name = "umap_pearson_loess",
                             feats_to_use = "hvf_loess",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 
#### Pearson #### 
combo_filt <- Giotto::runUMAP(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             dim_reduction_name = "pca_pearson_pearson",
                             name = "umap_pearson_pearson",
                             feats_to_use = "hvf_pearson",
                             dim_reduction_to_use = "pca",
                             dimensions_to_use = 1:20) 


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/05b_UMAP",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#---------------------------------- plotUMAP ----------------------------------#
##-------------------- Log-normalization by library size ---------------------##
################################################################################
fill_colors <- Giotto::getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 1234)

#### Covariance groups ####
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_default_groups",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería + Grupos de covarianza",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Default_Groups"))

#### Loess #### 
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_default_loess",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería + Regresion Loess",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Default_Loess"))

#### Pearson #### 
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_default_pearson",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería +  Residuos de Pearson",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Default_Pearson"))


################################################################################
#---------------------------------- plotUMAP ----------------------------------#
##--------------------- Pearson Residuals Normalization ----------------------##
################################################################################
#### Covariance groups ####
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_pearson_groups",
                cell_color = "Lesion.type",
                title = "Normalización Pearson + Grupos de covarianza",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Pearson_Groups"))

#### Loess #### 
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_pearson_loess",
                cell_color = "Lesion.type",
                title = "Normalizacion Pearson + Regresion Loess",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Pearson_Loess"))

#### Pearson #### 
Giotto::plotUMAP(combo_filt, 
                dim_reduction_name = "umap_pearson_pearson",
                cell_color = "Lesion.type",
                title = "Normalizacion Pearson + Residuos de Pearson",
                cell_color_code = fill_colors, axis_text = 12, axis_title = 12,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/UMAP/combo_filt",
                                  save_name = "Pearson_Pearson"))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#