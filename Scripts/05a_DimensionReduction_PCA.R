# ST Dimension Reduction - PCA
# Author: Natalia del Rey Díez
# Date: 08/04/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/04_HVG/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#--------------------- Log-normalization by library size ----------------------#
################################################################################
#### Covariance groups ####
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             name = "pca_default_groups",
                             feats_to_use = "hvf_groups",
                             center = TRUE,
                             scale_unit = TRUE,
                             ncp = 50) 

#### Loess #### 
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             name = "pca_default_loess",
                             feats_to_use = "hvf_loess",
                             center = TRUE,
                             scale_unit = TRUE,
                             ncp = 50) 
#### Pearson #### 
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "default",
                             reduction = "cells",
                             name = "pca_default_pearson",
                             feats_to_use = "hvf_pearson",
                             center = TRUE,
                             scale_unit = TRUE,
                             ncp = 50) 


################################################################################
#---------------------- Pearson Residuals Normalization -----------------------#
################################################################################
#### Covariance groups ####
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             name = "pca_pearson_groups",
                             feats_to_use = "hvf_groups",
                             center = FALSE,
                             scale_unit = FALSE,
                             ncp = 50) 

#### Loess #### 
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             name = "pca_pearson_loess",
                             feats_to_use = "hvf_loess",
                             center = FALSE,
                             scale_unit = FALSE,
                             ncp = 50) 
#### Pearson #### 
combo_filt <- Giotto::runPCA(combo_filt, 
                             expression_values = "pearson",
                             reduction = "cells",
                             name = "pca_pearson_pearson",
                             feats_to_use = "hvf_pearson",
                             center = FALSE,
                             scale_unit = FALSE,
                             ncp = 50) 


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/05a_PCA",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#---------------------------------- plotPCA -----------------------------------#
##-------------------- Log-normalization by library size ---------------------##
################################################################################
fill_colors <- Giotto::getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 1234)

#### Covariance groups ####
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_default_groups",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería + Grupos de covarianza",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Default_Groups"))

#### Loess #### 
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_default_loess",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería + Regresion Loess",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Default_Loess"))

#### Pearson #### 
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_default_pearson",
                cell_color = "Lesion.type",
                title = "Log-Normalización por tamaño de librería +  Residuos de Pearson",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Default_Pearson"))


################################################################################
#---------------------------------- plotPCA -----------------------------------#
##--------------------- Pearson Residuals Normalization ----------------------##
################################################################################
#### Covariance groups ####
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_pearson_groups",
                cell_color = "Lesion.type",
                title = "Normalización Pearson + Grupos de covarianza",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Pearson_Groups"))

#### Loess #### 
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_pearson_loess",
                cell_color = "Lesion.type",
                title = "Normalizacion Pearson + Regresion Loess",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Pearson_Loess"))

#### Pearson #### 
Giotto::plotPCA(combo_filt, 
                dim_reduction_name = "pca_pearson_pearson",
                cell_color = "Lesion.type",
                title = "Normalizacion Pearson + Residuos de Pearson",
                cell_color_code = fill_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "RedDim/PCA/combo_filt",
                                  save_name = "Pearson_Pearson"))


################################################################################
#--------------------------------- screePlot ----------------------------------#
################################################################################
#### Covariance groups ####
Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_default_groups",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_groups",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Default_Groups"))

Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_pearson_groups",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_groups",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Pearson_Groups"))

#### Loess #### 
Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_default_loess",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_loess",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Default_Loess"))

Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_pearson_loess",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_loess",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Pearson_Loess"))

#### Pearson #### 
Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_default_pearson",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_pearson",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Default_Pearson"))

Giotto::screePlot(combo_filt,
                  expression_values = "default",
                  dim_reduction_name = "pca_pearson_pearson",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf_pearson",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "RedDim/PCA/screePlot",
                                    save_name = "Pearson_Pearson"))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#