# ST Feature Selection: Highly Variable Genes (HVG)
# Author: Natalia del Rey Díez
# Date: 02/04/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/03_Normalizacion/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#----------------------------- Feature Selection ------------------------------#
################################################################################
#### Covariance groups ####
combo_filt = calculateHVF(combo_filt,
                          expression_values = "default",
                          method = "cov_groups",
                          reverse_log_scale = FALSE,
                          expression_threshold = 0,
                          nr_expression_groups = 20,
                          zscore_threshold = 1.5,
                          HVFname = "hvf_groups",
                          verbose = TRUE,
                          save_plot = TRUE,
                          save_param = list(save_dir = "../Graficos",
                                            save_folder = "HVG_Default",
                                            save_name = "Groups"))

#### Loess #### 
combo_filt = calculateHVF(combo_filt,
                          expression_values = "default",
                          method = "cov_loess",
                          reverse_log_scale = FALSE,
                          expression_threshold = 1,
                          HVFname = "hvf_loess",
                          difference_in_cov = 0.1,
                          verbose = TRUE,
                          save_plot = TRUE,
                          save_param = list(save_dir = "../Graficos",
                                            save_folder = "HVG_Default",
                                            save_name = "Loess"))

#### Pearson #### 
combo_filt = calculateHVF(combo_filt,
                          expression_values = "default",
                          method = "var_p_resid",
                          reverse_log_scale = FALSE,
                          expression_threshold = 1,
                          HVFname = "hvf_pearson",
                          var_threshold = 1.5,
                          verbose = TRUE,
                          save_plot = TRUE,
                          save_param = list(save_dir = "../Graficos",
                                            save_folder = "HVG_Default",
                                            save_name = "Pearson"))


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/04_HVG",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#-------------------------------- Comparation ---------------------------------#
################################################################################
# Filtered
loess_f <- fDataDT(combo_filt)$feat_ID[fDataDT(combo_filt)$hvf_loess == "yes"]
pearson_f <- fDataDT(combo_filt)$feat_ID[fDataDT(combo_filt)$hvf_pearson == "yes"]
groups_f <- fDataDT(combo_filt)$feat_ID[fDataDT(combo_filt)$hvf_groups == "yes"]

genes_filt <- list(Loess = loess_f,
                   Pearson = pearson_f,
                   Groups = groups_f)

# Intersection
inter_filt <- Reduce(intersect, genes_filt)
cat("HVG identified by the three methods:", inter_filt)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#                         