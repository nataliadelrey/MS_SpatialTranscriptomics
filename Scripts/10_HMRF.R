# ST Spatial Domains
# Author: Natalia del Rey Díez
# Date: 11/05/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)


################################################################################
#--------------------------- Loading Giotto object ----------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/08_PatronesEspaciales/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#------------------------- Spatially Variable Genes ---------------------------#
################################################################################
combo_filt = Giotto::binSpect(combo_filt, 
                              bin_method = "kmeans",    
                              expression_values = "default",  
                              spatial_network_name = "Delaunay_network", 
                              kmeans_algo = "kmeans",   
                              do_fisher_test = TRUE,
                              adjust_method = "fdr", 
                              calc_hub = TRUE,                        
                              get_av_expr = TRUE,   
                              get_high_expr = TRUE,  
                              implementation = "data.table", 
                              group_size = "automatic",
                              return_gobject = TRUE)


################################################################################
#----------------------------- Output Directory -------------------------------#
################################################################################
dir_path <- "../RData/ST_GSE279181/09_HMRF"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)}


################################################################################
#------------------------------------ HMRF ------------------------------------#
################################################################################
# Initialization
hmrf_Top500 = Giotto::initHMRF_V2(gobject = combo_filt, 
                                  expression_values = "default",
                                  spatial_network_name = "Delaunay_network",
                                  use_spatial_genes = "binSpect",
                                  gene_list_from_top = 500,
                                  filter_method = "none",
                                  use_pca = FALSE,
                                  gene_samples = 500,
                                  gene_sampling_rate = 1,
                                  gene_sampling_seed = 1234,
                                  use_metagene = FALSE,
                                  cl.method = "km",
                                  hmrf_seed = 100,
                                  k = 9,
                                  zscore = "none",
                                  nstart = 1000)
save(hmrf_Top500, file = "../RData/ST_GSE279181/09_HMRF/hmrf_Top500b.rda")

# HMRF model
resHMRF_Top500 = Giotto::doHMRF_V2(hmrf_Top500, 
                                   betas=c(0, 5, 10))

save(resHMRF_Top500, file = "../RData/ST_GSE279181/09_HMRF/resHMRF_Top500.rda")

# Add HMRF Domain Type to metadata
combo_filt <- Giotto::addHMRF_V2(gobject = combo_filt, 
                                 HMRFoutput = resHMRF_Top500, 
                                 name = "hmrf_Top500")


################################################################################
#------------------------------ Beta Comparison -------------------------------#
################################################################################
cell_color <- Giotto::getDistinctColors(9)
names(cell_color) <- as.character(1:length(cell_color))
cell_color["7"] <- "#B3DE69"

# beta = 5
beta05 = Giotto::spatPlot2D(gobject = combo_filt,
                            cell_color = "hmrf_Top500 k=9 b=5.00",
                            cell_color_code = cell_color,
                            group_by = "list_ID", 
                            group_by_subset = c("CO37", "MS411", "MS549T"),
                            point_size = 1.2,        
                            cow_n_col = 1,
                            axis_text = 7,
                            title = "HMRF beta = 5.0",
                            save_plot = TRUE,
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "PatronesEspaciales/HMRF",
                                              save_name = "beta5",
                                              base_width = 3,   
                                              base_height = 9)) 
# beta = 15
beta15 = Giotto::spatPlot2D(gobject = combo_filt,
                            cell_color = "hmrf_Top500 k=9 b=15.00",
                            cell_color_code = cell_color,
                            group_by = "list_ID",              
                            group_by_subset = c("CO37", "MS411", "MS549T"),
                            point_size = 1.2,       
                            cow_n_col = 1,
                            axis_text = 7,
                            title = "HMRF beta = 15.0",
                            save_plot = TRUE,
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "PatronesEspaciales/HMRF",
                                              save_name = "beta15",
                                              base_width = 3,   
                                              base_height = 9)) 

beta15_all = Giotto::spatPlot2D(gobject = combo_filt,
                                cell_color = "hmrf_Top500 k=9 b=15.00",
                                cell_color_code = cell_color,
                                group_by = "list_ID",              
                                point_size = 1.2,       
                                cow_n_col = 6,
                                axis_text = 7,
                                title = "HMRF beta = 15.0",
                                save_plot = TRUE,
                                save_param = list(save_dir = "../Graficos",
                                                  save_folder = "PatronesEspaciales/HMRF",
                                                  save_name = "beta15_all",
                                                  base_width = 18,   
                                                  base_height = 10)) 

# beta = 25
beta25 = Giotto::spatPlot2D(gobject = combo_filt,
                            cell_color = "hmrf_Top500 k=9 b=25.00",
                            cell_color_code = cell_color,
                            group_by = "list_ID", 
                            group_by_subset = c("CO37", "MS411", "MS549T"),
                            point_size = 1.2,       
                            cow_n_col = 1,
                            axis_text = 7,
                            title = "HMRF beta = 25.0",
                            save_plot = TRUE,
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "PatronesEspaciales/HMRF",
                                              save_name = "beta25",
                                              base_width = 3,   
                                              base_height = 9))  
# beta = 45
beta45 = Giotto::spatPlot2D(gobject = combo_filt,
                            cell_color = "hmrf_Top500 k=9 b=45.00",
                            cell_color_code = cell_color,
                            group_by = "list_ID", 
                            group_by_subset = c("CO37", "MS411", "MS549T"),
                            point_size = 1.2,       
                            cow_n_col = 1,
                            axis_text = 7,
                            title = "HMRF beta = 45.0",
                            save_plot = TRUE,
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "PatronesEspaciales/HMRF",
                                              save_name = "beta45",
                                              base_width = 3,   
                                              base_height = 9)) 
 

################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
WM_CTRL = Giotto::spatPlot2D(gobject = combo_filt,
                             cell_color = "hmrf_Top500 k=9 b=15.00",
                             coord_fix_ratio = 1, 
                             cell_color_code = c("1" = "#E7298A",
                                                 "2" = "#FF7F00", 
                                                 "3" = "grey",
                                                 "4" = "grey",
                                                 "5" = "grey",
                                                 "6" = "grey",
                                                 "7" = "grey",
                                                 "8" = "grey",
                                                 "9" = "grey"), 
                             group_by = "list_ID",        
                             group_by_subset = c("CO37", "CO40", "MS377I", 
                                                 "MS377N", "MS497I", "MS549H"),
                             point_size = 1.8,       
                             cow_n_col = 6,
                             axis_text = 7,
                             title = "",
                             show_legend = FALSE,
                             save_plot = TRUE,
                             save_param = list(save_dir = "../Graficos",
                                               save_folder = "PatronesEspaciales/HMRF",
                                               save_name = "WM_CTRL",
                                               base_width = 25,   
                                               base_height = 5)) 

WM_EM = Giotto::spatPlot2D(gobject = combo_filt,
                           cell_color = "hmrf_Top500 k=9 b=15.00",
                           group_by = "list_ID",    
                           cell_color_code = c("1" = "grey",
                                               "2" = "grey", 
                                               "3" = "grey",
                                               "4" = "#377EB8",
                                               "5" = "grey",
                                               "6" = "grey",
                                               "7" = "grey",
                                               "8" = "grey",
                                               "9" = "grey"), 
                           group_by_subset = c("CO37", "CO40", "MS377I", 
                                               "MS377N", "MS497I", "MS549T"),
                           point_size = 1.8,       
                           cow_n_col = 6,
                           axis_text = 7,
                           title = "",
                           show_legend = FALSE,
                           save_plot = TRUE,
                           save_param = list(save_dir = "../Graficos",
                                             save_folder = "PatronesEspaciales/HMRF",
                                             save_name = "WM_EM",
                                             base_width = 25,   
                                             base_height = 5)) 

LC = Giotto::spatPlot2D(gobject = combo_filt,
                        cell_color = "hmrf_Top500 k=9 b=15.00",
                        group_by = "list_ID",    
                        cell_color_code = c("1" = "grey",
                                            "2" = "grey", 
                                            "3" = "#E41A1C",
                                            "4" = "grey",
                                            "5" = "grey",
                                            "6" = "#A65628",
                                            "7" = "grey",
                                            "8" = "#BC80BD",
                                            "9" = "grey"), 
                        group_by_subset = c("CO37", "CO40", "MS377I", 
                                            "MS377N", "MS497I", "MS549T"),
                        point_size = 1.8,       
                        cow_n_col = 6,
                        axis_text = 7,
                        title = "",
                        show_legend = FALSE,
                        save_plot = TRUE,
                        save_param = list(save_dir = "../Graficos",
                                          save_folder = "PatronesEspaciales/HMRF",
                                          save_name = "Nucleo",
                                          base_width = 25,   
                                          base_height = 5))  

LR = Giotto::spatPlot2D(gobject = combo_filt,
                        cell_color = "hmrf_Top500 k=9 b=15.00",
                        group_by = "list_ID",    
                        cell_color_code = c("1" = "grey",
                                            "2" = "grey", 
                                            "3" = "grey",
                                            "4" = "grey",
                                            "5" = "grey",
                                            "6" = "grey",
                                            "7" = "#B3DE69",
                                            "8" = "grey",
                                            "9" = "grey"), 
                        group_by_subset = c("CO37", "CO40", "MS377I", 
                                            "MS377N", "MS497I", "MS549T"),
                        point_size = 1.8,       
                        cow_n_col = 6,
                        axis_text = 7,
                        title = "",
                        show_legend = FALSE,
                        save_plot = TRUE,
                        save_param = list(save_dir = "../Graficos",
                                          save_folder = "PatronesEspaciales/HMRF",
                                          save_name = "Borde",
                                          base_width = 25,   
                                          base_height = 5))  

GM = Giotto::spatPlot2D(gobject = combo_filt,
                        cell_color = "hmrf_Top500 k=9 b=15.00",
                        group_by = "list_ID",    
                        cell_color_code = c("1" = "grey",
                                            "2" = "grey", 
                                            "3" = "grey",
                                            "4" = "grey",
                                            "5" = "#4DAF4A",
                                            "6" = "grey",
                                            "7" = "grey",
                                            "8" = "grey",
                                            "9" = "#FFED6F"), 
                        group_by_subset = c("CO37", "CO40", "MS377I", 
                                            "MS377N", "MS497I", "MS549T"),
                        point_size = 1.8,       
                        cow_n_col = 6,
                        axis_text = 7,
                        title = "",
                        show_legend = FALSE,
                        save_plot = TRUE,
                        save_param = list(save_dir = "../Graficos",
                                          save_folder = "PatronesEspaciales/HMRF",
                                          save_name = "GM",
                                          base_width = 25,   
                                          base_height = 5))  

################################################################################
#----------------------------- Domain annotation ------------------------------#
################################################################################
cell_metadata <- pDataDT(combo_filt)

cell_metadata$anotacion <- with(
  cell_metadata, 
  ifelse(`hmrf_Top500 k=9 b=15.00` %in% c(1, 2) & Condition == "Control", "SB Sana",
         ifelse(`hmrf_Top500 k=9 b=15.00` == 4 & Condition == "Control", "SB Sana",
                ifelse(`hmrf_Top500 k=9 b=15.00` == 4 & Condition == "MS", "SB Perilesional",
                       ifelse(`hmrf_Top500 k=9 b=15.00` %in% c(5, 9), "SG",  
                              ifelse(`hmrf_Top500 k=9 b=15.00` %in% c(3, 6, 8) & Condition == "MS", "NL",
                                            ifelse(`hmrf_Top500 k=9 b=15.00` == 7, "BL", NA)))))))

mobj = Giotto::createCellMetaObj(cell_metadata)
combo_filt <- Giotto::setCellMetadata(gobject = combo_filt, 
                                      x = mobj)

Giotto::spatPlot2D(gobject = combo_filt,
                   cell_color = "anotacion", 
                   group_by = "list_ID",
                   point_size = 1.3,  
                   coord_fix_ratio = 1,
                   cow_n_col = 6,
                   axis_text = 7,
                   show_legend = FALSE,
                   save_plot = TRUE,
                   save_param = list(save_dir = "../Graficos",
                                     save_folder = "PatronesEspaciales/HMRF",
                                     save_name = "anotacion",
                                     base_width = 18,   
                                     base_height = 14))


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/09_HMRF",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#                        