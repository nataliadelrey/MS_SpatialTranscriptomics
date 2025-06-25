# ST Spatial Patterns
# Author: Natalia del Rey Díez
# Date: 29/04/2024

################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)


################################################################################
#--------------------------- Loading Giotto object ----------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/07_Deconvolucion/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#------------------------------ Delaunay Network ------------------------------#
################################################################################
combo_filt = Giotto::createSpatialDelaunayNetwork(combo_filt,
                                                  name = "Delaunay_network",
                                                  method = "deldir",
                                                  maximum_distance = "auto",
                                                  minimum_k = 3)

Giotto::plotStatDelaunayNetwork(combo_filt, 
                                method = "deldir",
                                maximum_distance = "auto",
                                minimum_k = 3,
                                save_plot = TRUE,
                                save_param = list(save_dir = "../Graficos",
                                                  save_folder = "/PatronesEspaciales/",
                                                  save_name = "Delaunay"))


################################################################################
#----------------------------- Output Directory -------------------------------#
################################################################################
dir_path <- "../RData/ST_GSE279181/08_PatronesEspaciales"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)}


################################################################################
#------------------------- Spatially Variable Genes ---------------------------#
################################################################################
# kmeans
kmeans = Giotto::binSpect(combo_filt, 
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
                          group_size = "automatic")
save(kmeans, file = "../RData/ST_GSE279181/08_PatronesEspaciales/kmeans.rda")

# rank
rank = Giotto::binSpect(combo_filt,
                        bin_method = "rank",
                        expression_values = "default",
                        spatial_network_name = "Delaunay_network",
                        percentage_rank = 30,
                        do_fisher_test = TRUE,
                        adjust_method = "fdr",
                        calc_hub = TRUE,
                        get_av_expr = TRUE,
                        get_high_expr = TRUE,
                        implementation = "data.table",
                        group_size = "automatic")
save(rank, file = "../RData/ST_GSE279181/08_PatronesEspaciales/rank.rda")


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/08_PatronesEspaciales",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#------------------------------- SVG FDR < 0.01 -------------------------------#
################################################################################
sig_kmeans <- kmeans[adj.p.value < 0.01, feats]
cat("SVG FDR < 0.01 kmeans:", length(sig_kmeans))

sig_rank <- rank[adj.p.value < 0.01, feats]
cat("SVG FDR < 0.01 rank:", length(sig_rank))

common <- intersect(sig_kmeans, sig_rank)
cat("SVG identified by the two methods::", length(common))


################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
combo_filt = processExpression(combo_filt,
                               expression_values = "default",
                               scaleParam(),
                               name = "scaled")

Giotto::spatFeatPlot2D(combo_filt, 
                       group_by  = "list_ID",
                       expression_values = "scaled", 
                       gradient_midpoint = 0,
                       gradient_style = "divergent",
                       group_by_subset = c("CO37", "CO85", "MS197U", "MS377T", "MS549T","MS549H"),
                       feats = c("MBP", "GFAP", "IGHG4"),
                       point_size = 1.2,
                       cow_align = "hv",
                       cow_n_col = 6,
                       save_plot = TRUE,
                       save_param = list(save_dir = "../Graficos",
                                         save_folder = "PatronesEspaciales/",
                                         save_name = "SVG",
                                         base_aspect_ratio = 1,
                                         base_width = 22,   
                                         base_height = 10))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#