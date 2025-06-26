# ST Spatial Domains
# Author: Natalia del Rey Díez
# Date: 25/05/2024

################################################################################
#------------------------------ Install packages ------------------------------#
################################################################################
# EnhancedVolcano
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")}

# ggvenn
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  install.packages("ggvenn")}


################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)
library(ggplot2)
library(EnhancedVolcano)
library(ggvenn)


################################################################################
#--------------------------- Loading Giotto object ----------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/09_HMRF/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#--------------------- Sex differences in the lesion core ---------------------#
################################################################################
LC_MS = Giotto::subsetGiotto(combo_filt, 
                             cell_ids = pDataDT(combo_filt)[anotacion == "NL"]$cell_ID)


LC_F_VS_M = findMarkers(LC_MS,
                        expression_values = "default",
                        cluster_column = "Sex",
                        method = "scran")
LC_F_VS_M <- LC_F_VS_M[[1]]
DEG_LC = LC_F_VS_M[which(FDR < 0.05 & abs(logFC.Mujer) > 0.5)]
cat("LC CTRL - DEGs for FDR < 0.05 and logFC absolute > 0.5:", nrow(DEG_LC))

custom_colors <- ifelse(LC_F_VS_M$FDR < 0.05 & LC_F_VS_M$summary.logFC < -0.5, 
                        "#FFC98B",
                        ifelse(LC_F_VS_M$FDR < 0.05 & LC_F_VS_M$summary.logFC > 0.5, 
                               "#CCA7FF", "grey"))
names(custom_colors) <- LC_F_VS_M$feats
lc <- EnhancedVolcano::EnhancedVolcano(LC_F_VS_M,
                                       legendPosition = "none",
                                       lab = LC_F_VS_M$feats,                 
                                       x = "summary.logFC",                
                                       y = "FDR",                      
                                       title = "Núcleo lesión: Mujer VS Hombre",
                                       subtitle = NULL,
                                       pCutoff = 0.05,   
                                       colAlpha = 1,
                                       FCcutoff = 0.5,                       
                                       pointSize = 2.0,
                                       labSize = 3,
                                       selectLab = DEG_LC$feats,
                                       legendLabSize = 12,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       colCustom = custom_colors,
                                       widthConnectors = 0.5,
                                       boxedLabels = FALSE)

ggsave("../Graficos/DEA/NucleoLesion_MvsH.png", 
       plot = lc, 
       width = 14, 
       height = 9, 
       dpi = 300)


################################################################################
#------------------ Sex differences in healthy white matter -------------------#
################################################################################
WM_CTRL = Giotto::subsetGiotto(combo_filt, 
                               cell_ids = pDataDT(combo_filt)[anotacion == "SB Sana"]$cell_ID)


WM_F_VS_M = findMarkers(WM_CTRL,
                        expression_values = "default",
                        cluster_column = "Sex",
                        method = "scran")

DEG_WM = WM_F_VS_M[[1]][which(FDR < 0.05 & abs(logFC.Mujer) > 0.5)]
cat("WM CTRL - DEGs for FDR < 0.05 and logFC absolute > 0.5:", nrow(DEG_WM))


################################################################################
#--- Differences between perilesional white matter and lesion core in women ---#
################################################################################
W_MS = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(
                              combo_filt)[Sex == "Mujer" & Condition == "MS" &
                                            (anotacion == "SB Perilesional" | 
                                               anotacion == "NL")]$cell_ID)


W_WMP_VS_LC = findMarkers(W_MS,
                          expression_values = "default",
                          cluster_column = "anotacion",
                          method = "scran")
W_WMP_VS_LC <- W_WMP_VS_LC[[1]]
DEG_W = W_WMP_VS_LC[which(FDR < 0.05 & abs(logFC.SB.Perilesional) > 0.5)]
cat("WMP VS LC Female - DEGs for FDR < 0.05 and logFC absolute > 0.5:", nrow(DEG_W))


################################################################################
#---- Differences between perilesional white matter and lesion core in men ----#
################################################################################
M_MS = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(
                              combo_filt)[Sex == "Hombre" & Condition == "MS" &
                                            (anotacion == "SB Perilesional" | 
                                               anotacion == "NL")]$cell_ID)


M_WMP_VS_LC = findMarkers(M_MS,
                          expression_values = "default",
                          cluster_column = "anotacion",
                          method = "scran")

DEG_M = M_WMP_VS_LC[[1]][which(FDR < 0.05 & abs(logFC.SB.Perilesional) > 0.5)]
cat("WMP VS LC Male - DEGs for FDR < 0.05 and logFC absolute > 0.5:", nrow(DEG_M))


################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
# Differences between men and women in both lesion core and white matter
LC_WM <- list("Núcleo Lesión EM" = DEG_LC$feats,
              "Sustancia Blanca CTRL" = DEG_WM$feats)
LC_WM_Venn = ggvenn(LC_WM, 
                    fill_color = c("#B84136", "#EDE0E0"), 
                    stroke_size = 0.5, 
                    set_name_size = 4)
ggsave("../Graficos/DEA/NL_SB_Venn.png", 
       plot = LC_WM_Venn, 
       width = 5, 
       height = 3, 
       dpi = 300)

# Differences between lesion core and white matter in both men and women
W_M <- list("♀ - Up" = DEG_W$feats[which(DEG_W$summary.logFC > 0)],
            "♀ - Down" = DEG_W$feats[which(DEG_W$summary.logFC < 0)],
            "♂ - Down" = DEG_M$feats[which(DEG_M$summary.logFC < 0)],
            "♂ - Up" = DEG_M$feats[which(DEG_M$summary.logFC > 0)])
W_M_Venn = ggvenn(W_M, 
                  fill_color = c("#FF8E4B", "#FEE3C5", "#E6D5FF",  "#B57BFF"), 
                  stroke_size = 0.5, 
                  set_name_size = 4)

ggsave("../Graficos/DEA/Muj_Hom_Venn.png", 
       plot = W_M_Venn, 
       width = 10, 
       height = 6, 
       dpi = 300)

# Heatmap
MS = Giotto::subsetGiotto(combo_filt, 
                            cell_ids = pDataDT(
                              combo_filt)[Condition == "MS" & 
                                            (anotacion == "SB Perilesional" | 
                                               anotacion == "NL")]$cell_ID)
feats <- unique(c(setdiff(DEG_W$feats, DEG_M$feats), 
                  setdiff(DEG_M$feats, DEG_W$feats)))
                      
Giotto::plotMetaDataHeatmap(gobject = MS,
                            expression_values = "default",
                            selected_feats = feats,
                            metadata_cols = c("Sex", "anotacion"),
                            first_meta_col = "anotacion",
                            second_meta_col = "Sex",
                            show_values = "zscores",
                            custom_feat_order = sort(feats, decreasing = TRUE),
                            x_text_size = 12, 
                            show_plot = TRUE,
                            y_text_size = 11,
                            save_plot = TRUE, 
                            save_param = list(save_dir = "../Graficos",
                                              save_folder = "DEA/",
                                              save_name = "Heatmap", 
                                              base_width = 14,   
                                              base_height = 20))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#