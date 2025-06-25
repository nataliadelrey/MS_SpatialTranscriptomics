# ST Quality Control
# Author: Natalia del Rey Díez
# Date: 13/03/2024

################################################################################
#------------------------------ Install packages ------------------------------#
################################################################################
# dplyr
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")}

# ggbeeswarm
if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
  install.packages("ggbeeswarm")}


################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(patchwork)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
object_names <- list.files(path = "../RData/ST_GSE279181/01_GiottoObjects", 
                           full.names = FALSE)

giotto_objects <- list()

for (obj_name in object_names) {
  giotto_object <- Giotto::loadGiotto(
    path_to_folder = paste0("../RData/ST_GSE279181/01_GiottoObjects/", obj_name),
    reconnect_giottoImage = TRUE,
    init_gobject = TRUE,
    verbose = TRUE)
  giotto_objects[[obj_name]] <- giotto_object
}


################################################################################
#------------------- Number of spots and genes per sample ---------------------#
################################################################################
num_genes_spots <- data.frame(Muestra = names(giotto_objects),  
                              Genes = sapply(giotto_objects, 
                                             function(g) dim(g)[1]),
                              Spots = sapply(giotto_objects, 
                                             function(g) dim(g)[2]),
                              Lesion.type = sapply(giotto_objects, 
                                                   function(g) 
                                                     pDataDT(g)$Lesion.type[1]) 
)
num_genes_spots


################################################################################
#---------------------------- Spot quality metrics ----------------------------#
################################################################################
giotto_objects <- lapply(giotto_objects, function(g) {
  Giotto::addCellStatistics(g,
                            feat_type = "rna",
                            spat_unit = "cell",
                            expression_values = "raw")
})


################################################################################
#---------------------------- Gene quality metrics ----------------------------#
################################################################################
giotto_objects <- lapply(giotto_objects, function(g) {
  Giotto::addFeatStatistics(g,
                            feat_type = "rna",
                            spat_unit = "cell",
                            expression_values = "raw")
})


################################################################################
#--------------------- Percentage of mitochondrial genes ----------------------#
################################################################################
mito_genes <- lapply(giotto_objects, function(g) 
  rownames(g)[grep("^MT-", rownames(g))]
)

giotto_objects <- lapply(giotto_objects, function(g) {
  Giotto::addFeatsPerc(g,
                       expression_values = "raw",
                       feats = unique(unlist(mito_genes)),
                       vector_name = "mito_ratio")
})


################################################################################
#-------------------------------- Join objects --------------------------------#
################################################################################
combo <- Giotto::joinGiottoObjects(gobject_list = giotto_objects,
                                   gobject_names = names(giotto_objects),
                                   join_method = "shift", 
                                   x_padding = 1000) 


################################################################################
#------------------------ Display filter combinations -------------------------#
################################################################################
Giotto::filterCombinations(combo,
                           feat_type = "rna",
                           spat_unit = "cell",
                           expression_values = "raw",
                           expression_thresholds = 1,
                           feat_det_in_min_cells = c(3, 5, 10, 15),
                           min_det_feats_per_cell = c(100, 100, 150, 200),
                           save_plot = TRUE,
                           save_param = list(save_dir = "../Graficos",
                                             save_name = "filterCombinations"))


################################################################################
#--------------------------------- Filtering ----------------------------------#
################################################################################
# Tag
combo <- Giotto::filterGiotto(combo,
                             spat_unit = "cell",
                             feat_type = "rna",
                             expression_values = "raw",
                             expression_threshold = 1,
                             feat_det_in_min_cells = 10,
                             min_det_feats_per_cell = 100,
                             tag_cells = TRUE, 
                             tag_cell_name = "filter_100",
                             tag_feats = TRUE, 
                             tag_feats_name = "filter_10",
                             verbose = FALSE)

combo <- Giotto::filterGiotto(combo,
                             spat_unit = "cell",
                             feat_type = "rna",
                             expression_values = "raw",
                             expression_threshold = 1,
                             feat_det_in_min_cells = 10,
                             min_det_feats_per_cell = 200,
                             tag_cells = TRUE, 
                             tag_cell_name = "filter_200",
                             tag_feats = TRUE, 
                             tag_feats_name = "filter_10",
                             verbose = FALSE)
# Filter
combo_filt <- Giotto::filterGiotto(combo,
                                  spat_unit = "cell",
                                  feat_type = "rna",
                                  expression_values = "raw",
                                  expression_threshold = 1,
                                  feat_det_in_min_cells = 10,
                                  min_det_feats_per_cell = 100, 
                                  verbose = FALSE)


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
# Tagged
Giotto::saveGiotto(gobject = combo,
                   foldername = "combo",
                   dir = "../RData/ST_GSE279181/02_QC",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)
# Filtered
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/02_QC",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#----------------------------- Visualize results ------------------------------#
################################################################################
nr_feats <- Giotto::spatPlot2D(gobject = combo,
                              cell_color = "nr_feats",
                              color_as_factor = FALSE,
                              gradient_style = "divergent",
                              gradient_midpoint = 1000,
                              group_by = "list_ID", 
                              group_by_subset = c("MS197D", "MS229", "MS377I", 
                                                  "MS94", "MS497T", "MS549H"),
                              point_size = 1,
                              point_alpha = 0.7,
                              cow_n_col = 6,
                              axis_text = 7,
                              title = "Número de genes por spot") 

filter200 <- Giotto::spatPlot2D(gobject = combo,
                               cell_color = "filter_200",
                               group_by = "list_ID", 
                               group_by_subset = c("MS197D", "MS229", "MS377I", 
                                                   "MS94", "MS497T", "MS549H"),
                               point_size = 1,
                               point_alpha = 0.7,
                               cow_n_col = 6,
                               axis_text = 7,
                               title = "Spots menos de 200 genes") 

filter100 <- Giotto::spatPlot2D(gobject = combo,
                              cell_color = "filter_100",
                              group_by = "list_ID", 
                              group_by_subset = c("MS197D", "MS229", "MS377I", 
                                                  "MS94", "MS497T", "MS549H"),
                              point_size = 1,
                              point_alpha = 0.7,
                              cow_n_col = 6,
                              axis_text = 7,
                              title = "Spots menos de 100 genes") 

filter100_all <- Giotto::spatPlot2D(gobject = combo,
                                    cell_color = "filter_100",
                                    group_by = "list_ID", 
                                    point_size = 1,
                                    point_alpha = 0.7,
                                    cow_n_col = 6,
                                    axis_text = 7,
                                    title = "Spots menos de 100 genes") 

mito_ratio <- Giotto::spatPlot2D(gobject = combo,
                               cell_color = "mito_ratio",
                               color_as_factor = FALSE,
                               gradient_style = "divergent",
                               gradient_midpoint = 20,
                               group_by = "list_ID", 
                               group_by_subset = c("MS197D", "MS229", "MS377I", 
                                                   "MS94", "MS497T", "MS549H"),
                               point_size = 1,
                               point_alpha = 0.7,
                               cow_n_col = 6,
                               axis_text = 7,
                               title = "Porcentaje de genes mitocondriales") 

ggsave("../Graficos/nr_feats.png",
       plot = nr_feats,
       width = 25,
       height = 5,
       dpi = 300)

ggsave("../Graficos/filter200.png",
       plot = filter200,
       width = 25,
       height = 5,
       dpi = 300)

ggsave("../Graficos/filter100.png",
       plot = filter100,
       width = 25,
       height = 5,
       dpi = 300)

ggsave("../Graficos/filter100_all.png",
       plot = filter100_all,
       width = 18,
       height = 10,
       dpi = 300)

ggsave("../Graficos/mito_ratio.png",
       plot = mito_ratio,
       width = 25,
       height = 5,
       dpi = 300)


################################################################################
#----------------------------------- Plots ------------------------------------#
################################################################################
# Spot statistics
cell_stats <- Giotto::pDataDT(combo)
cell_stats <- cell_stats %>%
  arrange(Lesion.type, list_ID) %>%  
  mutate(list_ID = factor(list_ID, levels = unique(list_ID)))

# Gene statistics
sample_names <- names(giotto_objects)
feat_stats <- lapply(seq_along(giotto_objects), function(i) {
  obj <- giotto_objects[[i]]
  sample_name <- sample_names[i]
  
  feature_data <- fDataDT(obj) %>%
    mutate(list_ID = sample_name)
  return(feature_data)
}) %>% bind_rows()

cell_stats_unique <- cell_stats %>%
  select(list_ID, Lesion.type) %>%
  distinct()

feat_10 <- fDataDT(combo) %>%
  select(feat_ID, filter_10) %>%
  distinct()

feat_stats <- feat_stats %>%
  left_join(cell_stats_unique, by = "list_ID") %>%
  left_join(feat_10, by = "feat_ID") %>%
  arrange(Lesion.type, list_ID) %>%  
  mutate(list_ID = factor(list_ID, levels = unique(list_ID)))

# Colors for the plots
fill_colors <- Giotto::getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 1234)

# Style of plots
theme_common <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10))

# Genes per spot
pspots <- ggplot(cell_stats, aes(x = list_ID, 
                                 y = nr_feats, 
                                 fill = Lesion.type)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_quasirandom(data = subset(cell_stats, filter_100 == 0),
                   aes(color = Lesion.type), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +
  labs(x = "Muestra", 
       y = "Nº de genes detectados por spot", 
       title = "Genes detectados en cada muestra") +
  scale_fill_manual(values = fill_colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = fill_colors, 
                     guide = "none") +
  theme_common

# Library size
plibsiz <- ggplot(cell_stats, aes(x = list_ID, 
                                  y = total_expr, 
                                  fill = Lesion.type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion.type), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(x = "Muestra", 
       y = "Tamaño de librería", 
       title = "Tamaño de librería por muestra") +
  scale_fill_manual(values = fill_colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = fill_colors, 
                     guide = "none") +
  theme_common

# Mitochondrial counts
pmito <- ggplot(cell_stats, aes(x = list_ID, 
                                y = mito_ratio, 
                                fill = Lesion.type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion.type), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(x = "Muestra", 
       y = "% de genes mitocondriales por spot", 
       title = "Porcentaje de genes mitocondriales en cada muestra") +
  scale_fill_manual(values = fill_colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = fill_colors, 
                     guide = "none") +
  theme_common

# Spots per gene
pgenes <- ggplot(feat_stats, aes(x = list_ID, 
                                 y = log10(nr_cells), 
                                 fill = Lesion.type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(data = subset(feat_stats, filter_10 == 0),
                   aes(color = Lesion.type), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +
  labs(x = "Muestra", 
       y = "Nº de spots por gen (log10)", 
       title = "Spots en los que se detecta cada gen por muestra") +
  scale_fill_manual(values = fill_colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = fill_colors, 
                     guide = "none") +
  theme_common

# Combine all plots and save
QC_plot <- (pspots + plot_spacer() + 
              plibsiz + plot_layout(widths = c(1, 0.1, 1))) / 
  (pmito + plot_spacer() + pgenes + 
     plot_layout(widths = c(1, 0.1, 1))) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("../Graficos/QC_ST.png",
       plot = QC_plot,
       width = 11,
       height = 15,
       dpi = 300)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#