# ST Normalization
# Author: Natalia del Rey Díez
# Date: 27/03/2024

################################################################################
#------------------------------ Install packages ------------------------------#
################################################################################
# Matrix
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")}


################################################################################
#------------------------------ Loading packages ------------------------------#
################################################################################
library(Giotto)
library(ggplot2)
library(Matrix)
library(ggbeeswarm)
library(patchwork)


################################################################################
#--------------------------- Loading Giotto objects ---------------------------#
################################################################################
combo_filt <- Giotto::loadGiotto(path_to_folder = "../RData/ST_GSE279181/02_QC/combo_filt",
                                 reconnect_giottoImage = TRUE,
                                 init_gobject = TRUE,
                                 verbose = TRUE)


################################################################################
#--------------------------- Default normalization ----------------------------#
################################################################################
combo_filt = Giotto::processExpression(gobject = combo_filt, 
                                       param = normParam("default"), 
                                       name = "default",
                                       expression_values = "raw")


################################################################################
#--------------------------- Pearson normalization ----------------------------#
################################################################################
combo_filt = Giotto::processExpression(gobject = combo_filt, 
                                       param = normParam("pearson"), 
                                       name = "pearson",
                                       expression_values = "raw")


################################################################################
#------------------------- TF-IDF + L2 normalization --------------------------#
################################################################################
combo_filt = Giotto::processExpression(gobject = combo_filt, 
                                       param = list(normParam("tf-idf"),
                                                    normParam("l2")), 
                                       name = "tfidf",
                                       expression_values = "raw")


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = combo_filt,
                   foldername = "combo_filt",
                   dir = "../RData/ST_GSE279181/03_Normalizacion",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)


################################################################################
#----------------------------------- Plots ------------------------------------#
################################################################################
# Normalized matrix
default = Giotto::getExpression(combo_filt, 
                                values = "default", 
                                output = "matrix")

pearson = Giotto::getExpression(combo_filt, 
                                values = "pearson", 
                                output = "matrix")

tfidf = Giotto::getExpression(combo_filt, 
                              values = "tfidf", 
                              output = "matrix")

# Calculate library size
library_sizes <- data.frame(Normalization = rep(c("default", "pearson", "tfidf"), 
                                                each = ncol(default)),
                            LibrarySize = c(colSums(default),
                                            colSums(pearson),
                                            colSums(tfidf)),
  Sample = rep(pDataDT(combo_filt)$list_ID, 3),
  Lesion = rep(pDataDT(combo_filt)$Lesion.type, 3))

# Plot
# Colors for the plots
colors <- Giotto::getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 1234)

# Style of plots
theme_common <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10))

# Default
default <- ggplot(subset(library_sizes, Normalization == "default"), 
                  aes(x = Sample,  
                      y = LibrarySize,
                      fill = Lesion)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(title = "Log-Normalización por tamaño de librería",
       x = "Muestra",
       y = "Tamaño de librería") +
  scale_fill_manual(values = colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = colors, 
                     guide = "none") +
  theme_common

# Pearson
pearson <- ggplot(subset(library_sizes, Normalization == "pearson"), 
                  aes(x = Sample,  
                      y = LibrarySize,
                      fill = Lesion)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(title = "Normalización por residuos de Pearson",
       x = "Muestra",
       y = "Tamaño de librería") +
  scale_fill_manual(values = colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = colors, 
                     guide = "none") +
  theme_common

# TF-IDF
tfidf <- ggplot(subset(library_sizes, Normalization == "tfidf"),
                aes(x = Sample,  
                    y = LibrarySize,
                    fill = Lesion)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(title = "Normalización TF-IDF + L2",
       x = "Muestra",
       y = "Tamaño de librería") +
  scale_fill_manual(values = colors, 
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(values = colors, 
                     guide = "none") +
  theme_common

# Combine all plots and save
default <- default +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pearson <- pearson +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Norm_plot <- (default / pearson / tfidf) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


ggsave("../Graficos/Norm_ST.png",
       plot = Norm_plot,
       width = 10,
       height = 7,
       dpi = 300)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#