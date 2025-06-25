# snRNA-Seq
# Author: Natalia del Rey Díez
# Date: 20/04/2024

################################################################################
#------------------------------ Install packages ------------------------------#
################################################################################
# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}

# zellkonverter
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter")}

# SingleCellExperiment
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")}

# scater
if (!requireNamespace("scater", quietly = TRUE)) {
  BiocManager::install("scater")}

# scDblFinder
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")}

# tibble
if (!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble")}

# tidyr
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")}


################################################################################
#----------------------------- Loading packages -------------------------------#
################################################################################
library(zellkonverter)  
library(SingleCellExperiment) 
library(scater)
library(scDblFinder)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggbeeswarm)
library(patchwork)
library(Giotto)


################################################################################
#----------------------------- Create SCE object ------------------------------#
################################################################################
sn_atlas_sce = zellkonverter::readH5AD(file = "../snRNA_GSE279180/GSE279180_sn_atlas.h5ad",
                                       reader = "R")

# Rename the count array
assayNames(sn_atlas_sce)[assayNames(sn_atlas_sce) == "X"] <- "counts"


################################################################################
#---------------------------------- Save SCE ----------------------------------#
################################################################################
system("mkdir -p ../RData/snRNA_GSE279180/00_SCE")
save(sn_atlas_sce, file = "../RData/snRNA_GSE279180/00_SCE/sn_atlas_sce.rda")


################################################################################
#------------------------------ Quality Control -------------------------------#
##----------------------- Number of nuclei per sample ------------------------##
################################################################################
cell_counts <- table(colData(sn_atlas_sce)$sample_id, 
                     colData(sn_atlas_sce)$celltype)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Sample", "Cell_Type", "Cell_Count")

# Colors for the plots
cell_colors <- Giotto::getDistinctColors(9)

# Style of plots
theme_common <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10))

pcell <- ggplot(cell_counts_df, aes(x = Sample, y = Cell_Count, fill = Cell_Type)) +
  geom_bar(stat = "identity") + 
  xlab("Muestras") + 
  ylab("Nº de núcleos") +
  ggtitle("Número de núcleos de cada tipo celular por muestra") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cell_colors, 
                    guide = guide_legend(title = NULL)) +
  theme_common

ggsave("../Graficos/snRNA/sn_cell_type.png",
       plot = pcell,
       width = 8,
       height = 4,
       dpi = 300)


################################################################################
#------------------------------ Quality Control -------------------------------#
##--------------------------- Gene quality metrics ---------------------------##
################################################################################
rowData(sn_atlas_sce)$nCells <- Matrix::rowSums(counts(sn_atlas_sce)>0)

# Number of nuclei in which each gene is expressed per sample
count_matrix <- counts(sn_atlas_sce)  
sample_info <- colData(sn_atlas_sce)$sample_id  

expressed <- count_matrix > 0

calculate_sample_expression <- function(sample_name) {
  sample_cells <- which(sample_info == sample_name)
  rowSums(expressed[, sample_cells, drop = FALSE])
}

sample_counts <- lapply(unique(sample_info), calculate_sample_expression)
names(sample_counts) <- unique(sample_info)

nuc_per_gene <- as.data.frame(do.call(cbind, sample_counts))
rownames(nuc_per_gene) <- rownames(sn_atlas_sce)

sample_groups <- data.frame(
  sample = colData(sn_atlas_sce)$sample_id,
  Lesion.type = colData(sn_atlas_sce)$lesion_type) %>%
  distinct()

nuc_long <- nuc_per_gene %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, 
                      names_to = "sample", 
                      values_to = "nuclei_count") %>%
  left_join(sample_groups, by = "sample") 

# Colors for the plots
colors <- Giotto::getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 1234)

pgenes <- ggplot(nuc_long, aes(x = sample, 
                               y = log10(nuclei_count), 
                               fill = Lesion.type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_quasirandom(aes(color = Lesion.type), 
                   width = 0.1, 
                   alpha = 0.7, 
                   size = 0.5) +  
  labs(x = "Muestra", 
       y = "Nº de núcleos por gen", 
       title = "Núcleos en los que se detecta cada gen por muestra") +
  scale_fill_manual(values = colors, guide = guide_legend(title = NULL)) +
  scale_color_manual(values = colors, guide = "none") +
  theme_common

ggsave("../Graficos/snRNA/cell_per_gene.png",
       plot = pgenes,
       width = 10,
       height = 5,
       dpi = 300)


################################################################################
#------------------------------ Quality Control -------------------------------#
##-------------------------- Nuclei quality metrics --------------------------##
################################################################################
sn_atlas_sce <- scater::addPerCellQC(sn_atlas_sce)


################################################################################
#------------------------------ Quality Control -------------------------------#
##-------------------- Percentage of mitochondrial genes ---------------------##
################################################################################
mito_genes <- rownames(sn_atlas_sce)[grep("^MT-", rownames(sn_atlas_sce))]
sn_atlas_sce$mito_ratio <- 
  Matrix::colSums(counts(sn_atlas_sce)[mito_genes, ]) / sn_atlas_sce$sum


################################################################################
#------------------------------ Quality Control -------------------------------#
##--------------------------------- Doublets ---------------------------------##
################################################################################
set.seed(123)
sn_atlas_sce <- scDblFinder::scDblFinder(sn_atlas_sce, 
                                         verbose = FALSE, 
                                         samples = "sample_id")


################################################################################
#---------------------------------- Save SCE ----------------------------------#
################################################################################
system("mkdir -p ../RData/snRNA_GSE279180/01_QC")
save(sn_atlas_sce, file = "../RData/snRNA_GSE279180/01_QC/sn_atlas_sce.rda")


################################################################################
#------------------------------ Quality Control -------------------------------#
##---------------------------------- Plots -----------------------------------##
################################################################################
# Genes per nuclei
pnuclei <- scater::plotColData(sn_atlas_sce, 
                               y = "detected", 
                               x = "sample_id", 
                               colour_by = "lesion_type",
                               scattermore = TRUE,
                               point_size = 0.5) +
  geom_hline(yintercept = 199, linetype = "dashed", color = "black") +
  scale_color_manual(values = colors, 
                     guide = guide_legend(title = NULL)) +
  xlab("Muestras") + 
  ylab("Nº de genes expresados por núcleo") +
  ggtitle("Genes detectados por muestra") + 
  theme_common

# Library size
plibsize <- scater::plotColData(sn_atlas_sce, 
                                y = "sum", 
                                x = "sample_id", 
                                colour_by = "lesion_type",
                                scattermore = TRUE,
                                point_size = 0.5) +
  scale_color_manual(values = colors, 
                     guide = guide_legend(title = NULL)) +
  xlab("Muestras") + 
  ylab("Tamaño de librería") +
  ggtitle("Tamaño de librería por muestra") + 
  theme_common 

# Mitochondrial counts
pmito <- scater::plotColData(sn_atlas_sce, 
                             y = "mito_ratio", 
                             x = "sample_id", 
                             colour_by = "lesion_type",
                             scattermore = TRUE,
                             point_size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_color_manual(values = colors, 
                     guide = guide_legend(title = NULL)) +
  xlab("Muestras") + 
  ylab("Proporción de genes mitocondriales por núcleo") +
  ggtitle("Proporción de genes mitocondriales por muestra") + 
  theme_common


# Combine all plots and save
QC_plot <- (pnuclei + plot_spacer() + plibsize + plot_spacer() + pmito) +
  plot_layout(widths = c(1.2, 0.05, 1.2, 0.05, 1.5), guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("../Graficos/snRNA/QC_snRNA.png",
       plot = QC_plot,
       width = 12,
       height = 5,
       dpi = 300)


################################################################################
#---------------------------- Create Giotto object ----------------------------#
################################################################################
sn_expression = assays(sn_atlas_sce)$counts 
sn_metadata = SingleCellExperiment::colData(sn_atlas_sce)
sn_metadata = data.table::as.data.table(sn_metadata,
                                        keep.rownames = "cell_ID")

giotto_SC = Giotto::createGiottoObject(expression = sn_expression,
                                       cell_metadata = sn_metadata)


################################################################################
#------------------------------- Normalization --------------------------------#
################################################################################
giotto_SC = Giotto::normalizeGiotto(giotto_SC)


################################################################################
#------------------------------------- HVG-------------------------------------#
################################################################################
giotto_SC = calculateHVF(giotto_SC,
                         method = "cov_groups")

table(fDataDT(giotto_SC)$hvf)


################################################################################
#---------------------------- Dimension Reduction -----------------------------#
################################################################################
# PCA
giotto_SC <- Giotto::runPCA(giotto_SC, 
                            expression_values = "normalized",
                            reduction = "cells",
                            feats_to_use = "hvf",
                            center = TRUE,
                            scale_unit = TRUE,
                            ncp = 50) 
Giotto::plotPCA(giotto_SC, 
                cell_color = "celltype",
                cell_color_code = cell_colors,
                save_plot = TRUE,
                save_param = list(save_dir = "../Graficos",
                                  save_folder = "snRNA/",
                                  save_name = "PCA"),
                title = "PCA snRNA-Seq")

Giotto::screePlot(giotto_SC,
                  dim_reduction_name = "pca",
                  reduction = "cells",
                  ncp = 50,
                  feats_to_use = "hvf",
                  save_plot = TRUE,
                  save_param = list(save_dir = "../Graficos",
                                    save_folder = "snRNA/",
                                    save_name = "screePlot"),
                  title = "screePlot snRNA-Seq")


# UMAP
giotto_SC <- Giotto::runUMAP(giotto_SC, 
                             expression_values = "normalized",
                             reduction = "cells",
                             feats_to_use = "hvf",
                             dimensions_to_use = 1:10) 

Giotto::plotUMAP(giotto_SC, 
                 cell_color = "celltype",
                 title = "UMAP snRNA-Seq",
                 cell_color_code = cell_colors,
                 save_plot = TRUE,
                 save_param = list(save_dir = "../Graficos",
                                   save_folder = "snRNA/",
                                   save_name = "UMAP"))


################################################################################
#-------------------------------- Save Giotto ---------------------------------#
################################################################################
Giotto::saveGiotto(gobject = giotto_SC,
                   foldername = "giotto_SC",
                   dir = "../RData/snRNA_GSE279180/02_Processed",
                   method = "RDS",
                   overwrite = TRUE,
                   export_image = TRUE,
                   image_filetype = "PNG",
                   include_feat_coord = TRUE,
                   verbose = TRUE)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#