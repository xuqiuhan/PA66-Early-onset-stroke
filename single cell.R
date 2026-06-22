
# Load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
# Paths
path_gse174574 <- "E:/GSE174574/date"
path_gse247474 <- "E:/GSE247474/date"
# Read 10X data
dir_name <- c(
  list.files(path_gse174574),
  list.files(path_gse247474)
)
scRNAlist <- list()
for (i in seq_along(dir_name)) {
  full_path <- file.path("E:/衰老人脑单核/date1", dir_name[i])
  
  if (!dir.exists(full_path)) {
    message("Directory not found: ", full_path)
    next
  }
  
  counts <- Read10X(data.dir = full_path)
  scRNAlist[[i]] <- CreateSeuratObject(
    counts,
    project = dir_name[i],
    min.cells = 3,
    min.features = 300
  )
}
# Merge samples
scRNA_merge <- merge(scRNAlist[[1]], y = scRNAlist[-1])
scRNA_merge[["RNA"]] <- JoinLayers(scRNA_merge[["RNA"]])
# QC and filtering
scRNA_merge[["percent.mt"]] <-
  PercentageFeatureSet(scRNA_merge, pattern = "^mt-")
VlnPlot(scRNA_merge,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
scRNA_merge <- subset(
  scRNA_merge,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)
# Standard preprocessing
scRNA_merge <- NormalizeData(scRNA_merge)
scRNA_merge <- FindVariableFeatures(scRNA_merge, nfeatures = 2000)
scRNA_merge <- ScaleData(scRNA_merge)
scRNA_merge <- RunPCA(scRNA_merge, npcs = 30)
# Harmony integration
scRNA_harmony <- RunHarmony(
  scRNA_merge,
  group.by.vars = "orig.ident"
)
ElbowPlot(scRNA_harmony)
scRNA_harmony <- FindNeighbors(
  scRNA_harmony,
  reduction = "harmony",
  dims = 1:20
)
scRNA_harmony <- FindClusters(
  scRNA_harmony,
  resolution = 0.5
)
scRNA_harmony <- RunUMAP(
  scRNA_harmony,
  reduction = "harmony",
  dims = 1:20
)
# Annotate experimental groups
group <- substring(scRNA_harmony$orig.ident, 1, 5)
group <- dplyr::recode(
  group,
  MCAO1 = "MCAO",
  MCAO2 = "MCAO",
  MCAO3 = "MCAO",
  shame1 = "sham",
  shame2 = "sham",
  shame3 = "sham"
)
scRNA_harmony$group.ident <- group
# Astrocyte identification
FeaturePlot(scRNA_harmony,
            features = c("Aqp4", "Slc1a2", "Ephx2"))
ast_ident <- "2"  # adjust based on UMAP
scRNA_AC <- subset(scRNA_harmony, idents = ast_ident)
# Astrocyte subclustering
scRNA_AC <- NormalizeData(scRNA_AC)
scRNA_AC <- FindVariableFeatures(scRNA_AC, nfeatures = 2000)
scRNA_AC <- ScaleData(scRNA_AC)
scRNA_AC <- RunPCA(scRNA_AC, npcs = 20)
scRNA_AC_harmony <- RunHarmony(
  scRNA_AC,
  group.by.vars = "orig.ident"
)
scRNA_AC_harmony <- FindNeighbors(
  scRNA_AC_harmony,
  reduction = "harmony",
  dims = 1:10
)
scRNA_AC_harmony <- FindClusters(
  scRNA_AC_harmony,
  resolution = 0.1
)
scRNA_AC_harmony <- RunUMAP(
  scRNA_AC_harmony,
  reduction = "harmony",
  dims = 1:10
)
# Astrocyte marker expression
ASCmarker <- c(
  "Glul", "Gja1", "Cst3", "Cldn10", "Ephx2",
  "Sparc", "Slc6a9", "Atp1a2", "Slc6a11",
  "Slc7a10", "S100a10", "Cd44", "Serpina3n",
  "Gfap", "Vim"
)
DotPlot(
  scRNA_AC_harmony,
  features = ASCmarker,
  group.by = "RNA_snn_res.0.1"
) +
  theme_classic()
# Differential expression
ACmarker <- FindMarkers(
  scRNA_AC_harmony,
  ident.1 = "AC3",
  ident.2 = "AC1",
  group.by = "RNA_snn_res.0.1",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
ACmarker$gene <- rownames(ACmarker)
# Volcano plot
ggplot(ACmarker, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(
    aes(color = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE ~ "NS"
    )),
    size = 1.2
  ) +
  scale_color_manual(
    values = c("Up" = "#E64B35FF",
               "Down" = "#4DBBD5FF",
               "NS" = "grey80")
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted p"
  ) +
  theme_classic()

# Define color palette
cols <- c(
  "MCAO" = "#E64B35FF",
  "sham" = "#4DBBD5FF"
)

# Violin plot for Ephx2
VlnPlot(
  object = scRNA_harmony,
  features = "Ephx2",
  group.by = "group.ident",
  layer = "data",
  pt.size = 0.1,
  cols = cols
) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    y = "Ephx2 expression (log-normalized)",
    title = "Astrocyte Ephx2 expression"
  )
