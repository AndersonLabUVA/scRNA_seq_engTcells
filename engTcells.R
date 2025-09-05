tcells <- subset(engTcellsOvarianCancer, newLabels =="T cell")
# Transform the ADT data to a log scale
tcells@assays$ADT@data <- log1p(tcells@assays$ADT@data)
RidgePlot(tcells, group.by = "orig.ident", features=c('CD8a')) + coord_cartesian(xlim = c(0,5))
RidgePlot(tcells, group.by = "orig.ident", features=c('CD90/CD90.1'))

# Perform standard Seurat normalization steps

tcells <- NormalizeData(tcells, normalization.method = "LogNormalize") %>%
  FindVariableFeatures(nfeatures = 10000) %>% ScaleData() %>% RunPCA()

ElbowPlot(tcells, ndims = 50)
pdf(file = paste0(work_dir_figures, project_name, date_processed, fileAdd,
                  "tcells_SeuratObj_Merged_ElbowPlot.pdf"))
ElbowPlot(tcells, ndims = 50)
dev.off()

tcells <- FindNeighbors(tcells, dims = 1:20)
tcells <- FindClusters(tcells, resolution = 1)
tcells <- RunUMAP(tcells, dims = 1:20)
DimPlot(tcells, reduction = "umap",group.by='seurat_clusters',label=TRUE)

plot_cells_by_metadata(tcells, "tcells", "orig.ident", work_dir_figures, project_name, 
                       date_processed, fileAdd, reduction_type = "umap")

Idents(tcells) <- "seurat_clusters"
rna.markers<-FindAllMarkers(tcells, assay='RNA', min.pct = 0.25)
write.csv(rna.markers, file= paste0(work_dir_outputs, project_name,
                                    date_processed,"tcells_RNA_FindAllMarkers.csv"))

# Save seurat object
saveRDS(tcells, file = paste0(work_dir_rds, project_name, date_processed, fileAdd, "tcells.Rds"))

Idents(tcells) <- "orig.ident"
rna.markers<-FindAllMarkers(engTcellsOvarianCancer, assay='RNA', min.pct = 0.25)
write.csv(rna.markers, file= paste0(work_dir_outputs, project_name,                                    date_processed,"RNA_FindAllMarkers_tCells_by_Orig.Ident.csv"))

```

``` {r engTcellsOnly}
RidgePlot(tcells, features = "CD8a", group.by = "orig.ident") + ggplot2::theme(legend.position = "none")
pdf(file = paste0(work_dir_figures, project_name, date_processed, fileAdd, "tcells.RidgPlot.CD8a.pdf"), 
    height=10, width = 10)
print(RidgePlot(tcells, features = "CD8a", group.by = "orig.ident")) + ggplot2::theme(legend.position = "none")
dev.off()

RidgePlot(tcells, features = "CD90/CD90.1", group.by = "orig.ident") + ggplot2::theme(legend.position = "none")
pdf(file = paste0(work_dir_figures, project_name, date_processed, fileAdd, "tcells.RidgPlot.CD90.1.pdf"), 
    height=10, width = 10)
print(RidgePlot(tcells, features = "CD90/CD90.1", group.by = "orig.ident")) + ggplot2::theme(legend.position = "none")
dev.off()

adt_expression_cd8a <- tcells@assays$ADT@data["CD8a", ]
adt_expression_cd90 <- tcells@assays$ADT@data["CD90/CD90.1", ]

# Create a data frame
df <- data.frame(CD8Expression = adt_expression_cd8a, CD90Expression = adt_expression_cd90)
# Add a column for the library identity
df$Library <- tcells$orig.ident

# Calculate the percentage of cells within the quadrant for each library
df$InQuadrant <- ifelse(df$CD8Expression >= 0.5 & df$CD90Expression >= 2, 1, 0)
percentages <- aggregate(InQuadrant ~ Library, df, function(x) sum(x)/length(x)*100)
df <- merge(df, percentages, by = "Library")

# Print the percentages
print(percentages)

# Plot the data
p1 <- ggplot(df, aes(x = CD8Expression, y = CD90Expression)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("ADT Expression of Cd8a") +
  ylab("ADT Expression of CD90/CD90.1") +
  ggtitle("Scatter plot of CD8a and CD90/CD90.1 Expression for All Libraries") +
  geom_vline(xintercept = 0.5, linetype="dashed", color = "red") +  # Add vertical line at gene expression = 1
  geom_hline(yintercept = 2, linetype="dashed", color = "red")  # Add horizontal line at ADT expression = 1.5

# Combine the individual library plots and the all libraries plot
p2 <- ggplot(df, aes(x = CD8Expression, y = CD90Expression)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("ADT Expression of Cd8a") +
  ylab("ADT Expression of CD90/CD90.1") +
  ggtitle("Scatter plot of CD8a and CD90/CD90.1 Expression") +
  geom_vline(xintercept = 0.5, linetype="dashed", color = "red") +  # Add vertical line at gene expression = 1
  geom_hline(yintercept = 2, linetype="dashed", color = "red") +  # Add horizontal line at ADT expression = 1.5
  facet_wrap(~ Library) # Separate plot for each library

# Display the plots
gridExtra::grid.arrange(p1, p2, ncol = 1)

pdf(file = paste0(work_dir_figures, project_name, date_processed, fileAdd, "Cd8a_VS_Thy1.1_Gating.pdf"), height=10, width = 10)
gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

#from this, get just the High Thy1.1 and CD8a T-cells
Thy1CD8 <- (GetAssayData(tcells, assay = "ADT", layer = "data")["CD90/CD90.1", ] >= 2.0 & 
              GetAssayData(tcells, assay = "ADT", layer = "data")["CD8a", ] >= 0.5)

tcells <- AddMetaData(tcells,Thy1CD8, col.name = 'Thy1CD8')
tcellsExo <- subset(tcells,subset=Thy1CD8==TRUE)
table(tcellsExo$orig.ident)
Idents(tcellsExo) <- tcellsExo$orig.ident
PreTransferTcells <- tcellsExo[, tcellsExo$orig.ident %in% c("Pre-transfer eng. T cells")]

tcellsExowithoutPreTransfer <- tcellsExo[, tcellsExo$orig.ident %in% c("eng. T cells + Isotype", 
                                                                       "eng. T cells + PD1",
                                                                       "eng. T cells + PD1 and Lag3", 
                                                                       "eng. T cells + PD1 and Tim3",
                                                                       "eng. T cells + Triple AB treatment")]
table(tcellsExowithoutPreTransfer$orig.ident)

# Perform standard Seurat normalization steps
tcellsExowithoutPreTransfer <- NormalizeData(tcellsExowithoutPreTransfer, normalization.method = "LogNormalize") %>%
  FindVariableFeatures(nfeatures = 10000) %>% ScaleData() %>% RunPCA()

ElbowPlot(tcellsExowithoutPreTransfer)

tcellsExowithoutPreTransfer <- FindNeighbors(tcellsExowithoutPreTransfer, dims = 1:20)
tcellsExowithoutPreTransfer <- FindClusters(tcellsExowithoutPreTransfer, resolution = 1)
tcellsExowithoutPreTransfer <- RunUMAP(tcellsExowithoutPreTransfer, dims = 1:20)
DimPlot(tcellsExowithoutPreTransfer, reduction = "umap",group.by='seurat_clusters',label=TRUE)
DimPlot(tcellsExowithoutPreTransfer, reduction = "umap",group.by='orig.ident',label=FALSE)

plot_cells_by_metadata(tcellsExowithoutPreTransfer, "tcellsExowithoutPreTransfer","seurat_clusters", 
                       work_dir_figures, project_name, date_processed, fileAdd)
plot_cells_by_metadata(tcellsExowithoutPreTransfer, "tcellsExowithoutPreTransfer","orig.ident",
                       work_dir_figures, project_name, date_processed, fileAdd)

DefaultAssay(tcellsExowithoutPreTransfer) <- "RNA"
rna.markers<-FindAllMarkers(tcellsExowithoutPreTransfer, assay='RNA', only.pos = TRUE, min.pct = 0.25)

write.csv(rna.markers, file= paste0(work_dir_outputs, project_name,
                                    date_processed,"tcellsExowithoutPreTransfer_RNAFindAllMarkers.csv"))

rna.markers.pct10 <- FindAllMarkers(tcellsExowithoutPreTransfer, assay='RNA', min.pct = 0.10, test.use = "MAST")
write.csv(rna.markers.pct10, file= paste0(work_dir_outputs, project_name,
                                          date_processed,"tcellsExowithoutPreTransfer_rna.markers.pct10_MAST_test.csv"))

rna.markers.pct15 <- FindAllMarkers(tcellsExowithoutPreTransfer, assay='RNA', min.pct = 0.15, test.use = "MAST")
write.csv(rna.markers.pct15, file= paste0(work_dir_outputs, project_name,
                                          date_processed,"tcellsExowithoutPreTransfer_rna.markers.pct15_MAST_test.csv"))

rna.markers.pct20 <- FindAllMarkers(tcellsExowithoutPreTransfer, assay='RNA', min.pct = 0.20, test.use = "MAST")
write.csv(rna.markers.pct20, file= paste0(work_dir_outputs, project_name,
                                          date_processed,"tcellsExowithoutPreTransfer_rna.markers.pct20_MAST_test.csv"))

rna.markers.pct25 <- FindAllMarkers(tcellsExowithoutPreTransfer, assay='RNA', min.pct = 0.25, test.use = "MAST")
write.csv(rna.markers.pct25, file= paste0(work_dir_outputs, project_name,
                                          date_processed,"tcellsExowithoutPreTransfer_rna.markers.pct25_MAST_test.csv"))

saveRDS(tcellsExowithoutPreTransfer, file = paste0(work_dir_rds, project_name, date_processed, 
                                                   fileAdd, "tcellsExowithoutPreTransfer.Rds"))
