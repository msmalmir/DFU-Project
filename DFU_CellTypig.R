#install.packages("devtools")
library(devtools)
#devtools::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)
# reading the data
#DFU_data <- read.csv(file = 'Gene-expression.csv')

#load the h5Seurat file 
#DFU_all <-  LoadH5Seurat("DFU-foot-24samples.h5seurat")
#DFU_all <-  LoadH5Seurat("Foot-all.h5seurat")
#DFU_all <- readRDS("DFU_all_60dims.rds")
#DFU_all <- DFU_foot_new_filter_celltyping
#rm(DFU_foot_new_filter_celltyping)
# save the new Seurat object
#saveRDS(DFU_all, file = "DFU_foot_24samples_V2.rds")

# from these analysis we decided to remove "Healthy3",  "Healing DFU5" and do the analysis again
# I am going to to the top and remove them from the Seurat object and getting the results again
#length(DFU_all@meta.data$nCount_RNA[DFU_all$origin1 == 'Healthy3'])
#length(DFU_all@meta.data$nCount_RNA[DFU_all$origin1 == 'Healing DFU5'])
#length(DFU_all@meta.data$nCount_RNA[DFU_all$origin1 == 'Non Healing1'])

# Create a logical vector to identify the cells to keep (those NOT in "Healthy3" or "Healing DFU5" or 'Non Healing1')
#keep_cells <- rownames(DFU_all@meta.data[!(DFU_all@meta.data$origin1 %in% c("Healthy3", "Healing DFU5", 'Non Healing1')), ])
#length(keep_cells)
# Subset the Seurat object using the list of cells to keep
#DFU_all <- subset(DFU_all, cells = keep_cells)
DFU_all_selected <- DFU_all
DFU_all_selected$scSorter <- DFU_all$scSorter_newGM
DFU_all_selected@meta.data <- DFU_all_selected@meta.data %>%
  dplyr::rename(sc2Assign = sc2Assign_newGM)
#saveRDS(DFU_all_selected, file = "DFU_foot_24samples_selected.rds")
# plot
DimPlot(DFU_all_selected, reduction = 'umap', label = TRUE,
        group.by = 'sc2Assign', pt.size = 0.4,
        label.size = 4)  
DimPlot(DFU_all_selected, reduction = 'umap', label = TRUE,
        group.by = 'scSorter', pt.size = 0.4,
        label.size = 4)  


# plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE,
        group.by = 'seurat_clusters', pt.size = 0.4,
         label.size = 4)  

unique()

# extracting the count matrix, scale matrix and norm matrix for downstream process
# we can consider HVG for downstream process Normalize the data (if not already normalized)
DefaultAssay(DFU_all) <- 'RNA'
DFU_all <- NormalizeData(DFU_all)
# Find highly variable features (genes)
DFU_all <- FindVariableFeatures(DFU_all, selection.method = "vst", nfeatures = 2000)
#
count.matrix <- as.data.frame(as.matrix(DFU_all@assays$RNA$counts)) 
# Step 1: Get the list of highly variable genes from Seurat object
hvg_genes <- VariableFeatures(DFU_all)
count.matrix <- count.matrix[hvg_genes, ]


norm.matrix <- as.data.frame( as.matrix(DFU_all@assays$RNA$data))
norm.matrix <- norm.matrix[hvg_genes, ]

#scaled.matrix <- as.data.frame(as.matrix(DFU_all@assays$integrated$scale.data)) 
#scaled.matrix <- scaled.matrix[rownames(scaled.matrix) %in% hvg_genes, ]






# performing the bicluster using both expression level and clinical info we got for each cell
install.packages("biclust")
install.packages("pheatmap")  # For heatmap visualization

# Load the necessary libraries
library(biclust)
library(pheatmap)  # Optional, for visualization

# Assume `count_matrix` is your count matrix with genes as rows and cells as columns
# And `clinical_groups` is a vector that contains clinical group information for each cell (column)

# Biclustering using the Plaid model (you can experiment with other methods)
#count.matrix <- as.matrix(count.matrix)
scaled.matrix <- as.matrix(scaled.matrix)
norm.matrix <- as.matrix(norm.matrix)
biclust_result <- biclust(norm.matrix, method = BCSpectral())

# View the result
summary(biclust_result)



# creating the gene marker set 
marker.genes <- list(
                    Fibro = union(
                      c("CFD", "APOD", "PLA2G2A", "DCN", "CXCL14", "PTGDS", "SFRP2", "COMP", "FBLN1", "GSN"),
                      c("APOD", "CFD", "DCN", "COL1A1", "COL1A2", "SFRP2", "COL3A1", "PTGDS", "CXCL14", "COMP")
                    ), 
                     SMC1 = union(
                       c("ACTA2", "TAGLN", "RGS5", "TPM2", "MYL9", "C11orf96", "MYH11", "NDUFA4L2", "CALD1", "NR2F2"),
                       c("ACTA2", "TAGLN", "RGS5", "TPM2", "MYL9", "C11orf96", "CALD1", "IGFBP7", "MYH11", "TPM1")
                     ), 
                     SMC2 = c("CENPF", "PTTG1", "H2AFZ", "TUBA1B", "TOP2A", "STMN1", "MKI67", "UBE2S", "TIMP1", "HIST1H4C"),
                     VasEndo = union(
                       c("IFI27", "ACKR1", "PECAM1", "CLDN5", "PLVAP", "VWF", "SOX18", "GNG11", "RAMP2", "AQP1"),
                       c("IFI27", "ACKR1", "CLDN5", "AQP1", "GNG11", "PECAM1", "SPARCL1", "PLVAP", "VWF", "RAMP2")
                     ), 
                     `T-Lympho` = T_Lympho <- union(
                       c("CD69", "IL32", "CXCR4", "CD52", "LTB", "KLRB1", "IL7R", "PTPRC", "DUSP2", "TRAC"),
                       c("LTB", "IL7R", "IL32", "TRAC", "CD69", "CD3D", "TRBC2", "CD52", "CD79B", "CXCR4")
                     ), 
                     DiffKera = union(
                       c("KRT1", "KRT10", "KRT2", "DMKN", "KRTDAP", "CALML5", "LGALS7B", "LY6D", "DSP", "PERP"),
                       c("KRT10", "KRT2", "KRT1", "DMKN", "KRTDAP", "CALML5", "LGALS7B", "LY6D", "DSP", "PERP")
                     ), 
                     BasalKera = union(
                       c("KRT14", "KRT16", "KRT6A", "KRT5", "KRT17", "KRT6B", "KRT6C", "S100A2", "S100A8", "S100A9"),
                       c("KRT14", "KRT5", "KRT16", "KRT6A", "S100A2", "KRT17", "KRT6B", "KRT6C", "COL17A1", "KRT15")
                     ),
                     `HE-Fibro` = union(
                       c("MMP1", "COL1A1", "ASPN", "POSTN", "COL3A1", "COL1A2", "COL12A1", "TNC", "LUM", "FN1"),
                       c("PLA2G2A", "MMP1", "CHI3L1", "TIMP1", "SFRP4", "FTH1", "FN1", "CHI3L2", "MT2A", "LUM")
                     ), 
                     `B-Lympho` = B_Lympho <- union(
                       c("CD37", "LTB", "MS4A1", "CD74", "CD52", "HLA-DRA", "HLA-DPB1", "CD79A", "HLA-DQA1", "IGHM"),
                       c("IGHM", "IGKC", "CD37", "MS4A1", "CD79A", "CD79B", "CD52", "CD74", "LTB", "CD74")
                     ), 
                     `M1-Macro` = union(
                       c("LYZ", "IL1B", "HLA-DRA", "HLA-DPB1", "CXCL8", "HLA-DRB1", "HLA-DPA1", "SRGN", "HLA-DQA1", "CD74"),
                       c("S100A9", "LYZ", "S100A8", "FCN1", "CTSS", "MNDA", "S100A12", "VCAN", "TYROBP", "AIF1")
                     ),
                     NKT = union(
                       c("GNLY", "NKG7", "CCL5", "CCL4", "XCL2", "DUSP2", "CD69", "XCL1", "GZMA", "CCL3"),
                       c("GNLY", "NKG7", "CCL5", "GZMB", "GZMA", "CST7", "FGFBP2", "PRF1", "CCL4", "GZMH")
                     ), 
                     `M2-Macro` = union(
                       c("RNASE1", "C1QA", "C1QB", "SELENOP", "C1QC", "FTL", "CD74", "CD14", "AIF1", "HLA-DRA"),
                       c("RNASE1", "C1QA", "C1QB", "SELENOP", "C1QC", "FTL", "CD74", "CD14", "AIF1", "HLA-DRA")
                     ), 
                     Schwann = c("MPZ", "S100B", "GPM6B", "NRXN1", "PLP1", "CDH19", "VWA1", "CRYAB", "PMP22", "MBP"),
                     Melano = Melano <- union(
                       c("DCT", "TYRP1", "PMEL", "MLANA", "QPCT", "CYB561A3", "MITF", "APOE", "MFSD12", "PLP1"),
                       c("DCT", "TYRP1", "PMEL", "MLANA", "GPM6B", "PLP1", "QPCT", "PMP22", "MFSD12", "CYB561A3")
                     ),
                     Mast = union(
                       c("TPSB2", "TPSAB1", "CTSG", "HPGD", "HPGDS", "CPA3", "GATA2", "MS4A2", "LGALS3", "SERPINB1"),
                       c("TPSB2", "TPSAB1", "CTSG", "HPGD", "HPGDS", "CPA3", "GATA2", "MS4A2", "LGALS3", "ANXA1")
                     )
                    , 
                     `Sweat/Seba` = c("DCD", "SCGB2A2", "MUCL1", "SCGB1B2P", "SCGB1D2", "PIP", "AZGP1", "KRT19", "KRT7", "AQP5"),
                     LymphEndo = union(
                       c("CCL21", "TFF3", "MMRN1", "CLDN5", "GNG11", "TFPI", "PPFIBP1", "CAVIN2", "PROX1", "LMO2"),
                       c("CCL21", "TFF3", "MMRN1", "CLDN5", "GNG11", "TFPI", "PPFIBP1", "CAVIN2", "LYVE1", "LMO2")
                     ), 
                     Plasma = union(
                       c("IGKC", "IGLC2", "IGHA1", "IGHG1", "IGLC3", "IGHG2", "IGHG3", "JCHAIN", "IGHG4", "MZB1"),
                       c("IGKC", "IGLC2", "IGHA1", "IGHG1", "IGLC3", "IGHG2", "IGHG3", "JCHAIN", "IGHG4", "IGLC7")
                     )
)

#c("MMP1" ,   "COL1A1",  "ASPN",    "POSTN",   "COL3A1",  "COL1A2",  "COL12A1", "TNC",     "LUM",     "FN1",
#  "PLA2G2A","CHI3L1",  "TIMP1",   "SFRP4" ,  "FTH1" ,   "CHI3L2" , "MT2A")
#c("CD37",     "LTB" ,     "MS4A1",    "CD74" ,    "CD52",     "HLA-DRA" , "HLA-DPB1" ,"CD79A" ,
#  "HLA-DQA1", "IGHM" ,  "IGKC" ,    "CD79B" )
# Creating list of gene vectors for each cell type

marker.genes <- list(
  SMC1 = c("ACTA2", "TAGLN", "RGS5", "TPM2", "MYL9", "C11orf96", "MYH11", "NDUFA4L2", "CALD1", "NR2F2"),
  VasEndo = c("IFI27", "ACKR1", "PECAM1", "CLDN5", "PLVAP", "VWF", "SOX18", "GNG11", "RAMP2", "AQP1"),
  Fibro = c("CFD", "APOD", "PLA2G2A", "DCN", "CXCL14", "PTGDS", "SFRP2", "COMP", "FBLN1", "GSN"),
  `T-Lympho` = c("CD69", "IL32", "CXCR4", "CD52", "LTB", "KLRB1", "IL7R", "PTPRC", "DUSP2", "TRAC"),
  DiffKera = c("KRT1", "KRT10", "KRT2", "DMKN", "KRTDAP", "CALML5", "LGALS7B", "LY6D", "DSP", "PERP"),
  `M1-Macro` = c("LYZ", "IL1B", "HLA-DRA", "HLA-DPB1", "CXCL8", "HLA-DRB1", "HLA-DPA1", "SRGN", "HLA-DQA1", "CD74"),
  `HE-Fibro` = #c("MMP1" ,"COL1A1","ASPN","POSTN","COL3A1","COL1A2","COL12A1","TNC","LUM","FN1","PLA2G2A","CHI3L1", "TIMP1","SFRP4" ,"FTH1" ,"CHI3L2","MT2A"), 
               c("MMP1", "COL1A1", "ASPN", "POSTN", "COL3A1", "COL1A2", "COL12A1", "TNC", "LUM", "FN1"),
  BasalKera = c("KRT14", "KRT16", "KRT6A", "KRT5", "KRT17", "KRT6B", "KRT6C", "S100A2", "S100A8", "S100A9"),
  `M2-Macro` = c("RNASE1", "C1QA", "C1QB", "SELENOP", "C1QC", "FTL", "CD74", "CD14", "AIF1", "HLA-DRA"),
  Mast = c("TPSB2", "TPSAB1", "CTSG", "HPGD", "HPGDS", "CPA3", "GATA2", "MS4A2", "LGALS3", "SERPINB1"),
  NKT = c("GNLY", "NKG7", "CCL5", "CCL4", "XCL2", "DUSP2", "CD69", "XCL1", "GZMA", "CCL3"),
  `B-Lympho` = #c("CD37","LTB" ,"MS4A1","CD74" ,"CD52","HLA-DRA","HLA-DPB1","CD79A","HLA-DQA1","IGHM","IGKC","CD79B" ),
              c("CD37", "LTB", "MS4A1", "CD74", "CD52", "HLA-DRA", "HLA-DPB1", "CD79A", "HLA-DQA1", "IGHM"),
  LymphEndo = c("CCL21", "TFF3", "MMRN1", "CLDN5", "GNG11", "TFPI", "PPFIBP1", "CAVIN2", "PROX1", "LMO2"),
  Melano = c("DCT", "TYRP1", "PMEL", "MLANA", "QPCT", "CYB561A3", "MITF", "APOE", "MFSD12", "PLP1"),
  SMC2 = c("CENPF", "PTTG1", "H2AFZ", "TUBA1B", "TOP2A", "STMN1", "MKI67", "UBE2S", "TIMP1", "HIST1H4C"),
  `Sweat/Seba` = c("DCD", "SCGB2A2", "MUCL1", "SCGB1B2P", "SCGB1D2", "PIP", "AZGP1", "KRT19", "KRT7", "AQP5"),
  Schwann = c("MPZ", "S100B", "GPM6B", "NRXN1", "PLP1", "CDH19", "VWA1", "CRYAB", "PMP22", "MBP"),
  Plasma = c("IGKC", "IGLC2", "IGHA1", "IGHG1", "IGLC3", "IGHG2", "IGHG3", "JCHAIN", "IGHG4", "IGLC7")
)

# Print the list of vectors
print(cell_type_genes)


unique(x = DFU_all$origin1)

DefaultAssay(DFU_all) <- "RNA"
FeaturePlot(DFU_all, 	 'PLA2G2A', reduction = 'umap')
# c('DCN', 'CHI3L1', 'MMP1', 'MMP3', 'TNFAIP6')
# c('CENPF', 'PTTG1', 'MKI67', 'TOP2A')
length(marker.genes)
names(marker.genes)

#reemove the genes not included in the count matrix
geneNames <- rownames(count.matrix)
marker.genes = lapply( marker.genes, FUN = function(x) {x[x %in% geneNames]})



#colorCodes <- c('#c75a7d',  '#aed6f1', '#007cba',  '#5ac7b9', 
#                 '#d4ac0d',  '#89e300', 'red' , '#e34e18',  )
colorCodes <- c('#800000',  '#9A6324','#808000', '#aed6f1','#000075',  '#3cb44b','#89e300','#5ac7b9','#c75a7d', '#007cba','red' ,
'#e6194B', '#f58231',  '#ffe119','#42d4f4', '#4363d8', '#911eb4', '#f032e6', 'gray')  

names(colorCodes) <- c("Fibro",      "SMC1",       "SMC2",      "VasEndo",    "T-Lympho",   "DiffKera",   "BasalKera",  "HE-Fibro",
  "B-Lympho",  "M1-Macro",   "NKT",        "M2-Macro",   "Schwann",    "Melano",     "Mast",
  "Sweat/Seba", "LymphEndo",  "Plasma" , 'unknown') 


# cell typing using ScType
start = Sys.time()
scType_result <- SC_ScType( scaled.matrix, marker.genes, DFU_all )
# Convert the 'clusters' column to character
scType_result$clusters <- as.character(scType_result$clusters)
# Create a named vector for cluster mapping
cluster_mapping <- setNames(scType_result$type, scType_result$clusters)
# Convert Seurat object cluster metadata to character
DFU_all$seurat_clusters <- as.character(DFU_all$seurat_clusters)
# Assign new cluster labels to Seurat object metadata
DFU_all$ScType <- unname(cluster_mapping[as.character(DFU_all$seurat_clusters)])
end = Sys.time()
print(end-start)

# plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE,
        group.by = 'ScType', pt.size = 0.4, cols= colorCodes,
        label.size = 4)  



# Method 2: WAffinity
# create WAffinity matrix.
start = Sys.time()
WAffinity_Score = SC_WAffinity( count.matrix, norm.matrix, marker.genes )
end = Sys.time()
print(end-start)
WAffinity_tabel <- WAffinity_Score[["score"]]
#weight_list <- WAffinity_Score[["weight"]]
DFU_all$WAffinity_newGM <- WAffinity_Score$predict

DimPlot(DFU_all, reduction = 'umap', label = TRUE, repel = TRUE, 
        group.by = 'WAffinity_newGM', pt.size = 0.5,cols= colorCodes,
        label.size = 4)


# Method 3: sc2Assign
start = Sys.time()
sc2Assign_results <- SC_sc2Assign(WAffinity_tabel, DFU_all, percentile = 0.25, reduction = 'umap')
misclassified_cells <- sc2Assign_results$misclassified_cells
reassigned_types <- sc2Assign_results$reassigned_types
# Update the reassigned cell types in the Seurat object outside the function
DFU_all$sc2Assign_newGM <- DFU_all$WAffinity_newGM
for (cell_name in misclassified_cells) {
  DFU_all@meta.data[cell_name, "sc2Assign_newGM"] <- reassigned_types[cell_name]
}
end = Sys.time()
print(end-start)

DimPlot(DFU_all, reduction = 'umap', label = TRUE, repel = TRUE, 
        group.by = 'sc2Assign_newGM', pt.size = 0.5, cols= colorCodes,
        label.size = 4)


# scSorter
start = Sys.time()
scSorter_result = SC_scSorter(count.matrix, marker.genes)
end = Sys.time()
print(end-start)
DFU_all$scSorter_newGM  <- scSorter_result$Pred_Type
DFU_all$scSorter_newGM[DFU_all$scSorter_newGM == 'Unknown'] <- 'unknown'
#plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE, repel = TRUE, 
        group.by = 'scSorter_newGM', pt.size = 0.5, cols= colorCodes,
        label.size = 4)



# SCINA
start = Sys.time()
SCINA_result <- SC_SCINA(count.matrix, marker.genes)
end = Sys.time()
print(end-start)
DFU_all$SCINA <- SCINA_result$cell_labels
#plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE, repel = TRUE, 
        group.by = 'SCINA', pt.size = 0.5, cols= colorCodes,
        label.size = 4)

#
# scCATCH
start = Sys.time()
scCATCH_result <- SC_scCATCH(count.matrix, marker.genes, DFU_all)
end = Sys.time()
print(end-start)
DFU_all$scCATCH <- scCATCH_result
#plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE, repel = TRUE, 
        group.by = 'scCATCH', pt.size = 0.5, cols= colorCodes,
        label.size = 4)






# as we do not have the complete scale matrix in the seurat object Karla to us and the ScType woring with scaled data 
# I think bad results we got from this method is due to the fact that we used highly variable genes filter scalsed data
# wich include only about 1000 genes. So we are going to create another seurat to have the scaled data and run the 
# ScType again

# Create Seurat object  and remove low detection genes
DFU_all_new <- CreateSeuratObject(counts = count.matrix, min.cells = 1, min.features = 1)
DFU_all_new <- NormalizeData(DFU_all_new)
DFU_all_new <- ScaleData(DFU_all_new, features = rownames(DFU_all_new))
#DFU_all_new <- RunPCA(DFU_all_new, features = rownames(DFU_all_new))
# DFU_all_new <- FindNeighbors(DFU_all_new, dims = 1:50)  
# DFU_all_new <- FindClusters(DFU_all_new, resolution = 1)    # the clusters will be used for ScType and scCATCH
#start = Sys.time()
#DFU_all_new <- RunMCA(DFU_all_new)
#end = Sys.time()
#print(end-start)
# DFU_all_new <- RunTSNE(DFU_all_new, dims = 1:50, method = "tsne")
# DFU_all_new <- RunUMAP(DFU_all_new, dims = 1:50, method = "umap", min.dist = 1, n.neighbors = 20)


scaled.matrix <- as.data.frame(as.matrix(DFU_all_new@assays$RNA$counts)) 
# load the clustering infromation
DFU_all_new$seurat_clusters <- DFU_all$seurat_clusters


# cell typing using ScType
start = Sys.time()
scType_result <- SC_ScType( scaled.matrix, marker.genes, DFU_all_new )
# Convert the 'clusters' column to character
scType_result$clusters <- as.character(scType_result$clusters)
# Create a named vector for cluster mapping
cluster_mapping <- setNames(scType_result$type, scType_result$clusters)
# Convert Seurat object cluster metadata to character
DFU_all_new$seurat_clusters <- as.character(DFU_all_new$seurat_clusters)
# Assign new cluster labels to Seurat object metadata
DFU_all_new$ScType <- unname(cluster_mapping[as.character(DFU_all_new$seurat_clusters)])
end = Sys.time()
print(end-start)

# upload the clustering results we got from the whole scaled matrix to DFU_all as we do not run umap in this seurat
DFU_all$ScType_com <- DFU_all_new$ScType

# plot
DimPlot(DFU_all, reduction = 'umap', label = TRUE,
        group.by = 'ScType', pt.size = 0.4, cols= colorCodes,
        label.size = 4)  











# combine different clinical group plots
library(patchwork)  # For combining plots
# Assuming you have already created the UMAP coordinates and the meta data

# Assuming you have already created the UMAP coordinates and the meta data

# Extract UMAP coordinates
umap_coords <- Embeddings(DFU_all, "umap")

# Create a data frame with UMAP coordinates, cell types, and clinical groups
umap_data <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Cell_Type = DFU_all@meta.data$sc2Assign_newGM,
  Clinical_Group = DFU_all@meta.data$origin
)

# Define color codes for cell types
colorCodes <- c('#800000',  '#9A6324','#808000', '#aed6f1','#000075',  '#3cb44b','#89e300','#5ac7b9','#c75a7d', 
                '#007cba','red', '#e6194B', '#f58231',  '#ffe119','#42d4f4', '#4363d8', '#911eb4', '#f032e6', 'gray')  

names(colorCodes) <- c("Fibro", "SMC1", "SMC2", "VasEndo", "T-Lympho", "DiffKera", "BasalKera", "HE-Fibro", 
                       "B-Lympho", "M1-Macro", "NKT", "M2-Macro", "Schwann", "Melano", "Mast", 
                       "Sweat/Seba", "LymphEndo", "Plasma", 'unknown')


# Filter for specific clinical groups in desired order
umap_data <- umap_data %>%
  filter(Clinical_Group %in% c("Healthy", "Healing DFU", "Non-healing DFU", "Diabetic")) %>%
  mutate(Clinical_Group = factor(Clinical_Group, levels = c("Healthy", "Healing DFU", "Non-healing DFU", "Diabetic")))

# Create a UMAP plot for each clinical group
plots <- umap_data %>%
  group_by(Clinical_Group) %>%
  do(
    plot = ggplot(., aes(x = UMAP_1, y = UMAP_2, color = Cell_Type)) +
      geom_point(size = 0.5) +
      labs(title = unique(.$Clinical_Group)) +
      scale_color_manual(values = colorCodes) +  # Use custom color mapping
      theme_minimal() +
      theme(
        legend.position = "none",  # Remove legend by default
        plot.title = element_text(size = 14, hjust = 0.5)
      )
  )

# Add the legend to the last plot only with larger dots
plots$plot[[length(plots$plot)]] <- plots$plot[[length(plots$plot)]] + 
  theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase legend marker size

# Combine the plots into a single plot for comparison, with all plots in one row
combined_plot <- wrap_plots(plots$plot, ncol = 4)

# Display the combined plot
print(combined_plot)











# staked plot
create_stacked_bar_plot <- function(data, cell_typing_method, method_column, cell_type_column, origin_column, plot_title = "Stacked Bar Plot") {
  
  # Filter data for the selected method and exclude "unknown" labels
  method_data <- data[data[[method_column]] != "unknown", ]
  
  # Convert the 'origin' column to a factor with the desired order
  method_data[[origin_column]] <- factor(method_data[[origin_column]], 
                                         levels = c("Healthy", "Healing DFU", "Non-healing DFU", "Diabetic"))
  
  # Create a count table for the selected method and categorize by origin
  method_counts <- table(method_data[[cell_type_column]], method_data[[origin_column]])
  
  # Convert the table to a data frame
  method_df <- as.data.frame(method_counts)
  
  # Rename columns for better clarity
  colnames(method_df) <- c("Cell_Type", "Origin", "Count")
  
  # Specify the desired order of cell types (update this list as needed)
  cell_type_order <- c("SMC1", "VasEndo", "Fibro", "T-Lympho", "DiffKera", "M1-Macro", 
                       "HE-Fibro", "BasalKera", "M2-Macro", "Mast", "NKT", "B-Lympho", 
                       "LymphEndo", "Melano", "SMC2", "Sweat/Seba", "Schwann", "Plasma")
  
  # Convert the Cell_Type column to a factor with the specified order
  method_df$Cell_Type <- factor(method_df$Cell_Type, levels = cell_type_order)
  
  # Create the stacked bar plot
  ggplot(method_df, aes(x = Cell_Type, y = Count, fill = Origin)) +
    geom_bar(stat = "identity", position = "fill") +  
    scale_y_continuous(labels = scales::percent_format()) +  
    labs(x = "Cell Type", y = "Percentage", title = plot_title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14),  
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    scale_fill_manual(values = c("Healthy" = "#00FF00", "Healing DFU" = "#FFA500", 
                                 "Non-healing DFU" = "#FF0000", "Diabetic" = "#800080"))  # Set custom colors
}

#scSorter 
create_stacked_bar_plot(data = DFU_all@meta.data, 
                        cell_typing_method = "scSorter", 
                        method_column = "scSorter", 
                        cell_type_column = "scSorter", 
                        origin_column = "origin", 
                        plot_title = "scSorter Stacked Bar Plot")

#ScType 
create_stacked_bar_plot(data = DFU_all@meta.data, 
                        cell_typing_method = "ScType", 
                        method_column = "ScType", 
                        cell_type_column = "ScType", 
                        origin_column = "origin", 
                        plot_title = "ScType Stacked Bar Plot")

#sc2Assign 
create_stacked_bar_plot(data = DFU_all@meta.data, 
                        cell_typing_method = "sc2Assign ", 
                        method_column = "sc2Assign ", 
                        cell_type_column = "sc2Assign ", 
                        origin_column = "origin", 
                        plot_title = "sc2Assign Stacked Bar Plot")


# Filter data for sc2Assign method
sc2Assign_data <- DFU_all@meta.data[DFU_all@meta.data$sc2Assign_newGM != "unknown", ]
# Convert the 'origin' column to a factor with the desired order
sc2Assign_data$origin <- factor(sc2Assign_data$origin, levels = c("Healthy", "Healing DFU", "Non-healing DFU", "Diabetic"))
# Create a count table for cell types identified by sc2Assign and categorized by origin
sc2Assign_counts <- table(sc2Assign_data$sc2Assign_newGM, sc2Assign_data$origin)
# Convert the table to a data frame
sc2Assign_df <- as.data.frame(sc2Assign_counts)
# Rename columns for better clarity
colnames(sc2Assign_df) <- c("Cell_Type", "Origin", "Count")
# Specify the desired order of cell types
cell_type_order <- c("SMC1", "VasEndo",'Fibro', "T-Lympho", "DiffKera", "M1-Macro", "HE-Fibro", "BasalKera", "M2-Macro", "Mast", "NKT", "B-Lympho", "LymphEndo", "Melano", "SMC2", "Sweat/Seba", "Schwann", "Plasma")
# Convert the Cell_Type column to a factor with the specified order
sc2Assign_df$Cell_Type <- factor(sc2Assign_df$Cell_Type, levels = cell_type_order)
# Create the stacked bar plot
library(ggplot2)
ggplot(sc2Assign_df, aes(x = Cell_Type, y = Count, fill = Origin)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_y_continuous(labels = scales::percent_format()) +  
  labs(x = "Cell Type", y = "Percentage", title = "sc2Assign_newGM Stacked Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 14),  
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  scale_fill_manual(values = c("Healthy" = "#00FF00", "Healing DFU" = "#FFA500", "Non-healing DFU" = "#FF0000", "Diabetic" = "#800080"))  # Set custom colors





# 1) can you find the p value for cell composition between each paired group? 
#We will determine which cell type we need to further study. The cell type with most significant p value will be examined first.
#They use Macrophage and fibroblast in the original paper. 

# Extract meta.data for processing
df <- DFU_all@meta.data

# Rename columns for easier use in the function
df <- df %>%
  dplyr::rename(Cell_Type = sc2Assign_newGM, Clinical_Group = origin, Sample_ID = origin1)

# 1. Calculate cell composition per patient (Sample_ID) for each clinical group and cell type
df_group_composition <- df %>%
  group_by(Clinical_Group, Sample_ID, Cell_Type) %>%
  summarise(Cell_Count = n()) %>%
  ungroup()

# df_group_composition will give us the cell counts of each cell type for each patient we want to add another column
# to have he persentage of cell type in each clinical group
# Count the number of cells per patient (origin1) for each cell type 
cell_counts <- df %>%
  group_by(Sample_ID, Cell_Type) %>%
  summarise(Cell_Count = n()) %>%
  ungroup()

# Calculate the total number of cells per cell type 
total_cells_per_group <- df %>%
  group_by(Sample_ID) %>%
  summarise(Total_Cells = n())

# Merge the total number of cells per cell type with the cell counts per patient
df_group_composition <- df_group_composition %>%
  left_join(total_cells_per_group, by = "Sample_ID") %>%
  mutate(Percentage = (Cell_Count / Total_Cells) * 100)


# Load the ggpubr package
library(ggpubr)
# list of the cell types in order 
# for sc2Assign 
custom_order <- c( "B-Lympho","BasalKera", "DiffKera","Fibro","HE-Fibro","LymphEndo","M1-Macro","M2-Macro",  
                   "Mast","NKT","Plasma","SMC1","SMC2","T-Lympho","VasEndo","Melano" ,   
                   "Schwann","Sweat/Seba" )
# for scSorter
custom_order <- c( "B-Lympho", "BasalKera", "DiffKera", "Fibro", "HE-Fibro", "LymphEndo", 
                   "M1-Macro", "M2-Macro", "Mast", "NKT", "Plasma", "SMC1", "SMC2", 
                   "T-Lympho", "VasEndo", "Melano", "Schwann", "Sweat/Seba", "unknown")

# cell counts functions

# 2. Function to perform Welch's t-test between two clinical groups for a specific cell type using the cell counts
perform_welch_test <- function(df_group_composition, cell_type, group1, group2) {
  # Filter data for the cell type and the two clinical groups
  data_subset <- df_group_composition %>% filter(Cell_Type == cell_type & Clinical_Group %in% c(group1, group2))
  
  # Ensure both groups have enough observations
  if (length(unique(data_subset$Clinical_Group)) == 2 && all(table(data_subset$Clinical_Group) > 1)) {
    # Perform two-sided Welch's t-test on the cell composition (Cell_Count) for the two groups
    t_test_result <- t.test(Cell_Count ~ Clinical_Group, data = data_subset, var.equal = FALSE)
    return(t_test_result$p.value)
  } else {
    return(NA)  # Return NA if there aren't enough observations in both groups
  }
}


# Function to create a box plot for cell composition of different cell types  
create_plot <- function(cell_type, df_patient_composition, group_combination) {
  ggplot(df_patient_composition %>% filter(Cell_Type == cell_type), 
         aes(x = factor(Clinical_Group, levels = c("Healthy", "Diabetic", "Healing DFU", "Non-healing DFU")), 
             y = Cell_Count, fill = Clinical_Group)) +  # Set the order of x-axis categories
    geom_boxplot(outlier.shape = NA) +
    #stat_summary(fun = mean, geom = "point", shape = 8, size = 3, color = "black") +
    geom_jitter(width = 0.2) +
    labs(title = cell_type, x = NULL, y = "Cell Count") +
    scale_fill_manual(values = c("Healthy" = "#00FF00", "Healing DFU" = "#FFA500", 
                                 "Non-healing DFU" = "#FF0000", "Diabetic" = "#800080")) +
    scale_x_discrete(labels = c("Healthy" = "Healthy", "Diabetic" = "Diabetic", 
                                "Healing DFU" = "Healers", "Non-healing DFU" = "Non-healers")) +  # Set custom labels
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14), 
          axis.title = element_text(size = 12), 
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.8),  # Rotate x-axis text
          axis.text.y = element_text(size = 10)) +  # Adjust y-axis text size
    stat_compare_means(aes(group = Clinical_Group), 
                       comparisons = group_combination, 
                       method = "t.test", label = "p.signif", hide.ns = TRUE,  vjust = 0.5)
}

# in the followingcode I am going to create a box plot for all 6 possible combination of the clinical groups instead of doing them 
# separately

# Get all possible combinations of groups without repetition
group_combinations <- combn(clinical_groups, 2, simplify = FALSE)

# creating the matrix and boxplot for the Healing DFU with respect to other groups
cell_types <- unique(df_group_composition$Cell_Type)
# Dataframe to store the results of Welch's t-test
results <- data.frame(Cell_Type = character(), Group1 = character(), Group2 = character(), P_value = numeric())
# Loop through each cell type and group comparison to perform Welch's t-test
for (cell_type in cell_types) {
  for (comparison in group_combinations) {
    p_value <- perform_welch_test(df_group_composition, cell_type, comparison[1], comparison[2])
    results <- rbind(results, data.frame(Cell_Type = cell_type, Group1 = comparison[1], Group2 = comparison[2], P_value = p_value))
  }
}
# View results of the Welch's test
print(results)
#write.csv(results,file = 'sc2Assign_newGM_cellularEnrch_2sidedWelch.csv',row.names =TRUE )
# Reorder the Cell_Type column in df_group_composition based on custom order
df_group_composition$Cell_Type <- factor(df_group_composition$Cell_Type, levels = custom_order)
# Create a list to store all the plots
plot_list <- list()
# Loop through each cell type and generate the corresponding plot
for (cell_type in custom_order) {  # Loop through custom_order to ensure consistent ordering
  if (cell_type %in% df_group_composition$Cell_Type) {  # Only create plots for cell types present in the data
    plot_list[[cell_type]] <- create_plot(cell_type, df_group_composition, group_combinations)
  }
}
# Arrange all plots into a single figure
combined_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 4)
# Display the combined plot
print(combined_plot)



# percentage of cell types each patient compose 

# doing it at once for all possible combination of clinical groups
perform_welch_test <- function(df_group_composition, cell_type, group1, group2) {
  # Filter data for the cell type and the two clinical groups
  data_subset <- df_group_composition %>% filter(Cell_Type == cell_type & Clinical_Group %in% c(group1, group2))
  
  # Ensure both groups have enough observations
  if (length(unique(data_subset$Clinical_Group)) == 2 && all(table(data_subset$Clinical_Group) > 1)) {
    # Perform two-sided Welch's t-test on the cell type composition (Percentage) for the two groups
    t_test_result <- t.test(Percentage ~ Clinical_Group, data = data_subset, var.equal = FALSE)
    return(t_test_result$p.value)
  } else {
    return(NA)  # Return NA if there aren't enough observations in both groups
  }
}


# Function to create a box plot for cell type composition of patients 
create_plot <- function(cell_type, df_patient_composition, group_combination) {
  ggplot(df_patient_composition %>% filter(Cell_Type == cell_type), 
         aes(x = factor(Clinical_Group, levels = c("Healthy", "Diabetic", "Healing DFU", "Non-healing DFU")), 
             y = Percentage, fill = Clinical_Group)) +  # Set the order of x-axis categories
    geom_boxplot(outlier.shape = NA) +
    #stat_summary(fun = mean, geom = "point", shape = 8, size = 3, color = "black") +
    geom_jitter(width = 0.2) +
    labs(title = cell_type, x = NULL, y = "% cell type composition of patients ") +
    scale_fill_manual(values = c("Healthy" = "#00FF00", "Healing DFU" = "#FFA500", 
                                 "Non-healing DFU" = "#FF0000", "Diabetic" = "#800080")) +
    scale_x_discrete(labels = c("Healthy" = "Healthy", "Diabetic" = "Diabetic", 
                                "Healing DFU" = "Healers", "Non-healing DFU" = "Non-healers")) +  # Set custom labels
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14), 
          axis.title = element_text(size = 12), 
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.8),  # Rotate x-axis text
          axis.text.y = element_text(size = 10)) +  # Adjust y-axis text size
    stat_compare_means(aes(group = Clinical_Group), 
                       comparisons = group_combination, 
                       method = "t.test", label = "p.signif", hide.ns = TRUE,  vjust = 0.5)
}

# in the followingcode I am going to create a box plot for all 6 possible combination of the clinical groups instead of doing them 
# separately

# Get all possible combinations of groups without repetition
group_combinations <- combn(clinical_groups, 2, simplify = FALSE)

# creating the matrix and boxplot for the Healing DFU with respect to other groups
cell_types <- unique(df_group_composition$Cell_Type)
# Dataframe to store the results of Welch's t-test
results <- data.frame(Cell_Type = character(), Group1 = character(), Group2 = character(), P_value = numeric())

# Loop through each cell type and group comparison to perform Welch's t-test
for (cell_type in cell_types) {
  for (comparison in group_combinations) {
    p_value <- perform_welch_test(df_group_composition, cell_type, comparison[1], comparison[2])
    results <- rbind(results, data.frame(Cell_Type = cell_type, Group1 = comparison[1], Group2 = comparison[2], P_value = p_value))
  }
}
# View results of the Welch's test
print(results)
#write.csv(results,file = 'sc2Assign_newGM_celltype_compos_2sidedWelch.csv',row.names =TRUE )
# Reorder the Cell_Type column in df_group_composition based on custom order
df_group_composition$Cell_Type <- factor(df_group_composition$Cell_Type, levels = custom_order)
# Create a list to store all the plots
plot_list <- list()
# Loop through each cell type and generate the corresponding plot
for (cell_type in custom_order) {  # Loop through custom_order to ensure consistent ordering
  if (cell_type %in% df_group_composition$Cell_Type) {  # Only create plots for cell types present in the data
    plot_list[[cell_type]] <- create_plot(cell_type, df_group_composition, group_combinations)
  }
}
# Arrange all plots into a single figure
combined_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 4)
# Display the combined plot
print(combined_plot)






# create two extra boxplots for M1/M2 and M1/HE-Fibro

# Step 1: Calculate the Ratios
# Load tidyr if not already loaded
library(tidyr)
library(dplyr)

# Filter M1-Macro and M2-Macro cell types
df_m1_m2_he <- df_group_composition %>%
  filter(Cell_Type %in% c("M1-Macro", "M2-Macro", 'HE-Fibro'))

# Separate M1-Macro and M2-Macro into different dataframes
df_m1 <- df_m1_m2_he %>%
  filter(Cell_Type == "M1-Macro") %>%
  select(Sample_ID, Clinical_Group, M1_Macro_Count = Cell_Count)

df_m2 <- df_m1_m2_he %>%
  filter(Cell_Type == "M2-Macro") %>%
  select(Sample_ID, M2_Macro_Count = Cell_Count)

df_he <- df_m1_m2_he %>%
  filter(Cell_Type == "HE-Fibro") %>%
  select(Sample_ID, HE_Fibro_Count = Cell_Count)

# Merge the two dataframes by Sample_ID
df_m1_m2_merged <- df_m1 %>%
  left_join(df_m2, by = "Sample_ID")

df_m1_he_merged <- df_m1 %>%
  left_join(df_he, by = "Sample_ID")

# Calculate the ratio M1/M2
# Replace NA values with 0 in the entire dataframe
df_m1_m2_merged[is.na(df_m1_m2_merged)] <- 0
df_m1_m2_merged <- df_m1_m2_merged %>%
  mutate(Ratio_M1_M2 = M1_Macro_Count / (M2_Macro_Count + M1_Macro_Count))

df_m1_he_merged <- df_m1_he_merged %>%
  mutate(Ratio_M1_HE = M1_Macro_Count / HE_Fibro_Count)


# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Function to create a boxplot with t-tests
create_boxplot_with_ttest <- function(df, ratio_column, y_label, group_comparisons) {
  ggplot(df, aes(x = factor(Clinical_Group, levels = c("Healthy", "Diabetic", "Healing DFU", "Non-healing DFU")), 
                 y = !!sym(ratio_column), fill = Clinical_Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(title = paste("Boxplot for", ratio_column), y = y_label, x = NULL) +
    scale_fill_manual(values = c("Healthy" = "#00FF00", "Healing DFU" = "#FFA500", 
                                 "Non-healing DFU" = "#FF0000", "Diabetic" = "#800080")) +
    scale_x_discrete(labels = c("Healthy" = "Healthy", "Diabetic" = "Diabetic", 
                                "Healing DFU" = "Healers", "Non-healing DFU" = "Non-healers")) +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14), 
          axis.title = element_text(size = 12), 
          axis.text.x = element_text(size = 10, angle = 45, vjust = 0.8),  
          axis.text.y = element_text(size = 10)) +
    stat_compare_means(aes(group = Clinical_Group), comparisons = group_comparisons, 
                       method = "t.test", label = "p.signif", hide.ns = TRUE,  vjust = 0.5)
}


# Create the boxplots for each condition
boxplot_m1_m2_hl <- create_boxplot_with_ttest(df_m1_m2_merged, "Ratio_M1_M2", "%M1 of Total Macrophages", group_comparisons_HL)
boxplot_m1_he_hl <- create_boxplot_with_ttest(df_m1_he_merged, "Ratio_M1_HE", "Ratio of M1-Macro to HE-Fibro", group_comparisons_HL)

boxplot_m1_m2_nhl <- create_boxplot_with_ttest(df_m1_m2_merged, "Ratio_M1_M2", "%M1 of Total Macrophages", group_comparisons_NHL)
boxplot_m1_he_nhl <- create_boxplot_with_ttest(df_m1_he_merged, "Ratio_M1_HE", "Ratio of M1-Macro to HE-Fibro", group_comparisons_NHL)

boxplot_m1_m2_db <- create_boxplot_with_ttest(df_m1_m2_merged, "Ratio_M1_M2", "%M1 of Total Macrophages", group_comparisons_DB)
boxplot_m1_he_db <- create_boxplot_with_ttest(df_m1_he_merged, "Ratio_M1_HE", "Ratio of M1-Macro to HE-Fibro", group_comparisons_DB)

# Stack all the boxplots vertically into one plot
combined_plots <- ggarrange(boxplot_m1_m2_hl, boxplot_m1_he_hl,
                            boxplot_m1_m2_nhl, boxplot_m1_he_nhl,
                            boxplot_m1_m2_db, boxplot_m1_he_db,
                            ncol = 2, nrow = 3)

# Print the combined plot
print(combined_plots)

# Define the clinical groups
clinical_groups <- c("Healthy", "Diabetic", "Healing DFU", "Non-healing DFU")

# Get all possible combinations of groups without repetition
group_combinations <- combn(clinical_groups, 2, simplify = FALSE)

#M1 to Total Macro
# Function to perform Welch's t-test between two clinical groups for M1/total macrophage ratio
perform_welch_test <- function(df_group_composition, group1, group2) {
  # Filter data for the two clinical groups
  data_subset <- df_group_composition %>%
    filter(Clinical_Group %in% c(group1, group2))
  
  # Ensure both groups have enough observations
  if (length(unique(data_subset$Clinical_Group)) == 2 && all(table(data_subset$Clinical_Group) > 1)) {
    # Perform two-sided Welch's t-test on the ratio
    t_test_result <- t.test(Ratio_M1_M2 ~ Clinical_Group, data = data_subset, var.equal = FALSE)
    return(t_test_result$p.value)
  } else {
    return(NA)  # Return NA if there aren't enough observations in both groups
  }
}


#  Create a dataframe to store the Welch's t-test results
results <- data.frame(Group1 = character(), Group2 = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each pair of groups and perform the Welch's t-test
for (comparison in group_combinations) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  # Perform Welch's t-test and store the p-value
  p_value <- perform_welch_test(df_m1_m2_merged, group1, group2)
  
  # Append the result to the dataframe
  results <- rbind(results, data.frame(Group1 = group1, Group2 = group2, P_value = p_value))
}

# View the results matrix
print(results)

#write.csv(results,file = 'sc2Assign_M1toTotalMacro_2sidedWelch_DB.csv',row.names =TRUE )


# M1 to HE-Fibro
perform_welch_test <- function(df_group_composition, group1, group2) {
  # Filter data for the two clinical groups
  data_subset <- df_group_composition %>%
    filter(Clinical_Group %in% c(group1, group2))
  
  # Ensure both groups have enough observations
  if (length(unique(data_subset$Clinical_Group)) == 2 && all(table(data_subset$Clinical_Group) > 1)) {
    # Perform two-sided Welch's t-test on the ratio
    t_test_result <- t.test(Ratio_M1_HE ~ Clinical_Group, data = data_subset, var.equal = FALSE)
    return(t_test_result$p.value)
  } else {
    return(NA)  # Return NA if there aren't enough observations in both groups
  }
}

#  Create a dataframe to store the Welch's t-test results
results <- data.frame(Group1 = character(), Group2 = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each pair of groups and perform the Welch's t-test
for (comparison in group_combinations) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  # Perform Welch's t-test and store the p-value
  p_value <- perform_welch_test(df_m1_he_merged, group1, group2)
  
  # Append the result to the dataframe
  results <- rbind(results, data.frame(Group1 = group1, Group2 = group2, P_value = p_value))
}

# View the results matrix
print(results)

#write.csv(results,file = 'sc2Assign_M1toHEFibro_2sidedWelch_DB.csv',row.names =TRUE )




# quality control of the cell typing (Dr. Jin suggest to create a matrix of all patient on the row and the cell
# types on the column and put the percentage of the cells recognized as each cell type is related to each 
# patient and do clustering based on this matrix)
library(tidyr)
unique(DFU_all$origin1)

# Extract metadata from Seurat object (patient information and cell types)
df_meta <- DFU_all@meta.data

# Choose the cell typing method (sc2Assign or scSorter)
df_meta <- df_meta %>%
  dplyr::select(origin1, sc2Assign_newGM)  # Select patient and cell type columns

# 1. Count the number of cells per patient (origin1) for each cell type (sc2Assign or scSorter)
cell_counts <- df_meta %>%
  group_by(origin1, sc2Assign_newGM) %>%
  summarise(Cell_Count = n()) %>%
  ungroup()

# 2. Calculate the total number of cells per cell type (sc2Assign or scSorter)
total_cells_per_celltype <- df_meta %>%
  group_by(sc2Assign_newGM) %>%
  summarise(Total_Cells = n())

# 3. Merge the total number of cells per cell type with the cell counts per patient
cell_percentage <- cell_counts %>%
  left_join(total_cells_per_celltype, by = "sc2Assign_newGM") %>%
  mutate(Percentage = (Cell_Count / Total_Cells) * 100)

# 4. Pivot the data to create a matrix with patients (origin1) as rows and cell types (sc2Assign or scSorter) as columns
percentage_matrix <- cell_percentage %>%
  pivot_wider(names_from = sc2Assign_newGM, values_from = Percentage, values_fill = list(Percentage = 0))

# 5. Aggregate data so that each patient (origin1) appears only once with all cell types summarized
percentage_matrix <- percentage_matrix %>%
  group_by(origin1) %>%
  summarise(across(everything(), sum))
sum(percentage_matrix$Melano)

percentage_matrix <- percentage_matrix %>% select(-Cell_Count)
percentage_matrix <- percentage_matrix %>% select(-Total_Cells)
# 5. Set rownames as the patient IDs (origin1)
Rownames <- percentage_matrix$origin1
percentage_matrix <- percentage_matrix %>% select(-origin1)
rownames(percentage_matrix) <- Rownames

# 6. View the matrix (patients x cell types matrix with percentages)
print(percentage_matrix)

# heat map
# Install pheatmap if not already installed
#install.packages("pheatmap")

# Load pheatmap
library(pheatmap)

# Plot heatmap with pheatmap
pheatmap(t(percentage_matrix), 
         cluster_rows = TRUE,   # Cluster rows (optional)
         cluster_cols = TRUE,   # Cluster columns (optional)
         display_numbers = TRUE, # Optional: to display the actual values in the cells
         fontsize_row = 10,     # Adjust the font size for row labels
         fontsize_col = 10)     # Adjust the font size for column labels




# performing the percentage matrix but not for cell types this time for each groups(means that I want to see how much % of each 
# clinical group are related to each cell type)
library(tidyr)
unique(DFU_all$origin1)

# Extract metadata from Seurat object (patient information and cell types)
df_meta <- DFU_all@meta.data

# Choose the cell typing method (sc2Assign or scSorter)
df_meta <- df_meta %>%
  dplyr::select(origin1, scSorter)  # Select patient and cell type columns

# 1. Count the number of cells per patient (origin1) for each cell type (sc2Assign or scSorter)
cell_counts <- df_meta %>%
  group_by(origin1, scSorter) %>%
  summarise(Cell_Count = n()) %>%
  ungroup()

# 2. Calculate the total number of cells per cell type (sc2Assign or scSorter)
total_cells_per_group <- df_meta %>%
  group_by(origin1) %>%
  summarise(Total_Cells = n())

# 3. Merge the total number of cells per cell type with the cell counts per patient
cell_percentage_for_groups <- cell_counts %>%
  left_join(total_cells_per_group, by = "origin1") %>%
  mutate(Percentage = (Cell_Count / Total_Cells) * 100)

# 4. Pivot the data to create a matrix with patients (origin1) as rows and cell types (sc2Assign or scSorter) as columns
percentage_matrix <- cell_percentage %>%
  pivot_wider(names_from = origin1, values_from = Percentage, values_fill = list(Percentage = 0))

# 5. Aggregate data so that each patient (origin1) appears only once with all cell types summarized
percentage_matrix <- percentage_matrix %>%
  group_by(sc2Assign_HVG) %>%
  summarise(across(everything(), sum))
sum(percentage_matrix$Melano)

percentage_matrix <- percentage_matrix %>% select(-Cell_Count)
percentage_matrix <- percentage_matrix %>% select(-Total_Cells)
# 5. Set rownames as the patient IDs (origin1)
Rownames <- percentage_matrix$sc2Assign_HVG
percentage_matrix <- percentage_matrix %>% select(-sc2Assign_HVG)
rownames(percentage_matrix) <- Rownames

# 6. View the matrix (patients x cell types matrix with percentages)
print(percentage_matrix)

# heat map
# Install pheatmap if not already installed
#install.packages("pheatmap")

# Load pheatmap
library(pheatmap)

# Plot heatmap with pheatmap
pheatmap(t(percentage_matrix), 
         cluster_rows = TRUE,   # Cluster rows (optional)
         cluster_cols = TRUE,   # Cluster columns (optional)
         display_numbers = TRUE, # Optional: to display the actual values in the cells
         fontsize_row = 10,     # Adjust the font size for row labels
         fontsize_col = 10)     # Adjust the font size for column labels


