library(Seurat)
library(ggplot2)
library(cowplot)
library(knitr)
library(dplyr)


#input E8 skin
setwd("/home/miller/Desktop/ZF_scAnalysis/chicken_Seurat_analysis/MillerCK_pooledRNAE8_ens106")
getwd()

filt.matrix <- Read10X_h5("./outs/filtered_feature_bc_matrix.h5",use.names = T)
E8  <- CreateSeuratObject(counts = filt.matrix, project = "E8",)
E8$culture <- "E8"

#input E65 skin
setwd("/home/miller/Desktop/ZF_scAnalysis/chicken_Seurat_analysis/MillerCK_pooledRNAE65_ens106")
getwd()

filt.matrix <- Read10X_h5("./outs/filtered_feature_bc_matrix.h5",use.names = T)
E65  <- CreateSeuratObject(counts = filt.matrix, project = "E8",)
E65$culture <- "E65"


#merge
chicken.merged <- merge(E8, y = E65, add.cell.ids = c("E8","E65"), project = "Chicken reconstitution")
chicken.merged

table(chicken.merged$orig.ident)


#grep all the MT genes
#Calculate percentage of mitochondrial on ribosomal counts
# Percent of mitochondrial counts
grep(pattern ='(^ND1$|^MT-|^ND3$|^ND4$|^ND4L$|^ND5$|^ND6$|^CYTB$|^COX3$|^ATP6$|^ATP8$)',rownames(chicken.merged@assays$RNA@counts),value = TRUE)

chicken.merged[["percent.mt"]] <- PercentageFeatureSet(chicken.merged, pattern = '(^ND1$|^MT-|^ND3$|^ND4$|^ND4L$|^ND5$|^ND6$|^CYTB$|^COX3$|^ATP6$|^ATP8$)')
str(chicken.merged@meta.data)

# Percent of ribosomal protein
grep(pattern ='(^RPL|^RPS|^MRP)',rownames(chicken.merged@assays$RNA@counts),value = TRUE)

chicken.merged[["percent.ribo"]] <- PercentageFeatureSet(chicken.merged, pattern = '(^RPL|^RPS|^MRP)')
str(chicken.merged@meta.data)


# Visualize QC metrics as a violin plot
VlnPlot(chicken.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2, pt.size = 0)

# Get cell IDs for cells with percent.mt > 25
chicken.merged_ids <- rownames(chicken.merged@meta.data[which(chicken.merged@meta.data$percent.mt >25), ])

head(chicken.merged_ids)


## filter them out
length(colnames(chicken.merged))
chicken.merged <- chicken.merged[,!colnames(chicken.merged) %in% chicken.merged_ids]
length(colnames(chicken.merged))

#Filtering cells based on number of genes and transcripts detected
#Remove cells with too few gene detected or with too many UMI counts
#Set low and hight thresholds on the number of detected genes

RNA.max <- round(mean(chicken.merged$nFeature_RNA) + 2 * sd(chicken.merged$nFeature_RNA), digits = -2)
RNA.min <- round(mean(chicken.merged$nFeature_RNA) - 2 * sd(chicken.merged$nFeature_RNA), digits = -2)

# Set minimum parameters to 0 if negative value
if (RNA.min < 0){
  RNA.min <- 0
} else {
  RNA.min <- RNA.min
}

# Set hight threshold on the number of transcripts
Cell.QC.Stat <- chicken.merged@meta.data
max.nCount_RNA.thr <- median(Cell.QC.Stat$nCount_RNA) + 3*mad(Cell.QC.Stat$nCount_RNA)

# Filter cells base on both metrics
chicken.merged_subset <- subset(chicken.merged, subset = nFeature_RNA < RNA.max & 
                                  nFeature_RNA > RNA.min &
                                  nCount_RNA < max.nCount_RNA.thr)


VlnPlot(chicken.merged_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2, pt.size = 0)


saveRDS(chicken.merged_subset, file = "Gary_merged_E8_E65_filtered_seuratobj")

chicken.merged
chicken.merged_subset
rm(chicken.merged)
rm(E65)
rm(E8)


###Use Scrublet to detect obvious doublets
#https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html#Use_Scrublet_to_detect_obvious_doublets

library(RColorBrewer)
library(reticulate)

#Export filtered matrix to input matrix format for scrublet
library(DropletUtils)

write10xCounts(x = chicken.merged_subset@assays$RNA@counts, path = '/home/miller/Desktop/Gary_E8_8h_43h_CK_integration/soupX_seuratQC_mergedGary')

# Check the current Python version
# !!!scrublet has to be installed and used in the same Python version!!!
reticulate::py_config()
# mypath: /home/miller/.local/share/r-miniconda/envs/r-reticulate/bin

# Loads Python Shell
repl_python()

#(Python) Export raw count matrix as input to Scrublet
import scrublet as scr
import scipy.io
import numpy as np
import os

#(Python) Load raw counts matrix and gene list
input_dir = '/home/miller/Desktop/Gary_E8_8h_43h_CK_integration'
counts_matrix = scipy.io.mmread(input_dir + '/soupX_seuratQC_mergedGary/matrix.mtx').T.tocsc()

#(Python) Initialize Scrublet object
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1)

#(Python) Run the default pipeline
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

# Import scrublet's doublet score
chicken.merged_subset@meta.data$Doubletscore <- py$doublet_scores

# Plot doublet score
ggplot(chicken.merged_subset@meta.data, aes(x = Doubletscore, stat(ndensity))) +
  geom_histogram(bins = 200, colour ="lightgrey")+
  geom_vline(xintercept = 0.2, colour = "red", linetype = 2) # Manually set threshold

# Manually set threshold at doublet score to 0.2
chicken.merged_subset@meta.data$Predicted_doublets <- ifelse(py$doublet_scores > 0.2, "Doublet","Singlet" )
table(chicken.merged_subset@meta.data$Predicted_doublets)
head(chicken.merged_subset@meta.data)

#remove doublet
chicken.merged_subset_scrublet <- subset(chicken.merged_subset, subset = Predicted_doublets=='Singlet')
head(chicken.merged_subset_scrublet@meta.data)

saveRDS(chicken.merged_subset_scrublet, file = "Gary_merged_E8_E65_filtered_scrublet_seuratobj.RDS")

rm(chicken.merged_subset)
rm(Cell.QC.Stat)
rm(filt.matrix)
rm(E65)

#Perform integration with SCTransform-normalized datasets
chicken.merged_subset_scrublet
table(chicken.merged_subset_scrublet@meta.data$culture)

# split the dataset into a list of two seurat objects (stim and CTRL)
sct.list <- SplitObject(chicken.merged_subset_scrublet, split.by = "culture")


E65 <- sct.list[["E65"]]
E65

E8 <- sct.list[["E8"]]
E8


E65 <- SCTransform(E65, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)

E8 <- SCTransform(E8, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)


rm(chicken.merged_subset_scrublet)
rm(sct.list)


ifnb.list <- list(E65=E65, E8 = E8)

rm(E65)
rm(E8)

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
head(features)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                  anchor.features = features)
integrated.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

setwd("/home/miller/Desktop/Gary_E8_8h_43h_CK_integration")

saveRDS(integrated.sct, file = "integrated_E65.E8_filtered_SCT_noPCAyet_seuratobj.RDS")
integrated.sct <- readRDS("integrated_E65.E8_filtered_SCT_noPCAyet_seuratobj.RDS")

integrated.sct <- RunPCA(integrated.sct, verbose = FALSE) 

ElbowPlot(integrated.sct)


integrated.sct <- RunUMAP(integrated.sct, reduction = "pca", dims = 1:12, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:12)

integrated.sct.cluster <- FindClusters(integrated.sct,resolution = 0.3)

DimPlot(integrated.sct.cluster, reduction = "umap", split.by = "culture")
DimPlot(integrated.sct.cluster, reduction = "umap", group.by = c("culture", "seurat_annotations"))

rm(integrated.sct)

saveRDS(integrated.sct.cluster, file = "Gary_SCTintegrated_E65.E8_res0.3_seuratobj.RDS")
integrated.sct.cluster <- readRDS("Gary_SCTintegrated_E65.E8_res0.3_seuratobj.RDS")



###Identify differential expressed genes across conditions
#First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column.
integrated.sct.cluster
head(integrated.sct.cluster@meta.data)

integrated.sct.cluster$cluster.culture <- paste(integrated.sct.cluster$seurat_clusters, integrated.sct.cluster$culture,
                                                sep = "_")
Idents(integrated.sct.cluster) <- "cluster.culture"

integrated.sct.cluster <- PrepSCTFindMarkers(integrated.sct.cluster)

integrated.sct.cluster

saveRDS(integrated.sct.cluster, file = "Gary_SCTintegratedFindmarkers_E65.E8_res0.3_seuratobj.RDS")

DimPlot(integrated.sct.cluster, 
        reduction = 'umap',split.by = "culture")

head(integrated.sct.cluster@meta.data)
DefaultAssay(integrated.sct.cluster) <- "SCT"
Idents(integrated.sct.cluster) <- "seurat_clusters"

#export seurat obj to loupe
remotes::install_github("10xGenomics/loupeR")
loupeR::setup()
library(loupeR)

create_loupe_from_seurat(
  integrated.sct.cluster,
  output_dir = "/home/miller/Desktop/Gary_E8_8h_43h_CK_integration/",
  output_name = "Gary_SCTintegrate_E65.E8_chicken_r0.3",
  dedup_clusters = FALSE,
  executable_path = NULL,
  force = T
)
