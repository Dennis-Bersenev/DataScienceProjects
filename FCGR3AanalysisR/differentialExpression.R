library(dplyr)
library(Seurat)
library(patchwork)

# Going through the Seurat workflow

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Welch's t-test for differential FCGR3A expression in different cell types ---------------------------------------------------------------------------
gene = pbmc@assays$RNA@var.features
FCGR3AMono_cells <- subset(pbmc, idents = "FCGR3A+ Mono")
nk_cells <- subset(pbmc, idents = "NK")

x_all = pbmc@assays $RNA@data[gene,] 
x_mono = FCGR3AMono_cells@assays$RNA@data[gene,]  
x_nk <- nk_cells@assays$RNA@data[gene,]

x_all = as.matrix(x_all)
x_all = t(x_all)

x_mono = as.matrix(x_mono)
x_mono = t(x_mono)

x_nk = as.matrix(x_nk)
x_nk = t(x_nk)


FCGR3A_all = unname(x_all[,"FCGR3A"])
FCGR3A_mono = unname(x_mono[,"FCGR3A"])
FCGR3A_nk = unname(x_nk[,"FCGR3A"])

t1 = t.test(FCGR3A_mono, FCGR3A_nk)
# t1

# Violin plot of FCGR3A across cell types
VlnPlot(pbmc, features = c("FCGR3A"))

# Entropy of raw counts of FCGR3A in FCGR3A monocytes and across all cells ---------------------------------------------------------------------------
#use the empirical distributions for p(x) from count data
#eg. the number of times the same umi appears/the total counts = p(umi)
#umis = genes
#across all cells you get a total count for each umi and all umis

raw_counts<-pbmc[["RNA"]]@counts
rx<-as.matrix(raw_counts)
rx<-t(rx)
#Column containing the raw counts of FCGR3A across all cells
FCGR3A_col_raw = unname(rx[,"FCGR3A"]) 
p_rx = table(unlist(FCGR3A_col_raw)) / nrow(rx) #the distribution across all 

raw_counts_mono <- FCGR3AMono_cells[["RNA"]]@counts
rx_mono<-as.matrix(raw_counts_mono)
rx_mono<-t(rx_mono)
#Column containing the raw counts of FCGR3A across fcgrea monocytes
FCGR3A_col_raw_mono = unname(rx_mono[,"FCGR3A"]) 
p_rx_mono = table(unlist(FCGR3A_col_raw_mono)) / nrow(rx_mono) #the distribution across mono

#Funct to calculate entropies given the distribution
entropy <- function(distro) {
  vec <- as.data.frame(distro)[,2]
  vec<-vec[vec>0]
  -sum(vec * log2(vec))
}

res_all<-entropy(p_rx)
res_mono<-entropy(p_rx_mono)
