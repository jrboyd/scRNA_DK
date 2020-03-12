library(Seurat)
library(magrittr)
source("../SF_AutoImmune_ssv/functions_setup.R")

data_dirs = dir("~/Dimitry_scRNA/Krementsov_10Xgenomics_iLabs_14530/", pattern = "^DK", full.names = TRUE)
names(data_dirs) = basename(data_dirs)

#load data
data_loaded = lapply(data_dirs, function(d){
    message(d)
    Seurat::Read10X(file.path(d, "outs/filtered_feature_bc_matrix"))
})
seurat_objs = lapply(data_loaded, function(d){
    CreateSeuratObject(d)
})             

#combine samples and normalize
dksc = merge(seurat_objs[[1]], seurat_objs[-1], add.cell.ids = names(seurat_objs))
ids = rownames(dksc@meta.data)
meta = strsplit(ids, "_") %>% unlist %>% matrix(., ncol = 5, byrow = TRUE)
colnames(meta) = c("source", "sampleId", "treatment", "rep", "barcode")
for(i in seq_len(ncol(meta))){
    dksc = AddMetaData(dksc, metadata = meta[,i], col.name = colnames(meta)[i])
}

dksc = NormalizeData(dksc, verbose = FALSE)
dksc = FindVariableFeatures(dksc, selection.method = "vst", 
                           nfeatures = 2000, verbose = FALSE)
dksc <- ScaleData(dksc, verbose = FALSE)
dksc <- RunPCA(dksc, npcs = 30, verbose = FALSE)
dksc <- RunUMAP(dksc, reduction = "pca", dims = 1:30)

#cell cycle
mm_cc.genes = lapply(Seurat::cc.genes.updated.2019, gene_name_gs2mm)
s.genes <- mm_cc.genes$s.genes
g2m.genes <- mm_cc.genes$g2m.genes
dksc <- CellCycleScoring(dksc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


dksc.split = SplitObject(dksc, "sampleId")

#clustering
dksc <- FindNeighbors(dksc, dims = 1:10)
dksc <- FindClusters(dksc, resolution = 0.5)



sapply(dksc.split,ncol)

#normalize and preprocess
dksc.norm = pbmcapply::pbmclapply(dksc.split, mc.cores = 10, function(obj){
    obj = NormalizeData(obj, verbose = FALSE)
    obj = FindVariableFeatures(obj, selection.method = "vst", 
                         nfeatures = 2000, verbose = FALSE)
    obj
})

dksc.anchor <- FindIntegrationAnchors(object.list = dksc.split, dims = 1:30, k.filter = 100)
dksc.integrated = IntegrateData(dksc.anchor)
DefaultAssay(dksc.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
dksc.integrated <- ScaleData(dksc.integrated, verbose = FALSE)
dksc.integrated <- RunPCA(dksc.integrated, npcs = 30, verbose = FALSE)
dksc.integrated <- RunUMAP(dksc.integrated, reduction = "pca", dims = 1:30)

#cell cycle
mm_cc.genes = lapply(Seurat::cc.genes.updated.2019, gene_name_gs2mm)
s.genes <- mm_cc.genes$s.genes
g2m.genes <- mm_cc.genes$g2m.genes
dksc.integrated <- CellCycleScoring(dksc.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#clustering
# DefaultAssay(dksc.integrated) = "RNA"
dksc.integrated <- FindNeighbors(dksc.integrated, dims = 1:10)
dksc.integrated <- FindClusters(dksc.integrated, resolution = 0.5)


DefaultAssay(dksc.integrated) = "RNA"

ts = format(Sys.time(), "%m%d%y")

saveRDS(dksc, file = paste0("datasets/DKSC.combined.", ts, ".Rds"))
saveRDS(dksc.integrated, file = paste0("datasets/DKSC.integrated.", ts, ".Rds"))
