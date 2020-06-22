library(Seurat)
library(data.table)
library(ggplot2)
library(magrittr)
library(seqsetvis)
library(openxlsx)

source("../SF_AutoImmune_ssv/functions_setup.R")
source("../SF_AutoImmune_ssv/functions_hres.R")

dksc = readRDS("datasets/DKSC.integrated.Rds")
to_comb = list(body = c(2,4,8,3,9), 
               fin = 7, 
               tail = c(0,5,6))
meta_dt = get_meta_dt(dksc, to_combine = to_comb)
ggplot(meta_dt, aes(x = UMAP_1, y = UMAP_2, color = meta_cluster)) +
    geom_point()
toadd = meta_dt$meta_cluster
names(toadd)= meta_dt$id
dksc = AddMetaData(dksc,toadd, col.name = "meta_cluster")
dksc@meta.data
?FindMarkers
mark_res = lapply(names(to_comb), function(g){
    FindMarkers(dksc, group.by = "meta_cluster", ident.1 = g)    
})
names(mark_res) = names(to_comb)
mark_dt = lapply(mark_res, as.data.table, keep.rownames = T) %>% rbindlist(., idcol = "region")
setnames(mark_dt, "rn", "gene_name")
FindMarkers(dksc, group.by = "meta_cluster", ident.1 = "fin")
FindMarkers(dksc, group.by = "meta_cluster", ident.1 = "tail")
