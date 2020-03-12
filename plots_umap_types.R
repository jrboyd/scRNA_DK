library(Seurat)
source("../SF_AutoImmune_ssv/functions_setup.R")
dksc = readRDS("datasets/DKSC.combined.Rds")
meta_dt = get_meta_dt(dksc)
dksc.integrated = readRDS("datasets/DKSC.integrated.Rds")
meta_dt.integrated = get_meta_dt(dksc.integrated)

mt_genes = rownames(dksc)[grepl("mt-", rownames(dksc))]


## combined umaps
p1 <- DimPlot(dksc, reduction = "umap", group.by = "sampleId")
p2 <- DimPlot(dksc, reduction = "umap", group.by = "treatment", 
              repel = TRUE)
p3 <- DimPlot(dksc, reduction = "umap", group.by = "rep",
              repel = TRUE)
p4 = DimPlot(dksc, reduction = "umap", group.by = "Phase",
             repel = TRUE)
p5 = DimPlot(dksc, reduction = "umap", group.by = "seurat_clusters",
             repel = TRUE) +
    labs(color = "Cluster")

p6 = ggplot(meta_dt, aes(x = UMAP_1, y = UMAP_2, color = log10(nCount_RNA))) + 
    geom_point(size = .4) +
    scale_color_viridis_c() + theme(legend.title = element_text(size = 10, angle = 90)) + 
    guides(color = guide_colorbar(title.position = "left"))

mt_rna_dt = get_rna_dt(dksc, mt_genes)[, .(mt_average = mean(expression)), .(id)]
mt_rna_dt = merge(mt_rna_dt, meta_dt[, .(id, source, sampleId, treatment, rep, UMAP_1, UMAP_2, seurat_clusters)], by = "id")
p7 = ggplot(mt_rna_dt, aes(x = UMAP_1, y = UMAP_2, color = mt_average)) + 
    geom_point(size = .4) +
    scale_color_viridis_c(option = "B")+
    labs(color = "mitochondrial average") +
    theme(legend.title = element_text(size = 10, angle = 90)) + 
    guides(color = guide_colorbar(title.position = "left"))

pg = cowplot::plot_grid(p1 + labs(title = "Unanchored", color = "sampleId"), 
                        p2 + labs(color = "treatment"), 
                        p3 + labs(color = "rep"),
                        p4 + labs(color = "Cell cycle"),
                        p5,
                        p6, 
                        p7,
                        ncol = 4)


pg
ggsave("combined_umap.pdf", pg, width = 7*2, height = 6.1)

## integrated umaps

p1 <- DimPlot(dksc.integrated, reduction = "umap", group.by = "sampleId")
p2 <- DimPlot(dksc.integrated, reduction = "umap", group.by = "treatment", 
              repel = TRUE)
p3 <- DimPlot(dksc.integrated, reduction = "umap", group.by = "rep",
              repel = TRUE)
p4 = DimPlot(dksc.integrated, reduction = "umap", group.by = "Phase",
             repel = TRUE)
p5 = DimPlot(dksc.integrated, reduction = "umap", group.by = "seurat_clusters",
             repel = TRUE) +
    labs(color = "Cluster")

p6 = ggplot(meta_dt.integrated, aes(x = UMAP_1, y = UMAP_2, color = log10(nCount_RNA))) + 
    geom_point(size = .4) +
    scale_color_viridis_c() + theme(legend.title = element_text(size = 10, angle = 90)) + 
    guides(color = guide_colorbar(title.position = "left"))

mt_rna_dt.integrated = get_rna_dt(dksc.integrated, mt_genes)[, .(mt_average = mean(expression)), .(id)]
mt_rna_dt.integrated = merge(mt_rna_dt.integrated, meta_dt.integrated[, .(id, source, sampleId, treatment, rep, UMAP_1, UMAP_2, seurat_clusters)], by = "id")
p7 = ggplot(mt_rna_dt.integrated, aes(x = UMAP_1, y = UMAP_2, color = mt_average)) + 
    geom_point(size = .4) +
    scale_color_viridis_c(option = "B") +
    labs(color = "mitochondrial average") +
    theme(legend.title = element_text(size = 10, angle = 90)) + 
    guides(color = guide_colorbar(title.position = "left"))

pg = cowplot::plot_grid(p1 + labs(title = "Anchored", color = "sampleId"), 
                        p2 + labs(color = "treatment"), 
                        p3 + labs(color = "rep"),
                        p4 + labs(color = "Cell cycle"),
                        p5,
                        p6, 
                        p7,
                        ncol = 4)
gc()
pg
ggsave("anchored_umap.pdf", pg, width = 7*2, height = 6.1)

## composition
p1 = DimPlot(dksc, reduction = "umap", group.by = "seurat_clusters",
             repel = TRUE, label = TRUE) + NoLegend() +
    labs(title = "Unanchored")

p2 = DimPlot(dksc.integrated, reduction = "umap", group.by = "seurat_clusters",
             repel = TRUE, label = TRUE) + NoLegend() +
    labs(title = "Anchored")

cnt_dt = meta_dt[, .N,.(sampleId, treatment, seurat_clusters)]
cnt_dt[, fraction := N / sum(N), .(sampleId)]

total_dt = cnt_dt[, .(total = sum(N)), .(seurat_clusters)]
setkey(total_dt, seurat_clusters)
lev = levels(total_dt$seurat_clusters)
levels(cnt_dt$seurat_clusters) = paste0(lev, " (", total_dt[.(lev)]$total, ")")

p1a = ggplot(cnt_dt, aes(x = sampleId, y = fraction, fill = treatment)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~seurat_clusters) + 
    theme(legend.position = "bottom")

p1b = ggplot(cnt_dt, aes(x = sampleId, y = N, fill = treatment)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~seurat_clusters) + 
    theme(legend.position = "bottom")




cnt_dt.integrated = meta_dt.integrated[, .N,.(sampleId, treatment, seurat_clusters)]
cnt_dt.integrated[, fraction := N / sum(N), .(sampleId)]

total_dt = cnt_dt.integrated[, .(total = sum(N)), .(seurat_clusters)]
setkey(total_dt, seurat_clusters)
lev = levels(total_dt$seurat_clusters)
levels(cnt_dt.integrated$seurat_clusters) = paste0(lev, " (", total_dt[.(lev)]$total, ")")

p2a = ggplot(cnt_dt.integrated, aes(x = sampleId, y = fraction, fill = treatment)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~seurat_clusters) +
    theme(legend.position =  "bottom")

p2b = ggplot(cnt_dt.integrated, aes(x = sampleId, y = N, fill = treatment)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~seurat_clusters) + 
    theme(legend.position = "bottom")


pg = cowplot::plot_grid(p1, p2, p1b, p2b, p1a, p2a, rel_heights = c(1, 1.3, 1.3), ncol = 2)
ggsave("cluster_composition.pdf", pg, width = 5.9, height = 10)
