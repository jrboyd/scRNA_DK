library(Seurat)
library(data.table)
library(ggplot2)
library(magrittr)
library(seqsetvis)
library(openxlsx)

source("../SF_AutoImmune_ssv/functions_setup.R")
source("../SF_AutoImmune_ssv/functions_hres.R")

dksc = readRDS("datasets/DKSC.integrated.Rds")
gls = readRDS("gene_references/AmitColonna_gene_lists.Rds")

mi_ref = "gene_references/microarray DEG list for Joe.xlsx"
sn = openxlsx::getSheetNames(mi_ref)
mi_dat = lapply(sn, function(x){
    as.data.table(read.xlsx(mi_ref, x)    )
})
names(mi_dat) = sn


mi_gls_up = lapply(mi_dat, function(x){
    gl = x[`P-val` <= .05 & Fold.Change > 0]$Gene.Symbol
    strsplit(gl, "; ") %>% unlist
})
mi_gls_down = lapply(mi_dat, function(x){
    gl = x[`P-val` <= .05 & Fold.Change < 0]$Gene.Symbol
    strsplit(gl, "; ") %>% unlist
})
names(mi_gls_up) = paste(names(mi_gls_up), "- up in KO")
names(mi_gls_down) = paste(names(mi_gls_down), "- down in KO")

to_score = list(
    up_tight = gls$up_tight,
    down_tight = gls$down_tight,
    up_loose = gls$up_loose,
    down_loose = gls$down_loose,
    `DEG with batch - down in KO` = mi_gls_down$`DEG with batch - down in KO`,
    `DEG without batch - down in KO` = mi_gls_down$`DEG without batch - down in KO`,
    `DEG with batch - up in KO` = mi_gls_up$`DEG with batch - up in KO`,
    `DEG without batch - up in KO` = mi_gls_up$`DEG without batch - up in KO`
)
lengths(to_score)
lapply(to_score, function(x)intersect(x, rownames(dksc))) %>% lengths

dksc.scored = AddModuleScore(dksc, to_score, name = paste0(names(to_score), "_"))
meta_dt = get_meta_dt(dksc.scored)
meta_dt$treatment = factor(meta_dt$treatment, levels = c("W", "K"))
td = colnames(meta_dt)[seq(ncol(meta_dt)-2-length(to_score), ncol(meta_dt)-3)]
td
i = 1
plots_score = lapply(seq_along(td), function(i){
    ggplot(meta_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = td[i])) +
        geom_point(size = .6) +
        labs(title = td[i], color = "Module Score") +
        scale_color_viridis_c(option = "B") +
        facet_wrap(~treatment) +
        theme(legend.position = "bottom", panel.background = element_rect(fill= "gray30"))
    
})
# cowplot::plot_grid(plotlist = plots_score[1:2])
# cowplot::plot_grid(plotlist = plots_score[3:4])
# cowplot::plot_grid(plotlist = plots_score[5:6])

pg = cowplot::plot_grid(plotlist = plots_score[1:8], ncol = 2)
ggsave("module_scores.png", pg, width = 9.3, height = 4*4)

dksc.markers = FindAllMarkers(dksc)
fwrite(dksc.markers, "markers_scCluster.csv")

to_heatmap = list(
    tight = c(gls$up_tight, gls$down_tight),
    loose = c(gls$up_loose, gls$down_loose),
    `DEG_with_batch` = c(mi_gls_down$`DEG with batch - down in KO`, mi_gls_up$`DEG with batch - up in KO`),
    `DEG_without_batch` = c(mi_gls_down$`DEG without batch - down in KO`, mi_gls_up$`DEG without batch - up in KO`),
    # scCluster_markers_all = dksc.markers$gene %>% unique,
    scCluster_markers_2fc = subset(dksc.markers, abs(avg_logFC) > 1)$gene %>% unique
)
lengths(to_heatmap)
setwd("heatmaps")
for(i in seq_along(to_heatmap)[5:6]){
    nam = names(to_heatmap)[i]
    x = to_heatmap[[i]]
    message(nam)
    rna_dt = get_rna_dt(dksc, x)
    rna_mat = my_dcast(rna_dt, value.var = "expression", colnames.var = "id", rownames.var = "gene_name")
    set.seed(0)
    hres = double_clustering(rna_mat)
    col_clust = cutree(as.hclust(hres$colDendrogram), 50)
    
    
    rna_dt = merge(rna_dt, data.table(id = names(col_clust), cluster_id = col_clust), by = "id")
    rna_dt = relevel_clusters(rna_dt, hres, row.var = "gene_name", column.var = "id")
    
    rna_dt[]
    
    agg_dt = rna_dt[, .(expression = mean(expression), fraction_K = sum(grepl("_K_", id)) / .N), .(gene_name, cluster_id)]
    agg_dt[, .(cluster_id, fraction_K)] %>% unique
    c_order = rna_dt[order(id)]$cluster_id %>% unique
    
    agg_dt$cluster_id = factor(agg_dt$cluster_id, levels = c_order)
    
    ggplot(rna_dt, aes(x = id, y = gene_name, fill= expression)) + geom_raster() +
        theme(panel.background = element_blank(),
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank()) +
        scale_fill_viridis_c() +
        coord_cartesian(expand = FALSE)
    
    p = ggplot(agg_dt, aes(x = cluster_id, y = gene_name, fill= expression, hidden = fraction_K)) + geom_raster() +
        theme(panel.background = element_blank())+
        #axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank()) +
        scale_fill_viridis_c() +
        coord_cartesian(expand = FALSE) +
        labs(title = nam, y = "")
    
    htmlwidgets::saveWidget(plotly::ggplotly(p), 
                            selfcontained = TRUE, 
                            file = paste0("heatmap_", nam, ".html"))
}
