#load and process data with monocle
#should be run from monocle directory
stopifnot(basename(getwd()) == "monocle")
library(monocle3)
library(Seurat)
library(ggplot2)
library(data.table)
library(magrittr)
if(!require(ssvRecipes)){
    devtools::install_github("jrboyd/ssvRecipes")
}

source("functions_monocle.R")
# source("functions_setup.R")
source("Step0_setup.R")

seu = readRDS(seurat_obj_rds)

#if seurat doesn't have cell cycle info
# seu = CellCycleScoring(seu, 
                 # s.features = gene_name_hs2mm(Seurat::cc.genes.updated.2019$s.genes), 
                 # g2m.features = gene_name_hs2mm(Seurat::cc.genes.updated.2019$g2m.genes))

if(!file.exists(monocle_obj_rds)){
    if(length(data_dirs) < 1) stop("No data directories located. Check `data_path` and `pattern`.")
    monocle_loaded = lapply(data_dirs, function(d){
        #load data
        mon = monocle3::load_cellranger_data(pipestance_path = d, genome = GEN)
        #creates cell barcodes compatible with Seurat
        prefix = dir2prefix[d]
        colnames(mon) = paste0(prefix, "_", sub("-.+", "", colnames(mon)))
        mon
    })
    total_loaded = sum(sapply(monocle_loaded, ncol))
    #combine all monocle objects
    if(length(monocle_loaded) > 1){
        mon = cbind(monocle_loaded[[1]], monocle_loaded[[2]])    
    }else{
        mon = monocle_loaded
    }
    if(length(monocle_loaded) > 2){
        for(i in seq(3, length(monocle_loaded))){
            mon = cbind(mon, monocle_loaded[[i]])
        }
    }
    #check we combined everything
    stopifnot(ncol(mon) == total_loaded)
    
    saveRDS(mon, file = monocle_obj_rds)
}else{
    mon = readRDS(monocle_obj_rds)
}

#see how seurat and monocle match
seqsetvis::ssvFeatureVenn(list(monocle = colnames(mon), seurat = colnames(seu))) +
    labs(title = "intersection between seurat and monocle barcodes")
n_before = ncol(mon)
mon = mon[, colnames(seu)]
n_after = ncol(mon)

message(n_before - n_after , " cells removed based on Seurat object. ", n_after, " cells in processed Monocle object.")


run_process = TRUE
if(file.exists(monocle_processed_rds)){
    run_process = shiny_overwrite()
}

if(run_process){
    #if you want to bring over cluster or other info from Seurat, do it here
    pData(mon)$seurat_names = rename_clust[seu$seurat_clusters[rownames(pData(mon))]]
    pData(mon)$seurat_clusters = seu$seurat_clusters[rownames(pData(mon))]
    pData(mon)$genotype = seu$treatment[rownames(pData(mon))]
    pData(mon)$rep = seu$rep[rownames(pData(mon))]
    pData(mon)$Phase = seu$Phase[rownames(pData(mon))]
    pData(mon)$G2M.Score = seu$G2M.Score[rownames(pData(mon))]
    pData(mon)$S.Score = seu$S.Score[rownames(pData(mon))]
    
    mon = mon[, colnames(seu)]
    
    table(pData(mon)$seurat_clusters)
    
    
    mon <- preprocess_cds(mon, num_dim = num_pca)
    
    red_line_y = .05
    plot_pc_variance_explained(mon) +
        annotate("line", x= c(0, num_pca), y = rep(red_line_y, 2), color = "red")
    
    
    if(align_apply){
        mon <- align_cds(mon, alignment_group = "genotype")    
    }
    
    mon = reduce_dimension(mon, reduction_method = "UMAP", cores = 20)
    mon = cluster_cells(mon)
    #minimal_branch_len is quite important and controls the detail of branches
    mon <- learn_graph(mon, learn_graph_control = list(minimal_branch_len = minimal_branch_len), use_partition = TRUE)
    plot_cells(mon)
    saveRDS(mon, monocle_processed_rds)
}else{
    mon = readRDS(monocle_processed_rds)
}
## Step 6: Order cells
colData(mon)$seurat_clusters
#random plots
monocle3::plot_cells(mon, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE, group_label_size = 8) +
    scale_color_manual(values = get_clusters_colors())
monocle3::plot_cells(mon, color_cells_by = "genotype")
monocle3::plot_cells(mon, color_cells_by = "genotype", show_trajectory_graph = FALSE) +
    scale_color_manual(values = c(wt = "black", df4 = "red")) +
    guides(color = guide_legend())

#there are some annoying monocle3 plotting quirks, this is my way around them
plot_data = my_plot_cells(mon, return_data = TRUE)
p_dt = as.data.table(plot_data$data_df)
setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_clusters)]
#combining p_dt with lab_dt allows labels to avoid points
lab_dt = rbind(lab_dt, p_dt[, .(seurat_clusters = "", UMAP_1, UMAP_2) ])
p_seurat = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
    geom_point(size = .3) +
    ggrepel::geom_text_repel(data = lab_dt, aes(label = seurat_clusters), color = "black", size = 5) +
    scale_color_manual(values = get_clusters_colors()) +
    theme_classic()
plot(p_seurat)
