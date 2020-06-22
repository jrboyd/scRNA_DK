#select specific cells and also root nodes.
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
#this monocle object will be created if needed
source("Step0_setup.R")

stopifnot(file.exists(monocle_processed_rds))
mon = readRDS(monocle_processed_rds)

#patch to bring in cell clusters
stopifnot(file.exists(seurat_obj_rds))
seu = readRDS(seurat_obj_rds)
if(is.null(pData(mon)$Phase)) pData(mon)$Phase = seu$Phase[rownames(pData(mon))]
if(is.null(pData(mon)$S.Score)) pData(mon)$S.Score = seu$S.Score[rownames(pData(mon))]
if(is.null(pData(mon)$G2M.Score)) pData(mon)$G2M.Score = seu$G2M.Score[rownames(pData(mon))]

plots = list(
    seurat_cluster_choose_plot(mon, show.legend = TRUE) + labs(title = "Seurat clusters"),
    seurat_genotype_choose_plot(mon, show.legend = TRUE) + labs(title = "Genotype"),
    partition_choose_plot(mon, show.legend = TRUE) + labs(title = "Monocle partition"),
    seurat_cell_cycle_choose_plot(mon, show.legend = TRUE)  + labs(title = "Cell cycle")
)


pg = cowplot::plot_grid(plotlist = plots)
plot_file = sub(".Rds", ".png", sub("monocle_processed", "plots", monocle_processed_rds))
ggsave(plot_file, pg, width = 8*1.5, height = 6*1.5)

monocle3::partitions(mon)

pg_score = cowplot::plot_grid(ncol = 1,
    seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "S.Score", show.trajectory = FALSE)  + 
        labs(title = "S") +
        facet_wrap(~genotype),
    seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "G2M.Score", show.trajectory = FALSE)  + 
        labs(title = "G2M") +
        facet_wrap(~genotype)
)
ggsave("monocle_cell_cycle_scores.png", pg_score, width = 6.5*1.5, height = 6*1.5)

#this is where the anlaysis really starts
if(!file.exists(monocle_selected_rds)){
    #select a region of interest (or assign sel_branch = mon to anlayze everything)
    sel_branch = my_choose_cells(mon, plot_FUN = seurat_cluster_choose_plot)
    sel_branch = my_choose_cells(mon, plot_FUN = seurat_genotype_choose_plot)
    #select root node
    sel_branch <- order_cells(sel_branch)
    #if you like your results, rename this file so you can repeat them.
    saveRDS(sel_branch, file = monocle_selected_rds)
    
}else{
    sel_branch = readRDS(monocle_selected_rds)
}

plot_cells(sel_branch, color_cells_by = "pseudotime")
