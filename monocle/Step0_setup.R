#this seurat object should already exist
seurat_obj_rds = "data_alyssa/DKSCintegrated7dimsQC.rds"
#this monocle object will be created if needed


#if you want to start fresh, alter this path or move existing
output_root_dir = "output"
monocle_obj_rds = file.path(output_root_dir, "monocle_loaded.Rds")

## Step 1 parameters
GEN = "mm10"
#monocle has an alignment procedure analogous to integration in Seurat
align_apply = FALSE
#number of PC to use for reduction
#I have it at 5 here because beyond that PCs explain less that 5% of the variance
num_pca = 5
#increase to remove small branches and decrease to see more
minimal_branch_len = 8 

out_dir = file.path(output_root_dir,
                    paste0(
                        ifelse(align_apply, "align", "unalign"), ".", 
                        num_pca, "_PC.",
                        minimal_branch_len, "_branchLen")
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

monocle_processed_rds = file.path(out_dir,
                                  paste0("monocle_processed.Rds"))

cluster_identities = read.table("cluster_idents.txt", sep = "\t", header = TRUE)
signature_genes = strsplit(cluster_identities$signature_genes, ",")

get_cluster_rename = function(){
    rename_clust = cluster_identities$bio_name
    names(rename_clust) = cluster_identities$seurat_clusters
    rename_clust
}

get_clusters_colors = function(){
    safeBrew = seqsetvis::safeBrew
    col.blue = safeBrew(7, "blues")[2:5+2]
    names(col.blue) = c("main", 11, 10, 7)
    col.purp = safeBrew(5, "BuPu")[c(3, 5)]
    names(col.purp) = c(4, 13)
    col.oranges = safeBrew(6, "Oranges")[c(3, 6)]
    names(col.oranges) = c(6, 15)
    col.reds = safeBrew(5, "Reds")[2:5]
    names(col.reds) = c(8, 9, 12, 14)
    col.all = c(col.blue, col.purp, col.oranges, col.reds)
    
    names(col.all) = get_cluster_rename()[names(col.all)]
    col.all
}

data_path = "~/Dimitry_scRNA/Krementsov_10Xgenomics_iLabs_14530/"
pattern = "^DK_[0-9]_[WK]_rep[0-9]$"
data_dirs = dir(data_path, pattern = pattern, full.names = TRUE)
data_dirs = normalizePath(data_dirs)

dir2prefix = basename(data_dirs); names(dir2prefix) = data_dirs


# seu = readRDS(seurat_obj_rds)
# library(patchwork)
# vln_plots = lapply(seq_len(nrow(cluster_identities)), function(i){
#     name = paste(cluster_identities$seurat_clusters[i], ":", cluster_identities$bio_name[i])
#     p = Seurat::VlnPlot(seu, signature_genes[[i]])
#     p + plot_annotation(title = name)
# })
# names(vln_plots) = cluster_identities$bio_name
# pdf("violins.pdf", width = 10, height = 7)
# vln_plots
# dev.off()

library(shiny)
library(miniUI)
library(ggplot2)

shiny_overwrite <- function(msg = "Processed data found, overwrite?", yes_prompt = "Yes, overwrite.", no_prompt = "No, use existing.") {
    
    ui <- miniPage(
        gadgetTitleBar(msg, left = NULL),
        radioButtons("overwrite", label = c(""), 
                     selected = "no",
                     choiceValues = c("yes", "no"), 
                     choiceNames = c(yes_prompt, no_prompt))
    )
    
    server <- function(input, output, session) {
        # Handle the Done button being pressed.
        observeEvent(input$done, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(input$overwrite == "yes")
        })
    }
    runGadget(ui, server, viewer = dialogViewer("overwrite"))    
}
