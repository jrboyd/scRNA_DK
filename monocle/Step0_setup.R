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
#this seurat object should already exist
seurat_obj_rds = "data_alyssa/DKSCintegrated7dimsQC.rds"
#this monocle object will be created if needed


#if you want to start fresh, alter this path or move existing
output_root_dir = "output_tmp"
monocle_obj_rds = file.path(output_root_dir, "monocle_loaded.Rds")

prefix_processed = "processed"
prefix_subset = "subset_branch"

## Step 1 parameters
GEN = "mm10"
#monocle has an alignment procedure analogous to integration in Seurat
align_apply = FALSE
#number of PC to use for reduction
#I have it at 5 here because beyond that PCs explain less that 5% of the variance
num_pca = 5
#increase to remove small branches and decrease to see more
minimal_branch_len = 6 

out_dir = file.path(output_root_dir,
                    paste0(
                        prefix_processed,
                        ".",
                        ifelse(align_apply, "align", "unalign"), ".", 
                        num_pca, "_PC.",
                        minimal_branch_len, "_branchLen")
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

monocle_processed_rds = file.path(out_dir,
                                  paste0("monocle_processed.Rds"))

cluster_identities = read.table("cluster_idents.txt", sep = "\t", header = TRUE)
signature_genes = strsplit(cluster_identities$signature_genes, ",")

rename_clust = cluster_identities$bio_name
names(rename_clust) = cluster_identities$seurat_clusters

get_cluster_rename = function(){
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

shiny_choose_branch <- function(output_path = "output", 
                                processed_pattern = paste0("^", prefix_processed), 
                                anlaysis_pattern = paste0("^", prefix_subset),
                                allow_new = TRUE) {
    proc_dirs = dir(output_path, pattern = processed_pattern)
    
    add_prefix = function(x){paste0(prefix_subset, ".", x)}
    rm_prefix = function(x){sub(paste0(prefix_subset, "\\."), "", x)}
    
    ui <- miniPage(
        gadgetTitleBar("Select processed and subset", left = NULL, right = NULL),
        selectInput("select_processed", label = "Select Processed", choices = proc_dirs),
        tags$br(),
        uiOutput("dynamic_subset"),
        actionButton("btnExisting", label = "Use Existing"),
        if(allow_new){
            tags$span(
                tags$h3("OR"),
                textInput("txt_new_subset", label = "New subset name"),
                actionButton("btnNew", label = "Use New")
            )
        }else{
            tags$span()
        }
        
    )
    
    server <- function(input, output, session) {
        
        output$dynamic_subset = renderUI({
            proc_dir = file.path(output_path, input$select_processed)
            available_subsets = dir(proc_dir, pattern = anlaysis_pattern)
            selectInput("select_subset", label = "Select Subset", choices = rm_prefix(available_subsets))
        })
        
        observeEvent(input$btnExisting, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, input$select_processed, add_prefix(input$select_subset)))
        })
        
        # Handle the Done button being pressed.
        observeEvent(input$btnNew, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, input$select_processed, add_prefix(input$txt_new_subset)))
        })
        
        
    }
    runGadget(ui, server, viewer = dialogViewer("subset"))    
}


