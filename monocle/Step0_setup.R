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

#don't mess with this please
prefix_processed = "processed"
prefix_subset = "subset_branch"
file_loaded = "monocle_loaded.Rds"
file_processed = "monocle_processed.Rds"
file_subset = "monocle_selected.Rds"

#this seurat object should already exist
# seurat_obj_rds = "data_alyssa/DKSCintegrated7dimsQC.rds"
seurat_obj_rds = "DKSCintegrated7dimsQC.rds"
stopifnot(file.exists(seurat_obj_rds))
#this monocle object will be created if needed


#if you want to start fresh, alter this path or move existing
output_root_dir = "output"
dir.create(output_root_dir, showWarnings = FALSE, recursive = TRUE)
monocle_obj_rds = file.path(output_root_dir, file_loaded)



## Step 1 parameters
GEN = "mm10"
#monocle has an alignment procedure analogous to integration in Seurat
align_apply = FALSE
#number of PC to use for reduction
#I have it at 5 here because beyond that PCs explain less that 5% of the variance
num_pca = 5
#increase to remove small branches and decrease to see more
minimal_branch_len = 6 

proc_desc = paste0(ifelse(align_apply, "align", "unalign"), ".", 
                   num_pca, "_PC.",
                   minimal_branch_len, "_branchLen")

get_proc_dir = function(name, out_dir = output_root_dir, prefix = prefix_processed, suffix = proc_desc){
    file.path(out_dir,
              paste0(
                  prefix,
                  ".",
                  name,
                  ".", 
                  suffix)
    )
}

default_proc_dir = get_proc_dir("full")

dir.create(default_proc_dir, showWarnings = FALSE, recursive = TRUE)

default_monocle_processed_rds = file.path(default_proc_dir,
                                  paste0(file_processed))

cluster_identities = read.table("cluster_idents.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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
    
    add_proc_prefix = function(x){paste0(prefix_processed, ".", x)}
    rm_proc_prefix = function(x){sub(paste0(prefix_processed, "\\."), "", x)}
    
    
    add_subset_prefix = function(x){paste0(prefix_subset, ".", x)}
    rm_subset_prefix = function(x){sub(paste0(prefix_subset, "\\."), "", x)}
    
    ui <- miniPage(
        gadgetTitleBar("Select processed and subset", left = NULL, right = NULL),
        selectInput("select_processed", label = "Select Processed", choices = rm_proc_prefix(proc_dirs)),
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
            proc_dir = file.path(output_path, add_proc_prefix(input$select_processed))
            available_subsets = dir(proc_dir, pattern = anlaysis_pattern)
            selectInput("select_subset", label = "Select Subset", choices = rm_subset_prefix(available_subsets))
        })
        
        observeEvent(input$btnExisting, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, add_proc_prefix(input$select_processed), add_subset_prefix(input$select_subset)))
        })
        
        # Handle the Done button being pressed.
        observeEvent(input$btnNew, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, add_proc_prefix(input$select_processed), add_subset_prefix(input$txt_new_subset)))
        })
        
        
    }
    runGadget(ui, server, viewer = dialogViewer("subset"))    
}

shiny_new_name = function(starting_value = "", output_path = "output", processed_pattern = paste0("^", prefix_processed)){
    stopifnot(dir.exists(output_path))
    proc_dirs = dir(output_path, pattern = processed_pattern)
    proc_names = sub(paste0(prefix_processed, "."), "", proc_dirs)
    
    add_proc_prefix = function(x){paste0(prefix_processed, ".", x)}
    rm_proc_prefix = function(x){sub(paste0(prefix_processed, "\\."), "", x)}
    
    ui <- miniPage(
        shinyjs::useShinyjs(),
        gadgetTitleBar("Create new processed dataset", left = NULL, right = NULL),
        tags$br(),
        textInput("txt_new_processed", label = "New processed name", value = starting_value),
        shinyjs::hidden(tags$h3("dataset must not exist, change the name.", style="color:red;", id = "warn_exist")),
        actionButton("btnNew", label = "Use New")
        
    )
    
    server <- function(input, output, session) {
        
        observeEvent({
            input$txt_new_processed
        }, {
            if(input$txt_new_processed %in% proc_names){
                shinyjs::disable("btnNew")
                shinyjs::show("warn_exist")
                
            }else{
                shinyjs::enable("btnNew")
                shinyjs::hide("warn_exist")
            }
        })
        # Handle the Done button being pressed.
        observeEvent(input$btnNew, {
            if(!input$txt_new_processed %in% proc_names){
                stopApp(input$txt_new_processed)
            }
            # Return the brushed points. See ?shiny::brushedPoints.
            
        })
        
        
    }
    runGadget(ui, server, viewer = dialogViewer("new name"))    
}

processed_plots = function(mon, out_dir){
    #there are some annoying monocle3 plotting quirks, this is my way around them
    plot_data = my_plot_cells(mon, return_data = TRUE)
    p_dt = as.data.table(plot_data$data_df)
    seg_dt = as.data.table(plot_data$edge_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_names, seurat_clusters)]
    #combining p_dt with lab_dt allows labels to avoid points
    lab_dt = rbind(lab_dt, p_dt[, .(seurat_names = "", seurat_clusters = "", UMAP_1, UMAP_2) ])
    
    p_seurat = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = seurat_names)) +
        geom_point(size = .6) +
        geom_segment(data = seg_dt, 
                     aes(x = source_prin_graph_dim_1, 
                         y = source_prin_graph_dim_2, 
                         xend = target_prin_graph_dim_1, 
                         yend = target_prin_graph_dim_2), 
                     color= "black", size = 1) +
        ggrepel::geom_text_repel(data = lab_dt, aes(label = seurat_names), color = "black", size = 5) +
        theme_classic() +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        labs(title = "Monocle processed data")
    plot(p_seurat)
    ggsave(file.path(out_dir, "plot_processed_seurat_clusters.png"), p_seurat, 
           width = 6, height = 5)
    ggsave(file.path(out_dir, "plot_processed_seurat_clusters.pdf"), p_seurat, 
           width = 6, height = 5)
    
    
    
    plots = list(
        seurat_cluster_choose_plot(mon, show.legend = TRUE, color_by = "seurat_names") + labs(title = "Seurat clusters"),
        seurat_genotype_choose_plot(mon, show.legend = TRUE, col_scale = c("W" = "black", "K" = "red")) + labs(title = "Genotype"),
        partition_choose_plot(mon, show.legend = TRUE) + labs(title = "Monocle partition"),
        seurat_cell_cycle_choose_plot(mon, show.legend = TRUE)  + labs(title = "Cell cycle")
    )
    pg = cowplot::plot_grid(plotlist = plots)
    
    ggsave(file.path(out_dir, "plot_procesed_overview.png"), pg, width = 8*1.5, height = 6*1.5)
    ggsave(file.path(out_dir, "plot_procesed_overview.pdf"), pg, width = 8*1.5, height = 6*1.5)
    
    pg_score = cowplot::plot_grid(ncol = 1,
                                  seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "S.Score", show.trajectory = FALSE)  + 
                                      labs(title = "S") +
                                      facet_wrap(~genotype) +
                                      theme(panel.background = element_rect(fill = "gray30")),
                                  seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "G2M.Score", show.trajectory = FALSE)  + 
                                      labs(title = "G2M") +
                                      facet_wrap(~genotype)+
                                      theme(panel.background = element_rect(fill = "gray30"))
    )
    ggsave(file.path(out_dir, "plot_processed_cell_cycle_scores.png"), pg_score, width = 6.5*1.5, height = 6*1.5)
    ggsave(file.path(out_dir, "plot_processed_cell_cycle_scores.pdf"), pg_score, width = 6.5*1.5, height = 6*1.5)
    
    red_line_y = .05
    p_var = plot_pc_variance_explained(mon) +
        annotate("line", x= c(0, num_pca), y = rep(red_line_y, 2), color = "red")
    
    ggsave(file.path(out_dir, "plot_processed_PCvariance.png"), p_var, width = 4, height = 4)
    ggsave(file.path(out_dir, "plot_processed_PCvariance.pdf"), p_var, width = 4, height = 4)
}
