#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(shinycssloaders)
library(tippy)
library(Seurat)

source("tippy.R")

gene_name_gs2mm = function(x){
    x = tolower(x)
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x = sub("(?<=[0-9])rik$", "Rik", x, perl = TRUE)
    is_rik = grepl("Rik$", x)
    substr(x[is_rik], 8, 8) = toupper(substr(x[is_rik], 8, 8))
    x
}

js <- '$("#selDataSource :input").each(function() {
    $(this).attr("id", "radio_" + $(this).val());
});'

js2 <- '$("#selDataSource span").each(function() {
    $(this).attr("id", "radioLabel_" + $(this).text());
});'
js2a <- '$("#selDataSource span").each(function() {
    $(this).parents()[1].id = "radioItem_" + $(this).text();
});'

js3 <- '$("#txtGene").parents()[1].id = "txtGene-container"'
js_selGeneLookup <- '$("#selGeneLookup").parents()[1].id = "selGeneLookup-container"'
get_meta_dt = function(seurat_obj){
    meta_dt = as.data.table(seurat_obj@meta.data)
    meta_dt$id = seurat_obj@meta.data %>% rownames
    
    umap_dt = as.data.table(seurat_obj@reductions$umap@cell.embeddings, keep.rownames = TRUE) %>% setnames(., "rn", "id")
    meta_dt =merge(meta_dt, umap_dt, by = "id")
    meta_dt$treatment = factor(meta_dt$treatment, levels = c("W", 'K'))
    meta_dt
}

get_rna_dt = function(seurat_obj, sel_genes = NULL){
    if(is.null(sel_genes)){
        rna_dt = as.data.frame(seurat_obj@assays$RNA)
    }else{
        rna_dt = as.data.frame(seurat_obj@assays$RNA[sel_genes, ])
        
    }
    rna_dt$gene_name = rownames(rna_dt)
    rna_dt = melt(as.data.table(rna_dt), variable.name = "id", value.name = "expression", id.vars = "gene_name")
    rna_dt
}



pbmc_sources_dir = "/home/user/Documents/R_workspace/Bcell_IKdf4"
if(!dir.exists(pbmc_sources_dir)){
    pbmc_sources_dir = "/slipstream/home/joeboyd/R/scRNA_DK/datasets/"  
}
if(!dir.exists(pbmc_sources_dir)) stop("bad dir")
# 
remove_toy = TRUE
pbmc_sources = list(
    "toy" = "DKSC.toy.Rds",
    "anchored" = "DKSC.integrated.Rds",
    "no anchor" = "DKSC.combined.Rds")


pbmc_loaded = lapply(pbmc_sources, function(x)NULL)
pbmc_sources = lapply(pbmc_sources, function(x)file.path(pbmc_sources_dir, x))

if(!file.exists(pbmc_sources$toy)){
    pbmc_small = readRDS(pbmc_sources$anchored)
    DefaultAssay(pbmc_small) = "RNA"
    k = grepl("^Cd", rownames(pbmc_small))
    pbmc_small = pbmc_small[k,]
    k = sample(seq(ncol(pbmc_small), 500))
    pbmc_small = pbmc_small[k,]
    saveRDS(pbmc_small, pbmc_sources$toy)
}



k = sapply(pbmc_sources, file.exists)
if(!all(k)) warning("some data sources not found, will be ignored.")
pbmc_loaded = pbmc_loaded[k]
pbmc_sources = pbmc_sources[k]

if(remove_toy){
    pbmc_loaded$toy = NULL
    pbmc_sources$toy = NULL
}else{
    pbmc_loaded$toy = readRDS(pbmc_sources$toy)    
}


stopifnot(names(pbmc_loaded) == names(pbmc_sources))

# Define UI for application that draws a histogram
ui <- fluidPage(
    shinyjs::useShinyjs(),
    actionButton("btnSelectUMAP", "Select UMAP"),
    checkboxInput("checkSplitGenotype", "Split By Genotype", value = TRUE),
    tippy_this("checkSplitGenotype", "If checked, all plots will facet data based on genotype."),
    # Application title
    titlePanel("scRNAseq of wt and ko microglial cells"),
    
    tabsetPanel(id = "tabset1", 
                # Sidebar with a slider input for number of bins
                tabPanel("Gene UMAP Query", id = "tag_goi",
                         sidebarLayout(
                             sidebarPanel(
                                 selectizeInput(inputId = "txtGene", label = "Select Genes", choices = NULL, multiple = FALSE),
                                 tags$script(js3),
                                 tippy_this(
                                     "txtGene-container",
                                     tooltip = "Delete and start typing gene names to see matching genes."
                                 )
                             ),
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                                 withSpinner(plotOutput("plotGOI", width = "550px", height = "310px")),
                                 tippy_this(
                                     "plotGOI",
                                     tooltip = paste("scRNAseq UMAP facetted by mutation status of cells with",
                                                     "expression values of the selected gene mapped to color.", 
                                                     "Select a different UMAP or gene on the left.",
                                                     "Right-click and save or copy image to export plot.")
                                 )
                             )
                         )
                )
    ),
    tippy_this("tabset1", "Gene UMAP Query")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    shinyjs::hide(id = "msgGoodModule")
    shinyjs::disable(id = "btnViewModuleDetail")
    pbmc_data = reactiveVal(NULL)
    active_source = reactiveVal(names(pbmc_loaded)[1])
    
    observeEvent(input$btnSelectUMAP, {
        showModal(modalDialog(
            radioButtons(inputId = "selDataSource", label = "Select UMAP", choices = names(pbmc_loaded), selected = active_source()),
            tags$script(js),
            tags$script(js2),
            tags$script(js2a),
            tippy_datasets(),
            actionButton("btnLoadModal", "Load"),
            title = "Select a different UMAP",
            "These can take a minute or two to load.",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    
    observeEvent({
        input$btnLoadModal
    }, {
        removeModal()
        sel = input$selDataSource
        stopifnot(sel %in% names(pbmc_loaded))
        active_source(sel)
        
    })
    
    observeEvent({
        active_source()
    }, {
        sel = active_source()
        if(is.null(pbmc_loaded[[sel]])){
            showNotification(paste0("loading ", sel, " UMAP..."), duration = NULL, id = "load_note")
            pbmc_loaded[[sel]] = readRDS(pbmc_sources[[sel]])
            Seurat::DefaultAssay(pbmc_loaded[[sel]]) = "RNA"
            removeNotification(id = "load_note")
        }
        
        pbmc_data(pbmc_loaded[[sel]])
        #invalidate derived data
        meta_dtR(NULL)
        module_dtR(NULL)
        
        all_genes = sort(rownames(pbmc_data()@assays$RNA))
        start = min(which(grepl("^Cd", all_genes)))
        
        all_genes = all_genes[c(seq(start, length(all_genes)), seq(1, start - 1))]
        def = "Cd84"
        if(def %in% all_genes){
            updateSelectizeInput(session, 'txtGene', choices = all_genes, selected = def, server = TRUE)    
        }else{
            updateSelectizeInput(session, 'txtGene', choices = all_genes, server = TRUE)    
        }
        
    })
    
    output$plotGOI <- renderPlot({
        goi = input$txtGene
        req(goi)
        
        x = pbmc_data()
        req(x)
        req(goi != '')
        meta_dt = get_meta_dt(x)
        meta_dt$treatment = factor(meta_dt$treatment, levels = c("W", "K"))
        rna_dt = get_rna_dt(x, goi)
        rna_dt = merge(meta_dt, rna_dt, by = "id")
        rna_dt = rna_dt[order(expression, decreasing = FALSE)]
        p = ggplot(rna_dt, aes(x = UMAP_1, y = UMAP_2, color = expression)) + 
            geom_point(size= .1) +
            theme_cowplot() +
            theme(panel.background = element_rect(fill = "gray70"), 
                  strip.background = element_blank()) +
            scale_color_viridis_c() +
            labs(title = goi, color = "expression")    
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~treatment)    
        }
        p
    })
    
    loaded_coarse = reactiveVal()
    loaded_fine = reactiveVal()
    meta_dtR = reactiveVal()
    
    observe({
        x = req(pbmc_data())
        module_resolution = req(input$selModuleResolution)
        
        meta_dt = get_meta_dt(x)
        
        meta_dt$treatment = factor(meta_dt$treatment, levels = c("W", "K"))
        
        if(module_resolution == "coarse"){
            cn = colnames(coarse_dt)[grepl("coarse_", colnames(coarse_dt))]
            meta_dt = merge(meta_dt, coarse_dt[, c("id", cn), with = FALSE], by = "id")    
        }else if(module_resolution == "fine"){
            cn = colnames(fine_dt)[grepl("fine_", colnames(fine_dt))]
            meta_dt = merge(meta_dt, fine_dt[, c("id", cn), with = FALSE], by = "id")    
        }
        meta_dt[["seurat_clustersAndGenotype"]] = factor(paste(meta_dt[["seurat_clusters"]], meta_dt[["treatment"]]))
        meta_dtR(meta_dt)
    })
    
    
    
    output$clustersAvailable = renderUI({
        cluster_variable = req(input$selClusterSet)
        meta_dt = req(meta_dtR())
        selectInput("selClusterIDs", label = "Select Clusters", choices = levels(meta_dt[[cluster_variable]]), 
                    selected = levels(meta_dt[[cluster_variable]])[1:3],
                    multiple = TRUE)
    })
    
    module_dtR = reactiveVal()
    
    output$coarseTable = DT::renderDT({
        # radioButtons("selClusterMethod",
        #              label = "Method",
        #              choices = c("identity", "variance", "average", "max", "min", "difference"))
        module_resolution = req(input$selModuleResolution)
        meta_dt = req(meta_dtR())
        sel_clust = req(input$selClusterIDs)
        cluster_variable = req(input$selClusterSet)
        sort_method = req(input$selClusterMethod)
        # browser()
        k = as.character(meta_dt[[cluster_variable]]) %in% sel_clust
        if(!any(k)) req(NULL)
        
        if(module_resolution == "coarse"){
            cn = colnames(coarse_dt)[grepl("coarse_", colnames(coarse_dt))]
        }else if(module_resolution == "fine"){
            cn = colnames(fine_dt)[grepl("fine_", colnames(fine_dt))]
        }
        
        suppressWarnings({
            melted = melt(meta_dt[k, ][, c(cluster_variable, cn), with = FALSE], id.vars = cluster_variable, variable.name = "module", measure.vars = cn)
            module_dt = melted[, .(mean_score = mean(value)), c(cluster_variable, "module")]    
        })
        
        if(sort_method == "identity"){
            lev = module_dt[[cluster_variable]] %>% unique %>% sort %>% as.character()
            module_dt = dcast(module_dt, paste0("module~", cluster_variable), value.var = "mean_score")
        }else if(sort_method == "variance"){
            lev = "variance"
            module_dt = module_dt[, .(variance = var(mean_score)), .(module)]
            
        }else if(sort_method == "average"){
            lev = "average"
            module_dt = module_dt[, .(average = mean(mean_score)), .(module)]
            
        }else if(sort_method == "max"){
            lev = "maximum"
            module_dt = module_dt[, .(maximum = max(mean_score)), .(module)]
            
        }else if(sort_method == "min"){
            lev = "minimum"
            module_dt = module_dt[, .(minimum = min(mean_score)), .(module)]
            
        }else if(sort_method == "difference"){
            lev = (module_dt[[cluster_variable]] %>% unique %>% as.character())
            if(length(lev) < 2) showNotification("select 2 clusters for difference", type = "warning")
            if(length(lev) > 2) showNotification("only first 2 clusters used for difference", type = "warning")
            lev = lev[1:2]
            module_dt = dcast(module_dt[get(cluster_variable) %in% lev], paste0("module~", cluster_variable), value.var = "mean_score")
            module_dt[, difference := get(lev[1]) - get(lev[2])]
            module_dt = module_dt[, c(1,4,2,3)]
            lev = c(lev, "difference")
            
        }else{
            stop("Unrecognized selClusterMethod: ", sort_method)
        }
        
        # add extra info
        if(module_resolution == "coarse"){
            module_dt = merge(module_dt, coarse_descriptions_dt, by = "module")        
        }else if(module_resolution == "fine"){
            module_dt = merge(module_dt, fine_descriptions_dt, by = "module")        
        }
        # sort
        module_dt = module_dt[order(module_dt[, 2][[1]], decreasing = TRUE),]
        module_dtR(module_dt)
        DT::datatable(module_dt, selection = "single") %>%
            DT::formatRound(lev, 3)
    })
    
    output$moduleLookup = renderUI({
        goi = req(input$selGeneLookup)
        sel = module_genes_dt[gene_name == goi,]
        stopifnot(nrow(sel) == 1)
        
        if(input$selModuleResolution == "coarse"){
            tags$span(paste0(goi, ": ", paste0("coarse_", sel$Coarse.module)),
                      style="color:blue; font-weight:bold")    
        }else{
            tags$span(paste0(goi, ": ", paste0("fine_", sel$Fine.module)),
                      style="color:blue; font-weight:bold")    
        }
    })
    
    observeEvent({
        input$coarseTable_rows_selected
    }, {
        sel = req(input$coarseTable_rows_selected)
        sel_dt = req(module_dtR())[sel,]
        stopifnot(nrow(sel_dt) == 1)
        updateTextInput(session = session, "txtModuleQuery", value = sel_dt$module)
        # showNotification(paste("selected", sel_dt$module))
    })
    
    isValidModule = reactiveVal(FALSE)
    
    observeEvent({
        input$txtModuleQuery
        meta_dtR()
    },{
        meta_dt = req(meta_dtR())
        module_query = req(input$txtModuleQuery)
        if(module_query %in% colnames(meta_dt)){
            shinyjs::hide(id = "msgBadModule")
            shinyjs::show(id = "msgGoodModule")
            shinyjs::enable(id = "btnViewModuleDetail")
            isValidModule(TRUE)
        }else{
            shinyjs::hide(id = "msgGoodModule")
            shinyjs::show(id = "msgBadModule")
            shinyjs::disable(id = "btnViewModuleDetail")
            isValidModule(FALSE)
        }
    })
    
    output$plotMiniClusters = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        
        k = as.character(meta_dt[[sel_set]]) %in% sel_id
        if(!any(k)) req(NULL)
        p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = sel_set)) + 
            geom_point(data = meta_dt[!k,], color = "gray10", size = .1) +
            geom_point(data = meta_dt[k,], size = .2) +
            labs(color = '', x= "", y = "") +
            guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
            theme_cowplot() +
            theme(legend.position = "bottom", 
                  axis.text = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.line = element_blank())
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~treatment) 
        }
        p
    })
    
    output$plotMiniModule = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        mod_query = input$txtModuleQuery
        if(is.null(mod_query)) req(NULL)
        if(isValidModule()){
            p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = mod_query)) + 
                geom_point(data = meta_dt, size = .2) +
                labs(color = '', x= "", y = "") +
                scale_color_viridis_c(option = "magma") +
                theme_cowplot() +
                theme(legend.position = "bottom", 
                      axis.text = element_blank(), 
                      axis.ticks = element_blank(), 
                      axis.line = element_blank(), 
                      legend.text = element_text(size = 6))
            if(input$checkSplitGenotype){
                p = p + facet_wrap(~treatment) 
            }
        }else{
            p = ggplot() + 
                annotate("text", x=0, y = 0, label = 'waiting for valid module') + 
                theme_void()
        }
        
        # ggsave("tmp.png", p)
        p
    })
    
    geneDetailDT = reactiveVal()
    
    output$dtGeneDetail = DT::renderDataTable({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_res = req(input$selModuleResolution)
        sel_mod = req(input$txtModuleQuery)
        if(sel_res == "fine"){
            gois = fine_genes_dt[module == sel_mod]$gene_name
        }else{
            gois = coarse_genes_dt[module == sel_mod]$gene_name
        }
        missed = setdiff(gois, rownames(pbmc_data()))
        rna_dt = get_rna_dt(pbmc_data(), intersect(gois, rownames(pbmc_data())))
        rna_dt = rna_dt[, .(mean_expression = mean(expression)), .(gene_name)]
        rna_dt = rna_dt[order(mean_expression, decreasing = TRUE)]
        rna_dt = rbind(rna_dt, data.table(gene_name = missed, mean_expression = NA))
        geneDetailDT(rna_dt)
        DT::datatable(rna_dt, selection = list(mode = 'single', selected = c(1), target = 'row')) %>%
            DT::formatRound("mean_expression", 3)
    })
    
    output$plotModuleDetail = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        mod_query = input$txtModuleQuery
        # k = as.character(meta_dt[[sel_set]]) %in% sel_id
        p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = mod_query)) + 
            # geom_point(data = meta_dt[!k,], color = "gray10", size = .1) +
            # geom_point(data = meta_dt[k,], size = .2) +
            geom_point(data = meta_dt, size = .2) +
            labs(color = '', x= "", y = "") +
            scale_color_viridis_c(option = "magma") +
            # guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
            theme_cowplot() +
            theme(legend.position = "bottom", 
                  axis.text = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.line = element_blank(),
                  legend.text = element_text(size = 6)) +
            labs(title = mod_query)
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~treatment) 
        }
        p
    })
    
    output$plotGeneDetail = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        sel_dt = req(geneDetailDT())
        if(nrow(sel_dt) == 0) req(NULL)
        
        sel_gene = sel_dt$gene_name[req(input$dtGeneDetail_rows_selected)]
        if(sel_gene %in% rownames(pbmc_data())){
            rna_dt = get_rna_dt(pbmc_data(), sel_gene)
            rna_dt = merge(meta_dt[, .(id, UMAP_1, UMAP_2, treatment)], rna_dt, by = "id")
            p = ggplot(data = rna_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = "expression")) + 
                geom_point(size = .2) +
                labs(color = '', x= "", y = "") +
                scale_color_viridis_c() +
                # guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
                theme_cowplot() +
                theme(legend.position = "bottom", 
                      axis.text = element_blank(), 
                      axis.ticks = element_blank(), 
                      axis.line = element_blank()) +
                labs(title = sel_gene)
            if(input$checkSplitGenotype){
                p = p + facet_wrap(~treatment) 
            }
        }else{
            p = ggplot() + annotate("text", x = .5, y = .5, label = "gene not found in scRNA data") +
                theme_void()
        }
        p
    })
    
    observeEvent({
        input$btnModalDetailDone
    }, {
        removeModal()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

