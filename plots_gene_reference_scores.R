library(Seurat)
library(data.table)
library(ggplot2)
library(magrittr)
library(seqsetvis)
source("../SF_AutoImmune_ssv/functions_setup.R")

dksc = readRDS("datasets/DKSC.integrated.Rds")
gls = readRDS("gene_references/AmitColonna_gene_lists.Rds")

lengths(gls)
lapply(gls, function(x)intersect(x, rownames(dksc))) %>% lengths

dksc = AddModuleScore(dksc, gls)

meta_dt = get_meta_dt(dksc)
setnames(meta_dt, paste0("Cluster", seq_along(gls)), names(gls))

theme_set(cowplot::theme_cowplot())

id_vars = c("id", "treatment", "UMAP_1", "UMAP_2", "seurat_clusters")

mod_dt = melt(meta_dt[, c(id_vars, names(gls)), with = FALSE], id.vars = id_vars)


mod_dt$treatment = factor(mod_dt$treatment, levels = c("W", "K"))
p_modules = ggplot(mod_dt[grepl("tight|loose", variable)], aes(x = UMAP_1, y=  UMAP_2, color = value)) + 
    geom_point(size = .2) +
    scale_color_viridis_c() +
    coord_fixed() +
    facet_grid(variable~treatment)

contig_dt = readRDS("contig_dt.Rds")
p_agreement = ggplot(melt(contig_dt, variable.name = "Colonna", value.name = "count", id.vars = "Amit"), 
                     aes(x = Amit, y = Colonna, label = count)) +
    geom_text() + 
    theme_classic() +
    labs(title = "Amit and Colonna DAM gene list agreement")

pg = cowplot::plot_grid(p_agreement + 
                       coord_fixed() +
                       theme(title = element_text(size = 10)), 
                   p_modules +
                       theme(strip.text = element_text(size = 8),
                             axis.text = element_text(size = 8),
                             axis.title = element_text(size = 8)), rel_widths = c(1, 1.3))
ggsave("AmitColonna_module_score")
