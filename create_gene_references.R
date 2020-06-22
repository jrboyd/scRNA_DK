library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(seqsetvis)
amit_dt = read.xlsx("gene_references/Amit 2017 Alzheimers DAM profile FDR 0.01.xlsx") %>% as.data.table()
colonna_dt = read.xlsx("gene_references/Colonna 2020 DAM Microglial transcripts.xlsx") %>% as.data.table()
amit_dt[, avg_logFC := as.numeric(`Fold-change.(DAM.to.homeostatic.microglia)`)]


hist(amit_dt$`Fold-change.(DAM.to.homeostatic.microglia)` %>% as.numeric())

ssvFeatureVenn(list(Amit = amit_dt[avg_logFC > 0]$Gene.name, 
                    Colonna = colonna_dt[avg_logFC > 0]$gene))

ssvFeatureVenn(list(Amit = amit_dt[avg_logFC < 0]$Gene.name, 
                    Colonna = colonna_dt[avg_logFC < 0]$gene))

amit_dt$source = "Amit"

colonna_dt$source = "Colonna"
amit_dt[, direction := ifelse(avg_logFC > 0, "up", "down")]
colonna_dt[, direction := ifelse(avg_logFC > 0, "up", "down")]

amit_dt = amit_dt[!is.na(direction)]

comb_dt = rbind(amit_dt[, .(gene_name = Gene.name, source, direction)], 
                colonna_dt[, .(gene_name = gene, source, direction)])
comb_dt[, .N, .(source, direction)]

contig_dt = dcast(comb_dt, gene_name~source, value.var = "direction", fill = "notDE")
contig_dt = contig_dt[, .N, .(Amit, Colonna)] %>% dcast(., Amit~Colonna, value.var = "N", fill = 0)
saveRDS(contig_dt, file = "contig_dt.Rds")
 

up_amit = comb_dt[source == "Amit" & direction == "up"]$gene_name %>% unique
up_colonna = comb_dt[source == "Colonna" & direction == "up"]$gene_name %>% unique

down_amit = comb_dt[source == "Amit" & direction == "down"]$gene_name %>% unique
down_colonna = comb_dt[source == "Colonna" & direction == "down"]$gene_name %>% unique


up_either = union(up_amit, up_colonna)
down_either = union(down_amit, down_colonna)

up_loose = setdiff(up_either, down_either)
down_loose = setdiff(down_either, up_either)

up_tight = intersect(up_amit, up_colonna)
down_tight = intersect(down_amit, down_colonna)

saveRDS(list(up_amit = up_amit, 
             up_colonna = up_colonna, 
             down_amit = down_amit, 
             down_colonna  = down_colonna,
             down_loose = down_loose,
             down_tight = down_tight,
             up_loose = up_loose,
             up_tight = up_tight), 
        file = "gene_references/AmitColonna_gene_lists.Rds")
