
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

pdf("DESeq2_analysis_no_interaction.pdf")


samples <- read_tsv("sample_info.txt")
samples$condition <- paste(samples$siRNA, samples$TimePoint, sep = "_")
samples$condition <- relevel(factor(samples$condition), ref = "TRF2_0")

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

vsd <- vst(ddsTxi, blind=TRUE)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$siRNA, vsd$TimePoint, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("siRNA", "TimePoint"))

ddsTxi <- DESeq(ddsTxi)

# Test for genes differentially expressed between timePoints in the TRF2 shRNA data
res48h <- results(ddsTxi, contrast = c("condition", "TRF2_48", "TRF2_0"))
res96h <- results(ddsTxi, contrast = c("condition", "TRF2_96", "TRF2_0"))

# Test for genes differentially expressed betweent timePoints in the control shRNA data
res48hControl <- results(ddsTxi, contrast = c("condition", "control_48", "control_0"))
res96hControl <- results(ddsTxi, contrast = c("condition", "control_96", "control_0"))

# Explicitly test for genes with a low LFC in the control timecourse
res48hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_48"))

res96hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_96"))


plotMA(res48h, main = "res48h")
plotMA(res96h, main = "res96h")

plotMA(res48hControl, main = "res48hControl")
plotMA(res96hControl, main = "res96hControl")

plotMA(res48hControl_lowLFC)
plotMA(res96hControl_lowLFC)


selected_genes_48h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res48h[which(res48h$padj < 0.1),])))

selected_genes_96h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res96h[which(res96h$padj < 0.1),])))

selected_genes_48h_control <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res48h[which(res48hControl$padj < 0.1),])))

selected_genes_96h_control <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res96h[which(res96hControl$padj < 0.1),])))

selected_genes_48h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res48h[which(res48hControl_lowLFC$padj < 0.1),])))

selected_genes_96h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res96h[which(res96hControl_lowLFC$padj < 0.1),])))

# Include only the genes that also have a low LFC in the control samples
selected_genes_48h <- selected_genes_48h[(names(selected_genes_48h) %in% names(selected_genes_48h_control_low_LFC))]
selected_genes_96h <- selected_genes_96h[(names(selected_genes_96h) %in% names(selected_genes_96h_control_low_LFC))]


gene_counts <- data.frame()

for (gene in names(selected_genes_96h)) {
  
  gene_id <- selected_genes_96h[gene]$gene_id
  symbol <- selected_genes_96h[gene]$symbol
  gene_counts <- rbind(gene_counts, cbind(plotCounts(ddsTxi, gene_id,
                                                     intgroup = c("siRNA", "TimePoint"),
                                                     returnData = TRUE), gene_id = gene_id, symbol = symbol))
  
}

gene_counts$sample <- rownames(gene_counts)

p <- gene_counts %>%
  separate(sample, into = c("shRNA", "timepoint", "replicate")) %>%
  select(-one_of(c("shRNA", "timepoint"))) %>%
  spread("TimePoint", "count") %>%
  group_by(gene_id, siRNA) %>%
  summarise(`relative 0` = 1,
            `relative 48` = mean(`48`) / mean(`0`),
            `relative 96` = mean(`96`) / mean(`0`)) %>%
  gather(timePoint, relative_counts, `relative 0`, `relative 48`, `relative 96`) %>%
  ggplot(aes(x = timePoint, y = relative_counts,
             group = interaction(gene_id, siRNA), color = siRNA)) +
  geom_line(alpha = I(1/5)) + scale_y_log10() + facet_wrap(~ siRNA)
print(p)


gene_counts <- data.frame()

for (gene in names(selected_genes_48h)) {
  
  gene_id <- selected_genes_48h[gene]$gene_id
  symbol <- selected_genes_48h[gene]$symbol
  gene_counts <- rbind(gene_counts, cbind(plotCounts(ddsTxi, gene_id,
                                                     intgroup = c("siRNA", "TimePoint"),
                                                     returnData = TRUE), gene_id = gene_id, symbol = symbol))
  
}

gene_counts$sample <- rownames(gene_counts)

p <- gene_counts %>%
  separate(sample, into = c("shRNA", "timepoint", "replicate")) %>%
  select(-one_of(c("shRNA", "timepoint"))) %>%
  spread("TimePoint", "count") %>%
  group_by(gene_id, siRNA) %>%
  summarise(`relative 0` = 1,
            `relative 48` = mean(`48`) / mean(`0`),
            `relative 96` = mean(`96`) / mean(`0`)) %>%
  gather(timePoint, relative_counts, `relative 0`, `relative 48`, `relative 96`) %>%
  ggplot(aes(x = timePoint, y = relative_counts,
             group = interaction(gene_id, siRNA), color = siRNA)) +
  geom_line(alpha = I(1/5)) + scale_y_log10() + facet_wrap(~ siRNA)
print(p)


# KEGG analysis

library("KEGGprofile")

library("topGO")

library("gplots")

sel_48h_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_48h)$entrezid),
                                             species = "hsa", download_latest = TRUE)

sel_96h_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_96h)$entrezid),
                                             species = "hsa", download_latest = TRUE)

sel_48hControl_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_48h_control)$entrezid),
                                             species = "hsa", download_latest = TRUE)

sel_96hControl_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_96h_control)$entrezid),
                                             species = "hsa", download_latest = TRUE)

textplot(sel_96h_KEGG_result$stastic)

textplot(sel_96hControl_KEGG_result$stastic)


# GO analysis


teste <- rownames(res96h)

names(teste) <- rownames(res96h)

teste[names(teste) %in% names(selected_genes_96h)] <- "selected"
teste[(names(teste) %in% names(selected_genes_96h)) == FALSE] <- "not_selected"
teste <- as.factor(teste)

sampleGOdata <- new("topGOdata", description = "teste", ontology = "BP",
                    allGenes = teste, geneSel = function(x) { (x == "selected") },
                    nodeSize = 10,
                    annot = annFUN.org,
                    mapping = "org.Hs.eg.db", ID = "ensembl")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher",
                   topNodes = 10)

textplot(allRes)


sampleGOdata <- new("topGOdata", description = "teste", ontology = "MF",
                    allGenes = teste, geneSel = function(x) { (x == "selected") },
                    nodeSize = 10,
                    annot = annFUN.org,
                    mapping = "org.Hs.eg.db", ID = "ensembl")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher",
                   topNodes = 10)

textplot(allRes)



dev.off()


