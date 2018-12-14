
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")
library("ggrepel")
library("gridExtra")

library("genefilter")
library("geneplotter")
library("VennDiagram")


setwd("/mnt/TRF2_siRNA/")

pdf("DESeq2_analysis_no_interaction_results/expression.pdf", width = 10)


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

#plotMA(res48hControl, main = "res48hControl")
#plotMA(res96hControl, main = "res96hControl")

plotMA(res48hControl_lowLFC, main = "res48hControl_lowLFC")
plotMA(res96hControl_lowLFC, main = "res96hControl_lowLFC")


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

# Write out diff-exp gene table
write.table(res48h[names(selected_genes_48h),] %>%
              data.frame(geneID = rownames(.),
                         geneName = selected_genes_48h[rownames(.)]$gene_name,
                         .) %>%
              arrange(padj),
            file = "DESeq2_analysis_no_interaction_results/selected_genes48h.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(res96h[names(selected_genes_96h),] %>%
              data.frame(geneID = rownames(.),
                         geneName = selected_genes_96h[rownames(.)]$gene_name,
                         .) %>%
              arrange(padj),
            file = "DESeq2_analysis_no_interaction_results/selected_genes96h.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

grid.newpage()

# Venn diagram with overlap between 48h and 96h genes
venn.plot <- draw.pairwise.venn(area1 = length(selected_genes_48h),
                                area2 = length(selected_genes_96h),
                                cross.area = length(intersect(names(selected_genes_48h),
                                                       names(selected_genes_96h))),
                              category = c("48h_genes",
                                           "96h_genes"),
                              euler.d = TRUE,
                              col = c("red", "green"))
grid.draw(venn.plot)


# Get plots showing changes in relative gene counts
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
  geom_line(alpha = I(1/5)) + scale_y_log10() + facet_wrap(~ siRNA) +
  ggtitle("selected_genes96h")
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
  geom_line(alpha = I(1/5)) + scale_y_log10() + facet_wrap(~ siRNA) +
  ggtitle("selected_genes48h")
print(p)


# Volcano plots
res96h_frame <- as.data.frame(res96h)
res96h_frame$id <- rownames(res96h)
volcano1 <- ggplot(res96h_frame %>% dplyr::filter(is.na(padj) == FALSE),
                   aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_label_repel(data = res96h_frame %>%
                     mutate(final_set = id %in% names(selected_genes_96h)) %>%
                     arrange(padj) %>%
                     head(10), mapping = aes(label = id, colour = final_set)) +
  labs(title = "design = ~ condition",
       subtitle = "contrast = c(\"condition\", \"TRF2_96\", \"TRF2_0\")")

# Plot volcano for the control non-changing genes at 96h
res96hControl_lowLFC_frame <- as.data.frame(res96hControl_lowLFC)
res96hControl_lowLFC_frame$id <- rownames(res96hControl_lowLFC)
volcanoControl <- ggplot(res96hControl_lowLFC_frame %>% dplyr::filter(is.na(padj) == FALSE),
                         aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_label_repel(data = res96hControl_lowLFC_frame %>%
                     mutate(final_set = id %in% names(selected_genes_96h)) %>%
                     arrange(padj) %>%
                     head(10), mapping = aes(label = id, colour = final_set)) +
  labs(title = "design = ~ condition",
       subtitle = "lfcThreshold = log2(1.25), altHypothesis = \"lessAbs\", contrast = c(\"condition\", \"TRF2_96\", \"TRF2_0\")")

grid.arrange(volcano1 + scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 20)),
             volcanoControl + scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 20)),
             ncol = 2)


dev.off()

# KEGG analysis

library("KEGGprofile")

library("topGO")

library("gplots")

sel_48h_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_48h)$entrezid),
                                             species = "hsa", download_latest = TRUE)

sel_96h_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_96h)$entrezid),
                                             species = "hsa", download_latest = TRUE)

#sel_48hControl_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_48h_control)$entrezid),
#                                             species = "hsa", download_latest = TRUE)
#
#sel_96hControl_KEGG_result <- find_enriched_pathway(unlist(values(selected_genes_96h_control)$entrezid),
#                                             species = "hsa", download_latest = TRUE)

write.table(sel_96h_KEGG_result$stastic %>% arrange(pvalueAdj),
            file = "DESeq2_analysis_no_interaction_results/top_KEGG_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#textplot(sel_96hControl_KEGG_result$stastic)


# GO analysis

# Trying to more explicitly compare genes with similar distribution of expression
overallBaseMean <- as.matrix(res96h[, "baseMean", drop = F])

sig_idx <- match(names(selected_genes_96h), rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

multidensity( list( 
  all= log2(res96h[is.na(res96h$padj) == FALSE ,"baseMean"]) ,
  foreground =log2(res96h[names(selected_genes_96h), "baseMean"]), 
  background =log2(res96h[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(names(selected_genes_96h),  backG) 
inSelection =  geneIDs %in% names(selected_genes_96h) 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tgd.BP <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
            annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.BP <- runTest(tgd.BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.BP <- runTest(tgd.BP, algorithm = "classic", statistic = "Fisher" )

allRes.BP <- GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP, 
                   Fisher.classic = resultTopGO.classic.BP,
                   orderBy = "Fisher.classic", topNodes = 200)

pdf("DESeq2_analysis_no_interaction_results/top_BP_GO_results.pdf", width = 12)
showSigOfNodes(tgd.BP,
               score(resultTopGO.classic.BP),
               firstSigNodes = sum(allRes.BP$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()


write.table(allRes.BP %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_BP_GO_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Quick look up to see if "response to dsRNA" terms are overpresented
GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP,
         Fisher.classic = resultTopGO.classic.BP,
         orderBy = "Fisher.classic", topNodes = 875) %>%
  filter(GO.ID == "GO:0043331")

# Quick look up to see if "cellular response to dsRNA" terms are overpresented
GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP,
         Fisher.classic = resultTopGO.classic.BP,
         orderBy = "Fisher.classic", topNodes = 875) %>%
  filter(GO.ID == "GO:0071359")


# Quick look up to see if "dsRNA processing" terms are overpresented
GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP,
         Fisher.classic = resultTopGO.classic.BP,
         orderBy = "Fisher.classic", topNodes = 875) %>%
  filter(GO.ID == "GO:0031050")

# Same for Molecular Function ontology
tgd.MF <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.MF <- runTest(tgd.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.MF <- runTest(tgd.MF, algorithm = "classic", statistic = "Fisher" )


allRes.MF <- GenTable(tgd.MF, Fisher.elim = resultTopGO.elim.MF, 
                      Fisher.classic = resultTopGO.classic.MF,
                      orderBy = "Fisher.classic", topNodes = 200)

pdf("DESeq2_analysis_no_interaction_results/top_MF_GO_results.pdf", width = 12)
showSigOfNodes(tgd.MF,
               score(resultTopGO.classic.MF),
               firstSigNodes = sum(allRes.MF$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()

write.table(allRes.MF %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_MF_GO_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# For upregulated genes only
overallBaseMean <- as.matrix(res96h[, "baseMean", drop = F])

sig_idx <- match(names(selected_genes_96h)[res96h[names(selected_genes_96h), "log2FoldChange"] > 0], rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

multidensity( list( 
  all= log2(res96h[is.na(res96h$padj) == FALSE ,"baseMean"]) ,
  foreground =log2(res96h[names(selected_genes_96h), "baseMean"]), 
  background =log2(res96h[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(names(selected_genes_96h),  backG) 
inSelection =  geneIDs %in% names(selected_genes_96h) 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tgd.BP <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.BP <- runTest(tgd.BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.BP <- runTest(tgd.BP, algorithm = "classic", statistic = "Fisher" )


allRes.BP <- GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP, 
                      Fisher.classic = resultTopGO.classic.BP,
                      orderBy = "Fisher.classic", topNodes = 200)

pdf("DESeq2_analysis_no_interaction_results/top_BP_GO_UP_results.pdf", width = 12)
showSigOfNodes(tgd.BP,
               score(resultTopGO.classic.BP),
               firstSigNodes = sum(allRes.BP$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()

write.table(allRes.BP %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_BP_GO_UP_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Same for Molecular Function ontology
tgd.MF <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.MF <- runTest(tgd.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.MF <- runTest(tgd.MF, algorithm = "classic", statistic = "Fisher" )


allRes.MF <- GenTable(tgd.MF, Fisher.elim = resultTopGO.elim.MF, 
                      Fisher.classic = resultTopGO.classic.MF,
                      orderBy = "Fisher.classic", topNodes = 200)

pdf("DESeq2_analysis_no_interaction_results/top_MF_GO_UP_results.pdf", width = 12)

showSigOfNodes(tgd.MF,
               score(resultTopGO.classic.MF),
               firstSigNodes = sum(allRes.MF$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()

write.table(allRes.MF %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_MF_GO_UP_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# For downregulated genes only
overallBaseMean <- as.matrix(res96h[, "baseMean", drop = F])

sig_idx <- match(names(selected_genes_96h)[res96h[names(selected_genes_96h), "log2FoldChange"] < 0], rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

multidensity( list( 
  all= log2(res96h[is.na(res96h$padj) == FALSE ,"baseMean"]) ,
  foreground =log2(res96h[names(selected_genes_96h), "baseMean"]), 
  background =log2(res96h[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(names(selected_genes_96h),  backG) 
inSelection =  geneIDs %in% names(selected_genes_96h) 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tgd.BP <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.BP <- runTest(tgd.BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.BP <- runTest(tgd.BP, algorithm = "classic", statistic = "Fisher" )

allRes.BP <- GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP, 
                      Fisher.classic = resultTopGO.classic.BP,
                      orderBy = "Fisher.classic", topNodes = 200)

write.table(allRes.BP %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_BP_GO_DOWN_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

pdf("DESeq2_analysis_no_interaction_results/top_BP_GO_DOWN_results.pdf", width = 12)
showSigOfNodes(tgd.BP,
               score(resultTopGO.classic.BP),
               firstSigNodes = sum(allRes.BP$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()

# Same for Molecular Function ontology
tgd.MF <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.MF <- runTest(tgd.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.MF <- runTest(tgd.MF, algorithm = "classic", statistic = "Fisher" )


allRes.MF <- GenTable(tgd.MF, Fisher.elim = resultTopGO.elim.MF, 
                      Fisher.classic = resultTopGO.classic.MF,
                      orderBy = "Fisher.classic", topNodes = 200)

write.table(allRes.MF %>% filter(Fisher.classic < 0.01),
            file = "DESeq2_analysis_no_interaction_results/top_MF_GO_DOWN_results.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

pdf("DESeq2_analysis_no_interaction_results/top_MF_GO_DOWN_results.pdf", width = 12)
showSigOfNodes(tgd.MF,
               score(resultTopGO.classic.MF),
               firstSigNodes = sum(allRes.MF$Fisher.classic < 0.01),
               useInfo = 'all')

dev.off()



# Check if the genes in Mukherjee et al show up in our differential expression analysis

Mukherjee_genes <- read_tsv("genes_from_Mukherjee_etal.txt")

sum(Mukherjee_genes$`Gene stable ID` %in% names(selected_genes_96h))
Mukherjee_genes[Mukherjee_genes$`Gene stable ID` %in% names(selected_genes_96h),]

res96h[Mukherjee_genes[Mukherjee_genes$`Gene stable ID` %in% rownames(res96h),]$`Gene stable ID`,]

# The Mukherjee TRF2-bound genes are not differentially-expressed in our dataset (with the exception of CDKN1A)
# In the paper they see the changes in HT1080 and MRC5 cells. Could it be because this is HeLa?
# Both cell types are ALT(-), like HeLa.

