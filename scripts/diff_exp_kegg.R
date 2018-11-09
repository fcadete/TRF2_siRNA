
library("DESeq2")
library("Repitools")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

library("KEGGprofile")

library("topGO")

library("gplots")

pdf("diff_expressed_kegg_analysis.pdf")

samples <- read_tsv("sample_info.txt")
samples$TimePoint <- as.character(samples$TimePoint)

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ siRNA + TimePoint + siRNA:TimePoint)


ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi,
               contrast = c("siRNA", "control", "TRF2"))


#Names of the genes
siRNATRF2.TimePoint48.results <- results(ddsTxi, name = "siRNATRF2.TimePoint48")
unique(ensembldb::genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_48h <- ensembldb::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))
seqlevels(selected_genes_48h) <- paste0("chr", seqlevels(selected_genes_48h))

# Number of diff-expressed genes in the siRNATRF2.TimePoint96 interaction
siRNATRF2.TimePoint96.results <- results(ddsTxi, name = "siRNATRF2.TimePoint96")
unique(ensembldb::genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_96h <- ensembldb::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))
seqlevels(selected_genes_96h) <- paste0("chr", seqlevels(selected_genes_96h))

all_genes <- ensembldb::genes(EnsDb.Hsapiens.v86)
seqlevels(all_genes) <- paste0("chr", seqlevels(all_genes))

expressed_genes <- ensembldb::genes(EnsDb.Hsapiens.v86,
                         filter = GeneIdFilter(rownames( siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$baseMean > 20),])))
seqlevels(expressed_genes) <- paste0("chr", seqlevels(expressed_genes))



sel_96h_KEGG_result <- find_enriched_pathway(values(selected_genes_96h)$entrezid,
                                             species = "hsa", download_latest = TRUE)

textplot(sel_96h_KEGG_result$stastic)

gene_counts <- data.frame()

for (pathway in rownames(sel_96h_KEGG_result$stastic)) {

    for (gene in sel_96h_KEGG_result$detail[[pathway]]) {

       gene_id <- selected_genes_96h[values(selected_genes_96h)$entrezid == gene]$gene_id
       symbol <- selected_genes_96h[values(selected_genes_96h)$entrezid == gene]$symbol
       p <- plotCounts(ddsTxi, gene_id,
                       intgroup = c("siRNA", "TimePoint"),
                       returnData = TRUE) %>%
              ggplot(aes(x = as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
               geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA)) +
               scale_colour_discrete("shRNA") + ggtitle(symbol)
       print(p)

       gene_counts <- rbind(gene_counts, cbind(plotCounts(ddsTxi, gene_id,
                       intgroup = c("siRNA", "TimePoint"),
                       returnData = TRUE), gene_id = gene_id, symbol = symbol, pathway = pathway))

   }
}

for (this_pathway in unique(gene_counts$pathway)) {
   p <- ggplot(filter(gene_counts, pathway == this_pathway),
               aes(x = as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
     geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA)) +
     scale_colour_discrete("shRNA") + facet_wrap(~ symbol, scales = "free_y") +
     ggtitle(sel_96h_KEGG_result$stastic[this_pathway, "Pathway_Name"])
print(p)

}



padj_96h <- siRNATRF2.TimePoint96.results$padj
padj_96h[is.na(padj_96h)] <- 1
names(padj_96h) <- rownames(siRNATRF2.TimePoint96.results)

sampleGOdata <- new("topGOdata", description = "teste", ontology = "BP",
                    allGenes = padj_96h, geneSel = function(x) { (x < 0.1) },
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
                    allGenes = padj_96h, geneSel = function(x) { (x < 0.1) },
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


sampleGOdata <- new("topGOdata", description = "teste", ontology = "CC",
                    allGenes = padj_96h, geneSel = function(x) { (x < 0.1) },
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

