source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

library("tximport")
biocLite("EnsDb.Hsapiens.v86")
library("EnsDb.Hsapiens.v86")

install.packages("tidyverse")
library("tidyverse")

samples <- read_tsv("sample_info.txt")
samples$TimePoint <- as.character(samples$TimePoint)

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ siRNA + TimePoint + siRNA:TimePoint)

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
res <- results(ddsTxi,
               contrast = c("siRNA", "control", "TRF2"))

TRF2_counts <- plotCounts(ddsTxi, "ENSG00000132604",
                          intgroup = c("siRNA", "TimePoint"),
                          returnData = TRUE)

ggplot(TRF2_counts, aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))



TRF1_counts <- plotCounts(ddsTxi, "ENSG00000147601",
                          intgroup = c("siRNA", "TimePoint"),
                          returnData = TRUE)
ggplot(TRF1_counts, aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))

res <- results(ddsTxi, name="siRNATRF2.TimePoint96")
plotCounts(ddsTxi, which.min(res$padj),
           intgroup = c("siRNA", "TimePoint"),
           returnData = TRUE) %>%
  ggplot(aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))

plotCounts(ddsTxi, "ENSG00000101680",
           intgroup = c("siRNA", "TimePoint"),
           returnData = TRUE) %>%
  ggplot(aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))


# Number of diff-expressed genes between TRF2 and control shRNA
with(results(ddsTxi, name = "siRNA_TRF2_vs_control"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNA_TRF2_vs_control"), main = "siRNA_TRF2_vs_control")

# Number of diff-expressed genes in the siRNATRF2.TimePoint48 interaction
with(results(ddsTxi, name = "siRNATRF2.TimePoint48"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNATRF2.TimePoint48"), main = "siRNATRF2.TimePoint48")

#Names of the genes
siRNATRF2.TimePoint48.results <- results(ddsTxi, name = "siRNATRF2.TimePoint48")
unique(genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_48h <- genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))
selected_genes_48h$chr_length <- seqlengths(selected_genes_48h)[as.character(seqnames(selected_genes_48h))]
selected_genes_48h <- selected_genes_48h %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))
ggplot(selected_genes_48h, aes(y=distance_to_nearest_end, x = "1")) + geom_jitter()

# Number of diff-expressed genes in the siRNATRF2.TimePoint96 interaction
with(results(ddsTxi, name = "siRNATRF2.TimePoint96"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNATRF2.TimePoint96"), main = "siRNATRF2.TimePoint96")

siRNATRF2.TimePoint96.results <- results(ddsTxi, name = "siRNATRF2.TimePoint96")
unique(genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_96h <- genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))
selected_genes_96h$chr_length <- seqlengths(selected_genes_96h)[as.character(seqnames(selected_genes_96h))]
selected_genes_96h <- selected_genes_96h %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))
ggplot(selected_genes_96h, aes(y=distance_to_nearest_end, x = "1")) + geom_jitter() + scale_y_log10()

all_genes <- genes(EnsDb.Hsapiens.v86)
all_genes$chr_length <- seqlengths(all_genes)[as.character(seqnames(all_genes))]
all_genes <- all_genes %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))

ggplot(data = rbind(data.frame(all_genes, type = "all genes"),
                    data.frame(selected_genes_48h, type = "48h"),
                    data.frame(selected_genes_96h, type = "96h")),
       mapping = aes(y = distance_to_nearest_end, x = type)) +
  geom_boxplot()

ggplot() +
  geom_boxplot(data = all_genes, mapping = aes(y=distance_to_nearest_end, x = "1")) +
  geom_jitter(data = selected_genes_48h, mapping = aes(y=distance_to_nearest_end, x = "1"), colour = "red") +
  geom_jitter(data = selected_genes_96h, mapping = aes(y=distance_to_nearest_end, x = "1"), colour = "blue") +
  scale_y_log10()
