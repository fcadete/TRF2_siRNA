
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

pdf("DESeq2_analysis.pdf")

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

p <- ggplot(TRF2_counts, aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))
print(p)


TRF1_counts <- plotCounts(ddsTxi, "ENSG00000147601",
                          intgroup = c("siRNA", "TimePoint"),
                          returnData = TRUE)
p <- ggplot(TRF1_counts, aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))
print(p)

res <- results(ddsTxi, name="siRNATRF2.TimePoint96")

p <- plotCounts(ddsTxi, which.min(res$padj),
           intgroup = c("siRNA", "TimePoint"),
           returnData = TRUE) %>%
  ggplot(aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))
print(p)

p <- plotCounts(ddsTxi, "ENSG00000101680",
           intgroup = c("siRNA", "TimePoint"),
           returnData = TRUE) %>%
  ggplot(aes(x= as.numeric(as.character(TimePoint)), y = count, colour = siRNA)) +
  geom_point() + geom_smooth() + scale_x_continuous("Time Point (hrs)") + scale_y_continuous(limits = c(0, NA))
print(p)

# Number of diff-expressed genes between TRF2 and control shRNA
with(results(ddsTxi, name = "siRNA_TRF2_vs_control"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNA_TRF2_vs_control"), main = "siRNA_TRF2_vs_control")

# Number of diff-expressed genes in the siRNATRF2.TimePoint48 interaction
with(results(ddsTxi, name = "siRNATRF2.TimePoint48"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNATRF2.TimePoint48"), main = "siRNATRF2.TimePoint48")

#Names of the genes
siRNATRF2.TimePoint48.results <- results(ddsTxi, name = "siRNATRF2.TimePoint48")
unique(GenomicFeatures::genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_48h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))
selected_genes_48h$chr_length <- seqlengths(selected_genes_48h)[as.character(seqnames(selected_genes_48h))]
selected_genes_48h <- selected_genes_48h %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))
p <- ggplot(selected_genes_48h, aes(y=distance_to_nearest_end, x = "1")) + geom_jitter()
print(p)

# Number of diff-expressed genes in the siRNATRF2.TimePoint96 interaction
with(results(ddsTxi, name = "siRNATRF2.TimePoint96"), sum(padj < 0.1, na.rm=TRUE))
plotMA(results(ddsTxi, name = "siRNATRF2.TimePoint96"), main = "siRNATRF2.TimePoint96")

siRNATRF2.TimePoint96.results <- results(ddsTxi, name = "siRNATRF2.TimePoint96")
unique(GenomicFeatures::genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_96h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))
selected_genes_96h$chr_length <- seqlengths(selected_genes_96h)[as.character(seqnames(selected_genes_96h))]
selected_genes_96h <- selected_genes_96h %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))
p <- ggplot(selected_genes_96h, aes(y=distance_to_nearest_end, x = "1")) + geom_jitter() + scale_y_log10()
print(p)

all_genes <- GenomicFeatures::genes(EnsDb.Hsapiens.v86)
all_genes$chr_length <- seqlengths(all_genes)[as.character(seqnames(all_genes))]
all_genes <- all_genes %>% as.data.frame() %>% group_by(gene_id) %>% mutate(distance_to_nearest_end = min(c(start, chr_length - end)))

p <- ggplot(data = rbind(data.frame(all_genes, type = "all genes"),
                    data.frame(selected_genes_48h, type = "48h"),
                    data.frame(selected_genes_96h, type = "96h")),
       mapping = aes(y = distance_to_nearest_end, x = type)) +
  geom_boxplot()
print(p)

p <- ggplot() +
  geom_boxplot(data = all_genes, mapping = aes(y=distance_to_nearest_end, x = "1")) +
  geom_jitter(data = selected_genes_48h, mapping = aes(y=distance_to_nearest_end, x = "1"), colour = "red") +
  geom_jitter(data = selected_genes_96h, mapping = aes(y=distance_to_nearest_end, x = "1"), colour = "blue") +
  scale_y_log10()
print(p)



selected_genes_48h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))

selected_genes_96h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))

gene_counts <- data.frame()


for (gene in names(selected_genes_96h)) {

       gene_id <- selected_genes_96h[gene]$gene_id
       symbol <- selected_genes_96h[gene]$symbol
       gene_counts <- rbind(gene_counts, cbind(plotCounts(ddsTxi, gene_id,
                       intgroup = c("siRNA", "TimePoint"),
                       returnData = TRUE), gene_id = gene_id, symbol = symbol))

   }



p <- ggplot(gene_counts,
            aes(x = as.numeric(as.character(TimePoint)), y = count,
                group = interaction(gene_id, siRNA), color = siRNA)) +
       geom_smooth(se = FALSE) + scale_y_log10()

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

p <- gene_counts %>%
   separate(sample, into = c("shRNA", "timepoint", "replicate")) %>%
   select(-one_of(c("shRNA", "timepoint"))) %>%
   spread("TimePoint", "count") %>%
   group_by(gene_id, siRNA) %>%
   summarise(`relative 0` = 1,
             `relative 48` = mean(`48`) / mean(`0`),
             `relative 96` = mean(`96`) / mean(`0`)) %>%
   gather(timePoint, relative_counts, `relative 0`, `relative 48`, `relative 96`) %>%
   spread(siRNA, relative_counts) %>%
   group_by(gene_id, timePoint) %>%
   summarise(fold_change = log2(TRF2 / control)) %>%
   ggplot(aes(x = timePoint, y = fold_change,
              group = gene_id)) +
   geom_line(alpha = I(1/5))
print(p)


gene_counts <- data.frame()


for (gene in names(selected_genes_48h)) {

       gene_id <- selected_genes_48h[gene]$gene_id
       symbol <- selected_genes_48h[gene]$symbol
       gene_counts <- rbind(gene_counts, cbind(plotCounts(ddsTxi, gene_id,
                       intgroup = c("siRNA", "TimePoint"),
                       returnData = TRUE), gene_id = gene_id, symbol = symbol))

   }



p <- ggplot(gene_counts,
            aes(x = as.numeric(as.character(TimePoint)), y = count,
                group = interaction(gene_id, siRNA), color = siRNA)) +
       geom_smooth(se = FALSE) + scale_y_log10()

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

p <- gene_counts %>%
   separate(sample, into = c("shRNA", "timepoint", "replicate")) %>%
   select(-one_of(c("shRNA", "timepoint"))) %>%
   spread("TimePoint", "count") %>%
   group_by(gene_id, siRNA) %>%
   summarise(`relative 0` = 1,
             `relative 48` = mean(`48`) / mean(`0`),
             `relative 96` = mean(`96`) / mean(`0`)) %>%
   gather(timePoint, relative_counts, `relative 0`, `relative 48`, `relative 96`) %>%
   spread(siRNA, relative_counts) %>%
   group_by(gene_id, timePoint) %>%
   summarise(fold_change = log2(TRF2 / control)) %>%
   ggplot(aes(x = timePoint, y = fold_change,
              group = gene_id)) +
   geom_line(alpha = I(1/5))
print(p)




dev.off()


