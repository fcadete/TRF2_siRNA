
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")
library("ggrepel")
library("gridExtra")

library("genefilter")
library("geneplotter")
library("VennDiagram")

library("Rsamtools")

#setwd("/mnt/TRF2_siRNA/")

pdf("DESeq2_analysis_no_interaction_results/TTAGGG_repeat_analysis.pdf", width = 10)


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

ddsTxi <- DESeq(ddsTxi)

# Test for genes differentially expressed between timePoints in the TRF2 shRNA data
res48h <- results(ddsTxi, contrast = c("condition", "TRF2_48", "TRF2_0"))
res96h <- results(ddsTxi, contrast = c("condition", "TRF2_96", "TRF2_0"))

# Explicitly test for genes with a low LFC in the control timecourse
res48hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_48"))

res96hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_96"))



selected_genes_48h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res48h[which(res48h$padj < 0.1),])))

selected_genes_96h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res96h[which(res96h$padj < 0.1),])))

selected_genes_48h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res48h[which(res48hControl_lowLFC$padj < 0.1),])))

selected_genes_96h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res96h[which(res96hControl_lowLFC$padj < 0.1),])))

# Include only the genes that also have a low LFC in the control samples
selected_genes_48h <- selected_genes_48h[(names(selected_genes_48h) %in% names(selected_genes_48h_control_low_LFC))]
selected_genes_96h <- selected_genes_96h[(names(selected_genes_96h) %in% names(selected_genes_96h_control_low_LFC))]

Dna <- getGenomeFaFile(EnsDb.Hsapiens.v86)

genes_in_data <- GenomicFeatures::genes(EnsDb.Hsapiens.v86, filter = GeneIdFilter(rownames(res96h)))

genes_in_data <- genes_in_data[seqnames(genes_in_data) %in% seqnames(seqinfo(Dna))]

geneSeqs <- getSeq(Dna, genes_in_data)

gene_telrep_matches <- vmatchPattern("TTAGGG", geneSeqs)

gene_telrep_number_matches <- unlist(lapply(gene_telrep_matches, length))

names(gene_telrep_number_matches) <- names(genes_in_data)


gene_telreps <- data_frame(gene_ID = names(genes_in_data),
                           telrep_matches = gene_telrep_number_matches,
                           gene_width = width(genes_in_data))


ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches)) + geom_point()

ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches)) + geom_point() + scale_x_log10()

ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches)) + geom_point() + scale_x_log10() + scale_y_log10()


# so the number of TTAGG repeats is very correlated with the length of the gene, which makes absolute sense

# We want the genes that deviate from this! There are some outliers in the scatterplot.

# First let's try to normalise just by getting the repeat to length ratio

ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches / gene_width)) + geom_point()
ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches / gene_width)) + geom_point() + scale_y_log10()
ggplot(gene_telreps, aes(x = gene_width, y = telrep_matches / gene_width)) + geom_point() + scale_x_log10() + scale_y_log10()


# Let's try to use linear modelling

mod <- lm(telrep_matches ~ gene_width, data = gene_telreps)

gene_telreps <- modelr::add_predictions(gene_telreps, mod) %>%
                modelr::add_residuals(mod)


ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = telrep_matches)) + geom_line(aes(y=pred), colour = "red", size = 1)

ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = resid))

gene_telreps$res96h <- gene_telreps$gene_ID %in% names(selected_genes_96h)


ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = telrep_matches, colour = res96h)) + geom_line(aes(y=pred), colour = "red", size = 1)

ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = resid, colour = res96h))

# And now let's try a Poisson GLM

poisson_mod <- glm(telrep_matches ~ gene_width, poisson, data = gene_telreps)

gene_telreps <- modelr::add_predictions(gene_telreps, poisson_mod) %>%
                modelr::add_residuals(poisson_mod)


ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = telrep_matches)) + geom_line(aes(y=pred), colour = "red", size = 1)

ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = resid))

gene_telreps$res96h <- gene_telreps$gene_ID %in% names(selected_genes_96h)


ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = telrep_matches, colour = res96h)) + geom_line(aes(y=pred), colour = "red", size = 1)

ggplot(gene_telreps, aes(x = gene_width)) + geom_point(aes( y = resid, colour = res96h))


dev.off()

