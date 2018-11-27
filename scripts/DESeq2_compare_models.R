
library(DESeq2)
library("tximport")
library("EnsDb.Hsapiens.v86")
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(VennDiagram)

pdf("DESeq2_compared_models.R", width = 12)

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


selected_genes_48h <- rownames(res48h[which(res48h$padj < 0.1),])

selected_genes_96h <- rownames(res96h[which(res96h$padj < 0.1),])

selected_genes_48h_control_low_LFC <- rownames(res48h[which(res48hControl_lowLFC$padj < 0.1),])

selected_genes_96h_control_low_LFC <- rownames(res96h[which(res96hControl_lowLFC$padj < 0.1),])

# Include only the genes that also have a low LFC in the control samples
selected_genes_48h <- selected_genes_48h[selected_genes_48h %in% selected_genes_48h_control_low_LFC]
selected_genes_96h <- selected_genes_96h[selected_genes_96h %in% selected_genes_96h_control_low_LFC]


# Plot Volcano for 96h time-point against 0h.
res96h_frame <- as.data.frame(res96h)
res96h_frame$id <- rownames(res96h)
volcano1 <- ggplot(res96h_frame %>% dplyr::filter(is.na(padj) == FALSE),
       aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_label_repel(data = res96h_frame %>%
                     mutate(final_set = id %in% selected_genes_96h) %>%
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
                     mutate(final_set = id %in% selected_genes_96h) %>%
                     arrange(padj) %>%
                     head(10), mapping = aes(label = id, colour = final_set)) +
  labs(title = "design = ~ condition",
       subtitle = "lfcThreshold = log2(1.25), altHypothesis = \"lessAbs\", contrast = c(\"condition\", \"TRF2_96\", \"TRF2_0\")")



# Do the same now but not including the control samples (so the dispersion estimates are different)
samples <- read_tsv("sample_info.txt")
samples$condition <- paste(samples$siRNA, samples$TimePoint, sep = "_")

samples <- samples %>% dplyr::filter(siRNA == "TRF2")

samples$condition <- relevel(factor(samples$condition), ref = "TRF2_0")

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxiNoControl <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

ddsTxiNoControl <- DESeq(ddsTxiNoControl)

# Test for genes differentially expressed between timePoints in the TRF2 shRNA data
res48hNoControl <- results(ddsTxiNoControl, contrast = c("condition", "TRF2_48", "TRF2_0"))
res96hNoControl <- results(ddsTxiNoControl, contrast = c("condition", "TRF2_96", "TRF2_0"))

# Do the comparison JUST with the control samples
samples <- read_tsv("sample_info.txt")
samples$condition <- paste(samples$siRNA, samples$TimePoint, sep = "_")

samples <- samples %>% dplyr::filter(siRNA == "control")

samples$condition <- relevel(factor(samples$condition), ref = "control_0")

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxiJustControl <- DESeqDataSetFromTximport(txi,
                                              colData = samples,
                                              design = ~ condition)

ddsTxiJustControl <- DESeq(ddsTxiJustControl)

# Test for genes differentially expressed between timePoints in the TRF2 shRNA data
res48hJustControl <- results(ddsTxiJustControl,
                             lfcThreshold = log2(1.25),
                             altHypothesis = "lessAbs",
                             contrast = c("condition", "control_48", "control_0"))
res96hJustControl <- results(ddsTxiJustControl,
                             lfcThreshold = log2(1.25),
                             altHypothesis = "lessAbs",
                             contrast = c("condition", "control_96", "control_0"))




selected_genes_48h_no_control <- rownames(res48hNoControl[which(res48hNoControl$padj < 0.1),])

selected_genes_96h_no_control <- rownames(res96hNoControl[which(res96hNoControl$padj < 0.1),])

selected_genes_48h_just_control_low_LFC <- rownames(res48hNoControl[which(res48hJustControl$padj < 0.1),])

selected_genes_96h_just_control_low_LFC <- rownames(res96hNoControl[which(res96hJustControl$padj < 0.1),])

# Include only the genes that also have a low LFC in the control samples
selected_genes_48h_no_control <- selected_genes_48h_no_control[selected_genes_48h_no_control %in% selected_genes_48h_just_control_low_LFC]
selected_genes_96h_no_control <- selected_genes_96h_no_control[selected_genes_96h_no_control %in% selected_genes_96h_just_control_low_LFC]


# Plot Volcano for 96h time-point against 0h.
res96hNoControl_frame <- as.data.frame(res96hNoControl)
res96hNoControl_frame$id <- rownames(res96hNoControl)
volcano2 <- ggplot(res96hNoControl_frame %>% dplyr::filter(is.na(padj) == FALSE),
       aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_label_repel(data = res96hNoControl_frame %>%
                     mutate(final_set = id %in% selected_genes_96h) %>%
                     arrange(padj) %>%
                     head(10), mapping = aes(label = id, colour = final_set)) +
  labs(title = "design = ~ condition, no control samples",
       subtitle = "contrast = c(\"condition\", \"TRF2_96\", \"TRF2_0\")")


grid.arrange(volcano1 + scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 20)),
             volcano2 + scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 20)),
             ncol = 2)


qqplot(res96h$stat, res96hNoControl$stat)
abline(a = 0, b = 1)


# Get the results for the interaction now

samples <- read_tsv("sample_info.txt")
samples$condition <- paste(samples$siRNA, samples$TimePoint, sep = "_")
samples$condition <- relevel(factor(samples$condition), ref = "TRF2_0")
samples$TimePoint <- as.character(samples$TimePoint)

files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxiInteraction <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ siRNA + TimePoint + siRNA:TimePoint)

ddsTxiInteraction <- DESeq(ddsTxiInteraction)

res96hInteraction <- results(ddsTxiInteraction, name = c("siRNATRF2.TimePoint96"))

selected_genes_96h_interaction <- rownames(res96hInteraction[which(res96hInteraction$padj < 0.1),])

selected_genes_96h_interaction <- selected_genes_96h_interaction[selected_genes_96h_interaction %in% selected_genes_96h_control_low_LFC]


res96hInteraction_frame <- as.data.frame(res96hInteraction)
res96hInteraction_frame$id <- rownames(res96hInteraction)
volcano3 <- ggplot(res96hInteraction_frame %>% dplyr::filter(is.na(padj) == FALSE),
                   aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_label_repel(data = res96hInteraction_frame %>%
                     mutate(final_set = id %in% selected_genes_96h) %>%
                     arrange(padj) %>%
                     head(10), mapping = aes(label = id, colour = final_set)) +
  labs(title = "design = ~ siRNA + TimePoint + siRNA:TimePoint",
       subtitle = "name = c(\"siRNATRF2.TimePoint96\")")


grid.arrange(volcano1 + scale_x_continuous(limits = c(-2.5, 2.5)),
             volcano2 + scale_x_continuous(limits = c(-2.5, 2.5)),
             volcano3 + scale_x_continuous(limits = c(-2.5, 2.5)),
             volcanoControl + scale_x_continuous(limits = c(-2.5, 2.5)),
             ncol = 4)

qqplot(res96h$stat, res96hInteraction$stat)
abline(a = 0, b = 1)


venn.plot <- draw.triple.venn(area1 = length(selected_genes_96h),
                              area2 = length(selected_genes_96h_no_control),
                              area3 = length(selected_genes_96h_interaction),
                              n12 = length(intersect(selected_genes_96h,
                                                     selected_genes_96h_no_control)),
                              n23 = length(intersect(selected_genes_96h_no_control,
                                                     selected_genes_96h_interaction)),
                              n13 = length(intersect(selected_genes_96h,
                                                     selected_genes_96h_interaction)),
                              n123 = length(Reduce(intersect,
                                                   list(selected_genes_96h,
                                                        selected_genes_96h_no_control,
                                                        selected_genes_96h_interaction))),
                              category = c("Full_samples",
                                           "Restricted_samples",
                                           "Interaction"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

grid.newpage()
grid.draw(venn.plot)

dev.off()
