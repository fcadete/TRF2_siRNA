
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

setwd("/mnt/TRF2_siRNA/")

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

gene_telreps <- data_frame(gene_ID = names(genes_in_data),
                           gene_width = width(genes_in_data))

gene_telreps$res96h <- gene_telreps$gene_ID %in% names(selected_genes_96h)


geneSeqs <- getSeq(Dna, genes_in_data)

for (number_reps in 1:5) {

  gene_telrep_matches <- vmatchPattern(paste(rep("TTAGGG", number_reps), collapse = ""), geneSeqs)

  gene_telrep_number_matches <- unlist(lapply(gene_telrep_matches, length))

  gene_telreps[[paste0("telrep_matches_", number_reps)]] <- gene_telrep_number_matches

}


ggplot(gene_telreps %>%
         gather(key = "number_reps", value = "telrep_matches",
                telrep_matches_1:telrep_matches_5),
       aes(x = gene_width, y = telrep_matches)) +
  geom_point() +
  facet_wrap(~ number_reps, scales = "free_y")

ggplot(gene_telreps %>%
         gather(key = "number_reps", value = "telrep_matches",
                telrep_matches_1:telrep_matches_5),
       aes(x = gene_width, y = telrep_matches)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ number_reps, scales = "free_y")

ggplot(gene_telreps %>%
         gather(key = "number_reps", value = "telrep_matches",
                telrep_matches_1:telrep_matches_5),
       aes(x = gene_width, y = telrep_matches)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ number_reps, scales = "free_y")


# so the number of TTAGG repeats is very correlated with the length of the gene, which makes absolute sense
# At least for the one tandem repeat

# We want the genes that deviate from this! There are some outliers in the scatterplot.

# First let's try to normalise just by getting the repeat to length ratio

ggplot(gene_telreps %>%
         gather(key = "number_reps", value = "telrep_matches",
                telrep_matches_1:telrep_matches_5),
       aes(x = gene_width, y = telrep_matches / gene_width)) +
  geom_point() +
  facet_wrap(~ number_reps, scales = "free_y")

ggplot(gene_telreps %>%
         gather(key = "number_reps", value = "telrep_matches",
                telrep_matches_1:telrep_matches_5),
       aes(x = gene_width, y = telrep_matches / gene_width)) +
  geom_point() +
  scale_y_log10() +
  facet_wrap(~ number_reps, scales = "free_y")

ggplot(mapping = aes(x = gene_width, y = telrep_matches / gene_width, colour = res96h)) +
  geom_point(data = gene_telreps %>%
               gather(key = "number_reps", value = "telrep_matches",
                      telrep_matches_1:telrep_matches_5) %>% filter(res96h == FALSE)) +
  geom_point(data = gene_telreps %>%
               gather(key = "number_reps", value = "telrep_matches",
                      telrep_matches_1:telrep_matches_5) %>% filter(res96h == TRUE)) +
  scale_y_log10() +
  facet_wrap(~ number_reps, scales = "free_y")

ggplot(mapping = aes(x = gene_width, y = telrep_matches / gene_width, colour = res96h)) +
  geom_point(data = gene_telreps %>%
               gather(key = "number_reps", value = "telrep_matches",
                      telrep_matches_1:telrep_matches_5) %>% filter(res96h == FALSE)) +
  geom_point(data = gene_telreps %>%
               gather(key = "number_reps", value = "telrep_matches",
                      telrep_matches_1:telrep_matches_5) %>% filter(res96h == TRUE)) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ number_reps, scales = "free_y")


# Let's try to use linear modelling for the one repeat (the other don't show such a linear correlation)

mod <- lm(telrep_matches_1 ~ gene_width, data = gene_telreps)

gene_telreps <- modelr::add_predictions(gene_telreps, mod) %>%
                modelr::add_residuals(mod)


ggplot(gene_telreps, aes(x = gene_width)) +
  geom_point(aes(y = telrep_matches_1)) +
  geom_line(aes(y=pred), colour = "red", size = 1)

ggplot(gene_telreps, aes(x = gene_width)) +
  geom_point(aes(y = resid))

ggplot(mapping = aes(x = gene_width)) +
  geom_point(data = filter(gene_telreps, res96h == FALSE), mapping = aes(y = telrep_matches_1, colour = res96h)) +
  geom_point(data = filter(gene_telreps, res96h == TRUE), mapping = aes(y = telrep_matches_1, colour = res96h)) +
  geom_line(data = gene_telreps, mapping = aes(y=pred), colour = "red", size = 1)

ggplot(mapping = aes(x = gene_width)) +
  geom_point(data = filter(gene_telreps, res96h == FALSE), mapping = aes(y = resid, colour = res96h)) +
  geom_point(data = filter(gene_telreps, res96h == TRUE), mapping = aes(y = resid, colour = res96h))

ggplot(mapping = aes(x = gene_width)) +
  geom_point(data = filter(gene_telreps, res96h == FALSE), mapping = aes(y = resid, colour = res96h)) +
  geom_point(data = filter(gene_telreps, res96h == TRUE), mapping = aes(y = resid, colour = res96h)) +
  scale_x_log10()


# Now do this just for promoter sequences.

promoterSeqs <- getSeq(Dna, promoters(genes_in_data[seqnames(genes_in_data) != "MT"],
                                      downstream = 1000))

promoter_telreps <- data_frame(gene_ID = names(genes_in_data[seqnames(genes_in_data) != "MT"]),
                           gene_width = width(genes_in_data[seqnames(genes_in_data) != "MT"]))

promoter_telreps$res96h <- promoter_telreps$gene_ID %in% names(selected_genes_96h)

for (number_reps in 1:5) {
  
  promoter_telrep_matches <- vmatchPattern(paste(rep("TTAGGG", number_reps), collapse = ""), promoterSeqs)
  
  promoter_telrep_number_matches <- unlist(lapply(promoter_telrep_matches, length))
  
  promoter_telreps[[paste0("promoter_telrep_matches_", number_reps)]] <- promoter_telrep_number_matches
  
}

ggplot(promoter_telreps %>%
         gather(key = "number_reps", value = "promoter_telrep_matches",
                promoter_telrep_matches_1:promoter_telrep_matches_5),
       aes(x = gene_width, y = promoter_telrep_matches)) +
  geom_point() +
  facet_wrap(~ number_reps, scales = "free_y")


ggplot(promoter_telreps %>%
         gather(key = "number_reps", value = "promoter_telrep_matches",
                promoter_telrep_matches_1:promoter_telrep_matches_5),
       aes(x = gene_width, y = promoter_telrep_matches)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ number_reps, scales = "free_y")

# No seeming linear correlation in promoters
# Seems appropriate, as the size of the promoter is of an arbitrary size


ggplot(mapping = aes(x = gene_width, y = promoter_telrep_matches, colour = res96h)) +
  geom_point(data = promoter_telreps %>%
               gather(key = "number_reps", value = "promoter_telrep_matches",
                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(res96h == FALSE)) +
  geom_point(data = promoter_telreps %>%
               gather(key = "number_reps", value = "promoter_telrep_matches",
                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(res96h == TRUE)) +
  facet_wrap(~ number_reps, scales = "free_y")


ggplot(mapping = aes(x = gene_width, y = promoter_telrep_matches, colour = res96h)) +
  geom_point(data = promoter_telreps %>%
               gather(key = "number_reps", value = "promoter_telrep_matches",
                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(res96h == FALSE)) +
  geom_point(data = promoter_telreps %>%
               gather(key = "number_reps", value = "promoter_telrep_matches",
                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(res96h == TRUE)) +
  scale_x_log10() +
  facet_wrap(~ number_reps, scales = "free_y")

# No enrichment at all of diff-expressed genes among those with more TTAGGG matches

dev.off()

