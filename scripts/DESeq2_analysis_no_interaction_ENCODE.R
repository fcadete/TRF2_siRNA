
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

pdf("DESeq2_analysis_no_interaction_ENCODE.pdf", width = 12)


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


seqlevels(selected_genes_48h) <- paste0("chr", seqlevels(selected_genes_48h))
seqlevels(selected_genes_96h) <- paste0("chr", seqlevels(selected_genes_96h))

all_genes <- genes(EnsDb.Hsapiens.v86)
seqlevels(all_genes) <- paste0("chr", seqlevels(all_genes))

expressed_genes <- genes(EnsDb.Hsapiens.v86,
                         filter = GeneIdFilter(rownames(
                                res96h[which(res96h$baseMean > 20),])))
seqlevels(expressed_genes) <- paste0("chr", seqlevels(expressed_genes))


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

hela_encode_metadata <- read_tsv("metadata.tsv", col_names = TRUE)

interesting_data <- hela_encode_metadata %>%
                      filter(`Output type` == "optimal idr thresholded peaks" | `Output type` == "replicated peaks",
                             `File format` == "bed narrowPeak",
                             `File Status` == "released",
                             `Assembly` == "GRCh38")

overlap_proportions <- as.data.frame(do.call(rbind,
     lapply(interesting_data$`File accession`, function(file_accession) {

       print(file_accession)

       this_bed <- rtracklayer::import(paste0("encode_files/", file_accession, ".bed.gz"),
                                       format = "BED",
                                       extraCols = extraCols_narrowPeak)

       all_genes_bed_overlaps <- countOverlaps(all_genes, this_bed) > 0
       all_genes_bed_overlap_proportion <- sum(all_genes_bed_overlaps == TRUE) /
                                            length(all_genes_bed_overlaps)

       expressed_genes_bed_overlaps <- countOverlaps(expressed_genes, this_bed) > 0
       expressed_genes_bed_overlap_proportion <- sum(expressed_genes_bed_overlaps == TRUE) /
                                                  length(expressed_genes_bed_overlaps)

       selected_genes_48h_bed_overlaps <- countOverlaps(selected_genes_48h, this_bed) > 0
       selected_genes_48h_bed_overlap_proportion <- sum(selected_genes_48h_bed_overlaps == TRUE) /
                                                     length(selected_genes_48h_bed_overlaps)

       selected_genes_96h_bed_overlaps <- countOverlaps(selected_genes_96h, this_bed) > 0
       selected_genes_96h_bed_overlap_proportion <- sum(selected_genes_96h_bed_overlaps == TRUE) /
                                                     length(selected_genes_96h_bed_overlaps)

       all_genes_bed_overlaps_promoters <- countOverlaps(promoters(all_genes), this_bed) > 0
       all_genes_bed_overlap_promoter_proportion <- sum(all_genes_bed_overlaps_promoters == TRUE) /
                                                     length(all_genes_bed_overlaps_promoters)
       
       expressed_genes_bed_overlaps_promoters <- countOverlaps(promoters(expressed_genes), this_bed) > 0
       expressed_genes_bed_overlap_promoter_proportion <- sum(expressed_genes_bed_overlaps_promoters == TRUE) /
                                                           length(expressed_genes_bed_overlaps_promoters)
       
       selected_genes_48h_bed_overlaps_promoters <- countOverlaps(promoters(selected_genes_48h), this_bed) > 0
       selected_genes_48h_bed_overlap_promoter_proportion <- sum(selected_genes_48h_bed_overlaps_promoters == TRUE) /
                                                              length(selected_genes_48h_bed_overlaps_promoters)
       
       selected_genes_96h_bed_overlaps_promoters <- countOverlaps(promoters(selected_genes_96h), this_bed) > 0
       selected_genes_96h_bed_overlap_promoter_proportion <- sum(selected_genes_96h_bed_overlaps_promoters == TRUE) /
                                                              length(selected_genes_96h_bed_overlaps_promoters)


       return(rbind(data.frame(`File accession` = file_accession, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_48h_bed_overlap_proportion, time_point = "48h",
                               all_genes_overlap = sum(all_genes_bed_overlaps == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps == FALSE),
                               selected_genes_overlap = sum(selected_genes_48h_bed_overlaps == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_48h_bed_overlaps == FALSE)),
                    data.frame(`File accession` = file_accession, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_96h_bed_overlap_proportion, time_point = "96h",
                               all_genes_overlap = sum(all_genes_bed_overlaps == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps == FALSE),
                               selected_genes_overlap = sum(selected_genes_96h_bed_overlaps == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_96h_bed_overlaps == FALSE)),

                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_48h_bed_overlap_promoter_proportion, time_point = "48h",
                               all_genes_overlap = sum(all_genes_bed_overlaps_promoters == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps_promoters == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps_promoters == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps_promoters == FALSE),
                               selected_genes_overlap = sum(selected_genes_48h_bed_overlaps_promoters == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_48h_bed_overlaps_promoters == FALSE)),
                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_96h_bed_overlap_promoter_proportion, time_point = "96h",
                               all_genes_overlap = sum(all_genes_bed_overlaps_promoters == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps_promoters == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps_promoters == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps_promoters == FALSE),
                               selected_genes_overlap = sum(selected_genes_96h_bed_overlaps_promoters == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_96h_bed_overlaps_promoters == FALSE))))
  }
)), stringsAsFactors = FALSE)

overlap_proportions <- left_join(overlap_proportions,
                                 select(interesting_data, `File accession`,
                                        `Experiment target`, `Experiment accession`),
                                 by = c(File.accession = "File accession"))


overlap_proportions_hg38 <- overlap_proportions %>%
                             mutate(all_genes = as.numeric(all_genes),
                                    expressed_genes = as.numeric(expressed_genes),
                                    selected = as.numeric(selected))


drip_seq_overlaps <- data.frame()

for (DRIP_seq in c("DRIP_seq_1", "DRIP_seq_2")) {

      drip_seq_peaks <- rtracklayer::import(paste0("drip_seq_raw_data/", DRIP_seq, "/", DRIP_seq, "_peaks.narrowPeak"))
      
       all_genes_bed_overlaps <- countOverlaps(all_genes, drip_seq_peaks) > 0
       all_genes_bed_overlap_proportion <- sum(all_genes_bed_overlaps == TRUE) /
                                            length(all_genes_bed_overlaps)

       expressed_genes_bed_overlaps <- countOverlaps(expressed_genes, drip_seq_peaks) > 0
       expressed_genes_bed_overlap_proportion <- sum(expressed_genes_bed_overlaps == TRUE) /
                                                  length(expressed_genes_bed_overlaps)

       selected_genes_48h_bed_overlaps <- countOverlaps(selected_genes_48h, drip_seq_peaks) > 0
       selected_genes_48h_bed_overlap_proportion <- sum(selected_genes_48h_bed_overlaps == TRUE) /
                                                     length(selected_genes_48h_bed_overlaps)

       selected_genes_96h_bed_overlaps <- countOverlaps(selected_genes_96h, drip_seq_peaks) > 0
       selected_genes_96h_bed_overlap_proportion <- sum(selected_genes_96h_bed_overlaps == TRUE) /
                                                     length(selected_genes_96h_bed_overlaps)

       all_genes_bed_overlaps_promoters <- countOverlaps(promoters(all_genes), drip_seq_peaks) > 0
       all_genes_bed_overlap_promoter_proportion <- sum(all_genes_bed_overlaps_promoters == TRUE) /
                                                     length(all_genes_bed_overlaps_promoters)
       
       expressed_genes_bed_overlaps_promoters <- countOverlaps(promoters(expressed_genes), drip_seq_peaks) > 0
       expressed_genes_bed_overlap_promoter_proportion <- sum(expressed_genes_bed_overlaps_promoters == TRUE) /
                                                           length(expressed_genes_bed_overlaps_promoters)
       
       selected_genes_48h_bed_overlaps_promoters <- countOverlaps(promoters(selected_genes_48h), drip_seq_peaks) > 0
       selected_genes_48h_bed_overlap_promoter_proportion <- sum(selected_genes_48h_bed_overlaps_promoters == TRUE) /
                                                              length(selected_genes_48h_bed_overlaps_promoters)
       
       selected_genes_96h_bed_overlaps_promoters <- countOverlaps(promoters(selected_genes_96h), drip_seq_peaks) > 0
       selected_genes_96h_bed_overlap_promoter_proportion <- sum(selected_genes_96h_bed_overlaps_promoters == TRUE) /
                                                              length(selected_genes_96h_bed_overlaps_promoters)

       drip_seq_overlaps <- rbind(drip_seq_overlaps,
                   data.frame(`File accession` = DRIP_seq, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_48h_bed_overlap_proportion, time_point = "48h",
                               all_genes_overlap = sum(all_genes_bed_overlaps == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps == FALSE),
                               selected_genes_overlap = sum(selected_genes_48h_bed_overlaps == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_48h_bed_overlaps == FALSE)),
                    data.frame(`File accession` = DRIP_seq, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_96h_bed_overlap_proportion, time_point = "96h",
                               all_genes_overlap = sum(all_genes_bed_overlaps == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps == FALSE),
                               selected_genes_overlap = sum(selected_genes_96h_bed_overlaps == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_96h_bed_overlaps == FALSE)),

                    data.frame(`File accession` = DRIP_seq, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_48h_bed_overlap_promoter_proportion, time_point = "48h",
                               all_genes_overlap = sum(all_genes_bed_overlaps_promoters == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps_promoters == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps_promoters == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps_promoters == FALSE),
                               selected_genes_overlap = sum(selected_genes_48h_bed_overlaps_promoters == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_48h_bed_overlaps_promoters == FALSE)),
                    data.frame(`File accession` = DRIP_seq, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_96h_bed_overlap_promoter_proportion, time_point = "96h",
                               all_genes_overlap = sum(all_genes_bed_overlaps_promoters == TRUE),
                               all_genes_non_overlap = sum(all_genes_bed_overlaps_promoters == FALSE),
                               expressed_genes_overlap = sum(expressed_genes_bed_overlaps_promoters == TRUE),
                               expressed_genes_non_overlap = sum(expressed_genes_bed_overlaps_promoters == FALSE),
                               selected_genes_overlap = sum(selected_genes_96h_bed_overlaps_promoters == TRUE),
                               selected_genes_non_overlap = sum(selected_genes_96h_bed_overlaps_promoters == FALSE)))

}


drip_seq_overlaps$`Experiment target` <- drip_seq_overlaps$File.accession
drip_seq_overlaps$`Experiment accession` <- drip_seq_overlaps$File.accession

overlap_proportions <- rbind(overlap_proportions, drip_seq_overlaps)
overlap_proportions_hg38 <- rbind(overlap_proportions_hg38, drip_seq_overlaps)


overlap_proportions_all_genes_tests <- rbind(full_join(overlap_proportions %>%
                  gather(key = "selection", value = "overlap", all_genes_overlap) %>%
                  select(File.accession, group, selection, time_point, overlap) %>%
                  mutate(selection = "all_genes"),
                overlap_proportions %>%
                  gather(key = "selection", value = "non_overlap", all_genes_non_overlap) %>%
                  select(File.accession, group, selection, time_point, non_overlap) %>%
                  mutate(selection = "all_genes")),
      
      full_join(overlap_proportions %>%
                  gather(key = "selection", value = "overlap", selected_genes_overlap) %>%
                  select(File.accession, group, selection, time_point, overlap) %>%
                  mutate(selection = "selected_genes"),
                overlap_proportions %>%
                  gather(key = "selection", value = "non_overlap", selected_genes_non_overlap) %>%
                  select(File.accession, group, selection, time_point, non_overlap) %>%
                  mutate(selection = "selected_genes"))) %>%
  group_by(File.accession, group, time_point) %>%
  summarise(test_matrix = list(matrix(c(overlap, non_overlap), nrow=2)),
            test_pvalue = fisher.test(test_matrix[[1]])$p.value)
overlap_proportions_all_genes_tests <- left_join(overlap_proportions_all_genes_tests,
                                 select(interesting_data, `File accession`,
                                        `Experiment target`, `Experiment accession`),
                                 by = c(File.accession = "File accession"))



overlap_proportions_expressed_genes_tests <- rbind(full_join(overlap_proportions %>%
                  gather(key = "selection", value = "overlap", expressed_genes_overlap) %>%
                  select(File.accession, group, selection, time_point, overlap) %>%
                  mutate(selection = "expressed_genes"),
                overlap_proportions %>%
                  gather(key = "selection", value = "non_overlap", expressed_genes_non_overlap) %>%
                  select(File.accession, group, selection, time_point, non_overlap) %>%
                  mutate(selection = "expressed_genes")),
      
      full_join(overlap_proportions %>%
                  gather(key = "selection", value = "overlap", selected_genes_overlap) %>%
                  select(File.accession, group, selection, time_point, overlap) %>%
                  mutate(selection = "selected_genes"),
                overlap_proportions %>%
                  gather(key = "selection", value = "non_overlap", selected_genes_non_overlap) %>%
                  select(File.accession, group, selection, time_point, non_overlap) %>%
                  mutate(selection = "selected_genes"))) %>%
  group_by(File.accession, group, time_point) %>%
  summarise(test_matrix = list(matrix(c(overlap, non_overlap), nrow=2)),
            test_pvalue = fisher.test(test_matrix[[1]])$p.value)
overlap_proportions_expressed_genes_tests <- left_join(overlap_proportions_expressed_genes_tests,
                                 select(interesting_data, `File accession`,
                                        `Experiment target`, `Experiment accession`),
                                 by = c(File.accession = "File accession"))


p <- ggplot(overlap_proportions_hg38,
            aes(x = `Experiment target`, y = time_point)) +
      geom_tile(mapping = aes(fill = log2(selected / all_genes))) +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = overlap_proportions_all_genes_tests,
                mapping = aes(label = ifelse(test_pvalue < 0.01, "*", " ")))
print(p)
p <- ggplot(overlap_proportions_hg38,
            aes(x = `Experiment target`, y = time_point)) +
      geom_tile(mapping = aes(fill = log2(selected / expressed_genes))) +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = overlap_proportions_expressed_genes_tests,
                mapping = aes(label = ifelse(test_pvalue < 0.01, "*", " ")))

print(p)



wig_files_GRCh38 <- hela_encode_metadata %>%
                       filter(`File format`== "bigWig",
                              `Output type` == "fold change over control",
                              `Assembly` == "GRCh38",
                              `File Status` == "released",
                              `Biological replicate(s)` == "1, 2")


signal_windows <- mclapply(wig_files_GRCh38$`File accession`, function(file_accession) {

       print(file_accession)

       this_wig <- rtracklayer::import(paste0("encode_files/", file_accession, ".bigWig"),
                                       format = "bigWig", as = "RleList")

       selected_genes_48h_TSS <- promoters(selected_genes_48h, downstream = 1000)
       selected_genes_48h_views <- Views(this_wig,
                                         GRangesList(lapply(names(this_wig),
                                          function(x) selected_genes_48h_TSS[seqnames(selected_genes_48h_TSS) == x])))
       selected_genes_48h_views <- selected_genes_48h_views[unlist(lapply(selected_genes_48h_views, length)) > 0]
       selected_genes_48h_matrix <- do.call(rbind, lapply(selected_genes_48h_views, as.matrix))
       selected_genes_48h_frame <- data.frame(reshape2::melt(selected_genes_48h_matrix),
                                              group = "48h_genes", file_accession)

       selected_genes_96h_TSS <- promoters(selected_genes_96h, downstream = 1000)
       selected_genes_96h_views <- Views(this_wig,
                                         GRangesList(lapply(names(this_wig),
                                          function(x) selected_genes_96h_TSS[seqnames(selected_genes_96h_TSS) == x])))
       selected_genes_96h_views <- selected_genes_96h_views[unlist(lapply(selected_genes_96h_views, length)) > 0]
       selected_genes_96h_matrix <- do.call(rbind, lapply(selected_genes_96h_views, as.matrix))
       selected_genes_96h_frame <- data.frame(reshape2::melt(selected_genes_96h_matrix),
                                              group = "96h_genes", file_accession)

       expressed_genes_TSS <- promoters(expressed_genes, downstream = 1000)
       expressed_genes_views <- Views(this_wig,
                                      GRangesList(lapply(names(this_wig),
                                       function(x) expressed_genes_TSS[seqnames(expressed_genes_TSS) == x])))
       expressed_genes_views <- expressed_genes_views[unlist(lapply(expressed_genes_views, length)) > 0]
       expressed_genes_matrix <- do.call(rbind, lapply(expressed_genes_views, as.matrix))
       expressed_genes_frame <- data.frame(reshape2::melt(expressed_genes_matrix),
                                           group = "Expressed", file_accession)

       return(rbind(selected_genes_48h_frame, selected_genes_96h_frame, expressed_genes_frame))

  }, mc.cores = 3

)

signal_windows <- do.call(rbind, signal_windows)

signal_windows <- left_join(signal_windows,
                            select(hela_encode_metadata, `File accession`, `Experiment target`),
                            by = c("file_accession" = "File accession"))

signal_windows <- signal_windows %>% ungroup() %>% mutate(diff_expressed = ifelse(group != "Expressed", TRUE, FALSE))

save(signal_windows, file="signal_windows.Rdata")


drip_seq_frame <- data.frame()

for (DRIP_seq in c("DRIP_seq_1", "DRIP_seq_2")) {

       drip_seq_signal <- rtracklayer::import(paste0("drip_seq_raw_data/", DRIP_seq, "/", DRIP_seq, "_logLR.bdg"),
                                              format = "bedGraph")
       drip_seq_signal <- coverage(drip_seq_signal, weight=drip_seq_signal$score)

       selected_genes_48h_TSS <- promoters(selected_genes_48h, downstream = 1000)
       selected_genes_48h_views <- Views(drip_seq_signal,
                                         GRangesList(lapply(names(drip_seq_signal),
                                          function(x) selected_genes_48h_TSS[seqnames(selected_genes_48h_TSS) == x])))
       selected_genes_48h_views <- selected_genes_48h_views[unlist(lapply(selected_genes_48h_views, length)) > 0]
       selected_genes_48h_matrix <- do.call(rbind, lapply(selected_genes_48h_views, as.matrix))
       selected_genes_48h_frame <- data.frame(reshape2::melt(selected_genes_48h_matrix),
                                              group = "48h_genes", file_accession = DRIP_seq)

       selected_genes_96h_TSS <- promoters(selected_genes_96h, downstream = 1000)
       selected_genes_96h_views <- Views(drip_seq_signal,
                                         GRangesList(lapply(names(drip_seq_signal),
                                          function(x) selected_genes_96h_TSS[seqnames(selected_genes_96h_TSS) == x])))
       selected_genes_96h_views <- selected_genes_96h_views[unlist(lapply(selected_genes_96h_views, length)) > 0]
       selected_genes_96h_matrix <- do.call(rbind, lapply(selected_genes_96h_views, as.matrix))
       selected_genes_96h_frame <- data.frame(reshape2::melt(selected_genes_96h_matrix),
                                              group = "96h_genes", file_accession = DRIP_seq)

       expressed_genes_TSS <- promoters(expressed_genes, downstream = 1000)
       expressed_genes_views <- Views(drip_seq_signal,
                                      GRangesList(lapply(names(drip_seq_signal),
                                       function(x) expressed_genes_TSS[seqnames(expressed_genes_TSS) == x])))
       expressed_genes_views <- expressed_genes_views[unlist(lapply(expressed_genes_views, length)) > 0]
       expressed_genes_matrix <- do.call(rbind, lapply(expressed_genes_views, as.matrix))
       expressed_genes_frame <- data.frame(reshape2::melt(expressed_genes_matrix),
                                           group = "Expressed", file_accession = DRIP_seq)

       drip_seq_frame <- rbind(drip_seq_frame,
                               selected_genes_48h_frame,
                               selected_genes_96h_frame,
                               expressed_genes_frame)

}      

drip_seq_frame$`Experiment target` <- drip_seq_frame$file_accession
drip_seq_frame <- drip_seq_frame %>% ungroup() %>% mutate(diff_expressed = ifelse(group != "Expressed", TRUE, FALSE))

signal_windows <- rbind(signal_windows, drip_seq_frame)

p <- signal_windows %>%
       group_by(group, `Experiment target`, Var2) %>%
       summarise(mean_per_pos = mean(value)) %>%
       ggplot(aes(x = Var2 - 2000, y = mean_per_pos, colour = group)) +
         geom_line() + facet_wrap(~ `Experiment target`, scales = "free_y") +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

p <- signal_windows %>%
       group_by(diff_expressed, `Experiment target`, Var2) %>%
       summarise(mean_per_pos = mean(value)) %>%
       ggplot(aes(x = Var2 - 2000, y = mean_per_pos, colour = diff_expressed)) +
         geom_line() + facet_wrap(~ `Experiment target`, scales = "free_y") +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)



dev.off()


