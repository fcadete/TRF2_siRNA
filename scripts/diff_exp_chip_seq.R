
library("DESeq2")
library("Repitools")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

library("parallel")

pdf("diff_expressed_encode_analysis.pdf", width = 10)

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
unique(genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_48h <- genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint48.results[which(siRNATRF2.TimePoint48.results$padj < 0.1),])))
seqlevels(selected_genes_48h) <- paste0("chr", seqlevels(selected_genes_48h))

# Number of diff-expressed genes in the siRNATRF2.TimePoint96 interaction
siRNATRF2.TimePoint96.results <- results(ddsTxi, name = "siRNATRF2.TimePoint96")
unique(genes(EnsDb.Hsapiens.v86,
             filter = GeneIdFilter(rownames(
               siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))$symbol)

# Get and plot the distance that these genes have to the nearest end of their respective chromosomes
selected_genes_96h <- genes(EnsDb.Hsapiens.v86,
                        filter = GeneIdFilter(rownames(
                          siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$padj < 0.1),])))
seqlevels(selected_genes_96h) <- paste0("chr", seqlevels(selected_genes_96h))

all_genes <- genes(EnsDb.Hsapiens.v86)
seqlevels(all_genes) <- paste0("chr", seqlevels(all_genes))

expressed_genes <- genes(EnsDb.Hsapiens.v86,
                         filter = GeneIdFilter(rownames(
                                siRNATRF2.TimePoint96.results[which(siRNATRF2.TimePoint96.results$baseMean > 20),])))
seqlevels(expressed_genes) <- paste0("chr", seqlevels(expressed_genes))

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

hela_encode_metadata <- read_tsv("metadata.tsv", col_names = TRUE)

interesting_data <- hela_encode_metadata %>%
                      filter(`Output type` == "optimal idr thresholded peaks",
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
                               selected = selected_genes_48h_bed_overlap_proportion, time_point = "48h"),
                    data.frame(`File accession` = file_accession, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_96h_bed_overlap_proportion, time_point = "96h"),

                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_48h_bed_overlap_promoter_proportion, time_point = "48h"),
                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_96h_bed_overlap_promoter_proportion, time_point = "96h")))
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

p <- ggplot(overlap_proportions_hg38,
            aes(x = `Experiment target`, y = time_point, fill = log2(selected / all_genes))) +
      geom_tile() +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
p <- ggplot(overlap_proportions_hg38,
            aes(x = `Experiment target`, y = time_point, fill = log2(selected / expressed_genes))) +
      geom_tile() +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)


chain <- rtracklayer::import.chain("hg19ToHg38.over.chain")

interesting_data <- hela_encode_metadata %>%
                      filter(`Output type` == "optimal idr thresholded peaks",
                             `File format` == "bed narrowPeak",
                             `File Status` == "released",
                             `Assembly` == "hg19")

overlap_proportions <- as.data.frame(do.call(rbind,
     lapply(interesting_data$`File accession`, function(file_accession) {

       print(file_accession)

       this_bed <- rtracklayer::import(paste0("encode_files/", file_accession, ".bed.gz"),
                                       format = "BED",
                                       extraCols = extraCols_narrowPeak)

       this_bed <- unlist(rtracklayer::liftOver(this_bed, chain))

       all_genes_bed_overlaps <- countOverlaps(all_genes, this_bed) > 0
       all_genes_bed_overlap_proportion <- sum(all_genes_bed_overlaps == TRUE) / length(all_genes_bed_overlaps)

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
                               selected = selected_genes_48h_bed_overlap_proportion, time_point = "48h"),
                    data.frame(`File accession` = file_accession, group = "all",
                               all_genes = all_genes_bed_overlap_proportion,
                               expressed_genes = expressed_genes_bed_overlap_proportion,
                               selected = selected_genes_96h_bed_overlap_proportion, time_point = "96h"),

                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_48h_bed_overlap_promoter_proportion, time_point = "48h"),
                    data.frame(`File accession` = file_accession, group = "promoters",
                               all_genes = all_genes_bed_overlap_promoter_proportion,
                               expressed_genes = expressed_genes_bed_overlap_promoter_proportion,
                               selected = selected_genes_96h_bed_overlap_promoter_proportion, time_point = "96h")))
  }
)), stringsAsFactors = FALSE)

overlap_proportions <- left_join(overlap_proportions,
                                 select(interesting_data, `File accession`,
                                        `Experiment target`, `Experiment accession`),
                                 by = c(File.accession = "File accession"))

overlap_proportions_hg19 <- overlap_proportions %>%
                             mutate(all_genes = as.numeric(all_genes),
                                    expressed_genes = as.numeric(expressed_genes),
                                    selected = as.numeric(selected))

p <- ggplot(overlap_proportions_hg19,
            aes(x = `Experiment target`, y = time_point, fill = log2(selected / all_genes))) +
      geom_tile() +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
p <- ggplot(overlap_proportions_hg19,
            aes(x = `Experiment target`, y = time_point, fill = log2(selected / expressed_genes))) +
      geom_tile() +
      coord_equal() +
      scale_fill_gradient2(low = "darkblue", high = "gold", mid = "white") +
      facet_grid(group ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
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

       selected_genes_48h_TSS <- promoters(selected_genes_48h)
       selected_genes_48h_views <- Views(this_wig,
                                         GRangesList(lapply(names(this_wig),
                                          function(x) selected_genes_48h_TSS[seqnames(selected_genes_48h_TSS) == x])))
       selected_genes_48h_views <- selected_genes_48h_views[unlist(lapply(selected_genes_48h_views, length)) > 0]
       selected_genes_48h_matrix <- do.call(rbind, lapply(selected_genes_48h_views, as.matrix))
       selected_genes_48h_frame <- data.frame(reshape2::melt(selected_genes_48h_matrix),
                                              group = "48h_genes", file_accession)

       selected_genes_96h_TSS <- promoters(selected_genes_96h)
       selected_genes_96h_views <- Views(this_wig,
                                         GRangesList(lapply(names(this_wig),
                                          function(x) selected_genes_96h_TSS[seqnames(selected_genes_96h_TSS) == x])))
       selected_genes_96h_views <- selected_genes_96h_views[unlist(lapply(selected_genes_96h_views, length)) > 0]
       selected_genes_96h_matrix <- do.call(rbind, lapply(selected_genes_96h_views, as.matrix))
       selected_genes_96h_frame <- data.frame(reshape2::melt(selected_genes_96h_matrix),
                                              group = "96h_genes", file_accession)

       expressed_genes_TSS <- promoters(expressed_genes)
       expressed_genes_views <- Views(this_wig,
                                      GRangesList(lapply(names(this_wig),
                                       function(x) expressed_genes_TSS[seqnames(expressed_genes_TSS) == x])))
       expressed_genes_views <- expressed_genes_views[unlist(lapply(expressed_genes_views, length)) > 0]
       expressed_genes_matrix <- do.call(rbind, lapply(expressed_genes_views, as.matrix))
       expressed_genes_frame <- data.frame(reshape2::melt(expressed_genes_matrix),
                                           group = "Expressed", file_accession)

       return(rbind(selected_genes_48h_frame, selected_genes_96h_frame, expressed_genes_frame))

  }, mc.cores = 12

)

signal_windows <- do.call(rbind, signal_windows)

signal_windows <- left_join(signal_windows,
                            select(hela_encode_metadata, `File accession`, `Experiment target`),
                            by = c("file_accession" = "File accession"))

signal_windows <- signal_windows %>% ungroup() %>% mutate(diff_expressed = ifelse(group != "Expressed", TRUE, FALSE))

save(signal_windows, file="signal_windows.Rdata")


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

