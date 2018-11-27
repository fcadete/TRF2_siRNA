
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")

pdf("DESeq2_analysis_lncRNA.pdf")


samples <- read_tsv("sample_info.txt")
samples$TimePoint <- as.character(samples$TimePoint)

read_table <- data.frame()
for (sample in samples$Filename) {

file_name <- paste0("star_on_hg38_output/", sample, "ReadsPerGene.out.tab")

this_sample <- read_tsv(file_name, col_names = c("gene", "unstranded", "sense", "antisense"))

read_table <- rbind(read_table, data.frame(sample = sample, this_sample))

}

read_table <- read_table %>%
                 dplyr::select(sample, gene, sense) %>%
                 dplyr::filter(!grepl("N_", gene))


read_matrix <- read_table %>%
                 tidyr::spread(sample, sense)

rownames(read_matrix) <- read_matrix$gene

read_matrix <- read_matrix %>% dplyr::select(-gene)

ddsLncRNA <- DESeqDataSetFromMatrix(read_matrix,
                                   colData = samples,
                                   design = ~ siRNA + TimePoint + siRNA:TimePoint)

vsd <- vst(ddsLncRNA, blind=TRUE)

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



read_table_gencode <- data.frame()
for (sample in samples$Filename) {

file_name <- paste0("star_on_hg38_output_gencode/", sample, "ReadsPerGene.out.tab")

this_sample <- read_tsv(file_name, col_names = c("gene", "unstranded", "sense", "antisense"))

read_table_gencode <- rbind(read_table_gencode, data.frame(sample = sample, this_sample))

}

read_table_gencode <- read_table_gencode %>%
                 dplyr::select(sample, gene, sense) %>%
                 dplyr::filter(!grepl("N_", gene))


p <- ggplot(read_table_gencode,
       aes(x=sense+1, colour = sample, linetype = gene %in% unique(read_table$gene))) +
   geom_density() +
   scale_x_log10()
print(p)


p <- ggplot(read_table_gencode %>% filter(sample %in% c("C_NT_1", "C_NT_2", "C_NT_3")),
       aes(x=sense+1, colour = sample, linetype = gene %in% unique(read_table$gene))) +
   geom_density() +
   scale_x_log10()
print(p)

dev.off()


