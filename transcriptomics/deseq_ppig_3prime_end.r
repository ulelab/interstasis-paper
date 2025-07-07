rm(list = ls())

library(tidyverse)
library(DESeq2)
library(ggsci)

# Read data -------------------------------------------------------------

gene_info <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details_exons.txt")

germs_clusters <- read_tsv("/camp/lab/ulej/home/users/farawar/GASR/germs/data/germs_CDS_umap_clusters.tsv.gz") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "region"), convert = T, remove = F, sep = ":|-|@") %>%
  left_join(gene_info)

counts_table <- read_tsv("/camp/lab/ulej/home/users/farawar/GASR/rnaseq/ePB_ER_Timecourse_1/nf-core-results/star_salmon/salmon.merged.gene_counts.tsv")

counts_matrix_raw <- counts_table %>%
  dplyr::select(-gene_name) %>%
  filter(gene_id %in% gene_info$gene_id) %>% #Taking forward protein coding transcripts only (this is 3' end seq, after all)
  column_to_rownames("gene_id") %>%
  as.matrix()

rscm <-counts_matrix_raw %>% rowSums()

counts_matrix_raw[rscm > quantile(rscm, 0.5),] ->
  counts_matrix

sample_annotation <- data.frame(sample = colnames(counts_matrix)) %>%
  mutate(timepoint = word(sample, 1, sep = "_") %>% fct_relevel("h0", "h4", "h8", "h12"),
         fraction = word(sample, 2, sep = "_") %>% fct_relevel("N")) %>%
  column_to_rownames("sample")

# Fractionation efficiency ------------------------------------------------

unfiltered_counts_matrix <- counts_table %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

unfiltered_dds <- DESeqDataSetFromMatrix(countData = round(unfiltered_counts_matrix),
                                         colData = sample_annotation,
                                         design = ~ fraction * timepoint)

fractionation_genes <- c("NEAT1", "RPL12")

gene_counts_table <- lapply(fractionation_genes, function(x) {
  gene_id <- counts_table$gene_id[counts_table$gene_name == x]
  counts_table <- plotCounts(unfiltered_dds, gene = gene_id, intgroup = c("fraction", "timepoint"), returnData = T) %>%
    rownames_to_column("sample_name")
  return(counts_table)
}) %>%
  set_names(fractionation_genes) %>%
  bind_rows(.id = "gene_name") %>%
  mutate(gene_name = fct_relevel(gene_name, "MALAT1", "NEAT1"),
         fraction = case_when(fraction == "N" ~ "Nuc", T ~ "Cyto") %>% fct_relevel("Nuc"))

fractionation_genes_plot <- ggplot(gene_counts_table, aes(x = timepoint, y = count, fill = fraction, group = sample_name)) +
  facet_wrap(. ~ gene_name, scales = "free_y") +
  theme_classic() +
  scale_y_log10() +
  # expand_limits(y = 0) +
  scale_fill_d3() +
  labs(x = "", y = "normalised counts", fill = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black")) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21)

fractionation_genes_plot

# ggsave(filename = "GASR/export_reporter/quantseq/neve_timecourse_1/plots/fractionation_genes_counts.pdf", 
#        plot = fractionation_genes_plot,
#        device = "pdf", units = "in",
#        height = 3, width = 6)


# Example genes -----------------------------------------------------------

example_genes <- c("EIF3A", "BRD4", "PSAP", "HNRNPDL")

example_gene_counts_table <- lapply(example_genes, function(x) {
  gene_id <- counts_table$gene_id[counts_table$gene_name == x]
  counts_table <- plotCounts(unfiltered_dds, gene = gene_id, intgroup = c("fraction", "timepoint"), returnData = T) %>%
    rownames_to_column("sample_name")
  return(counts_table)
}) %>%
  set_names(example_genes) %>%
  bind_rows(.id = "gene_name") %>%
  mutate(gene_name = gene_name %>% fct_relevel("PSAP", "HNRNPDL"),
    fraction = case_when(fraction == "N" ~ "Nuc", T ~ "Cyto") %>% fct_relevel("Nuc"))

example_genes_plot <- ggplot(example_gene_counts_table, aes(x = timepoint, y = count, fill = fraction, group = sample_name)) +
  facet_wrap(. ~ gene_name, scales = "free_y", ncol = 2) +
  theme_classic() +
  # scale_y_log10() +
  # expand_limits(y = 0) +
  scale_fill_d3() +
  labs(x = "", y = "normalised counts", fill = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black")) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21)

example_genes_plot

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/example_gene_counts.pdf", example_genes_plot, width = 5, height = 5)


# individual timepoint DESeq -------------------------------------------------------------

timepoints <- c("h0", "h4", "h8", "h12")

timepoint_df_list_shrunk <- lapply(timepoints, function(my_timepoint) {
  
  # my_timepoint <- 'h0'
  
  timepoint_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,grepl(my_timepoint, colnames(counts_matrix))]),
                                          colData = sample_annotation %>% filter(timepoint == my_timepoint),
                                          design = ~ fraction)
  
  timepoint_deseq <- DESeq(timepoint_dds)
  
  timepoint_shrunk <- lfcShrink(timepoint_deseq, coef = 2, type = "apeglm")
  
  timepoint_shrunk_df <-  as.data.frame(timepoint_shrunk) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_info) %>%
    mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
           gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
           ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
           crich = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
    ) %>%
    mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                  T ~ "Not AG-rich") %>%
             fct_relevel("Not AG-rich"))
  
  return(timepoint_shrunk_df)
  
})

timepoint_df_list_shrunk %>%
  set_names(timepoints) %>%
  bind_rows(.id = "timepoint") ->
  bound_timepoints_df_shrunk

write_tsv(bound_timepoints_df_shrunk %>% select(-c(aagaa, gagga, ggagc, crich, germs_type)), "GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_each_timepoint.tsv")

# Interaction DESeq -------------------------------------------------------------------

interaction_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix)[,grepl("h0|h12", colnames(counts_matrix))],
                                          colData = sample_annotation %>% filter(timepoint %in% c("h0", "h12")),
                                          design = ~ fraction * timepoint)

interaction_deseq <- DESeq(interaction_dds)

interaction_shrunk <- lfcShrink(interaction_deseq, coef = 4, type = "apeglm")

interaction_shrunk_df <- as.data.frame(interaction_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info)

write_tsv(interaction_shrunk_df, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_interaction_df.tsv")


