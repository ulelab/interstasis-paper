rm(list = ls())

library(tidyverse)
library(DESeq2)
library(ggsci)
library(patchwork)

# Read data -------------------------------------------------------------

gene_info <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt")

germs_clusters <- read_tsv("/camp/lab/ulej/home/users/farawar/GASR/germs/data/germs_CDS_umap_clusters.tsv.gz") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "region"), convert = T, remove = F, sep = ":|-|@") %>%
  left_join(gene_info) %>%
  mutate(cluster_width = end - start)

counts_table <- read_tsv("/camp/lab/ulej/home/shared/interstasis/rnaseq/clkint3_ncf/latest_rnaseq_version/nf-core/star_salmon/salmon.merged.gene_counts.tsv")

# counts_table %>% select(-control_1_n) -> counts_table

counts_matrix <- counts_table %>%
  dplyr::select(-gene_name) %>%
  filter(gene_id %in% gene_info$gene_id) %>% #Taking forward protein coding transcripts only (this is 3' end seq, after all)
  mutate(rowsum = rowSums(across(where(is.numeric)))) %>%
  filter(rowsum > quantile(rowsum, 0.5)) %>%
  dplyr::select(-rowsum) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# rowsum_cutoff <- quantile(rowSums(counts_matrix), 0.25)
# 
# counts_matrix <- counts_matrix[rowSums(counts_matrix) >= rowsum_cutoff,]

# colnames(counts_matrix) <- toupper(colnames(counts_matrix)) %>% str_replace("CONTROL", "Control")

sample_annotation <- data.frame(sample = colnames(counts_matrix)) %>%
  mutate(condition = word(sample, 1, sep = "_") %>% fct_relevel("DMSO"),
         fraction = word(sample, 2, sep = "_") %>% fct_relevel("N")) %>%
  column_to_rownames("sample")

# Fractionation efficiency ------------------------------------------------

unfiltered_counts_matrix <- counts_table %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# colnames(unfiltered_counts_matrix) <- toupper(colnames(unfiltered_counts_matrix)) %>% str_replace("CONTROL", "Control")

unfiltered_dds <- DESeqDataSetFromMatrix(countData = round(unfiltered_counts_matrix),
                                         colData = sample_annotation,
                                         design = ~ fraction * condition)

fractionation_genes <- c("MALAT1", "NEAT1", "GAPDH", "RPL12")

gene_counts_table <- lapply(fractionation_genes, function(x) {
  gene_id <- counts_table$gene_id[counts_table$gene_name == x]
  counts_table <- plotCounts(unfiltered_dds, gene = gene_id, intgroup = c("fraction", "condition"), returnData = T) %>%
    rownames_to_column("sample_name")
  return(counts_table)
}) %>%
  set_names(fractionation_genes) %>%
  bind_rows(.id = "gene_name") %>%
  mutate(gene_name = fct_relevel(gene_name, "MALAT1", "NEAT1"),
         fraction = case_when(fraction == "N" ~ "Nuc", T ~ "Cyto") %>% fct_relevel("Nuc"))

fractionation_genes_plot <- ggplot(gene_counts_table, aes(x = condition, y = count, fill = fraction, group = sample_name)) +
  facet_wrap(. ~ gene_name, scales = "free_y") +
  theme_classic() +
  scale_y_log10() +
  # expand_limits(y = 0) +
  scale_fill_d3() +
  labs(x = "", y = "normalised counts", fill = "") +
  theme(axis.text.x = element_text(size = 10, color = "black")) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21)

fractionation_genes_plot

# ggsave(filename = "GASR/rnaseq/LCD_CheapSeq_1/plots/fractionation_genes_counts.pdf", 
#        plot = fractionation_genes_plot,
#        device = "pdf", units = "in",
#        height = 4, width = 6)

# Control DESeq -------------------------------------------------------------

dmso_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,grepl("DMSO", colnames(counts_matrix))]),
                                      colData = sample_annotation %>% filter(condition == "DMSO"),
                                      design = ~ fraction)

dmso_deseq <- DESeq(dmso_dds)

dmso_shrunk <- lfcShrink(dmso_deseq, coef = 2, type = "apeglm")

dmso_shrunk_df <-  as.data.frame(dmso_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

ggplot(dmso_shrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)),
       aes(x = log2FoldChange, fill = germs_type)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-4, 4)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama()

# CLKINT3 DESeq -------------------------------------------------------------

clkint3_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,grepl("CLKINT3", colnames(counts_matrix))]),
                                      colData = sample_annotation %>% filter(condition == "CLKINT3"),
                                      design = ~ fraction)

clkint3_deseq <- DESeq(clkint3_dds)

clkint3_shrunk <- lfcShrink(clkint3_deseq, coef = 2, type = "apeglm")

clkint3_shrunk_df <-  as.data.frame(clkint3_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

ggplot(clkint3_shrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)),
       aes(x = log2FoldChange, fill = germs_type)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-4, 4)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama()

# All in one df ------------------------------------------------

all_df <- bind_rows(list("DMSO" = dmso_shrunk_df,
                         "CLKINT3" = clkint3_shrunk_df),
                    .id = "condition") %>%
  mutate(condition = fct_relevel(condition, "DMSO"))

nuclear_cytoplasmic_boxplots <- ggplot(all_df %>%
                                         filter(baseMean > quantile(baseMean, 0.5)),
                                       aes(x = germs_type, y = -log2FoldChange, fill = condition)) +
  geom_boxplot(alpha = 0.75, width = 0.4, position = position_dodge(width = 0.5)) +
  # geom_violin(alpha = 0.35, width = 0.5, position = position_dodge(width = 0.5), show.legend = F) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_fill_aaas() +
  labs(y = "log2 fold change\n(nuc / cyto)") +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"))

nuclear_cytoplasmic_boxplots

# ggsave(filename = "GASR/export_reporter/quantseq/short_vs_long_1/plots/nuclear_cytoplasmic_boxplots.pdf",
#        plot = nuclear_cytoplasmic_boxplots,
#        device = "pdf", units = "in",
#        height = 4, width = 4)
# 
# lm(data = short_and_long_df %>%
#      filter(baseMean > quantile(baseMean, 0.5)),
#    formula = log2FoldChange ~ germs_type * condition) %>%
#   summary()


# Nuclear  ----------------------------------------------------------------

nuc_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,grepl("_N_", colnames(counts_matrix))]),
                                   colData = sample_annotation %>% filter(fraction == "N"),
                                   design = ~ condition)

nuc_deseq <- DESeq(nuc_dds)

nuc_shrunk <- lfcShrink(nuc_deseq, coef = 2, type = "apeglm")

nuc_shrunk_df <-  as.data.frame(nuc_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

ggplot(nuc_shrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)),
       aes(x = log2FoldChange, fill = germs_type)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama()

# Cytoplasm  ----------------------------------------------------------------

cyto_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix[,grepl("_C_", colnames(counts_matrix))]),
                                  colData = sample_annotation %>% filter(fraction == "C"),
                                  design = ~ condition)

cyto_deseq <- DESeq(cyto_dds)

cyto_shrunk <- lfcShrink(cyto_deseq, coef = 2, type = "apeglm")

cyto_shrunk_df <-  as.data.frame(cyto_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

ggplot(cyto_shrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)),
       aes(x = log2FoldChange, fill = germs_type)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama()

# Interaction shrunk DESeq -------------------------------------------------------------------

interaction_dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix),
                                          colData = sample_annotation,
                                          design = ~ fraction * condition)

interaction_deseq <- DESeq(interaction_dds)

resultsNames(interaction_deseq)

interaction_shrunk <- lfcShrink(interaction_deseq, coef = 4, type = "apeglm")

interaction_shrunk_df <- as.data.frame(interaction_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

write_tsv(interaction_shrunk_df, "GASR/clk/tables/neve_quantseq_interaction_df.newmapping.tsv")

# Interaction plot --------------------------------------------------------

interaction_plot <- ggplot(interaction_shrunk_df %>% 
                                    filter(baseMean > quantile(baseMean, 0.5)), 
                                  aes(x = log2FoldChange, fill = germs_type)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-2, 2)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama() +
  labs(x = "log2 fold change (fraction * condition)", fill = "") +
  theme(plot.title = element_text(hjust = 0.5))

interaction_plot

# ggsave(filename = "GASR/rnaseq/LCD_CheapSeq_1/plots/merged_interaction_plots.pdf",
#        plot = merged_plot,
#        device = "pdf", units = "in",
#        height = 4, width = 12)

# Making target lists -----------------------------------------------------

interaction_shrunk_df %>%
  filter(padj < 0.05, 
         baseMean > quantile(baseMean, 0.5),
         ) %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  more_retained_in_nucleus

interaction_shrunk_df %>%
  filter(padj < 0.05, 
         baseMean > quantile(baseMean, 0.5),
         ) %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  more_exported_from_nucleus

interaction_shrunk_df %>%
  filter(baseMean > quantile(baseMean, 0.5),) %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  expressed_in_interaction_set

# Gene ontology -----------------------------------------------------------


library(topGO)
library(biomaRt)

all_genes <- gene_info$gene_name

#Get GO terms from Ensembl.
db <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "https://www.ensembl.org")
go_ids <- getBM(attributes = c('go_id', 'external_gene_name', 'namespace_1003'), 
                filters = 'external_gene_name', values = all_genes, mart = db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO <- unstack(go_ids[,c(1,2)])

go_term_table <- function(targets, background){
  
  target_genes <- targets %>%
    unique()
  
  target_genes_f <- target_genes[target_genes %in% go_ids[,2]]
  
  geneList <- factor(as.integer(background %in% target_genes_f))
  names(geneList) <- background
  
  big_go_table <- lapply(c("BP", "MF", "CC"), function(x){
    
    GOdata <- new('topGOdata', 
                  ontology = x, 
                  allGenes = geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = gene_2_GO)
    
    weight_fisher_result <- runTest(GOdata, algorithm='weight01', statistic='fisher')
    # classic_fisher_result <- runTest(GOdata, algorithm='classic', statistic='fisher')
    
    allGO <- usedGO(GOdata)
    
    all_res <- GenTable(GOdata,
                        # Fis = classic_fisher_result,
                        weightFisher = weight_fisher_result,
                        orderBy = 'weightFisher',
                        topNodes = length(allGO), 
                        numChar = 1000)
    
  }) %>%
    set_names(c("BP", "MF", "CC")) %>%
    bind_rows(.id = "go_type") %>%
    mutate(fold_enrichment = Significant / Expected,
           weightFisher = as.numeric(weightFisher)) %>%
    replace_na(list(weightFisher = 0))
}

retained_in_nucleus <- go_term_table(more_retained_in_nucleus, expressed_in_interaction_set)
exported_from_nucleus <- go_term_table(more_exported_from_nucleus, expressed_in_interaction_set)

nuclear_speck_genes <- go_ids %>%
  filter(go_id == "GO:0016607") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

upr_genes <- go_ids %>%
  filter(go_id == "GO:0051082") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

cytotrans_genes <- go_ids %>%
  filter(go_id == "GO:0002181") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

ggplot(interaction_shrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)), 
       aes(x = log2FoldChange, fill = gene_name %in% upr_genes)) + 
  geom_density(alpha = 0.5) + 
  # coord_cartesian(xlim = c(-2.5, 1.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama() +
  labs(x = "log2 fold change (fraction * timepoint", fill = "")

# What about expression ---------------------------------------------------

drug_shrunk <- lfcShrink(interaction_deseq, coef = 3, type = "apeglm")

drug_shrunk_df <- as.data.frame(drug_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))


drug_unshrunk <- results(interaction_deseq, name = resultsNames(interaction_deseq)[3])

drug_unshrunk_df <- as.data.frame(drug_unshrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info) %>%
  mutate(aagaa = gene_name %in% filter(germs_clusters, cluster == 2)$gene_name,
         gagga = gene_name %in% filter(germs_clusters, cluster == 1)$gene_name,
         ggagc = gene_name %in% filter(germs_clusters, cluster == 7)$gene_name,
         gc = gene_name %in% filter(germs_clusters, cluster == 4)$gene_name,
  ) %>%
  mutate(germs_type = case_when(aagaa | gagga ~ "AG-rich",
                                T ~ "Not AG-rich") %>%
           fct_relevel("Not AG-rich"))

ggplot(drug_unshrunk_df %>% 
         filter(baseMean > quantile(baseMean, 0.5)), 
       aes(x = log2FoldChange, fill = aagaa)) + 
  geom_density(alpha = 0.5) + 
  coord_cartesian(xlim = c(-3, 3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama() +
  labs(x = "log2 fold change (condition)", fill = "") +
  theme(plot.title = element_text(hjust = 0.5))

# Interaction effect based on initial ratio -------------------------------

# I am reversing the sign for the initial l2c values because I want them to be nuclear-cytoplasmic ratios.
interaction_with_initial_df <- inner_join(interaction_shrunk_df,
                                          dmso_shrunk_df %>%
                                            dplyr::select(gene_id, initial_l2fc = log2FoldChange)) %>%
  mutate(initial_l2fc = -initial_l2fc) %>%
  mutate(intial_localisation_class = cut_number(initial_l2fc, 6) %>%
           fct_relabel(~ str_remove_all(.x, pattern = fixed("(")) %>%
                         str_remove_all(pattern = fixed("-")) %>%
                         str_remove_all(pattern = fixed("]")) %>%
                         str_remove_all(pattern = fixed("[")) %>%
                         str_replace_all(",", " - "))) %>%
  drop_na()

interaction_plot_by_initial_l2fc <- ggplot(interaction_with_initial_df %>%
                                             filter(baseMean > quantile(baseMean, 0.5)), 
                                           aes(x = log2FoldChange, fill = germs_type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ intial_localisation_class) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(x = "log2 fold change (fraction * condition", fill = "") +
  scale_fill_jama()

interaction_plot_by_initial_l2fc

# ggsave(filename = "GASR/export_reporter/quantseq/short_vs_long_1/plots/interaction_density_plot_initial_l2fc.pdf",
#        plot = interaction_plot_by_initial_l2fc,
#        device = "pdf", units = "in",
#        height = 4, width = 8)

interaction_plot_filtered <- ggplot(interaction_with_initial_df %>%
                                      filter(baseMean > quantile(baseMean, 0.5),
                                             initial_l2fc < 1),
                                    aes(x = log2FoldChange, fill = germs_type)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(-2, 2)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_fill_jama() +
  labs(x = "log2 fold change (fraction * condition", fill = "")

interaction_plot_filtered

# ggsave(filename = "GASR/export_reporter/quantseq/short_vs_long_1/plots/interaction_density_plot_filtered.pdf",
#        plot = interaction_plot_filtered,
#        device = "pdf", units = "in",
#        height = 4, width = 6)

# Nice example genes ------------------------------------------------------

example_genes <- c("OGFR", "INO80", "CREBBP", "U2SURP", "SAFB2", "TOP2A")

example_gene_counts_table <- lapply(example_genes, function(x) {
  gene_id <- counts_table$gene_id[counts_table$gene_name == x]
  counts_table <- plotCounts(unfiltered_dds, gene = gene_id, intgroup = c("fraction", "condition"), returnData = T) %>%
    rownames_to_column("sample_name")
  return(counts_table)
}) %>%
  set_names(example_genes) %>%
  bind_rows(.id = "gene_name") %>%
  mutate(fraction = case_when(fraction == "N" ~ "Nuc", T ~ "Cyto") %>% fct_relevel("Nuc"))

example_genes_plot <- ggplot(example_gene_counts_table, aes(x = condition, y = count, fill = fraction, group = sample_name)) +
  facet_wrap(. ~ gene_name, scales = "free_y") +
  theme_classic() +
  scale_y_log10() +
  # expand_limits(y = 0) +
  scale_fill_d3() +
  labs(x = "", y = "normalised counts", fill = "") +
  theme(axis.text.x = element_text(size = 10, color = "black")) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21)

example_genes_plot

# ggsave(filename = "GASR/export_reporter/quantseq/short_vs_long_1/plots/example_genes_counts.pdf", 
#        plot = example_genes_plot,
#        device = "pdf", units = "in",
#        height = 4, width = 6)
# 
