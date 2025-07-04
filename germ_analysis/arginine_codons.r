# Setup -------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(Biostrings)
library(broom)
library(ggsci)
library(patchwork)

# Create codon and aa tables ----------------------------------------------

aa_df <- tibble(aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

nts <- c("A", "C", "G", "T")

codon_matrix <- do.call(expand.grid, rep(list(nts), 3))

codon_table <- tibble(codon = do.call(paste0, c(codon_matrix))) %>%
  mutate(aa = as.character(translate(DNAStringSet(codon), no.init.codon = T))) %>%
  group_by(aa) %>%
  mutate(codon_number = c(1:length(aa))) %>%
  ungroup() %>%
  arrange(aa)

# Load stuff --------------------------------------------------------------

txseq <- readDNAStringSet("GASR/lists/longest_gencode29.fa")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

aa_per_region_df <- read_tsv("GASR/sequence_multivalency/data/20220315_123win_98per_k5_10nei/data/percentile_98_windowsize_123_klen_5_le_regions_with_aa_props.txt.gz")

r_rich_df <- aa_per_region_df %>% filter(R >= 0.2)

r_rich_df %>%
  dplyr::select(c(1, 5:24)) %>%
  pivot_longer(cols = -1, names_to = "aa", values_to = "proportion") %>%
  ggplot(aes(x = aa, fill = aa, y = proportion)) +
  geom_boxplot(alpha = 0.5) + 
  theme_classic()

low_ent_positions <- r_rich_df %>%
  dplyr::select(region, transcript_id, protein_region_start, protein_region_end)


# Get codons and transcriptome wide codon usage --------------------------------------------------------------

#Get CDS, protein and codon sequences for target transcripts.
seq_df <- tibble(transcript_id = names(txseq),
                 transcript_seq = as.character(txseq)) %>%
  right_join(transcript_details, by = "transcript_id") %>%
  mutate(cds_seq = pmap(list(cds_start, cds_end, transcript_seq), 
                        function(cds_start, cds_end, transcript_seq){
                          str_sub(transcript_seq, cds_start, cds_end)
                        }) %>% 
           unlist()) %>%
  mutate(protein_seq = translate(cds_seq %>% 
                                   DNAStringSet()) %>%
           as.character()) %>%
  mutate(codons = purrr::map(cds_seq, function(cds_seq){
    substring(cds_seq, seq(1, nchar(cds_seq), 3), seq(3, nchar(cds_seq), 3))
  }))

transcriptome_codon_usage <- seq_df$codons %>% 
  unlist() %>%
  table() %>%
  as_tibble() %>%
  set_names(c("codon", "count")) %>%
  left_join(codon_table, by = "codon") %>%
  arrange(aa, codon) %>%
  group_by(aa) %>%
  mutate(proportion_usage_transcriptome = count / sum(count)) %>%
  dplyr::select(-c(count, codon_number))

# Chop to codons ----------------------------------------------------------

#Get the codons within the low entropy region
codons_in_leds <- low_ent_positions %>%
  left_join(seq_df %>%
              filter(transcript_id %in% r_rich_df$transcript_id) %>%
              dplyr::select(transcript_id, codons)) %>%
  mutate(codons_in_domain = pmap(list(protein_region_start, protein_region_end, codons),
                                 function(protein_region_start, protein_region_end, codons){
                                   codons[protein_region_start:protein_region_end]
                                 })) %>%
  mutate(codons_frequency = map(codons_in_domain, ~{
    tibble(codon = .x) %>%
      group_by(codon) %>%
      count() %>%
      right_join(codon_table, by = "codon")
  })) %>% 
  dplyr::select(-c(protein_region_start, protein_region_end, codons, codons_in_domain)) %>%
  unnest(codons_frequency) %>%
  replace_na(list(n = 0))

codon_usage_in_leds <- codons_in_leds %>%
  group_by(region, aa) %>%
  mutate(codon_usage = n/sum(n)) %>%
  ungroup %>%
  left_join(transcript_details)
  

#Get only R codons
r_codon_usage_in_leds <- codon_usage_in_leds %>%
  filter(aa == "R")

#Format for PCA
r_codons_matrix <- r_codon_usage_in_leds %>%
  as.data.frame() %>%
  dplyr::select(region, codon, codon_usage) %>%
  pivot_wider(names_from = codon, values_from = codon_usage) %>%
  column_to_rownames("region")

# Cor matrix -----------------------------------------------------------------

cor_df <- cor(r_codons_matrix, method = "spearman") %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  pivot_longer(cols = -1, names_to = "y", values_to = "cor") %>%
  mutate(cor = case_when(x == y ~ 0,
                         T ~ cor),
         x = fct_relevel(x, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"),
         y = fct_relevel(y, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC") %>% fct_rev())

cor_plot <- ggplot(cor_df, aes(x = x, y = y, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation")

ggsave("GASR/germs/plots/arginine/correlation_matrix.pdf",
       cor_plot,
       device = "pdf", units = "in",
       height = 3, width = 4)

# PCA ---------------------------------------------------------------------

pca_codons <- prcomp(r_codons_matrix, scale. = F)

pca_df <- as.data.frame(pca_codons$x) %>%
  rownames_to_column("region")

pca_df_classes <- pca_df %>%
  mutate(pca_class = cut_number(PC1, 3, labels = c("CG", "middle third", "AG"))) %>%
  left_join(r_codon_usage_in_leds) %>%
  mutate(codon = fct_relevel(codon, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"))

write_tsv(pca_df_classes, "GASR/germs/data/arginine_codon_usage_pca_clusters.tsv.gz")

pca_plot <- ggplot(pca_df_classes %>%
         dplyr::select(PC1, PC2, pca_class) %>% 
         distinct(),
       aes(x = PC1, y = PC2, fill = pca_class)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  theme_classic() +
  labs(fill = "",
       x = paste0("PC1:", round(summary(pca_codons)$importance[2, 1] * 100),"% of variance"),
       y = paste0("PC2:", round(summary(pca_codons)$importance[2, 2] * 100),"% of variance"),) +
  scale_fill_manual(values = c("#1b9e77", "white", "#d95f02")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave("GASR/germs/plots/arginine/pca_plot.pdf",
       pca_plot,
       device = "pdf", units = "in",
       height = 2, width = 5)

# Codon usage in clusters -------------------------------------------------------------

#Codon usage within the classes
codon_usage_boxplot <- pca_df_classes %>%
  filter(pca_class != "middle third") %>%
  ggplot(aes(x = factor(codon), fill = factor(pca_class), y = codon_usage)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "proportion of codons", fill = "codon bias") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

codon_usage_significance <- pca_df_classes %>%
  filter(pca_class != "middle third") %>%
  group_by(codon) %>%
  summarise(ttest = list(t.test(codon_usage ~ pca_class) %>% 
                           tidy)) %>%
  unnest(ttest) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

ggsave("GASR/germs/plots/arginine/codon_usage_boxplot.pdf",
       codon_usage_boxplot,
       device = "pdf", units = "in",
       height = 3, width = 5)


# AA representation -------------------------------------------------------

aa_prop_and_cluster <- r_rich_df %>%
  dplyr::select(c(1, 5:24)) %>%
  pivot_longer(cols = -1, names_to = "aa", values_to = "proportion") %>%
  left_join(
    pca_df_classes %>%
      dplyr::select(region, pca_class, PC1) %>%
      distinct()
  )

filtered_aa_prop_and_cluster <- aa_prop_and_cluster %>%
  filter(pca_class != "middle third") %>%
  filter(aa %in% c("A", "D", "E", "G", "K", "L", "P", "R", "S")) %>%
  mutate(aa = fct_relevel(aa, "R", "D", "E", "K", "S",))

aa_usage_significance <- filtered_aa_prop_and_cluster %>%
  group_by(aa) %>%
  summarise(ttest = list(t.test(proportion ~ pca_class) %>% 
    tidy)) %>%
  unnest(ttest) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

write_tsv(aa_usage_significance, "GASR/germs/plots/arginine/aa_usage_significance.tsv")

aa_usage_boxplot <- ggplot(filtered_aa_prop_and_cluster,
       aes(x = aa, fill = factor(pca_class), y = proportion)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.25) + 
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "proportion of amino acids", fill = "") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave("GASR/germs/plots/arginine/aa_usage_boxplot.pdf",
       aa_usage_boxplot,
       device = "pdf", units = "in",
       height = 3, width = 4)


aa_prop_and_cluster %>%
  mutate(charged = aa %in% c("D", "E", "K", "R")) %>%
  group_by(region, pca_class, PC1) %>%
  summarise(proportion_charged = sum(proportion * charged)) %>%
  ungroup() ->
  prop_charged_and_cluster

prop_charged_and_cluster %>%
  filter(pca_class != "middle third") %>%
  ggplot() +
  aes(x = proportion_charged, fill = pca_class) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "proportion charged\namino acids", fill = "") +
  geom_density(alpha = 0.5) ->
  proportion_charged_plot

prop_charged_and_cluster %>%
  filter(pca_class != "middle third") %>%
  summarise(ttest = t.test(proportion_charged ~ pca_class) %>% tidy()) %>%
  unnest(ttest) ->
  prop_charged_ttest

write_tsv(prop_charged_ttest, "GASR/germs/plots/arginine/prop_charged_ttest.tsv")

ggsave("GASR/germs/plots/arginine/proportion_charged_plot.pdf",
       proportion_charged_plot,
       device = "pdf", units = "in",
       height = 3, width = 3)

# topGO -------------------------------------------------------------------

library(topGO)
library(biomaRt)

pca_class_names <- c("AG", "CG")

all_genes <- unique(transcript_details$gene_name)

#Get GO terms from Ensembl.
all_db <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "https://www.ensembl.org")
all_go_ids <- getBM(attributes = c('go_id', 'external_gene_name', 'namespace_1003'), 
                filters = 'external_gene_name', values = all_genes, mart = all_db)
# build the gene 2 GO annotation list (needed to create topGO object)
all_gene_to_GO <- unstack(all_go_ids[,c(1,2)])

all_genes_go_enrichment <- lapply(pca_class_names, function(class_name){
  
  clust_genes <- unique(pca_df_classes$gene_name[pca_df_classes$pca_class == class_name])
  
  clust_genes_f <- clust_genes[clust_genes %in% all_go_ids[,2]]
  
  geneList <- factor(as.integer(all_genes %in% clust_genes_f))
  names(geneList) <- all_genes
  
  big_go_table <- lapply(c("BP", "MF", "CC"), function(x){
    
    GOdata <- new('topGOdata', 
                  ontology = x, 
                  allGenes = geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = all_gene_to_GO)
    
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
    mutate(fold_enrichment = Significant / Expected)
  
  return(big_go_table)
  
  }) %>%
  set_names(pca_class_names) %>%
  bind_rows(.id = "class")  %>%
  mutate(weightFisher = weightFisher %>% str_replace("< 1e-30", "1e-30") %>% as.numeric())

r_genes <- unique(r_rich_df$gene_name)

r_genes_go_enrichment <- lapply(pca_class_names, function(class_name){
  
  clust_genes <- unique(pca_df_classes$gene_name[pca_df_classes$pca_class == class_name])
  
  clust_genes_f <- clust_genes[clust_genes %in% all_go_ids[,2]]
  
  geneList <- factor(as.integer(r_genes %in% clust_genes_f))
  names(geneList) <- r_genes
  
  big_go_table <- lapply(c("BP", "MF", "CC"), function(x){
    
    GOdata <- new('topGOdata', 
                  ontology = x, 
                  allGenes = geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = all_gene_to_GO)
    
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
    mutate(fold_enrichment = Significant / Expected)
  
  return(big_go_table)
  
}) %>%
  set_names(pca_class_names) %>%
  bind_rows(.id = "class") %>%
  mutate(weightFisher = as.numeric(weightFisher))

selected_terms <- c("nuclear speck", "mRNA binding", "mRNA processing")

go_enrichment_plot <- all_genes_go_enrichment %>%
  mutate(class = class %>% str_replace("bottom third", "CG") %>% 
           str_replace("top third", "AG") %>% fct_relevel("CG")) %>%
  filter(Term %in% selected_terms) %>%
  ggplot(aes(x = Term, y = fold_enrichment, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "fold enrichment over expected", fill = "codon bias") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave("GASR/germs/plots/arginine/go_term_enrichment.pdf",
       go_enrichment_plot,
       device = "pdf", units = "in",
       height = 4, width = 4)


# Analysis of codon usage as a function of GO terms -----------------------

nuclear_speck_genes <- all_go_ids %>%
  filter(go_id == "GO:0016607") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

mrna_binding_genes <- all_go_ids %>%
  filter(go_id == "GO:0003729") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

rna_binding_genes <- all_go_ids %>%
  filter(go_id == "GO:0003723") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

mrna_processing_genes <- all_go_ids %>%
  filter(go_id == "GO:0006397") %>%
  dplyr::select(external_gene_name) %>%
  unlist(use.names = F)

go_annotated_codon_usage <- codon_usage_in_leds %>%
  mutate(nuclear_speck = case_when(gene_name %in% nuclear_speck_genes ~ "nuclear speck",
                                   T ~ "non-nuclear speck") %>%
           fct_relevel("non-nuclear speck"),
         mrna_binding = case_when(gene_name %in% mrna_binding_genes ~ "mRNA binding",
                                  T ~ "non-mRNA binding"),
         rna_binding = case_when(gene_name %in% rna_binding_genes ~ "RNA binding",
                                 T ~ "non-RNA binding"),
         mrna_processing = case_when(gene_name %in% mrna_processing_genes ~ "mRNA processing",
                                     T ~ "non-mRNA processing"),) %>%
  mutate(codon = fct_relevel(codon, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"))

speckle_codon_usage <- go_annotated_codon_usage %>% 
  drop_na() %>%
  group_by(nuclear_speck, codon, aa) %>%
  summarise(codon_usage = mean(codon_usage))

nuclear_speck_codon_plot <- ggplot(go_annotated_codon_usage %>%
         filter(aa == "R"),
       aes(x = codon, y = codon_usage, fill = factor(nuclear_speck))) +
  # facet_wrap(. ~ aa, scales = "free_x") +
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c(scales::muted("red3"), scales::muted("blue3"))) +
  labs(y = "proportion of codons", fill = "") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

rna_binding_codon_plot <- ggplot(go_annotated_codon_usage %>%
         filter(aa == "R"),
       aes(x = codon, y = codon_usage, fill = factor(rna_binding))) +
  # facet_wrap(. ~ aa, scales = "free_x") +
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c(scales::muted("deeppink"), scales::muted("chartreuse"))) +
  labs(y = "proportion of codons", fill = "") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave("GASR/germs/plots/arginine/go_term_codon_biases.pdf",
       nuclear_speck_codon_plot / rna_binding_codon_plot,
       device = "pdf", units = "in",
       height = 6, width = 6)




