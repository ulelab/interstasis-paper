rm(list = ls())

library(tidyverse)
library(topGO)
library(biomaRt)
library(factoextra)
library(cluster)

# Load data ---------------------------------------------------------------

gene_info <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt")
all_genes <- gene_info$gene_name

germs_clusters <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "region"), sep = ":|-|@", remove = F, convert = T) %>%
  left_join(gene_info)

# GO analysis  ------------------------------------------------------------

#Get GO terms from Ensembl.
db <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "https://www.ensembl.org")
go_ids <- getBM(attributes = c('go_id', 'external_gene_name', 'namespace_1003'), 
                filters = 'external_gene_name', values = all_genes, mart = db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO <- unstack(go_ids[,c(1,2)])


# GeRMS clusters ----------------------------------------------------------

all_clust_go_weight <- lapply(unique(germs_clusters$cluster),
                              function(clust_name){
                                # clust_name = 2
                                
                                clust_genes <- germs_clusters$gene_name[germs_clusters$cluster == clust_name] %>%
                                  unique()
                                
                                clust_genes_f <- clust_genes[clust_genes %in% go_ids[,2]]
                                
                                geneList <- factor(as.integer(all_genes %in% clust_genes_f))
                                names(geneList) <- all_genes
                                
                                big_go_table <- lapply(c("BP", "MF", "CC"), function(x){
                                  
                                  # x = "BP"
                                  
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
                                  mutate(fold_enrichment = Significant / Expected)
                                
                                return(big_go_table)
                                
                              }) %>%
  set_names(unique(germs_clusters$cluster) %>% paste0("cluster_", .)) %>%
  bind_rows(.id = "cluster")

all_clust_go_weight <- all_clust_go_weight %>%
  mutate(weightFisher = weightFisher %>% str_replace("< 1e-30", "1e-30") %>% as.numeric())

# all_clust_go_weight <- all_clust_go_weight %>% dplyr::select(-cluster_name)

all_clust_go_weight <- left_join(all_clust_go_weight,
                                 germs_clusters %>%
                                   mutate(cluster = paste0("cluster_", cluster)) %>%
                                   dplyr::select(cluster, cluster_name) %>%
                                   distinct())

write_tsv(all_clust_go_weight, file = "GASR/germs/data/germs_cds_clust_full_go_terms.tsv.gz")

# Analyse GO info ---------------------------------------------------------

all_clust_go_weight <- read_tsv(file = "GASR/germs/data/germs_cds_clust_full_go_terms.tsv.gz")

write_tsv(all_clust_go_weight %>%
            group_by(cluster_name) %>% 
            mutate(padj = p.adjust(weightFisher)) %>%
            ungroup() %>% filter(padj < 0.05) %>%
            dplyr::select(cluster_name, Term, go_type, fold_enrichment, padj, Significant, Expected), "GASR/germs/data/germs_cds_clust_significant_go_terms.tsv")

filtered_go <- all_clust_go_weight %>%
  filter(cluster != "cluster_0") %>%
  # replace_na(list(weightFisher = 1e-30)) %>%
  group_by(cluster, go_type) %>%
  mutate(padj = p.adjust(weightFisher, "BH")) %>%
  filter(padj < 0.05) %>%
  arrange(cluster, go_type, padj) %>%
  filter(Annotated > 10) %>%
  dplyr::slice(1:4)

unique_term_num <- unique(filtered_go$Term) %>% length()

unique_terms_all_clust <- filtered_go %>%
  ungroup() %>%
  dplyr::select(GO.ID, Term, go_type) %>%
  distinct() %>%
  left_join(all_clust_go_weight) %>%
  # mutate(Fis = Fis %>% as.numeric()) %>%
  # replace_na(list(Fis = 1e-30)) %>%
  mutate(padj = p.adjust(weightFisher, "BH")) %>%
  mutate(log_signed_p = -log2(padj),
         log_fold = log2(fold_enrichment)) %>%
  mutate(log_fold_capped = case_when(log_fold < 0 ~ 0,
                                     log_fold > 4 ~ 4,
                                     T ~ log_fold),
         fold_capped = case_when(fold_enrichment > 10 ~ 10, fold_enrichment < 1 ~ 1, T ~ fold_enrichment)) %>%
  mutate(log_signed_p_capped = case_when(log_signed_p >= 30 ~ 30,
                                         T ~ log_signed_p)) %>%
  filter(cluster != "cluster_0")

wide_by_clust <- unique_terms_all_clust %>%
  dplyr::select(cluster, Term, log_signed_p) %>%
  pivot_wider(names_from = "cluster",
              values_from = "log_signed_p") %>%
  column_to_rownames("Term")

cluster_of_clusters <- eclust(wide_by_clust %>% t(), "hclust", k.max = 5)
fviz_dend(cluster_of_clusters) 

cluster_of_terms <-  eclust(wide_by_clust, "hclust")
fviz_dend(cluster_of_terms)

cluster_order <- germs_clusters$cluster_name %>% unique() %>% str_subset("No", negate = T)
cluster_order <- cluster_order[c(2,1,7,8,3,4,6,5)]

# unique_terms_all_clust$Term <- unique_terms_all_clust$Term %>% fct_relevel(cluster_of_terms$labels[cluster_of_terms$order])
unique_terms_all_clust$Term <- unique_terms_all_clust$Term %>% fct_relevel(filtered_go$Term %>% unique())

unique_terms_all_clust$cluster_name <- unique_terms_all_clust$cluster_name %>% fct_relevel(cluster_order)

ggplot(unique_terms_all_clust, aes(x = cluster_name, y = Term, fill = log_signed_p_capped)) +
  geom_tile() +
  labs(x = "",
       y= "",
       fill = "-log2 p. adj") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "top") +
  scale_fill_viridis_c(option = "magma") +
  # scale_fill_gradient2(
  #   low = "blue",
  #   mid = "white",
  #   high = "red",
  #   midpoint = 0) +
  facet_grid(go_type~ ., scales = "free_y")

size_heatmap <- ggplot(unique_terms_all_clust, 
                       aes(x = cluster_name, y = Term, 
                           fill = log_fold_capped, size = log_signed_p_capped)) +
  geom_point(shape = 21, stroke = 0.6) +
  scale_radius(range = c(0.1, 10)) +
  labs(x = "",
       y= "",
       fill = "log2 fold enrichment",
       size = "-log2 p. adj.") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.box = "vertical", legend.margin = margin(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  facet_grid(go_type~ ., scales = "free_y")

ggsave(filename = "GASR/germs/plots/germs_CDS_gene_ontology_heatmap.pdf",
       size_heatmap,
       device = "pdf", units = "in",
       height = 10, width = 8.5)



# mini heatmap ------------------------------------------------------------

germs_clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df

my_favourite_things <- c(
  "nuclear speck",
  "nucleolus",
  "chromatin",
  "RNA binding",
  "RNA splicing",
  "mRNA processing",
  "histone binding",
  # "transcription coactivator activity",
  # "transcription corepressor activity",
  "positive regulation of transcription by RNA polymerase II",
  "SH3 domain binding"
)

all_clust_go_weight %>%
  mutate(cluster = cluster %>% str_remove("cluster_") %>% as.numeric()) %>%
  filter(Term %in% my_favourite_things) %>%
  mutate(Term = Term %>%
           str_replace("positive regulation of transcription by RNA polymerase II", "pos. reg. of txn")) %>%
  group_by(cluster_name) %>%
  mutate(padj = p.adjust(weightFisher, "BH")) %>%
  mutate(log_signed_p = -log2(padj),
         log_fold = log2(fold_enrichment)) %>%
  mutate(log_fold_capped = case_when(log_fold < 0 ~ 0,
                                     log_fold > 4 ~ 4,
                                     T ~ log_fold)) %>%
  mutate(log_signed_p_capped = case_when(log_signed_p >= 20 ~ 20,
                                         T ~ log_signed_p)) %>%
  left_join(cluster_rename_df) %>%
  filter(!(cluster_name_short %in% c("CUG-rep", "CU-rich", "None", "Pur. + C"))) %>%
  ggplot() +
       aes(x = cluster_name_short, y = Term, 
           fill = log_fold_capped, size = log_signed_p_capped) +
  geom_point(shape = 21, stroke = 0.6) +
  scale_radius(range = c(0.1, 10)) +
  labs(x = "",
       y= "",
       fill = "log2 fold enrichment",
       size = "-log2 p. adj.") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.9, hjust = 0.5, colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.box = "vertical", legend.margin = margin(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  facet_grid(go_type~ ., scales = "free_y") ->
  mini_gene_ontology_plot

ggsave(file = "GASR/germs/plots/mini_gene_ontology_plot.pdf", mini_gene_ontology_plot,
       width = 4, height = 5)
  

# midi heatmap ------------------------------------------------------------

germs_clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df

my_favourite_things_expanded <- c(
  # GA/AG
  "nuclear speck",
  "nucleolus",
  "RNA binding",
  "RNA splicing",
  "mRNA processing",
  "histone binding",
  "nucleosome assembly",
  "chromatin remodeling",
  "chromatin binding",
  # GC
  "anterior/posterior pattern specification",
  # "neuron fate commitment",
  "neuromuscular process controlling balance",
  # "transcription regulator complex",
  "chromatin",
  "DNA-binding transcription factor activity, RNA polymerase II-specific",
  "negative regulation of transcription by RNA polymerase II",
  "positive regulation of transcription by RNA polymerase II",
  # CAG
  # "nuclear matrix",
  # CC
  "actin binding",
  "SH3 domain binding",
  # CT
  "keratin filament"
)

all_clust_go_weight %>%
  filter(cluster_name != "No Cluster") %>%
  group_by(cluster) %>%
  mutate(padj = p.adjust(weightFisher, method = "BH")) %>%
  ungroup() %>%
  filter(Term %in% my_favourite_things_expanded) %>%
  mutate(Term = Term %>% fct_relevel(my_favourite_things_expanded),
         log_signed_p = -log2(padj),
         log_fold = log2(fold_enrichment)) %>%
  dplyr::select(-cluster) %>%
  mutate(log_fold_capped = case_when(log_fold < 0 ~ 0,
                                     log_fold > 4 ~ 4,
                                     T ~ log_fold),
         log_signed_p_capped = case_when(log_signed_p >= 20 ~ 20,
                                         log_signed_p <= -log2(0.05) ~ 0,
                                         T ~ log_signed_p)) %>%
  left_join(cluster_rename_df) %>%
  arrange(Term) ->
  go_plotting_df
  

ggplot(go_plotting_df %>%
         filter(!(cluster_name_short %in% c("CUG-rep", "Pur. + C")))) +
  aes(x = cluster_name_short, y = Term, 
      fill = log_fold_capped, size = log_signed_p_capped) +
  geom_point(shape = 21, stroke = 0.6) +
  scale_radius(range = c(0.1, 10)) +
  labs(x = "",
       y= "",
       fill = "log2 fold enrichment",
       size = "-log2 p. adj.") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.box = "vertical", legend.margin = margin(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  scale_fill_viridis_c(option = "magma", direction = -1, 
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", 
                                              ticks.linewidth = 0.5, 
                                              frame.linewidth = 0.5)) ->
  midi_gene_ontology_plot

ggsave(file = "GASR/germs/plots/midi_gene_ontology_plot.pdf", midi_gene_ontology_plot,
       width = 6, height = 6)
  
  
  
  
  
  
  
  

