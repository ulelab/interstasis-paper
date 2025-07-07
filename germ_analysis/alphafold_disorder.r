rm(list = ls())

library(tidyverse)

# Loading data ------------------------------------------------------------

gencode_to_swissprot <- read_tsv("GASR/lists/gencode.v29.metadata.SwissProt.txt",
                                 col_names = c("transcript_id", "uniprot_id", "id2"),
                                 col_type = "c")

af2_disorder <- read_tsv("General/alpha_fold_disorder/af2_human_disorder_pred.tsv.gz",
                         skip = 1,
                         col_names = c("name", "pos", "aa", "lddt", "disorder", "rsa", "ss", "disorder_25", "binding_25_0581"),
                         col_types = "cdcdddcdd")
  
transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

transcript_details_uniprot <- inner_join(gencode_to_swissprot, transcript_details)

transcript_af2_disorder <- af2_disorder %>%
  group_by(name) %>%
  nest(data = c(pos, aa, lddt, disorder, rsa, ss, disorder_25, binding_25_0581)) %>%
  ungroup() %>%
  mutate(uniprot_id = word(name, 2, sep = fixed("-"))) %>%
  dplyr::select(-name) %>%
  inner_join(transcript_details_uniprot)

#Does the length of the alpha fold vector correctly match the length of the CDS? If not, discard.
transcript_af2_disorder <- transcript_af2_disorder %>%
  ungroup() %>%
  mutate(disorder_nrow = map(data, ~{nrow(.x)}) %>% unlist()) %>%
  filter((disorder_nrow + 1) == (cds_length/3)) %>%
  dplyr::select(-disorder_nrow)

#GeRMS clusters
mv <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz",
                                 col_types = cols()) %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "bullshit"), sep = ":|-|@", remove = F, convert = T) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  mutate(end = end + 4)

# Get disorder values -----------------------------------------------------

median_disorder <- transcript_af2_disorder %>%
  unnest(data) %>% 
  ungroup() %>%
  dplyr::select(disorder_25) %>%
  unlist() %>%
  median()

median_binding <- transcript_af2_disorder %>% 
  unnest(data) %>% 
  ungroup() %>%
  dplyr::select(binding_25_0581) %>%
  unlist() %>%
  median()

median_lddt <- transcript_af2_disorder %>% 
  unnest(data) %>% 
  ungroup() %>%
  dplyr::select(lddt) %>%
  unlist() %>%
  median()

disorder_in_germs_regions <- transcript_af2_disorder %>%
  ungroup() %>%
  inner_join(mv) %>%
  filter(cluster != 0) %>%
  mutate(aa_start = case_when(start - cds_start < 0 ~ 1,
                              T ~ round((start - cds_start)/3) + 1),
         aa_end = case_when(cds_end - end <= 3 ~ (cds_length/3) - 1, #Stop codon should be excluded, hence the 3.
                              T ~ round((end - cds_start)/3))) %>%
  mutate(data_in_region = pmap(list(data, aa_start, aa_end),
                                   function(data, aa_start, aa_end) {
                                     data[aa_start:aa_end,]
                                   }),
         data_not_in_region = pmap(list(data, aa_start, aa_end),
                               function(data, aa_start, aa_end) {
                                 if(aa_start <= 30) {
                                   aa_start <- 1
                                 } else {
                                   aa_start = aa_start - 30
                                 }
                                 
                                 if(aa_end >= (nrow(data) - 30)) {
                                   aa_end <- nrow(data)
                                 } else {
                                   aa_end = aa_end + 30
                                 }
                                 
                                 data[-c(aa_start:aa_end),]
                                 
                               })) %>%
  mutate(median_disorder = map(data_in_region, ~{ median(.x$disorder_25, na.rm = T) }) %>% unlist()) %>%
  mutate(median_binding = map(data_in_region, ~{ median(.x$binding_25_0581, na.rm = T) }) %>% unlist()) %>%
  mutate(median_lddt = map(data_in_region, ~{ median(.x$lddt, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_disorder = map(data_in_region, ~{ mean(.x$disorder_25, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_binding = map(data_in_region, ~{ mean(.x$binding_25_0581, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_lddt = map(data_in_region, ~{ mean(.x$lddt, na.rm = T) }) %>% unlist())
  
disorder_in_germs_regions %>%
  dplyr::select(-data_in_region) %>%
  dplyr::rename(data_in_region = data_not_in_region) %>%
  mutate(median_disorder = map(data_in_region, ~{ median(.x$disorder_25, na.rm = T) }) %>% unlist()) %>%
  mutate(median_binding = map(data_in_region, ~{ median(.x$binding_25_0581, na.rm = T) }) %>% unlist()) %>%
  mutate(median_lddt = map(data_in_region, ~{ median(.x$lddt, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_disorder = map(data_in_region, ~{ mean(.x$disorder_25, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_binding = map(data_in_region, ~{ mean(.x$binding_25_0581, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_lddt = map(data_in_region, ~{ mean(.x$lddt, na.rm = T) }) %>% unlist()) %>%
  mutate(cluster = 0,
         cluster_name = "non-GeRM",
         cluster_colour = "#CCCCCC") ->
  control_disorder 

bind_rows(list(disorder_in_germs_regions %>% select(-data_not_in_region),
               control_disorder)) %>%
  arrange(cluster) ->
  combined_disorder

viewable <- combined_disorder %>%
  dplyr::select(gene_name, uniprot_id, cluster_name, aa_start, aa_end,
                median_disorder, median_binding, median_lddt,
                mean_disorder, mean_binding, mean_lddt)

# write_tsv(viewable, file = "GASR/germs/data/alphafold_disorder/predictions_within_germs_regions.txt")

combined_disorder %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("non-GeRM", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df

# colours_df <- combined_disorder %>%
#   # filter(cluster > 0) %>%
#   dplyr::select(cluster_name, cluster_colour, cluster) %>%
#   distinct() %>%
#   arrange(cluster)

combined_disorder <- combined_disorder %>%
  left_join(cluster_rename_df) %>%
  mutate(cluster_colour = factor(cluster_colour, levels = cluster_rename_df$cluster_colour),
         cluster_name = factor(cluster_name, levels = cluster_rename_df$cluster_name)) %>%
  arrange(cluster_name_short)

lddt_plot <- ggplot(combined_disorder,
       aes(x = mean_lddt , fill = factor(cluster_name))) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ factor(cluster_name), nrow = 4) +
  geom_vline(xintercept = median_lddt, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(x = "median alphafold pLDDT in GeRMS region",
       fill = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = cluster_rename_df$cluster_colour)


lddt_ecdf_plot <- ggplot(combined_disorder,
                    aes(x = mean_lddt , color = factor(cluster_name_short))) +
  stat_ecdf() +
  geom_vline(xintercept = median_lddt, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(x = "mean pLDDT in region (alphafold)",
       y = "cumulative density",
       color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        ) +
  scale_color_manual(values = cluster_rename_df$cluster_colour)

disorder_ecdf_plot <- ggplot(combined_disorder,
                         aes(x = mean_disorder , color = factor(cluster_name_short))) +
  stat_ecdf() +
  geom_vline(xintercept = median_disorder, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(x = "mean disorder in\nencoded region (alphafold)",
       y = "cumulative density",
       color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) +
  scale_color_manual(values = cluster_rename_df$cluster_colour)

ggsave("GASR/germs/plots/ecdf_alphafold_disorder_germs.pdf", disorder_ecdf_plot,
       device = "pdf", units = "in",
       height = 3, width = 4.5)

ggsave("GASR/germs/plots/ecdf_alphafold_lddt_germs.pdf", lddt_ecdf_plot,
       device = "pdf", units = "in",
       height = 3, width = 4.5)

# Stats -------------------------------------------------------------------

cluster_rename_df %>% 
  filter(new_order != 0) %>% 
  select(cluster_name_short) %>% 
  unlist(use.names = F) ->
  testing_clusters

testing_clusters %>%
  lapply(., function(x) {
    germ_ent <- combined_disorder$mean_disorder[combined_disorder$cluster_name_short == x]
    nongerm_ent <-combined_disorder$mean_disorder[combined_disorder$cluster_name_short == "non-GeRM"]
    
    return(t.test(nongerm_ent, germ_ent) %>% tidy())
  }) %>%
  set_names(testing_clusters) %>%
  bind_rows(.id = "cluster_name") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  disord_ttests

write_tsv(disord_ttests, "GASR/germs/data/stats/disorder_in_germ_regions.tsv")

testing_clusters %>%
  lapply(., function(x) {
    germ_ent <- combined_disorder$mean_lddt[combined_disorder$cluster_name_short == x]
    nongerm_ent <-combined_disorder$mean_lddt[combined_disorder$cluster_name_short == "non-GeRM"]
    
    return(t.test(nongerm_ent, germ_ent) %>% tidy())
  }) %>%
  set_names(testing_clusters) %>%
  bind_rows(.id = "cluster_name") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  lddt_ttests

write_tsv(lddt_ttests, "GASR/germs/data/stats/lddt_in_germ_regions.tsv")

# Plots ignoring classes --------------------------------------------------

combined_disorder %>%
  mutate(is_germ = case_when(new_order == 0 ~ "low GeRM",
                             T ~ "high GeRM") %>%
           fct_relevel("high GeRM")) ->
  combined_disorder

mergedcat_lddt_ecdf_plot <- ggplot(combined_disorder,
                         aes(x = mean_lddt , color = is_germ)) +
  stat_ecdf() +
  geom_vline(xintercept = median_lddt, linetype = "dashed", alpha = 0.8) +
  theme_classic() +
  labs(x = "mean pLDDT in region (alphafold)",
       y = "cumulative density",
       color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) +
  scale_color_manual(values = c("black", "grey70"))

mergedcat_disorder_ecdf_plot <- ggplot(combined_disorder,
                                   aes(x = mean_disorder , color = is_germ)) +
  stat_ecdf() +
  geom_vline(xintercept = median_disorder, linetype = "dashed", alpha = 0.8) +
  theme_classic() +
  labs(x = "mean solvent accessibility\nin region (alphafold)",
       y = "cumulative density",
       color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) +
  scale_color_manual(values = c("black", "grey70"))

ggsave("GASR/germs/plots/ecdf_alphafold_lddt_all_germs.pdf", mergedcat_lddt_ecdf_plot,
       device = "pdf", units = "in",
       height = 2, width = 3.5)

ggsave("GASR/germs/plots/ecdf_alphafold_disorder_all_germs.pdf", mergedcat_disorder_ecdf_plot,
       device = "pdf", units = "in",
       height = 2, width = 3.5)

combined_disorder %>%
  summarise(ttest = t.test(median_lddt ~ is_germ, data = .) %>% tidy()) %>%
  unnest(ttest) ->
  mergedcat_lddt_ttest

combined_disorder %>%
  summarise(ttest = t.test(median_lddt ~ is_germ, data = .) %>% tidy()) %>%
  unnest(ttest) ->
  mergedcat_disorder_ttest

write_tsv(mergedcat_lddt_ttest, "GASR/germs/plots/alphafold_lddt_all_germs_ttests.tsv")
write_tsv(mergedcat_disorder_ttest, "GASR/germs/plots/alphafold_disorder_all_germs_ttests.tsv")
