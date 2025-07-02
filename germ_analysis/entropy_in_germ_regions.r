rm(list = ls())

library(tidyverse)
library(broom)

# Loading data ------------------------------------------------------------

protein_entropy <- read_rds(file = "GASR/germs/data/all_protein_entropy.rds.gz")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

transcript_entropy <- protein_entropy %>%
  inner_join(transcript_details) %>%
  mutate(entropy_nrow = map(entropy_vector, 
                            ~{length(.x)}) %>%
           unlist()) %>%
  # Does the length of the entropy vector correctly match the length of the CDS? If not, discard.
  filter((entropy_nrow + 42) == (cds_length / 3)) %>% # Length + smoothing window (41) + 1 (stop codon)
  dplyr::select(-entropy_nrow)


#GeRM clusters
mv <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz",
               col_types = cols()) %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "bullshit"), sep = ":|-|@", remove = F, convert = T) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  mutate(end = end + 4)

# Get entropy values -----------------------------------------------------

median_entropy <- transcript_entropy %>%
  dplyr::select(entropy_vector) %>%
  unlist() %>%
  median()

p10_entropy <- transcript_entropy %>%
  dplyr::select(entropy_vector) %>%
  unlist() %>%
  quantile(0.1)

entropy_in_germs_regions <- transcript_entropy %>%
  ungroup() %>%
  inner_join(mv) %>%
  filter(cluster != 0) %>%
  mutate(aa_start = case_when(start - cds_start < 0 ~ 1,
                              T ~ round((start - cds_start)/3) + 1),
         aa_end = case_when(cds_end - end <= 3 ~ n_aa, 
                            # Subtract half of the smoothing window - seems fair to me
                            T ~ round((end - cds_start)/3) - 20)) %>%
  mutate(data_in_region = pmap(list(entropy_vector, aa_start, aa_end),
                               function(entropy_vector, aa_start, aa_end) {
                                 entropy_vector[aa_start:aa_end]
                               }),
         data_not_in_region = pmap(list(entropy_vector, aa_start, aa_end),
                                   function(entropy_vector, aa_start, aa_end) {
                                     if(aa_start <= 30) {
                                       aa_start <- 1
                                     } else {
                                       aa_start = aa_start - 30
                                     }
                                     
                                     if(aa_end >= (length(entropy_vector) - 30)) {
                                       aa_end <- length(entropy_vector)
                                     } else {
                                       aa_end = aa_end + 30
                                     }
                                     
                                     entropy_vector[-c(aa_start:aa_end)]
                                     
                                   })) %>%
  mutate(median_entropy = map(data_in_region, ~{ median(.x, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_entropy = map(data_in_region, ~{ mean(.x, na.rm = T) }) %>% unlist())

entropy_in_germs_regions %>%
  dplyr::select(-data_in_region) %>%
  dplyr::rename(data_in_region = data_not_in_region) %>%
  mutate(median_entropy = map(data_in_region, ~{ median(.x, na.rm = T) }) %>% unlist()) %>%
  mutate(mean_entropy = map(data_in_region, ~{ mean(.x, na.rm = T) }) %>% unlist()) %>%
  mutate(cluster = 0,
         cluster_name = "non-GeRM",
         cluster_colour = "#CCCCCC") ->
  control_entropy

bind_rows(list(entropy_in_germs_regions %>% select(-data_not_in_region),
               control_entropy)) %>%
  arrange(cluster) ->
  combined_entropy

viewable <- combined_entropy %>%
  dplyr::select(gene_name, cluster_name, aa_start, aa_end,
                median_entropy, mean_entropy)

# write_tsv(viewable, file = "GASR/germs/data/alphafold_disorder/predictions_within_germs_regions.txt")

combined_entropy %>%
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

combined_entropy <- combined_entropy %>%
  left_join(cluster_rename_df) %>%
  mutate(cluster_colour = factor(cluster_colour, levels = cluster_rename_df$cluster_colour),
         cluster_name = factor(cluster_name, levels = cluster_rename_df$cluster_name)) %>%
  arrange(cluster_name_short)

entropy_plot <- ggplot(combined_entropy,
                    aes(x = mean_entropy , fill = factor(cluster_name))) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ factor(cluster_name), nrow = 4) +
  geom_vline(xintercept = median_entropy, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(x = "mean aa entropy in GeRM region",
       fill = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = cluster_rename_df$cluster_colour)

mean_entropy_ecdf_plot <- ggplot(combined_entropy,
                         aes(x = mean_entropy , color = factor(cluster_name_short))) +
  stat_ecdf() +
  geom_vline(xintercept = median_entropy, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(xlim = c(0.4,0.9)) +
  theme_classic() +
  labs(x = "mean aa entropy in region",
       y = "cumulative density",
       color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) +
  scale_color_manual(values = cluster_rename_df$cluster_colour)

ggsave("GASR/germs/plots/ecdf_entropy_germs.pdf", mean_entropy_ecdf_plot,
       device = "pdf", units = "in",
       height = 3, width = 4.5)

ggplot(combined_entropy %>%
         mutate(is_germ = case_when(new_order == 0 ~ "low GeRM",
                                    T ~ "high GeRM") %>%
                  fct_relevel("high GeRM")),
       aes(x = mean_entropy, color = is_germ)) +
  stat_ecdf() +
  geom_vline(xintercept = median_entropy, linetype = "dashed", alpha = 0.8) +
  geom_vline(xintercept = p10_entropy, linetype = "dashed", alpha = 0.8, color = "red3") +
  coord_cartesian(xlim = c(0.4,0.9)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic() +
  labs(x = "mean aa entropy in region",
       y = "cumulative density",
       color = "") +
  scale_color_manual(values = c("black", "grey70")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) ->
  all_germ_entropy_edcf_plot

ggsave("GASR/germs/plots/ecdf_entropy_all_germs.pdf", all_germ_entropy_edcf_plot,
       device = "pdf", units = "in",
       height = 2, width = 3.5)

# Stats -------------------------------------------------------------------

cluster_rename_df %>% 
  filter(new_order != 0) %>% 
  select(cluster_name_short) %>% 
  unlist(use.names = F) ->
  testing_clusters

testing_clusters %>%
  lapply(., function(x) {
    germ_ent <- combined_entropy$mean_entropy[combined_entropy$cluster_name_short == x]
    nongerm_ent <- combined_entropy$mean_entropy[combined_entropy$cluster_name_short == "non-GeRM"]
    
    return(t.test(nongerm_ent, germ_ent) %>% tidy())
    }) %>%
  set_names(testing_clusters) %>%
  bind_rows(.id = "cluster_name") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  ent_ttests

write_tsv(ent_ttests, "GASR/germs/data/stats/entropy_in_germ_regions.tsv")

combined_entropy %>%
  mutate(is_germ = case_when(new_order == 0 ~ "low GeRM",
                             T ~ "high GeRM") %>%
           fct_relevel("high GeRM")) %>%
  summarise(ttest = t.test(mean_entropy ~ is_germ, data = .) %>% tidy()) %>%
  unnest(ttest) ->
  all_germ_entropy_ttest

write_tsv(all_germ_entropy_ttest, "GASR/germs/data/stats/entropy_in_all_germ_regions.tsv")

