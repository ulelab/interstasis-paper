rm(list = ls())

library(tidyverse)


# Load --------------------------------------------------------------------

clusters <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")
germ_scores <- read_tsv("GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")


clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  arrange(cluster) %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df
  

clusters %>%
  separate(peak_identifier, 
           into = c("transcript_id", "peak_start", 
                    "peak_end", "region"),
           sep = ":|-|@",
           remove = F,
           convert = T) %>%
  left_join(germ_scores) %>%
  filter(cluster > 0) %>%
  group_by(cluster_name, kmer) %>%
  summarise(sum_mv = sum(kmer_multivalency)) %>%
  group_by(cluster_name) %>%
  mutate(prop_mv = sum_mv/sum(sum_mv)) %>%
  group_by(cluster_name, kmer) %>%
  summarise(mean_prop_mv = mean(prop_mv)) %>%
  ungroup() %>%
  left_join(cluster_rename_df) ->
  kmer_mv_props_per_cluster

# clusters %>%
#   separate(peak_identifier, 
#            into = c("transcript_id", "peak_start", 
#                     "peak_end", "region"),
#            sep = ":|-|@",
#            remove = F,
#            convert = T) %>%
#   left_join(germ_scores) %>%
#   filter(cluster > 0) %>%
#   group_by(peak_identifier, cluster_name, kmer) %>%
#   summarise(sum_mv = sum(kmer_multivalency)) %>%
#   group_by(peak_identifier, cluster_name) %>%
#   mutate(prop_mv = sum_mv/sum(sum_mv)) ->
#   whaaat

write_tsv(kmer_mv_props_per_cluster, "GASR/germs/data/prop_mv_per_kmer_cds.tsv.gz")

n_kmers <- 3

kmer_mv_props_per_cluster %>%
  group_by(cluster_name) %>%
  mutate(rank = rank(-mean_prop_mv)) %>%
  ungroup() %>%
  group_by(kmer) %>%
  mutate(min_rank = min(rank)) %>%
  filter(min_rank <= n_kmers) %>%
  ungroup() %>%
  arrange(new_order, rank) %>%
  complete(cluster_name_short, kmer) %>% 
  filter(cluster_name_short != "None") %>%
  replace_na(list(mean_prop_mv = 0)) -> table_for_heatmap
  

table_for_heatmap %>%
  group_by(cluster_name_short) %>%
  filter(rank <= n_kmers) %>%
  ungroup() %>%
  arrange(cluster_name_short, rank) %>%
  select(kmer) %>%
  unique() %>%
  unlist(use.names = T) ->
  kmer_order

table_for_heatmap %>%
  select(cluster_name_short, kmer, mean_prop_mv) %>%
  mutate(mean_prop_mv = case_when(mean_prop_mv >= 0.08 ~ 0.08,
                                  T ~ mean_prop_mv)) %>%
# Cap at .08, then modify the figure legend in illustrator to show this
  mutate(kmer = kmer %>% fct_relevel(kmer_order) %>% fct_relevel("AGCAG", after = Inf)) %>%
  ggplot() +
  aes(x = cluster_name_short, y = kmer, fill = mean_prop_mv) +
  geom_tile(color = "black", linewidth = 0.5) +
  # scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = "red3",
                      breaks = c(0, .04, .08)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, vjust = 1, hjust = 1),
        legend.position = "top") +
  labs(x = "",
       y = "",
       fill = "proportion GeRM in cluster") ->
  heatmap

ggsave("GASR/germs/plots/heatmap_of_most_common_kmers_in_cds_clusters.pdf", heatmap,
       width = 3, height = 6)  



# Load --------------------------------------------------------------------

# Rerunning all with mouse shit

clusters <- read_tsv("GASR/germs/data/mouse/germs_CDS_umap_clusters.tsv.gz")
germ_scores <- read_tsv("GASR/germs/data/mouse/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")


clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  arrange(cluster) %>%
  mutate(cluster_name_short = c("None", "C-rich", "GC-rich",
                                "CUG-rep", "CAG-rep", "GA-rich")) %>%
  mutate(new_order = c(0, 4, 2, 3, 5, 1)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df


clusters %>%
  separate(peak_identifier, 
           into = c("transcript_id", "peak_start", 
                    "peak_end", "region"),
           sep = ":|-|@",
           remove = F,
           convert = T) %>%
  left_join(germ_scores) %>%
  filter(cluster > 0) %>%
  group_by(cluster_name, kmer) %>%
  summarise(sum_mv = sum(kmer_multivalency)) %>%
  group_by(cluster_name) %>%
  mutate(prop_mv = sum_mv/sum(sum_mv)) %>%
  group_by(cluster_name, kmer) %>%
  summarise(mean_prop_mv = mean(prop_mv)) %>%
  ungroup() %>%
  left_join(cluster_rename_df) ->
  kmer_mv_props_per_cluster

# clusters %>%
#   separate(peak_identifier, 
#            into = c("transcript_id", "peak_start", 
#                     "peak_end", "region"),
#            sep = ":|-|@",
#            remove = F,
#            convert = T) %>%
#   left_join(germ_scores) %>%
#   filter(cluster > 0) %>%
#   group_by(peak_identifier, cluster_name, kmer) %>%
#   summarise(sum_mv = sum(kmer_multivalency)) %>%
#   group_by(peak_identifier, cluster_name) %>%
#   mutate(prop_mv = sum_mv/sum(sum_mv)) ->
#   whaaat

write_tsv(kmer_mv_props_per_cluster, "GASR/germs/data/mouse/prop_mv_per_kmer_cds.tsv.gz")

n_kmers <- 3

kmer_mv_props_per_cluster %>%
  group_by(cluster_name) %>%
  mutate(rank = rank(-mean_prop_mv)) %>%
  ungroup() %>%
  group_by(kmer) %>%
  mutate(min_rank = min(rank)) %>%
  filter(min_rank <= n_kmers) %>%
  ungroup() %>%
  arrange(new_order, rank) %>%
  complete(cluster_name_short, kmer) %>% 
  filter(cluster_name_short != "None") %>%
  replace_na(list(mean_prop_mv = 0)) ->
  table_for_heatmap


table_for_heatmap %>%
  group_by(cluster_name_short) %>%
  filter(rank <= n_kmers) %>%
  ungroup() %>%
  arrange(cluster_name_short, rank) %>%
  select(kmer) %>%
  unique() %>%
  unlist(use.names = T) ->
  kmer_order

table_for_heatmap %>%
  select(cluster_name_short, kmer, mean_prop_mv) %>%
  mutate(mean_prop_mv = case_when(mean_prop_mv >= 0.08 ~ 0.08,
                                  T ~ mean_prop_mv)) %>%
  # Cap at .08, then modify the figure legend in illustrator to show this
  mutate(kmer = kmer %>% fct_relevel(kmer_order) %>% fct_relevel("AGCAG", after = Inf)) %>%
  ggplot() +
  aes(x = cluster_name_short, y = kmer, fill = mean_prop_mv) +
  geom_tile(color = "black", linewidth = 0.5) +
  # scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = "red3",
                      breaks = c(0, .04, .08), 
                      guide = guide_colorbar(frame.colour = "black", 
                                             ticks.colour = "black", 
                                             frame.linewidth = 0.5, ticks.linewidth = 0.25)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, vjust = 1, hjust = 1),
        legend.position = "top") +
  labs(x = "",
       y = "",
       fill = "proportion GeRM in cluster") ->
  heatmap

ggsave("GASR/germs/plots/mouse/heatmap_of_most_common_kmers_in_cds_clusters.pdf", heatmap,
       width = 3, height = 5)  
