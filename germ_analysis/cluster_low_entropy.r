library(tidyverse)
library(umap)
library(dbscan)

# Load --------------------------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

aa_per_region_df <- read_tsv("GASR/sequence_multivalency/data/20220315_123win_98per_k5_10nei/data/percentile_98_windowsize_123_klen_5_le_regions_with_aa_props.txt.gz")

aa_per_region_matrix <- aa_per_region_df[c(1,5:24)] %>%
  column_to_rownames("region") %>%
  as.matrix()

# Map ---------------------------------------------------------------------

region_umap <- umap(aa_per_region_matrix, n_components = 4, n_neighbors = 50, min_dist = 0.001)
region_umap_2d <- umap(aa_per_region_matrix, n_components = 2, n_neighbors = 50, min_dist = 0.001)

region_umap_df_2d <- region_umap_2d$layout %>% 
  as.data.frame() %>%
  rownames_to_column("peak_identifier")

region_umap_df <- region_umap$layout %>% 
  as.data.frame() %>%
  rownames_to_column("peak_identifier")

# Cluster -----------------------------------------------------------------

region_optics <- optics(region_umap_df %>%
                          dplyr::select(-region_region_id), 
                        minPts = nrow(region_umap_df)/100)

# What's the maximum reachability distance? Take the 99th percentile to deal with instances where there are a few very distant points.
maximum_reach <- quantile(region_optics$reachdist, .99)

sorted_reach <- region_optics$reachdist[is.finite(region_optics$reachdist)] %>% sort()

theshold_tests <- seq(0, maximum_reach, length.out = 500)

above_threshold <- lapply(theshold_tests,
                          function(x) {
                            thresh_diff <- sorted_reach - x
                            thresh_diff[thresh_diff > 0] %>% 
                              sum() %>%
                              return()
                          }) %>%
  unlist()

plot(above_threshold)

cluster_threshold <- KneeArrower::findCutoff(above_threshold, theshold_tests, method = "first", 0.15)$y

# Define the clusters using the threshold.
region_hdb <- extractDBSCAN(region_optics, cluster_threshold)
# region_hdb <- extractDBSCAN(region_optics, 0.42)

plot(region_hdb)

region_cluster_identity <- as.data.frame(aa_per_region_matrix[region_umap_df$peak_identifier,]) %>%
  rownames_to_column("peak_identifier") %>%
  mutate(cluster = region_hdb$cluster) %>%
  group_by(cluster) %>%
  mutate(group_size = dplyr::n()) %>%
  ungroup() %>%
  mutate(cluster = case_when((group_size < 30) ~ as.numeric(0),
                             (group_size >= 30) ~ as.numeric(cluster))) %>%
  select(-group_size)

region_cluster_identity_long <- region_cluster_identity %>%
  pivot_longer(cols = !c(peak_identifier, cluster), names_to = "kmer", values_to = "prop_total_mv")

cluster_averages <- region_cluster_identity_long %>%
  group_by(cluster, kmer) %>%
  summarise(average_mv = mean(prop_total_mv))


#We can only get 12 colours out of carto_pal
number_of_clusters <- region_cluster_identity$cluster %>% unique() %>% length()

if(number_of_clusters > 13) {
  cluster_colours <- c("#000000",
                       rcartocolor::carto_pal(12, "Prism"),
                       sample(
                         viridis(number_of_clusters - 13)
                       ))
} else {
  cluster_colours <- c("#000000",
                       sample(
                         rcartocolor::carto_pal(number_of_clusters - 1, "Prism")
                       ))
}

set.seed(100)

cluster_names <- cluster_averages %>% 
  arrange(cluster, desc(average_mv)) %>% 
  group_by(cluster) %>% dplyr::slice(1:3) %>%
  summarise(cluster_name = paste(kmer, collapse = "/")) %>%
  mutate(cluster_name = case_when(cluster == 0 ~ "No Cluster",
                                  T ~ cluster_name)) %>%
  mutate(cluster_name = case_when(duplicated(cluster_name) ~ paste0(cluster_name, " 2"),
                                  T ~ cluster_name)) %>%
  mutate(cluster_colour = cluster_colours)



named_region_umap_df <- region_umap_df %>%
  mutate(cluster = region_cluster_identity$cluster) %>%
  left_join(cluster_names, by = "cluster") %>%
  mutate(cluster_colour = cluster_colour %>% fct_relevel(cluster_names$cluster_colour)) %>%
  mutate(cluster_name = cluster_name %>% fct_relevel(cluster_names$cluster_name)) %>%
  dplyr::select(-c(V1, V2, V3, V4)) %>%
  left_join(region_umap_df_2d) %>%
  arrange(cluster)

toc()

# UMAP plot with clusters.
umap_plot <- ggplot(named_region_umap_df %>%
                      filter(cluster > 0), 
                    aes(x = V1, y = V2, fill = cluster_colour)) +  
  geom_point(size = 2, colour = alpha("black", 0.5), pch = 21) +
  scale_fill_identity(guide = "legend", 
                      labels = named_region_umap_df$cluster_name,
                      breaks = named_region_umap_df$cluster_colour) +
  guides(fill = guide_legend(ncol = 1, 
                             title = element_blank(),
                             override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
  )

region_hdb_df <- tibble(cluster = region_hdb$cluster[region_hdb$order],
                        reachability = region_hdb$reachdist[region_hdb$order],
                        order = c(1:length(region_hdb$reachdist))) %>%
  left_join(cluster_names, by = "cluster") %>%
  replace_na(list(cluster_names = cluster_names$cluster_name[1],
                  cluster_colour = cluster_names$cluster_colour[1]))

hdb_plot <- ggplot(region_hdb_df, aes(x = order, y = reachability, fill = factor(cluster_colour))) + 
  geom_col(width = 1.05) +
  scale_fill_identity(guide = "legend",
                      labels = region_hdb_df$cluster_name,
                      breaks = region_hdb_df$cluster_colour) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = cluster_threshold, linetype = "dashed", alpha = 0.5)

composite_plot <- (umap_plot / hdb_plot) + plot_layout(heights = c(2, 1), guides = "collect")


ggsave(filename = "GASR/germs/plots/low_entropy_umap_plot.pdf", 
       composite_plot,
       device = "pdf", units = "in",
       height = 6, width = 6)

write_tsv(named_region_umap_df, file = "GASR/germs/data/low_entropy_clusters.tsv.gz")

