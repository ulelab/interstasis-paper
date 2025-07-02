list.of.packages <- c("germs", "Biostrings", "tidyverse", "parallel", "tictoc",
                      "patchwork", "dbscan", "umap", "viridis") 

for(i in list.of.packages){
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}

rm(list = ls())
setwd("/camp/lab/ulej/home/users/farawar/")

# Parameters --------------------------------------------------------------

k_length <- 5
window_size <- 123
smoothing_size <- 123
lambda <- 1
scaling_function <- function(x) exp(lambda * x)
threshold_percentile <- 0.98

# Fetch transcripts -------------------------------------------------------

sequences <- readDNAStringSet("genomes/hs/fasta/longest_gencode29.fa") %>% as.character()

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_length - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

sequences <- sequences[match(transcript_details$transcript_id, names(sequences))] #Reorder sequences to match transcript detail table

# Make HDM and PDV --------------------------------------------------------

hdm <- create_hamming_distance_matrix(k_length, lambda = lambda, scale_fun = scaling_function)
pdv <- create_positional_distance_vector(window_size, k_length)

# Score multivalency ------------------------------------------------------

tic("Making the multivalency dataframe.")
all_kmer_multivalency <- mclapply(seq_along(sequences), function(i) {
  calculate_kmer_multivalencies_df(sequences[i],
                                   names(sequences)[i],
                                   k_length,
                                   window_size,
                                   smoothing_size,
                                   hdm,
                                   pdv)
}, mc.cores = 8) %>%
  bind_rows()
toc()


# Peak calling multivalency -----------------------------------------------

#Remove the 0'd out values padding all transcripts.
threshold <- quantile(all_kmer_multivalency$smoothed_kmer_multivalency[all_kmer_multivalency$smoothed_kmer_multivalency > 0], threshold_percentile)

call_multivalent_regions <- function(input_df, threshold){
  # input_df <- all_kmer_multivalency %>% filter(sequence_name == "ENST00000322723.8")
  # threshold <- threshold
  
  above_threshold <- input_df$smoothed_kmer_multivalency > threshold
  
  input_df$peak_start <- 0
  input_df$peak_end <- 0
  
  #If there are no multivalent positions, return an empty list.
  if(sum(above_threshold) == 0){
    return(input_df)
  }
  
  mv_rle <- rle(above_threshold)
  
  mv_end <- cumsum(mv_rle$lengths) #End positions for all runs.
  mv_start <- c(1, mv_end[-length(mv_end)] + 1)[mv_rle$values] #Start positions for runs above threshold.
  mv_end <- mv_end[mv_rle$values] #End positions for runs above threshold.
  
  # Adjust the positions based on the smoothing window.
  mv_start <- mv_start - ((smoothing_size - 1) / 2)
  mv_end <- mv_end + ((smoothing_size - 1) / 2)
  
  # Often, we end up finding multiple multivalent regions which, when considering the window size, are made up of
  # overlapping regions of the transcript. If these regions have an overlap of at least 2/3 the window size, they
  # are merged.
  
  mv_dist <- (mv_end[-length(mv_end)] + (window_size * (2/3))) - mv_start[-1] #Find the distance between clusters. Many
  
  olap <- mv_dist >= 0
  
  mv_start <- mv_start[!c(F, olap)] # Keep the first start of the overlaps
  mv_end <- mv_end[!c(olap, F)] # Keep the last end of the overlaps.
  
  #Annotated the peaks
  for(i in c(1:length(mv_start))) {
    input_df$peak_start[c(mv_start[i]:mv_end[i])] <- mv_start[i]
    input_df$peak_end[c(mv_start[i]:mv_end[i])] <- mv_end[i]
  }
  
  #Make a list of mv starts and ends.
  return(input_df)
  
}

tic("Calling multivalent peaks.")
all_kmer_peaks <- all_kmer_multivalency %>%
  group_split(sequence_name) %>%
  mclapply(., function(x) call_multivalent_regions(x, threshold), mc.cores = 8) %>%
  bind_rows()
toc()

# Annotating multivalency -------------------------------------------------

# I'm calling any peak that is at least 1/3 in a UTR a UTR peak.

tic("Annotating multivalent_peaks")
annotated_peaks <- all_kmer_peaks %>%
  dplyr::rename("transcript_id" = sequence_name) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  mutate(peak_width = peak_end - peak_start) %>%
  mutate(peak_location = case_when(peak_start == 0 ~ "None",
                                   peak_start + (peak_width / 3) < cds_start ~ "UTR5",
                                   peak_end - (peak_width / 3) > cds_end ~ "UTR3",
                                   T ~ "CDS")) %>%
  dplyr::select(transcript_id, kmer, kmer_multivalency, smoothed_kmer_multivalency, peak_start, peak_end, peak_location)
toc()

write_tsv(annotated_peaks, file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")
annotated_peaks <- read_tsv(file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")

# Preparing for projection and clustering --------------------------------------------------------------

tic("Creating matrix of kmer multivalency per peak")
peak_list_for_umap <- annotated_peaks %>%
  filter(peak_location != "None") %>%
  mutate(peak_identifier = paste0(transcript_id, ":", peak_start, "-", peak_end, "@", peak_location)) %>%
  split(f = as.factor(.$peak_location)) %>%
  lapply(., function(x) {
    x %>% 
      group_by(peak_identifier, kmer) %>%
      summarise(total_kmer_multivalency = sum(kmer_multivalency)) %>%
      group_by(peak_identifier) %>%
      mutate(prop_region_multivalency = total_kmer_multivalency / sum(total_kmer_multivalency)) %>%
      dplyr::select(peak_identifier, kmer, prop_region_multivalency) %>%
      pivot_wider(names_from = kmer, values_from = prop_region_multivalency) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("peak_identifier") %>%
      as.matrix() })
toc()

# We will encounter duplicated genes or repeated regions that have extremely high similarity to each other and screw up the clustering.
# To remove these, I calculate the similarity between all peaks using dist, then remove any peaks that show unusually high similarity.
# This is the slowest part of the entire analysis by far.

tic("Filter out repetitive regions")
filtered_peak_list_for_umap <- lapply(peak_list_for_umap, function(x) {
  
  peak_dist <- dist(x)
  
  # Add the smallest non-zero value to all distances in order to facilitate log transformation.
  smallest_dist <- unique(peak_dist) %>% sort() %>% .[2]
  logged_dist <- log2(peak_dist + smallest_dist)
  
  # 5 SD below mean are considered outliers here.
  ldist_vec <- as.vector(logged_dist)
  dist_cutoff <- mean(ldist_vec) - (5 * sd(ldist_vec))
  
  #Return the peak names
  blacklist <- which(as.matrix(logged_dist) < dist_cutoff, arr.ind = TRUE) %>% rownames() %>% unique()
  
  # x <- peak_list_for_umap$CDS

  return(x[!(rownames(x) %in% blacklist),])
  
})
toc()

# Projecting and clustering -----------------------------------------------

umap_cluster_results_list <- lapply(filtered_peak_list_for_umap, function(x) {
  # x <- filtered_peak_list_for_umap$UTR3
  
  # Experimenting with different UMAP parameters has very little effect on the final result.
  tic("Creating UMAP projections.")
  region_umap <- umap(x, n_components = 4, n_neighbors = 50, min_dist = 0.001)
  region_umap_low_dim <- umap(x)
  toc()
  
  region_umap_df <- region_umap$layout %>% 
    as.data.frame() %>%
    rownames_to_column("peak_identifier")
  
  region_umap_low_dim_df <- region_umap_low_dim$layout %>% 
    as.data.frame() %>%
    rownames_to_column("peak_identifier")
  
  tic("Clustering with OPTICS and DBSCAN")
  # Make initial optics.
  region_optics <- optics(region_umap_df %>%
                            dplyr::select(-peak_identifier), 
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
  
  # 0.4 seems to give pretty reasonable clusters. 
  # This is basically subjective, but we don't really care, as long as the clusters we make have predictive power.
  cluster_threshold <- KneeArrower::findCutoff(above_threshold, theshold_tests, method = "first", 0.4)$y
  
  # Define the clusters using the threshold.
  region_hdb <- extractDBSCAN(region_optics, cluster_threshold)
  
  region_cluster_identity <- as.data.frame(x[region_umap_df$peak_identifier,]) %>%
    rownames_to_column("peak_identifier") %>%
    mutate(cluster = region_hdb$cluster) %>%
    group_by(cluster) %>%
    mutate(group_size = dplyr::n()) %>%
    ungroup() %>%
    mutate(cluster = case_when((group_size < 20) ~ as.numeric(0),
                               (group_size >= 20) ~ as.numeric(cluster))) %>%
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
    left_join(region_umap_low_dim_df) %>%
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
  
  return(
    list("plot" = composite_plot,
         "data" = named_region_umap_df)
  )
  
})

# Write out plots and data ------------------------------------------------

ggsave(filename = "GASR/germs/plots/germs_CDS_umap_plot.pdf", 
       umap_cluster_results_list$CDS$plot,
       device = "pdf", units = "in",
       height = 6, width = 6)

ggsave(filename = "GASR/germs/plots/germs_UTR3_umap_plot.pdf", 
       umap_cluster_results_list$UTR3$plot,
       device = "pdf", units = "in",
       height = 6, width = 6)

ggsave(filename = "GASR/germs/plots/germs_UTR5_umap_plot.pdf", 
       umap_cluster_results_list$UTR5$plot,
       device = "pdf", units = "in",
       height = 6, width = 6)

write_tsv(umap_cluster_results_list$CDS$data, file = "GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")
write_tsv(umap_cluster_results_list$UTR3$data, file = "GASR/germs/data/germs_UTR3_umap_clusters.tsv.gz")
write_tsv(umap_cluster_results_list$UTR5$data, file = "GASR/germs/data/germs_UTR5_umap_clusters.tsv.gz")

# plotting CDS with new names ---------------------------------------------

just_cds <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")

just_cds %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df

just_cds %>%
  left_join(cluster_rename_df) %>%
  mutate(cluster_colour = factor(cluster_colour, levels = cluster_rename_df$cluster_colour)) %>%
  arrange(cluster_colour) ->
  just_cds_named


cds_umap_plot <- ggplot(just_cds_named %>%
                          filter(cluster > 0), 
                    aes(x = V1, y = V2, fill = cluster_colour)) +  
  geom_point(size = 3, colour = alpha("black", 1), pch = 21,) +
  scale_fill_identity(guide = "legend",
                      labels = just_cds_named$cluster_name_short,
                      breaks = just_cds_named$cluster_colour,
                      ) +
  guides(fill = guide_legend(ncol = 1, 
                             title = element_blank(),
                             override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
  )

ggsave(filename = "GASR/germs/plots/germs_CDS_umap_plot_ordered_named.pdf", cds_umap_plot,
       width = 6, height = 4.5)
