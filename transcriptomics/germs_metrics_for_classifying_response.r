rm(list = ls())

library(tidyverse)
library(patchwork)
library(broom)

# Load --------------------------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details_exons.txt", 
                               col_types = cols())

transcript_details_junctions <- read_rds("GASR/lists/longest_proteincoding_transcript_hs_junctionpositions.rds") %>%
  mutate(exon_density = number_all_exons / tx_length)

# Once this has been run once, we don't need to run it again.
# 
# germs <- read_tsv(file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")
# 
# # For normalising the multivalency
# minkm <- min(germs$kmer_multivalency)
# medkm <- median(germs$kmer_multivalency)
# 
# # Probably won't use smoothed values here, but good to keep it consistent. Ignore the 0s from the edges of transcripts.
# minskm <- min(germs$smoothed_kmer_multivalency[germs$smoothed_kmer_multivalency > 0])
# medskm <- median(germs$smoothed_kmer_multivalency[germs$smoothed_kmer_multivalency > 0])
# 
# # Normalise the multivalency
# germs <- germs %>%
#   mutate(med_norm =
#            (kmer_multivalency - minkm)/
#            (medkm - minkm),
#          med_norm_smooth =
#            (smoothed_kmer_multivalency - minskm)/
#            (medskm - minskm))
# 
# left_join(germs, transcript_details, by = "transcript_id") %>%
#   # dplyr::select(-c(peak_start, peak_end, peak_location)) %>%
#   group_by(transcript_id) %>%
#   mutate(start = row_number()) %>%
#   ungroup() %>%
#   mutate(end = start + 4,
#          region = case_when(start < cds_start ~ "UTR5", #Yes, I am aware that I am not taking k-length into account here.
#                             start > cds_end ~ "UTR3", #How would you prefer I classify boundary cases?
#                             T ~ "CDS")) %>%
#   group_by(region, kmer) %>%
#   mutate(kmer_percentile = rank(med_norm) / dplyr::n(),
#          scale_perc = scale(kmer_percentile),
#          scale_germ = scale(med_norm)) %>%
#   ungroup() ->
#   germs_regions
# 
# # Load proportion of multivalency per cluster -----------------------------
# 
# kmer_mv_prop <- read_tsv("GASR/germs/data/prop_mv_per_kmer_cds.tsv.gz") %>%
#   arrange(new_order)
# 
# kmer_mv_prop %>%
#   select(cluster_name_short, cluster_colour, new_order) %>%
#   distinct() %>%
#   mutate(cluster_name_short = cluster_name_short %>%
#            fct_relevel(kmer_mv_prop$cluster_name_short %>% unique())) ->
#   cluster_name_df
# 
# kmer_mv_prop %>%
#   arrange(new_order, desc(mean_prop_mv) ) %>%
#   mutate(cluster_name_short = cluster_name_short %>%
#            fct_relevel(kmer_mv_prop$cluster_name_short %>% unique())) %>%
#   group_by(cluster_name_short) %>%
#   mutate(rank = row_number()) %>%
#   mutate(cum_sum = cumsum(mean_prop_mv)) %>%
#   ungroup() ->
#   kmer_mv_prop_cumsum
# 
# ggplot(kmer_mv_prop_cumsum) +
#   aes(x = rank, y = cum_sum) +
#   geom_line() +
#   geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
#   facet_wrap(. ~ cluster_name_short, nrow = 2) +
#   theme_classic() +
#   theme(axis.text = element_text(color = 'black')) +
#   labs(x = "kmer rank", y = "cumulative sum of\nmultivalency proportion") ->
#   cumsum_plot
# 
# ggsave("GASR/germs/plots/cum_sum_of_mv_prop_in_cds_clusters.pdf", cumsum_plot,
#        width = 6, height = 3)
# 
# kmer_mv_prop_cumsum %>%
#   filter(cum_sum <= .5) ->
#   representative_kmers
# 
# 
# # Calculate multivalency sums for transcripts -----------------------------
# 
# # I change clusters to "cluster_x" based on their ordering (new order), as actual names will behave badly.
# representative_kmers %>%
#   select(kmer, new_order, mean_prop_mv) %>%
#   mutate(new_order = paste0("cluster_", new_order)) %>%
#   pivot_wider(names_from = "new_order", values_from = "mean_prop_mv") %>%
#   replace(is.na(.), 0) ->
#   representative_kmers_wide
# 
# germs_regions %>%
#   filter(region == "CDS") %>%
#   group_by(transcript_id, gene_name, kmer) %>%
#   summarise(sum_scale_germ = sum(scale_germ)) %>%
#   ungroup() %>%
#   inner_join(representative_kmers_wide) %>%
#   mutate_at(vars(contains("cluster")), ~ . * sum_scale_germ) %>%
#   group_by(transcript_id, gene_name) %>%
#   summarise_at(vars(contains("cluster")), sum) %>%
#   ungroup() %>%
#   left_join(transcript_details %>% select(transcript_id, gene_name, cds_length)) ->
#   scores_per_mv_cluster
# 
# scores_per_mv_cluster  %>%
#   pivot_longer(cols = contains("cluster"), names_to = "new_order", values_to = "total_mv_score") %>%
#   mutate(new_order = new_order %>% word(2, sep = "_") %>% as.numeric(),
#          # Normalising by cds length, but adjusting the scaling factor so that the numbers are not super small.
#          # This changes asbsolutely nothing, but it's easier to read the numbers.
#          length_norm_score = total_mv_score / (cds_length/median(transcript_details$cds_length))) %>%
#   left_join(cluster_name_df) ->
#   long_multivalency_scores_without_junctions
# 
# # Exon junction positions -------------------------------------------------
# 
# # This is kind of arbitrary size. 
# # Most correct way might be to have a decaying function with distance from EJC, but less intuitive to explain.
# junction_mask_size <- 100
# 
# transcript_details_junctions %>%
#   mutate(junction_mask = pmap(.l = list(junction_positions, tx_length),
#                               .f = function(junction_positions, tx_length) {
#                                 # junction_positions <- transcript_details_junctions$junction_positions[transcript_details_junctions$gene_name == "ZNRF4"]
#                                 # tx_length <- transcript_details_junctions$tx_length[transcript_details_junctions$gene_name == "ZNRF4"]
#                                 
#                                 junction_positions <- unlist(junction_positions) %>% unname()
#                                 
#                                 ones <- rep(1, tx_length - 4)
#                                 
#                                 if(length(junction_positions) == 0) {
#                                   return(ones)
#                                 }
#                                 
#                                 # EJC is ~20nt upstream of junction.
#                                 
#                                 junction_mask_starts <- junction_positions - ((junction_mask_size / 2) + 20)
#                                 junction_mask_ends <- junction_positions + ((junction_mask_size / 2) - 20)
#                                 
#                                 junction_mask_starts[junction_mask_starts < 1] <- 1
#                                 junction_mask_ends[junction_mask_ends > (tx_length - 4)] <- tx_length - 4
#                                 
#                                 lapply(c(1:length(junction_positions)), function(x) c(junction_mask_starts[x]:junction_mask_ends[x])) %>%
#                                   unlist() ->
#                                   junction_mask_indices
#                                 
#                                 ones[junction_mask_indices] = 0
#                                 
#                                 return(ones)
#                                 
#                               })) %>%
#   select(transcript_id, junction_mask) %>%
#   unnest(junction_mask = junction_mask) %>%
#   group_by(transcript_id) %>%
#   mutate(position = row_number()) %>%
#   ungroup() ->
#   junction_mask_df
# 
# # Approach 1: Masking junctions and scoring like before-------------------------------------------
# # 
# # germs_regions %>%
# #   group_by(transcript_id) %>%
# #   mutate(position = row_number()) %>%
# #   ungroup() %>%
# #   inner_join(junction_mask_df, by = c("transcript_id", "position")) %>%
# #   select(-position) %>%
# #   filter(region == "CDS") %>%
# #   group_by(transcript_id, gene_name, kmer) %>%
# #   summarise(sum_scale_germ = sum(scale_germ * junction_mask)) %>%
# #   ungroup() %>%
# #   inner_join(representative_kmers_wide) %>%
# #   mutate_at(vars(contains("cluster")), ~ . * sum_scale_germ) %>%
# #   group_by(transcript_id, gene_name) %>%
# #   summarise_at(vars(contains("cluster")), sum) %>%
# #   ungroup() %>%
# #   left_join(transcript_details %>% select(transcript_id, gene_name, cds_length)) ->
# #   scores_per_mv_cluster_junction_masked
# # 
# # scores_per_mv_cluster_junction_masked  %>%
# #   pivot_longer(cols = contains("cluster"), names_to = "new_order", values_to = "total_mv_score_junc") %>%
# #   mutate(new_order = new_order %>% word(2, sep = "_") %>% as.numeric(),
# #          # Normalising by cds length, but adjusting the scaling factor so that the numbers are not super small.
# #          # This changes asbsolutely nothing, but it's easier to read the numbers.
# #          length_norm_score_junc = total_mv_score_junc / (cds_length/median(transcript_details$cds_length))) %>%
# #   left_join(cluster_name_df) ->
# #   long_multivalency_scores_junc
# 
# # Exon junction positions part 2: percentage of multivalency -------------------------------------------------
# 
# germs_regions %>%
#   group_by(transcript_id) %>%
#   mutate(position = row_number()) %>%
#   ungroup() %>%
#   inner_join(junction_mask_df, by = c("transcript_id", "position")) %>%
#   select(-position) %>%
#   filter(region == "CDS") %>%
#   group_by(transcript_id, gene_name, kmer) %>%
#   summarise(sum_germ = sum(med_norm),
#             sum_germ_j = sum(med_norm * junction_mask)) %>%
#   ungroup() %>%
#   inner_join(representative_kmers_wide) %>%
#   mutate(across(contains("cluster"), 
#                 .fns = list(raw = ~ . * sum_germ), 
#                 .names = "{fn}:{col}"),
#          across(c(contains("cluster"), -contains("raw")), 
#                 .fns = list(junc = ~ . * sum_germ_j),
#                 .names = "{fn}:{col}")
#   ) %>%
#   select(-c(contains("cluster"), -contains("raw"), -contains("junc"),)) %>%
#   group_by(transcript_id, gene_name) %>%
#   summarise_at(vars(contains("cluster")), sum) %>%
#   ungroup() %>%
#   pivot_longer(cols = contains("cluster"), names_to = "cluster", values_to = "sum_mv") %>%
#   separate(cluster, into = c("score_type", "cluster"), sep = ":", convert = T) %>%
#   pivot_wider(names_from = score_type, values_from = sum_mv) %>%
#   mutate(percentage_multivalency = junc / raw) %>%
#   left_join(transcript_details %>% select(transcript_id, gene_name, cds_length)) %>%
#   mutate(new_order = cluster %>% word(2, sep = "_") %>% as.numeric()) %>%
#   left_join(cluster_name_df) ->
#   long_multivalency_scores_junc_percentage
# 
# 
# inner_join(long_multivalency_scores_without_junctions, long_multivalency_scores_junc_percentage) %>%
#   select(transcript_id, gene_name, cluster_name_short, cluster_colour, new_order, 
#          total_mv_score, length_norm_score, percentage_multivalency) ->
#   long_multivalency_scores_with_junction
# 
# write_tsv(long_multivalency_scores_with_junction, "GASR/germs/data/summed_multivalency_score_per_cluster_per_transcript.tsv.gz")

long_multivalency_scores <- read_tsv("GASR/germs/data/summed_multivalency_score_per_cluster_per_transcript.tsv.gz")


# Load DESeq tables and annotate with mv scores --------------------------------------------------------------------

er_interaction <- read_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_interaction_df.tsv")
er_cytoplasmic <- read_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_cytoplasmic_0vs12.tsv")
er_nuclear <- read_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_nuclear_0vs12.tsv")
er_timepoints <- read_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/neve_quantseq_each_timepoint.tsv")

er_interaction %>%
  filter(baseMean > quantile(baseMean, 0.25)) %>%
  left_join(long_multivalency_scores) %>%
  mutate(log2FoldChange = -log2FoldChange) %>% #We want positive scores to reflect being more nuclear after dox
  filter(new_order %in% c(1,2,3)) %>%
  # filter(new_order %in% c(6)) %>%
  group_by(transcript_id, gene_name, cds_length, log2FoldChange, padj) %>%
  summarise(max_score = max(total_mv_score),
    max_score_length = max(length_norm_score),
    ) %>%
  ungroup() %>%
  mutate(shifted_log_max_score = log2(max_score + abs(min(max_score)) + 0.1),
         scaled_log2_max_score = scale(shifted_log_max_score) %>% as.numeric(),
         shifted_log_max_score_length = log2(max_score_length + abs(min(max_score_length)) + 0.1),
         scaled_log2_max_score_length = scale(shifted_log_max_score_length) %>% as.numeric(),
         scaled_log2_cds_length = scale(log2(cds_length)) %>% as.numeric(),
        ) %>%
  dplyr::select(-c(shifted_log_max_score, shifted_log_max_score_length)) ->
  er_interaction_with_mv_scores

er_interaction %>%
  filter(baseMean > quantile(baseMean, 0.25)) %>%
  left_join(long_multivalency_scores) %>%
  mutate(log2FoldChange = -log2FoldChange) %>% #We want positive scores to reflect being more nuclear after dox
  filter(new_order %in% c(1,2,3)) %>%
  # filter(new_order %in% c(6)) %>%
  group_by(transcript_id, gene_name, cds_length, log2FoldChange, padj) %>%
  filter(total_mv_score == max(total_mv_score),) %>%
  ungroup() %>%
  select(transcript_id, gene_name, cds_length, log2FoldChange, padj, total_mv_score, percentage_multivalency) %>%
  mutate(shifted_log_max_score = log2(total_mv_score + abs(min(total_mv_score)) + 0.1),
         scaled_log2_max_score = scale(shifted_log_max_score) %>% as.numeric(),
         scaled_log2_cds_length = scale(log2(cds_length)) %>% as.numeric(),) %>%
  dplyr::select(-c(shifted_log_max_score)) ->
  er_interaction_with_mv_junc_perc

# Save summary table for publication
er_interaction_with_mv_junc_perc %>%
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%")),
         junction_class = case_when(percentage_multivalency > 0.5 ~ "distant from EJC",
                                    T ~ "close to EJC") %>%
           fct_relevel("close to EJC")) %>%
  select(transcript_id, gene_name, log2FoldChange, padj,
         purine_multivalency_score = total_mv_score,
         percentage_score_distant_from_junction = percentage_multivalency,
         purine_multivalency_class = score_class,
         purine_multivalency_location = junction_class) %>%
  write_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/publication_summary_data_interaction.tsv")

er_timepoints %>%
  group_by(transcript_id) %>%
  mutate(meanexpression = mean(baseMean)) %>%
  ungroup() %>%
  filter(meanexpression > quantile(meanexpression, 0.25)) %>%
  left_join(long_multivalency_scores) %>%
  # mutate(log2FoldChange = -log2FoldChange) %>%
  filter(new_order %in% c(1,2,3)) %>%
  group_by(timepoint, transcript_id, gene_name, cds_length, log2FoldChange, padj) %>%
  summarise(max_score = max(total_mv_score),
            max_score_length = max(length_norm_score)) %>%
  ungroup() %>%
  mutate(shifted_log_max_score = log2(max_score + abs(min(max_score)) + 0.1),
         scaled_log2_max_score = scale(shifted_log_max_score) %>% as.numeric(),
         shifted_log_max_score_length = log2(max_score_length + abs(min(max_score_length)) + 0.1),
         scaled_log2_max_score_length = scale(shifted_log_max_score_length) %>% as.numeric(),
         # POSITIVE = MORE NUCLEAR
         log2FoldChange = -log2FoldChange) ->
  er_timepoints_with_mv_scores

#Save publication stuff
er_timepoints_with_mv_scores %>%
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%")),
         timepoint = paste0(timepoint, "h dox") %>%
           str_sub(2)) %>%
  select(timepoint, transcript_id, gene_name, log2FoldChange, padj,
         purine_multivalency_score = max_score,
         purine_multivalency_class = score_class) %>%
  write_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/publication_summary_data_timepoints.tsv")

# Interaction of length and multivalency ----------------------------------

# POSITIVE L2FC IS NOW MORE NUCLEAR AFTER DOX

er_interaction_with_mv_scores %>% 
  mutate(capped_lfc = case_when(log2FoldChange > 1.5 ~ 1.5,
                                log2FoldChange < -1.5 ~ -1.5,
                                T ~ log2FoldChange)) %>%
  ggplot() + 
  aes(x = scaled_log2_max_score_length, y =  log2(cds_length), fill = capped_lfc) + 
  geom_point(shape = 21) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = log2(median(er_interaction_with_mv_scores$cds_length)), linetype = "dashed") +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "scaled log2 length-normalised\npurine multivalency",
       y = "log2 cds length",
       fill = "capped log2 fold change\n(dox * localisation)") ->
  interaction_scatter_plot_length_vs_lengthscaledmultivalency

er_interaction_with_mv_scores %>% 
  mutate(capped_lfc = case_when(log2FoldChange > 1.5 ~ 1.5,
                                log2FoldChange < -1.5 ~ -1.5,
                                T ~ log2FoldChange)) %>%
  ggplot() + 
  aes(x = scaled_log2_max_score, y = log2(cds_length), fill = capped_lfc) + 
  geom_point(shape = 21) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = log2(median(er_interaction_with_mv_scores$cds_length)), linetype = "dashed") +
  theme(axis.text = element_text(color = "black")) +
  coord_cartesian(c(-5,7)) +
  labs(x = "scaled log2\npurine multivalency",
       y = "log2 cds length",
       fill = "capped log2 fold change\n(dox * localisation)") ->
  interaction_scatter_plot_length_vs_unscaledmultivalency

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_scatter_plot_length_vs_lengthscaledmultivalency.pdf",
       interaction_scatter_plot_length_vs_lengthscaledmultivalency,
       width = 6, height = 4)

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_scatter_plot_length_vs_multivalency.pdf",
       interaction_scatter_plot_length_vs_unscaledmultivalency,
       width = 6, height = 4)

lm(log2FoldChange ~ scaled_log2_cds_length * scaled_log2_max_score_length, data = er_interaction_with_mv_scores) %>%
  tidy() -> er_int_length_mv_model

write_tsv(er_int_length_mv_model, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_int_with_length_mv_linearmodel.tsv")


# Interaction of multivalency and intron content --------------------------

er_interaction_with_mv_scores %>% 
  left_join(transcript_details_junctions) %>%
  mutate(capped_lfc = case_when(log2FoldChange > 1.5 ~ 1.5,
                                log2FoldChange < -1.5 ~ -1.5,
                                T ~ log2FoldChange)) %>%
  sample_frac(1) %>%
  ggplot() + 
  aes(x = scaled_log2_max_score_length, y =  log2(exon_density), fill = capped_lfc) + 
  # facet_wrap(. ~ (total_mv_score > quantile(er_interaction_with_mv_junc_perc$total_mv_score, .9))) +
  geom_point(shape = 21) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = log2(median(transcript_details_junctions$exon_density)), linetype = "dashed") +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "scaled log2 length-normalised\npurine multivalency",
       y = "log2 exon density",
       fill = "capped log2 fold change\n(dox * localisation)") ->
  exon_density_plot

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_scatter_plot_exondensity_vs_lengthscaledmultivalency.pdf",
       exon_density_plot,
       width = 6, height = 4)

lm(log2FoldChange ~ scale(log2(exon_density)) * scaled_log2_max_score_length, data = er_interaction_with_mv_scores %>% 
     left_join(transcript_details_junctions)) %>%
  tidy() ->
  exon_density_model

write_tsv(exon_density_model, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/exon_density_linearmodel.tsv")


# Plotting interaction density effect --------------------------------------------------------------------

er_interaction_with_mv_scores %>% 
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%"))) ->
  er_int_with_classes 

er_int_with_classes %>%
  ggplot() +
  aes(x = log2FoldChange, fill = score_class) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "log2 fold change\n(fraction * dox)",
       fill = "purine\nmultivalency") +
  scale_fill_manual(values = c("gray30", "orange2", "red3")) +
  coord_cartesian(xlim = c(-1.5, 3)) ->
  er_interaction_plot

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_density_plot_0vs12.pdf",
       er_interaction_plot,
       width = 4.5, height = 3)

list("bottom_vs_top" = er_int_with_classes %>%
  filter(score_class != "top 10%") %>%
  t.test(log2FoldChange ~ score_class, data = .) %>%
  tidy(),
  "bottom_vs_mid" = er_int_with_classes %>%
    filter(score_class != "top 2%") %>%
    t.test(log2FoldChange ~ score_class, data = .) %>%
    tidy(),
  "mid_vs_top" = er_int_with_classes %>%
    filter(score_class != "bottom 90%") %>%
    t.test(log2FoldChange ~ score_class, data = .) %>%
    tidy()) %>%
  bind_rows(.id = "comparison") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  er_int_ttests
  
write_tsv(er_int_with_classes, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_int_with_classes.tsv.gz")
write_tsv(er_int_ttests, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_int_with_classes_ttests.tsv")


# Junction multivalency analysis ------------------------------------------

er_interaction_with_mv_junc_perc %>% 
  drop_na() %>%
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%")),
         junction_class = case_when(percentage_multivalency > 0.5 ~ "distant from EJC",
                                    T ~ "close to EJC") %>%
           fct_relevel("close to EJC")) ->
  er_int_junc_with_classes 

er_int_junc_with_classes %>%
  ggplot() +
  aes(x = percentage_multivalency) +
  geom_density(alpha = 0.5)

er_int_junc_with_classes %>%
  count(score_class, junction_class)

er_int_junc_with_classes %>%
  ggplot() +
  aes(x = log2FoldChange, color = score_class, linetype = junction_class) +
  facet_wrap(. ~ score_class) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "log2 fold change\n(fraction * dox)",
       linetype = "GeRM location",
       color = "purine\nmultivalency"
       ) +
  scale_color_manual(values = c("gray30", "orange2", "red3")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  coord_cartesian(xlim = c(-1.5, 3)) ->
  er_int_with_junc_perc_plot_all_classes

er_int_junc_with_classes %>%
  filter(score_class == "top 2%") %>%
  ggplot() +
  aes(x = log2FoldChange, linetype = junction_class) +
  geom_density(color = "red3") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "log2 fold change\n(fraction * dox)",
       linetype = "GeRM location",
       color = "purine\nmultivalency"
  ) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  coord_cartesian(xlim = c(-1.5, 3)) ->
  er_int_with_junc_perc_plot_top2

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_density_plot_top2_ejc_distance.pdf",
       er_int_with_junc_perc_plot_top2,
       width = 4.5, height = 3)

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/interaction_density_plot_allclasses_ejc_distance.pdf",
       er_int_with_junc_perc_plot_all_classes,
       width = 8.5, height = 3)


er_int_junc_with_classes %>%
  select(score_class, junction_class, log2FoldChange) %>%
  nest(data = -score_class) %>%
  mutate(ttest = map(data, ~{t.test(log2FoldChange ~ junction_class, data = .x) %>% tidy() })) %>%
  unnest(ttest) %>%
  select(-data) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  er_int_with_junc_perc_ttests

er_int_junc_with_classes %>%
  filter(score_class == "top 2%") %>%
  ggplot() +
  aes(y = log2FoldChange, x = percentage_multivalency) +
  geom_point() +
  geom_smooth(method = "glm", se = F)

lm(log2FoldChange ~ percentage_multivalency, data = er_int_junc_with_classes %>%
  filter(score_class == "top 2%")) %>% summary()

write_tsv(er_int_with_junc_perc_ttests, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_int_with_junc_perc_ttests.tsv")

# Time course -------------------------------------------------------------

er_timepoints_with_mv_scores %>%
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%")),
         timepoint = timepoint %>%
           str_remove("h") %>%
           as.numeric()) ->
  er_timepoints_with_classes

er_timepoints_with_classes %>%
  ggplot() +
  aes(x = factor(timepoint), y = log2FoldChange, fill = score_class) +
  facet_wrap(. ~ score_class) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.5, width = 0.33, outlier.size = .5, show.legend = F) +
  coord_cartesian(ylim = c(-3.5, 4)) +
  theme_classic() +
  labs(fill = "purine\nmultivalency",
       x = "dox induction",
       y = "log2 nucleus / cytoplasm") +
  theme(axis.text = element_text(color = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("gray30", "orange2", "red3")) ->
  timecourse_boxplots

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/er_timecourse_boxplots.pdf",
       timecourse_boxplots, width = 6, height = 3)

  # My dplyr skills have failed me - resorting to foul base R (boomer shit). 
list("bottom_90" = pairwise.t.test(x = er_timepoints_with_classes$log2FoldChange[er_timepoints_with_classes$score_class == "bottom 90%"],
                g = er_timepoints_with_classes$timepoint[er_timepoints_with_classes$score_class == "bottom 90%"]) %>%
  tidy(),
  "top_10" = pairwise.t.test(x = er_timepoints_with_classes$log2FoldChange[er_timepoints_with_classes$score_class == "top 10%"],
                                g = er_timepoints_with_classes$timepoint[er_timepoints_with_classes$score_class == "top 10%"]) %>%
    tidy(),
  "top_2" = pairwise.t.test(x = er_timepoints_with_classes$log2FoldChange[er_timepoints_with_classes$score_class == "top 2%"],
                             g = er_timepoints_with_classes$timepoint[er_timepoints_with_classes$score_class == "top 2%"]) %>%
    tidy()) %>%
  bind_rows(.id = "score_class") %>%
  filter(group2 == "0") %>%
  mutate(padj = p.adjust(p.value)) ->
  comparisons_to_zero

write_tsv(er_timepoints_with_classes, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_timecourse_with_classes.tsv.gz")
write_tsv(comparisons_to_zero, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_timecourse_ttests_to_zero.tsv")

# Gene ontology -----------------------------------------------------------

library(topGO)
library(biomaRt)

all_genes <- transcript_details$gene_name

#Get GO terms from Ensembl.
db <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "https://www.ensembl.org")
go_ids <- getBM(attributes = c('go_id', 'external_gene_name', 'namespace_1003'), 
                filters = 'external_gene_name', values = all_genes, mart = db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO <- unstack(go_ids[,c(1,2)])

# Make the target and background lists
er_interaction_with_mv_scores %>%
  filter(padj < 0.05, 
         log2FoldChange > 0) %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  more_retained_in_nucleus

er_interaction_with_mv_scores %>%
  filter(padj < 0.05, 
         log2FoldChange < 0) %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  more_exported_from_nucleus

er_interaction_with_mv_scores %>%
  dplyr::select(gene_name) %>%
  unlist(use.names = F) ->
  expressed_in_interaction_set

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
           weightFisher = weightFisher %>% str_replace("< 1e-30", "1e-30") %>% as.numeric()) %>%
    replace_na(list(weightFisher = 0))
}

retained_in_nucleus <- go_term_table(more_retained_in_nucleus, all_genes)
exported_from_nucleus <- go_term_table(more_exported_from_nucleus, all_genes)

retained_in_nucleus$padj <- p.adjust(retained_in_nucleus$weightFisher, method = "BH") 

my_fav_terms <- c(
  "nucleolus",
  "RNA binding",
  "regulation of mRNA processing",
  "nuclear speck",
  "chromatin binding",
  "chromatin remodeling"
  )

term_levels <- c(
  "nucleolus",
  "RNA binding",
  "reg. of mRNA proc.",
  "nuclear speck",
  "chromatin binding",
  "chromatin remodeling"
)

retained_in_nucleus %>%
  dplyr::select(Term, padj, fold_enrichment) %>%
  filter(Term %in% my_fav_terms) %>%
  arrange(Term) %>%
  mutate(Term = Term %>% 
           str_replace("regulation of mRNA processing", "reg. of mRNA proc.") %>%
           fct_relevel(term_levels),
         neg_log_padj = -log2(padj)) %>% 
  ggplot() +
  aes(y = Term, x = fold_enrichment, fill = neg_log_padj, size = neg_log_padj) +
  geom_segment(aes(x = 0, xend = fold_enrichment, y = Term, yend = Term), linewidth = 1, show.legend = F) +
  geom_point(shape = 21) +
  scale_size_continuous(limits=c(0, 64), breaks = seq(0, 64, by = 16)) +
  theme_classic() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 64), breaks = seq(0, 64, by = 16)) +
  guides(fill = guide_legend(), size = guide_legend()) +
  labs(x = "fold enrichment", y = "", size = "-log2\np. adj", fill = "-log2\np. adj") +
  theme(axis.text = element_text(color = "black")) ->
  go_terms_plot

write_tsv(exported_from_nucleus, "GASR/export_reporter/quantseq/neve_timecourse_1/tables/gene_ontology_more_exported_from_nuc.tsv.gz")

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/goterms_more_exported.pdf",
       go_terms_plot, width = 4, height = 2)

# nuclear_speck_genes <- go_ids %>%
#   filter(go_id == "GO:0016607") %>%
#   dplyr::select(external_gene_name) %>%
#   unlist(use.names = F)
# 
# upr_genes <- go_ids %>%
#   filter(go_id == "GO:0051082") %>%
#   dplyr::select(external_gene_name) %>%
#   unlist(use.names = F)
# 
# cytotrans_genes <- go_ids %>%
#   filter(go_id == "GO:0002181") %>%
#   dplyr::select(external_gene_name) %>%
#   unlist(use.names = F)


# AA properties of significant genes --------------------------------------

library(Biostrings)
library(roll)
prots <- readAAStringSet("GASR/lists/longest_gencode29_prots.fa")

tibble(transcript_id = names(prots),
       protein = as.character(prots)) %>%
  inner_join(transcript_details) %>%
  mutate(aa_props = map(protein, 
                        ~{.x %>%
                            str_split("") %>% 
                            table() %>% 
                            as_tibble() %>% 
                            set_names("aa", "count")}),
         more_retained = case_when(gene_name %in% more_retained_in_nucleus ~ "retained",
                                   gene_name %in% more_exported_from_nucleus ~ "exported",
                                   T ~ "unchanged") %>%
           fct_relevel("unchanged")) %>%
  dplyr::select(-c(cds_start, cds_length, tx_length, cds_end, protein)) %>%
  unnest(aa_props) %>%
  filter(aa != "*") ->
  aa_prop_per_protein


aa_prop_per_protein %>%
  # Don't care about 'exported' genes (and they resemble unaffected ones here anyway).
  filter(more_retained != "exported") %>%
  group_by(gene_id, more_retained, aa) %>%
  summarise(count = sum(count)) %>%
  group_by(gene_id, more_retained) %>%
  mutate(prop = count/sum(count)) %>%
  ungroup() %>% 
  mutate(aa = aa %>%
           fct_relevel("K", "R", "E", "D", "S", "Q", "P", "A", "G", "L", "T", "N", "V", "W", "M", "Y")) ->
  aa_prop_comparison
  
aa_prop_comparison %>% filter(aa %in% c("K", "R", "E", "D", "S")) -> just_charged_aa

just_charged_aa %>%
  group_by(aa) %>%
  summarise(ttest = t.test(prop ~ more_retained) %>% tidy()) %>%
  unnest(ttest) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  charged_aa_prop_ttests

ggplot(just_charged_aa, aes(x = aa, y = prop, fill = more_retained)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), width = 0.2, outlier.shape = NA, show.legend = F) +
  scale_fill_manual(values = c("gray70", "red3")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "", y = "proportion of protein", fill = "") +
  coord_cartesian(ylim = c(0, .15)) ->
  charged_aa_prop_plot


ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/charged_aa_prop_plot.pdf", charged_aa_prop_plot,
       width = 5, height = 3)

write_tsv(charged_aa_prop_ttests, "GASR/export_reporter/quantseq/neve_timecourse_1/plots/charged_aa_prop_ttests.tsv")

tibble(transcript_id = names(prots),
       protein = as.character(prots)) %>%
  inner_join(transcript_details) %>%
  # sample_frac(0.02) %>%
  mutate(prop_charged = map(protein, 
                            ~{
                              aa_window_size = 40
                              charge_threshold = 0.4
                              
                              zeros <- rep(0, nchar(.x))
                              
                              one_locs <- .x %>% str_locate_all("K|E|D|R") %>% as.data.frame() %>% select(start) %>% unlist(use.names = F)
                              
                              zeros[one_locs] <- 1
                              
                              prop_charged <- roll_mean(zeros, aa_window_size)[-c(1:(aa_window_size-1))]
                              
                              return(sum(prop_charged > charge_threshold) / length(prop_charged))
                            }) %>% unlist(),
         
         more_retained = case_when(gene_name %in% more_retained_in_nucleus ~ "retained",
                                   gene_name %in% more_exported_from_nucleus ~ "exported",
                                   gene_name %in% er_interaction_with_mv_scores$gene_name ~ "unchanged",
                                   T ~ "no_data") %>%
           fct_relevel("unchanged", "retained")) %>%
  filter(more_retained != "no_data") %>%
  dplyr::select(-c(cds_start, cds_length, tx_length, cds_end, protein)) ->
  prop_charged

ggplot(prop_charged %>% filter(more_retained != "exported"), aes(x = more_retained, y = prop_charged, fill = more_retained)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.25, show.legend = F) +  
  scale_fill_manual(values = c("gray70", "red3")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "",
       y = "proportion of protein\n>40% charged",
       fill = "") ->
  proportion_highly_charged_boxplot

ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/proportion_highly_charged_boxplot.pdf", proportion_highly_charged_boxplot,
          width = 2.25, height = 3)

ggplot(prop_charged %>% filter(more_retained != "exported"), aes(x = prop_charged, fill = more_retained)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("gray70", "red3")) +
  theme_classic() +
  coord_cartesian(xlim = c(0, 0.6)) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "proportion of protein\nhighly charged",
       fill = "") ->
  proportion_highly_charged_densityplot

wilcox.test(prop_charged ~ more_retained, data = prop_charged %>% filter(more_retained != "exported")) %>% tidy() ->
  proportion_highly_charged_wilcox

# ggsave("GASR/export_reporter/quantseq/neve_timecourse_1/plots/proportion_highly_charged_plot.pdf", proportion_highly_charged_densityplot,
#           width = 4, height = 3)

write_tsv(proportion_highly_charged_wilcox, "GASR/export_reporter/quantseq/neve_timecourse_1/plots/proportion_highly_charged_wilcoxon_test.tsv")


# CLKi --------------------------------------------------------------------


clk_interaction <- read_tsv("GASR/clk/tables/neve_quantseq_interaction_df.newmapping.tsv")

clk_interaction %>%
  filter(baseMean > quantile(baseMean, 0.25)) %>%
  left_join(long_multivalency_scores) %>%
  mutate(log2FoldChange = -log2FoldChange) %>% #We want positive scores to reflect being more nuclear after dox
  filter(new_order %in% c(1,2,3)) %>%
  # filter(new_order %in% c(6)) %>%
  drop_na() %>%
  group_by(transcript_id, gene_name, cds_length, log2FoldChange, padj) %>%
  summarise(max_score = max(total_mv_score),
            max_score_length = max(length_norm_score)) %>%
  ungroup() %>%
  mutate(shifted_log_max_score = log2(max_score + abs(min(max_score)) + 0.1),
         scaled_log2_max_score = scale(shifted_log_max_score) %>% as.numeric(),
         shifted_log_max_score_length = log2(max_score_length + abs(min(max_score_length)) + 0.1),
         scaled_log2_max_score_length = scale(shifted_log_max_score_length) %>% as.numeric(),
         scaled_log2_cds_length = scale(log2(cds_length)) %>% as.numeric()) %>%
  dplyr::select(-c(shifted_log_max_score, shifted_log_max_score_length)) ->
  clk_interaction_with_mv_scores

#Save for publication

clk_interaction_with_mv_scores %>%
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%"))) %>%
  select(transcript_id, gene_name, log2FoldChange, padj,
         purine_multivalency_score = max_score,
         purine_multivalency_class = score_class) %>%
  write_tsv("GASR/clk/tables/clk_publication_summary_data_interaction.newmapping.tsv")


clk_interaction_with_mv_scores %>% 
  mutate(capped_lfc = case_when(log2FoldChange > 1.5 ~ 1.5,
                                log2FoldChange < -1.5 ~ -1.5,
                                T ~ log2FoldChange)) %>%
  ggplot() + 
  aes(x = scaled_log2_max_score_length, y = log2(cds_length), fill = capped_lfc) + 
  geom_point(shape = 21) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = log2(median(clk_interaction_with_mv_scores$cds_length)), linetype = "dashed") +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "scaled log2 length-normalised\npurine multivalency",
       y = "log2 cds length",
       fill = "capped log2 fold change\n(dox * localisation)") ->
  clk_interaction_scatter_plot_length_vs_lengthscaledmultivalency

clk_interaction_with_mv_scores %>% 
  mutate(capped_lfc = case_when(log2FoldChange > 1.5 ~ 1.5,
                                log2FoldChange < -1.5 ~ -1.5,
                                T ~ log2FoldChange)) %>%
  ggplot() + 
  aes(x = scaled_log2_max_score, y = log2(cds_length), fill = capped_lfc) + 
  geom_point(shape = 21) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic() +
  # coord_cartesian(xlim = c(-5,7)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = log2(median(clk_interaction_with_mv_scores$cds_length)), linetype = "dashed") +
  theme(axis.text = element_text(color = "black")) +
  coord_cartesian(c(-5,7)) +
  labs(x = "scaled log2\npurine multivalency",
       y = "scaled log2 cds length",
       fill = "capped log2 fold change\n(dox * localisation)") ->
  clk_interaction_scatter_plot_length_vs_unscaledmultivalency

lm(log2FoldChange ~ scaled_log2_cds_length * scaled_log2_max_score_length, data = clk_interaction_with_mv_scores) %>%
  tidy() -> clk_int_length_mv_model

write_tsv(clk_int_length_mv_model, "GASR/clk/plots/clk_int_with_length_mv_linearmodel.newmapping.tsv")

ggsave("GASR/clk/plots/interaction_scatter_plot_length_vs_lengthscaledmultivalency.newmapping.pdf",
       clk_interaction_scatter_plot_length_vs_lengthscaledmultivalency,
       width = 6, height = 4)

ggsave("GASR/clk/plots/interaction_scatter_plot_length_vs_multivalency.newmapping.pdf",
       clk_interaction_scatter_plot_length_vs_unscaledmultivalency,
       width = 6, height = 4)

# CLK - Plotting interaction density effect --------------------------------------------------------------------

clk_interaction_with_mv_scores %>% 
  mutate(score_class = cut(scaled_log2_max_score, 
                           quantile(scaled_log2_max_score, c(0, .9, .98, 1)), 
                           include.lowest = T,
                           labels = c("bottom 90%", "top 10%", "top 2%"))) ->
  clk_int_with_classes 

clk_int_with_classes %>%
  ggplot() +
  aes(x = log2FoldChange, fill = score_class) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "log2 fold change\n(fraction * drug)",
       fill = "purine\nmultivalency") +
  scale_fill_manual(values = c("#B4B4B4", "#AA94B1", "#B9529F")) +
  coord_cartesian(xlim = c(-1.5, 2.5)) ->
  clk_interaction_plot

ggsave("GASR/clk/plots/clk_interaction_plot.newmapping.pdf",
       clk_interaction_plot,
       width = 4.5, height = 3)

list("bottom_vs_top" = clk_int_with_classes %>%
       filter(score_class != "top 10%") %>%
       t.test(log2FoldChange ~ score_class, data = .) %>%
       tidy(),
     "bottom_vs_mid" = clk_int_with_classes %>%
       filter(score_class != "top 2%") %>%
       t.test(log2FoldChange ~ score_class, data = .) %>%
       tidy(),
     "mid_vs_top" = clk_int_with_classes %>%
       filter(score_class != "bottom 90%") %>%
       t.test(log2FoldChange ~ score_class, data = .) %>%
       tidy()) %>%
  bind_rows(.id = "comparison") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  clk_int_ttests

write_tsv(clk_int_ttests, "GASR/clk/plots/clk_int_ttests.newmapping.tsv")



# Correlations between CLK and PPIG ---------------------------------------

inner_join(clk_interaction_with_mv_scores %>%
  select(transcript_id, gene_name, clki_l2fc = log2FoldChange),
  er_interaction_with_mv_scores %>%
    select(transcript_id, gene_name, ppig_l2fc = log2FoldChange, scaled_log2_max_score)) ->
  merged_clk_er

lm(ppig_l2fc ~ clki_l2fc, data = merged_clk_er) %>% summary()
cor(merged_clk_er$clki_l2fc, merged_clk_er$ppig_l2fc)  

ggplot(merged_clk_er) +
  aes(x = ppig_l2fc, y = clki_l2fc, fill = scaled_log2_max_score) +
  geom_point(shape = 21, size = 2) +
  # scale_fill_gradient2(mid = "white", low = "#B4B4B4", high = "#B9529F", midpoint = 0,
  #                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_fill_gradientn(colours = c("#B4B4B4", "white", "#B9529F"), values = c(0, 0.66 , 1),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  # colorspace::scale_fill_continuous_diverging(palette = "Purple-Green", p1 = 0.75, p2 = 0.75) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-2.5, 3.5)) +
  annotate(geom = "text", x = -1.75, y = 2, label = "R = 0.44") +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "log2 fold change\n(fraction * PPIG expression)",
       y = "log2 fold change\n(fraction * CLK-IN-T3)",
       fill = "scaled\npurine\nmultivalency") ->
  correlation_between_datasets_plot

ggsave("GASR/clk/plots/clk_vs_ppig_correlation.newmapping.pdf",
       correlation_between_datasets_plot,
       width = 4.5, height = 3)
