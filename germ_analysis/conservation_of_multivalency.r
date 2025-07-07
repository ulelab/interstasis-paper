rm(list = ls())

library(tidyverse)
library(broom)
# library(ggbeeswarm)

window_size <- 123
k_length <- 5

# Load da shit ------------------------------------------------------------

# For normalisation purposes
germs <- read_tsv(file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")

# For normalising the multivalency
minkm <- min(germs$kmer_multivalency)
medkm <- median(germs$kmer_multivalency)

just_smoothed_mv <- germs %>%
  group_by(transcript_id) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  select(transcript_id, position, smoothed_kmer_multivalency) %>%
  filter(smoothed_kmer_multivalency > 0) %>%
  mutate(normalised_skm = (smoothed_kmer_multivalency - minkm)/
           (medkm - minkm))

# just_smoothed_mv <- just_smoothed_mv %>% dplyr::rename("position" = "transcript_position")

rm(germs)

mutated_multivalency <- read_tsv("GASR/germs/data/mutated_multivalency_correct.tsv.gz")

mutated_multivalency %>%
  # There are a few NA values around because they are too close to the start/end of the transcript - not sure how these got through.
  drop_na() %>%
  mutate(native = (native - minkm)/
           (medkm - minkm),
         mutated = (mutated - minkm)/
           (medkm - minkm),
         diff = native - mutated,
         ratio = native / mutated) ->
  mutated_multivalency

bullshit <- mutated_multivalency %>% group_by(transcript_id) %>% summarise(bullshit = "*" %in% aa) %>% filter(bullshit)

transcript_conservation <- read_tsv("genomes/conservation/conservation_in_transcriptome_coordinates.tsv.gz")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_length - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

mutated_multivalency_conservation <- mutated_multivalency %>%
  dplyr::rename("position" = "transcript_position") %>%
  mutate(codon2_position = case_when(position_in_codon == 1 ~ position + 1,
                                     position_in_codon == 3 ~ position - 1)) %>%
  left_join(transcript_conservation, by = c("transcript_id", "position")) %>%
  left_join(transcript_conservation %>% dplyr::rename(codon2_position = position, 
                                                      codon2_phylop_17_primates = phylop_17_primates,
                                                      codon2_phylop_470_mammals = phylop_470_mammals,
                                                      codon2_phylop_100_vertebrates = phylop_100_vertebrates), 
            by = c("transcript_id", "codon2_position")) %>%
  inner_join(transcript_details) %>%
  # Add information about local multivalency
  inner_join(just_smoothed_mv, by = c("transcript_id", "position")) %>%
  filter(!(transcript_id %in% bullshit$transcript_id)) %>%
  mutate(codon2_diff_17p = phylop_17_primates - codon2_phylop_17_primates,
         codon2_diff_100v = phylop_100_vertebrates - codon2_phylop_100_vertebrates,
         codon2_diff_470m = phylop_470_mammals - codon2_phylop_470_mammals,) %>%
  group_by(codon, position_in_codon) %>%
  # Scale by subtraction because phylop is log scaled.
  mutate(scaled_17p = phylop_17_primates - mean(phylop_17_primates),
         scaled_100v = phylop_100_vertebrates - mean(phylop_100_vertebrates),
         scaled_470m = phylop_470_mammals - mean(phylop_470_mammals),
         codon2_scaled_17p = codon2_diff_17p - mean(codon2_diff_17p),
         codon2_scaled_100v = codon2_diff_100v - mean(codon2_diff_100v),
         codon2_scaled_470m = codon2_diff_470m - mean(codon2_diff_470m),
         ) %>%
  ungroup() 

# Codons like to promote multivalency -------------------------------------

mutated_multivalency_conservation %>%
  group_by(aa, codon, position_in_codon) %>%
  summarise(mean_ratio = mean(ratio)) ->
  mean_ratios_per_codon

mutated_multivalency_conservation %>%
  count(aa, codon) %>%
  group_by(aa) %>%
  mutate(usage = n/sum(n)) %>%
  mutate(scaled_usage = scale(usage) %>% as.numeric(),
         otherscale_usage = usage/mean(usage)) ->
  usage_df

mean_ratios_per_codon %>%
  left_join(usage_df) %>% 
  ggplot() +
  aes(x = "", y = mean_ratio, fill = log2(otherscale_usage)) +
  geom_beeswarm(cex = 5, shape = 21, size = 3) +
  scale_fill_gradient2(midpoint = 0) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black")) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  labs(y = "ratio\nnative to mutated GeRM",
       x = "",
       fill = "log2\nrelative\ncodon\nusage")

mean_ratios_per_codon %>%
  left_join(usage_df) %>%
  ungroup() %>%
  reframe(model = lm(log2(mean_ratio) ~ log2(otherscale_usage)) %>% tidy()) %>%
  unnest(model) ->
  codon_usage_ratio_model

mean_ratios_per_codon %>%
  left_join(usage_df) %>% 
  ggplot() +
  aes(x = log2(otherscale_usage), y = log2(mean_ratio)) +
  geom_abline(slope = codon_usage_ratio_model$estimate[codon_usage_ratio_model$term == "log2(otherscale_usage)"],
              intercept = codon_usage_ratio_model$estimate[codon_usage_ratio_model$term == "(Intercept)"],
              linetype = "dashed", linewidth = 1) +
  geom_point(shape = 21, size = 3) +
  # scale_fill_gradient2(midpoint = 0) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = codon_usage_ratio_model$estimate[codon_usage_ratio_model$term == "(Intercept)"],
             linetype = "dashed", alpha = 0.7,
             color = "orangered3") +
  labs(y = "log2 ratio\nnative to mutated GeRM",
       x = "log2 relative codon usage") ->
  codon_usage_ratio_plot

write_tsv(codon_usage_ratio_model, "GASR/germs/plots/conservation_mutation/codon_usage_ratio_model.tsv")
ggsave("GASR/germs/plots/conservation_mutation/codon_usage_ratio_plot.pdf", codon_usage_ratio_plot,
       width = 3, height = 3)

library(ggrepel)

mutated_multivalency_conservation %>%
  group_by(aa, codon, position_in_codon) %>%
  summarise(mean_vert = mean(phylop_100_vertebrates),
            mean_mam = mean(phylop_470_mammals),
            mean_c2_vert = mean(codon2_phylop_100_vertebrates),
            mean_c2_mam = mean(codon2_phylop_470_mammals),
            mean_c2diff_vert = mean(codon2_diff_100v),
            mean_c2diff_mam = mean(codon2_diff_470m)) %>%
  ungroup() ->
  mean_conservations_per_codon_pos

ggplot(mean_conservations_per_codon_pos, aes(x = mean_c2_mam, y = mean_mam, fill = mean_c2diff_mam)) +
  geom_point(data = mean_conservations_per_codon_pos %>% filter(aa == "R"), size = 3, shape = 25) +
  geom_point(data = mean_conservations_per_codon_pos %>% filter(aa != "R"), size = 3, shape = 21) +
  scale_fill_gradient2(midpoint = median(mean_conservations_per_codon_pos$mean_c2diff_mam)) +
  geom_label_repel(mean_conservations_per_codon_pos %>% filter(aa == "R"),
                   mapping = aes(label = paste0(codon, ": pos. ", position_in_codon), x = mean_c2_mam, y = mean_mam), 
                   size = 3,
                   box.padding = 1, alpha = 0.5, max.overlaps = Inf) +
  theme_classic() + theme(axis.text = element_text(color = "black")) +
  coord_cartesian(xlim = c(3, 7), ylim = c(-5, 4)) +
  labs(x = "mean position 2 PhyloP\n(470 mammals)", y = "mean variable position Phylop\n(470 mammals)", fill = "mean PhyloP difference\nbetween variable position\nand position 2") ->
  codon_conservation_versus_2nd_pos

ggsave("GASR/germs/plots/conservation_mutation/codon_conservation_versus_2nd_pos.pdf", codon_conservation_versus_2nd_pos,
       width = 6, height = 4)


# Conservation across classes ---------------------------------------------


mutated_multivalency_conservation %>%
  mutate(ratio_class = cut(ratio, breaks = quantile(ratio, c(0, 0.5, 0.7, 0.9, 0.95, 1)), include_lowest = T,
                           labels = c("0 to 50th", "to 70th", "to 90th", "to 95th", "to 100th")),
         diff_class = cut(diff, breaks = quantile(diff, c(0, 0.2, 0.4, 0.6, 0.8, 0.95, 1)), include_lowest = T),
         mv_smv_ratio = native/normalised_skm) %>%
  drop_na() %>%
  group_by(codon, position_in_codon) %>%
  mutate(ratio_per_codon_class = cut(ratio, breaks = quantile(ratio, c(0, 0.5, 0.7, 0.9, 0.95, 1)), include_lowest = T,
                           labels = c("0 to 50th", "to 70th", "to 90th", "to 95th", "to 100th")),
         ratio_per_codon_class_highvsmid = cut(ratio, breaks = quantile(ratio, c(0, 0.3, 0.7, 0.9, 1)), include_lowest = T,
                                     labels = c("bottom 30%", "middle 40%", "70 to 90th", "top 10%")),
         mv_class = cut(native, quantile(native, c(0, 0.5, 0.9, 1)), include_lowest = T,
                        labels = c("low GeRM", "medium GeRM", "high GeRM")),
         smv_class = cut(normalised_skm, quantile(native, c(0, 0.5, 0.9, 1)), include_lowest = T,
                         labels = c("low GeRM", "medium GeRM", "high GeRM")),
         mv_smv_class = cut(normalised_skm, quantile(native, c(0, 0.33, 0.66, 1)), include_lowest = T,
                            labels = c("lower than local GeRM", "same local GeRM", "higher than local GeRM")),
    p100v_class = cut(scaled_100v, breaks = quantile(scaled_100v, c(0, 0.2, 0.4, 0.6, 0.8, 1)), include_lowest = T,
                           labels = c("0 to 20th", "to 40th", "to 60th", "to 80th", "to 100th")),
         p470m_class = cut(scaled_470m, breaks = quantile(scaled_470m, c(0, 0.2, 0.4, 0.6, 0.8, 1)), include_lowest = T,
                           labels = c("0 to 20th", "to 40th", "to 60th", "to 80th",  "to 100th")),) %>%
  ungroup() %>%
  drop_na() ->
  mmc_classes

ggplot(mmc_classes) +
  aes(x = ratio_per_codon_class, y = scaled_100v) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  coord_cartesian(ylim = c(-4,4)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised PhyloP across 100 vertebrates") ->
  n100v_by_ratioclass_plot
  
ggplot(mmc_classes) +
  aes(x = ratio_per_codon_class, y = scaled_470m) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  # facet_wrap(. ~ mv_class) +
  coord_cartesian(ylim = c(-5,7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised PhyloP across 470 mammals") ->
  n470m_by_ratioclass_plot

mmc_classes %>%
  reframe(blah = pairwise.t.test(x = phylop_100_vertebrates, g = ratio_per_codon_class, p.adjust.method = "BH") %>% tidy()) %>%
  unnest(everything()) ->
  n100v_by_ratioclass_ttests

mmc_classes %>%
  reframe(blah = pairwise.t.test(x = phylop_470_mammals, g = ratio_per_codon_class, p.adjust.method = "BH") %>% tidy()) %>%
  unnest(everything()) ->
  n470m_by_ratioclass_ttests

write_tsv(n100v_by_ratioclass_ttests, "GASR/germs/plots/conservation_mutation/n100v_by_ratioclass_ttests.tsv")
write_tsv(n470m_by_ratioclass_ttests, "GASR/germs/plots/conservation_mutation/n470m_by_ratioclass_ttests.tsv")

ggsave("GASR/germs/plots/conservation_mutation/n100v_by_ratioclass_plot.pdf", plot = ratio_by_norm_100v_plot,
       device = "pdf", units = "in",
       height = 5, width = 4)

ggsave("GASR/germs/plots/conservation_mutation/n470m_by_ratioclass_plot.pdf", plot = n470m_by_ratioclass_plot,
       device = "pdf", units = "in",
       height = 5, width = 4)


# Second position difference ----------------------------------------------

ggplot(mmc_classes) +
  aes(x = ratio_per_codon_class, y = scaled_470m) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ smv_class) +
  coord_cartesian(ylim = c(-5,7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised PhyloP across 470 mammals") ->
  n470m_by_ratioclass_bysmoothedmv_plot

ggplot(mmc_classes) +
  aes(x = smv_class, y = codon2_phylop_470_mammals) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  # facet_wrap(. ~ smv_class) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "local GeRM",
       y = "2nd pos PhyloP across 470 mammals") ->
  pos2_470m_bysmoothedmv_plot

ggplot(mmc_classes) +
  aes(x = smv_class, y = codon2_scaled_470m) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  # facet_wrap(. ~ smv_class) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(-10,10)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "local GeRM",
       y = "normalised 2nd position - mutatable position\n(PhyloP, 470 mammals)") ->
  n470m2ndposdiff_bysmoothedmv_plot

ggplot(mmc_classes) +
  aes(x = ratio_per_codon_class, y = codon2_scaled_470m) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ smv_class) +
  coord_cartesian(ylim = c(-10,10)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised 2nd position - mutatable position\n(PhyloP, 470 mammals)") ->
  n470m2ndposdiff_by_ratioclass_bysmoothedmv_plot

ggsave("GASR/germs/plots/conservation_mutation/n470m_by_ratioclass_bysmoothedmv_plot.pdf", n470m_by_ratioclass_bysmoothedmv_plot,
       width = 8, height = 4)

ggsave("GASR/germs/plots/conservation_mutation/pos2_470m_bysmoothedmv_plot.pdf", pos2_470m_bysmoothedmv_plot,
       width = 4, height = 4)

ggsave("GASR/germs/plots/conservation_mutation/n470m2ndposdiff_bysmoothedmv_plot.pdf", n470m2ndposdiff_bysmoothedmv_plot,
       width = 4, height = 4)

ggsave("GASR/germs/plots/conservation_mutation/n470m2ndposdiff_by_ratioclass_bysmoothedmv_plot.pdf", n470m2ndposdiff_by_ratioclass_bysmoothedmv_plot,
       width = 8, height = 4)

# Low entropy -------------------------------------------------------------

le <- read_tsv("GASR/germs/data/low_entropy_clusters.tsv.gz",
               col_types = cols()) %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end"), sep = ":|-", remove = F, convert = T) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  mutate(start = ((start * 3) - 3) + cds_start, # A region that starts at position 1 should end up with its start being equal to the cds_start.
         end = (end * 3) + cds_start + 2) # A region that ends on the last amino acid should end up with its end being equal to the cds_end - 3 (no stop codon in entropy) 

left_join(mmc_classes, 
          le[c("transcript_id", "start", "end", "cluster", "cluster_name")],
          multiple = "all") %>%
  mutate(is_le = case_when((position >= start) & (position <= end) ~ "LCD",
                           T ~ "non-LCD") %>%
           fct_relevel("non-LCD")) %>%
  mutate(cluster_type = case_when(is_le == "LCD" ~ cluster_name,
                                  T ~ "non-LCD") %>%
           fct_relevel("non-LCD")) %>%
  select(-c(start, end, cluster, cluster_name)) %>%
  distinct() ->
  mmc_le

ggplot(mmc_le) +
  aes(x = ratio_per_codon_class, y = codon2_scaled_100v) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ is_le) +
  coord_cartesian(ylim = c(-10,10)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised 2nd position - mutatable position\n(PhyloP, 100 vertebrates)") ->
  lcd_2ndposnormalised_conservation_by_ratio

ggsave("GASR/germs/plots/conservation_mutation/lcd_2ndposnormalised_conservation_by_ratio_wide.pdf", lcd_2ndposnormalised_conservation_by_ratio,
       width = 7.5, height = 3)

ggsave("GASR/germs/plots/conservation_mutation/lcd_2ndposnormalised_conservation_by_ratio_thin.pdf", lcd_2ndposnormalised_conservation_by_ratio,
       width = 5, height = 3)



ggplot(mmc_le) +
  aes(x = ratio_per_codon_class, y = phylop_100_vertebrates) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ is_le) +
  coord_cartesian(ylim = c(-5,7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "PhyloP, 100 vertebrates)") ->
  lcd_mutable_conservation_by_ratio_vertebrates

mmc_le %>%
  group_by(is_le) %>%
  reframe(ttests = pairwise.t.test(phylop_100_vertebrates, ratio_per_codon_class) %>% tidy()) %>%
  unnest(ttests) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  lcd_mutable_conservation_by_ratio_vertebrates_ttests

write_tsv(lcd_mutable_conservation_by_ratio_vertebrates_ttests, 
          "GASR/germs/plots/conservation_mutation/lcd_mutable_conservation_by_ratio_vertebrates_ttests.tsv")

ggsave("GASR/germs/plots/conservation_mutation/lcd_mutable_conservation_by_ratio_vertebrates.pdf", lcd_mutable_conservation_by_ratio_vertebrates,
       width = 5, height = 3)

ggplot(mmc_le) +
  aes(x = ratio_per_codon_class, y = phylop_470_mammals) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ is_le) +
  coord_cartesian(ylim = c(-5,7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "PhyloP, 470 mammals)")

ggplot(mmc_le) +
  aes(x = is_le, y = codon2_phylop_100_vertebrates) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  # facet_wrap(. ~ smv_class) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(-4, 12)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "local GeRM",
       y = "middle pos PhyloP, 100 vertebrates") ->
  lcd_codon2_phylop_100_vertebrates_plot

ggsave("GASR/germs/plots/conservation_mutation/lcd_codon2_phylop_100_vertebrates_plot.pdf", lcd_codon2_phylop_100_vertebrates_plot,
       width = 2, height = 3)

wilcox.test(codon2_phylop_100_vertebrates ~ is_le, data = mmc_le) %>% tidy() -> lcd_codon2_phylop_100_vertebrates_wilcoxon_test

write_tsv(lcd_codon2_phylop_100_vertebrates_wilcoxon_test, "GASR/germs/plots/conservation_mutation/lcd_codon2_phylop_100_vertebrates_wilcoxon_test.tsv")


ggplot(mmc_le) +
  aes(x = ratio_per_codon_class, y = codon2_scaled_470m) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.3) +
  facet_wrap(. ~ is_le) +
  coord_cartesian(ylim = c(-10,10)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "percentile of ratio\nnative to mutated GeRM",
       y = "normalised 2nd position - mutatable position\n(PhyloP, 470 mammals)") ->
  lcd_2ndposnormalised_conservation_by_ratio_mammals

ggsave("GASR/germs/plots/conservation_mutation/lcd_2ndposnormalised_mammalconservation_by_ratio_thin.pdf", lcd_2ndposnormalised_conservation_by_ratio_mammals,
       width = 5, height = 3)

mmc_le %>%
  group_by(is_le) %>%
  reframe(ttests = pairwise.t.test(codon2_scaled_100v, ratio_per_codon_class) %>% tidy()) %>%
  unnest(ttests) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  lcd_2ndposnormalised_conservation_by_ratio_vertebrates_ttests

write_tsv(lcd_2ndposnormalised_conservation_by_ratio_vertebrates_ttests, 
          "GASR/germs/plots/conservation_mutation/lcd_2ndposnormalised_conservation_by_ratio_vertebrates_ttests.tsv")

mmc_le %>%
  group_by(is_le) %>%
  reframe(ttests = pairwise.t.test(codon2_scaled_470m, ratio_per_codon_class) %>% tidy()) %>%
  unnest(ttests) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  lcd_2ndposnormalised_conservation_by_ratio_mammals_ttests

write_tsv(lcd_2ndposnormalised_conservation_by_ratio_mammals_ttests, 
          "GASR/germs/plots/conservation_mutation/lcd_2ndposnormalised_conservation_by_ratio_mammals_ttests.tsv")
