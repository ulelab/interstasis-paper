rm(list = ls())

library(tidyverse)
library(broom)

# Load --------------------------------------------------------------------

proteomics_H <- read_tsv("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/neve_reanalysed_trypsin_H.txt")

proteomics_H %>%
  rename(uniprot = Accession) %>%
  pivot_longer(cols = -1) %>%
  replace(is.na(.), 0) %>%
  mutate(name = name %>%
           str_replace_all(fixed("Abundances (Normalized): "), "") %>%
           str_replace_all(fixed(": "), "_") %>%
           str_replace_all(fixed(", Sample"), "") %>%
           str_replace_all(fixed("Heavy"), "labelled") %>%
           str_replace_all(fixed("Medium"), "labelled") %>%
           str_replace_all(fixed("Light"), "unlabelled")) ->
  proteomics

ppig_int <- read_tsv("/nemo/lab/ulej/home/users/farawar/GASR/export_reporter/quantseq/neve_timecourse_1/tables/er_int_with_classes.tsv.gz") %>%
  select(gene_name, purine_multivalency_class = score_class, purine_multivalency_score = scaled_log2_max_score)

tibble(sample = paste0("F", c(1:16)),
       timepoint = c("4hr", "8hr") %>% rep(each = 8),
       dox = c("dox", "control") %>% rep(each = 4) %>% rep(2),
       replicate = rep(c(1:4), 4)) ->
  sample_info

# correct swap (can be seen looking at mscar abundances)
sample_info$dox[12] <- "control"
sample_info$dox[16] <- "dox"
sample_info$dox[4] <- "control"
sample_info$dox[8] <- "dox"

# Get gene names for uniprot ----------------------------------------------

mscar_tibble <- tibble(uniprot = "Q00000", gene_name = "mScarlet")

uniprot_to_gene_name_without_mscar <- read_tsv("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/uniprot_to_gene_name.tsv")

list(mscar_tibble, uniprot_to_gene_name_without_mscar) %>%
  bind_rows() ->
  uniprot_to_gene_name

# Make big table ----------------------------------------------------------

proteomics %>%
  mutate(uniprot = uniprot %>% word(1, sep = "-")) %>%
  inner_join(uniprot_to_gene_name, relationship = "many-to-many") %>%
  select(-uniprot) %>%
  group_by(gene_name, name) %>%
  summarise(abundance = sum(value)) %>%
  ungroup() %>%
  separate(name, into = c("sample", "labelling")) %>%
  left_join(sample_info) ->
  annotated_proteomics

annotated_proteomics %>%
  pivot_wider(names_from = labelling, values_from = abundance) %>%
  mutate(total_abundance = (labelled + unlabelled),
         proportion_labelled = labelled / total_abundance) ->
  annotated_ratios

annotated_ratios %>%
  group_by(gene_name) %>%
  summarise(total_abundance = sum(unlabelled + labelled)) ->
  total_abundance

annotated_ratios %>%
  group_by(sample) %>%
  summarise(total_sample_abundance = sum(unlabelled + labelled)) %>%
  mutate(relative_sample_abundance = total_sample_abundance / mean(total_sample_abundance)) ->
  total_sample_abundance

annotated_ratios %>%
  left_join(total_sample_abundance) %>%
  mutate(total_abundance_normalised = total_abundance / relative_sample_abundance,
         unlabelled_normalised = unlabelled / relative_sample_abundance) ->
  annotated_ratios_normalised

annotated_ratios %>%
  drop_na() %>%
  filter(unlabelled + labelled > 0) %>%
  count(gene_name) %>%
  filter(n >= 16) %>%
  select(gene_name) %>% unlist(use.names = F) ->
  present_in_all

annotated_ratios_classes <- annotated_ratios %>% inner_join(ppig_int)
annotated_ratios_normalised_classes <- annotated_ratios_normalised %>% inner_join(ppig_int)

annotated_ratios_normalised_classes %>%
  filter(gene_name %in% present_in_all) %>% 
  select(sample, timepoint, dox, gene_name, purine_multivalency_score, purine_multivalency_class,proportion_labelled) %>%
  ggplot(aes(y = proportion_labelled, x = purine_multivalency_class, color = dox)) +
  geom_boxplot() +
  facet_wrap(. ~ timepoint) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_colour_manual(values = c("gray30", "purple4"))

# Just 8 hr ---------------------------------------------------------------

annotated_ratios_classes %>%
  filter(timepoint == "8hr",
         unlabelled + labelled > 0) %>%
  drop_na() %>%
  count(gene_name) %>%
  filter(n >= 8) %>%
  select(gene_name) %>% unlist(use.names = F) ->
  present_in_all_8

annotated_ratios_classes %>%
  filter(timepoint == "8hr") %>%
  filter(dox == "control",
         labelled > 0) %>%
  count(gene_name) %>%
  filter(n >= 4) %>%
  select(gene_name) %>%
  unlist(use.names = F) ->
  labelled_in_control

annotated_ratios_normalised_classes %>%
  filter(gene_name %in% present_in_all_8,
         gene_name %in% labelled_in_control,
         timepoint == "8hr") %>% 
  group_by(purine_multivalency_class, dox, replicate) %>%
  summarise(mean_ratio = mean(proportion_labelled, na.rm = T)) ->
  mean_ratios_df

annotated_ratios_normalised_classes %>%
  filter(gene_name %in% present_in_all_8,
         gene_name %in% labelled_in_control,
         timepoint == "8hr") %>% 
  select(sample, timepoint, dox, gene_name, purine_multivalency_score, purine_multivalency_class,proportion_labelled) %>%
  ggplot(aes(y = proportion_labelled, x = purine_multivalency_class, color = dox)) +
  # geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.5), outlier.shape = NA, alpha = 0.5, show.legend = F) +
  geom_point(aes(y = mean_ratio), data = mean_ratios_df, size = 3, shape = 21,
             position = position_dodge2(width = 0.5), show.legend = F) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_colour_manual(values = c("gray30", "purple4"))

annotated_ratios_classes %>%
  filter(
    gene_name %in% present_in_all_8,
         gene_name %in% labelled_in_control,
         timepoint == "8hr") %>% 
  select(-c(sample, unlabelled, labelled, timepoint, total_abundance)) %>%
  pivot_wider(names_from = dox, values_from = proportion_labelled) %>%
  mutate(change_in_proportion = dox - control,
         # Adding 1% of the control to deal with infinite sized change in labelling
         ratio_change_proportion = (dox + (control/100)) / control) ->
  change_in_ratios

change_in_ratios %>%
  group_by(purine_multivalency_class, replicate) %>%
  summarise(mean_ratio = mean(ratio_change_proportion, na.rm = T)) ->
  mean_ratio_change_df

pairwise.t.test(log2(mean_ratio_change_df$mean_ratio), mean_ratio_change_df$purine_multivalency_class, p.adjust.method = "BH") %>%
  tidy() ->
  change_in_silac_ratio_ga_multivalency_ttests

ggplot(change_in_ratios %>%
         group_by(gene_name, purine_multivalency_class) %>%
         summarise(ratio_change_proportion = mean(ratio_change_proportion))) +
  aes(y = ratio_change_proportion, x = purine_multivalency_class, fill = purine_multivalency_class) +
  # geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5, show.legend = F) +
  geom_point(aes(y = mean_ratio), data = mean_ratio_change_df, size = 3, shape = 21,
             position = position_dodge2(width = 0.5), show.legend = F) +
  theme_classic() +
  # scale_y_log10() +
  coord_cartesian(ylim = c(0, 1.6)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GA multivalency", y = "change in SILAC labelling\n(dox / control)", fill = "", color = "") +
  scale_fill_manual(values = c("#b4b4b4", "#aa94b1", "#b9529f")) ->
  change_in_silac_ratio_ga_multivalency_plot

change_in_ratios %>%
  group_by(gene_name, purine_multivalency_class) %>%
  summarise(ratio_change_proportion = mean(ratio_change_proportion)) ->
  ratio_mean_per_replicate

ggplot(ratio_mean_per_replicate) +
  aes(y = ratio_change_proportion, x = purine_multivalency_class, fill = purine_multivalency_class) +
  # geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.75, outlier.shape = NA, alpha = 0.5, show.legend = F) +
  geom_point(aes(y = mean_ratio), data = mean_ratio_change_df, size = 3, shape = 21,
             position = position_dodge2(width = 0.5), show.legend = F) +
  theme_classic() +
  # scale_y_log10() +
  coord_cartesian(ylim = c(0.1, 1.7)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  labs(x = "", y = "change in SILAC labelling\n(dox / control)", fill = "", color = "") +
  scale_fill_manual(values = c("#b4b4b4", "#aa94b1", "#b9529f")) ->
  change_in_silac_ratio_ga_multivalency_plot_dense

ggsave("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/change_in_silac_ratio_ga_multivalency_plot.pdf",
       change_in_silac_ratio_ga_multivalency_plot, width = 3, height = 3)

ggsave("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/change_in_silac_ratio_ga_multivalency_plot_dense.pdf",
       change_in_silac_ratio_ga_multivalency_plot_dense, width = 1.5, height = 3)

write_tsv(change_in_silac_ratio_ga_multivalency_ttests, "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/change_in_silac_ratio_ga_multivalency_ttests.tsv")

write_tsv(mean_ratio_change_df, "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/change_in_silac_ratio_ga_multivalency_groupmean_data.tsv")
write_tsv(ratio_mean_per_replicate,  "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/change_in_silac_ratio_ga_multivalency_replicatemean_data.tsv")

# Proportion of proteome ----------------------------------------------------

annotated_ratios_normalised_classes %>%
  filter(timepoint == "8hr",
         gene_name %in% present_in_all_8) %>%
  group_by(dox, replicate, purine_multivalency_class) %>%
  # summarise(summed_abundance = sum(total_abundance_normalised)) %>%
  summarise(summed_abundance = sum(labelled)) %>%
  group_by(dox, replicate) %>%
  mutate(proportion_proteome = summed_abundance / sum(summed_abundance)) ->
  summed_abundance_df

ggplot(summed_abundance_df, aes(x = purine_multivalency_class, y = proportion_proteome, fill = dox)) +
  geom_point(position = position_dodge2(width = 0.5), shape = 21, size = 3) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = c("#b4b4b4", "#b6549a")) +
  coord_cartesian(ylim = c(0.005, 1)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GA multivalency", y = "proportion of proteome\n(8hr SILAC labelling)", fill = "") ->
  proportion_of_proteome_ga_multivalency_plot

summed_abundance_df %>%
  group_by(purine_multivalency_class) %>%
  summarise(ttest = t.test(proportion_proteome ~ dox) %>% tidy()) %>%
  unnest(ttest) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  proportion_of_proteome_ga_multivalency_ttest

ggsave("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/proportion_of_proteome_ga_multivalency_plot.pdf",
       proportion_of_proteome_ga_multivalency_plot, width = 4, height = 3)

write_tsv(proportion_of_proteome_ga_multivalency_ttest, "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/proportion_of_proteome_ga_multivalency_ttest.tsv")

write_tsv(summed_abundance_df, "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/proportion_of_proteome_ga_multivalency_data.tsv")

# MCD proportions ---------------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())
aa_per_region_df <- read_tsv("GASR/sequence_multivalency/data/20220315_123win_98per_k5_10nei/data/percentile_98_windowsize_123_klen_5_le_regions_with_aa_props.txt.gz")

aa_per_region_df %>%
  mutate(net_charge = (K+R) - (E+D),
         total_charge = K+R+E+D,
         aa_len = protein_region_end - protein_region_start) %>% 
  filter(total_charge >= 0.4, R >= 0.2,  net_charge > 0, aa_len >= 50) %>%
  select(transcript_id) %>%
  distinct() %>%
  left_join(transcript_details) %>%
  mutate(ensg_id = word(gene_id, 1, sep = fixed("."))) ->
  mcds

annotated_ratios_normalised %>%
  filter(
    # timepoint == "8hr",
    gene_name %in% c("mScarlet", present_in_all_8)) %>%
  mutate(mcd = case_when(gene_name %in% mcds$gene_name ~ "R-MCD",
                         gene_name == "mScarlet" ~ "mScarlet",
                         T ~ "no MCD") %>%
           fct_relevel("no MCD", "R-MCD")) %>%
  group_by(dox, timepoint, replicate, mcd) %>%
  # summarise(summed_abundance = sum(total_abundance_normalised)) %>%
  summarise(summed_abundance = sum(unlabelled_normalised)) %>%
  group_by(dox, timepoint, replicate) %>%
  mutate(proportion_proteome = summed_abundance / sum(summed_abundance)) %>%
  filter(mcd != "no MCD") ->
  summed_abundance_mcd_df

summed_abundance_mcd_df %>%
  filter(dox == "dox", timepoint == "8hr") %>%
  ggplot(aes(x = mcd, y = proportion_proteome, fill = mcd)) +
  geom_point(position = position_dodge2(width = 0.5), shape = 21, size = 3, show.legend = F) +
  # scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = c("#b6549a", "#ed1c96")) +
  coord_cartesian(ylim = c(0, .015)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(y = "proportion of total proteome",
       x = "") ->
  mscarlet_vs_rmcd_plot

write_tsv(summed_abundance_mcd_df %>%
            filter(dox == "dox", timepoint == "8hr"),
          "/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/mscarlet_vs_rmcd.tsv")

ggsave("/nemo/lab/ulej/home/users/farawar/GASR/SILAC/neve_two/plots/mscarlet_vs_rmcd_plot.pdf",
       mscarlet_vs_rmcd_plot, width = 2, height = 3)


