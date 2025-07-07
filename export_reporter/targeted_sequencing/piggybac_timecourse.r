rm(list = ls())

library(tidyverse)
library(stringdist)
library(Biostrings)
library(scales)
library(furrr)

library(ggridges)
library(ggpubr)
library(patchwork)
library(broom)

# library(ggridges)
# library(apcluster)


# Barcodes ----------------------------------------------------------------

consensus_barcodes <-  read_tsv("GASR/export_reporter/nanopore/minimap2/processed_data/plasmid_consensus_barcodes.tsv.gz")

consensus_barcodes$rc_barcode <- consensus_barcodes$barcode %>% DNAStringSet() %>% reverseComplement() %>% as.character()

# Reads -------------------------------------------------------------------

segmented_filtered <- read_tsv("GASR/export_reporter/targeted_sequencing/data/third_experiment_barcodes_upis.txt.gz")
# segmented_filtered <- segmented_filtered %>% 
#   mutate(sample = paste0(experiment, "_", sample))

barcode_occurence <- table(segmented_filtered$plasmid_bc) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>%
  set_names(c("plasmid_bc", "frequency"))

# Inflection plot ---------------------------------------------------------

barcode_occurence %>%
  mutate(rn = row_number(),
         cs = cumsum(frequency)) %>%
  mutate(csp = cs / max(cs)) %>%
  ggplot(aes(x = rn, y = csp)) +
  geom_line() +
  geom_hline(yintercept = 0.97, linetype = "dashed") +
  # geom_hline(yintercept = nrow(segmented_filtered) * 0.25, linetype = "dashed") +
  geom_vline(xintercept = 1600, linetype = "dashed") +
  # geom_vline(xintercept = 2000, linetype = "dashed") +
  theme_classic() +
  scale_x_continuous(breaks = breaks_pretty(5),
                     labels = comma) +
  coord_cartesian(xlim = c(0, 10000)) +
  labs(x = "unique barcodes",
       y = "cumulative proportion")


# Annotate  ---------------------------------------------------------------

barcode_occurence[1:1600,] %>%
  mutate(matches = future_map(plasmid_bc, ~{
    matches <- which(stringdist(.x, consensus_barcodes$rc_barcode) < 4)
  })) %>%
  mutate(match_number = map(matches, ~{ length(.x) }) %>% unlist()) %>%
  filter(match_number > 0) %>%
  mutate(archi = future_map(matches, ~{
    # .x <- c(7766, 7767, 7768, 7769, 7770, 7771, 12543, 10)
    
    architectures <- table(paste0(consensus_barcodes$opt_pattern[.x], "@", consensus_barcodes$intron_pattern[.x])) %>%
      as.data.frame() %>%
      set_names(c("architecture", "frequency")) %>%
      arrange(desc(frequency))
    
    if(nrow(architectures) == 1) {
      return(architectures$architecture[1])
    }
    
    if(architectures$frequency[1] >= (architectures$frequency[2] * 2)) {
      return(architectures$architecture[1])
    } else {
      return("unclear")
    }
  })) %>%
  unnest(archi) %>%
  filter(archi != "unclear") %>%
  select(-c(matches, frequency)) %>%
  separate(archi, into = c("opt_pattern", "intron_pattern"), sep = "@") ->
  barcode_architectures

barcode_architectures %>%
  mutate(bc_number = row_number()) %>%
  mutate(my_fun_df = future_pmap(.l = list(opt_pattern = opt_pattern,
                                           intron_pattern = intron_pattern),
                                 .f = function(opt_pattern, intron_pattern){
                                   opt_pattern %>%
                                     str_split("-") %>%
                                     unlist() %>%
                                     str_detect("Opt") ->
                                     opt_presence
                                   
                                   intron_pattern %>%
                                     str_split("-") %>%
                                     unlist() %>% 
                                     as.numeric() %>% 
                                     as.logical() ->
                                     intron_presence
                                   
                                   tibble(position = c(1:8),
                                          exon = opt_presence,
                                          intron = intron_presence) 
                                 }
  )) %>%
  dplyr::select(plasmid_bc, bc_number, match_number, my_fun_df) %>%
  unnest(my_fun_df) ->
  new_barcode_to_struct
  
write_tsv(new_barcode_to_struct, "GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/piggybac_timecourse_barcode_architectures.tsv.gz")
new_barcode_to_struct <- read_tsv("GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/piggybac_timecourse_barcode_architectures.tsv.gz")

# Analyse -----------------------------------------------------------------

usable_counts_per_sample <- segmented_filtered %>% 
  filter(plasmid_bc %in% barcode_architectures$plasmid_bc) %>%
  count(experiment, sample, fraction, replicate)

count_usable_upis <- segmented_filtered %>%
  mutate(sample = paste0(experiment, "_", sample)) %>%
  # sample_frac(0.01) %>%
  filter(plasmid_bc %in% barcode_architectures$plasmid_bc) %>%
  count(experiment, sample, exp_bc, replicate, fraction, plasmid_bc) %>%
  group_by(sample) %>%
  mutate(cpm = (n/sum(n)) * 1e6) %>%
  ungroup() %>%
  mutate(bc_number = match(plasmid_bc, barcode_architectures$plasmid_bc)) %>%
  group_by(bc_number) %>%
  mutate(min_cpm = min(cpm),
         min_counts = min(n)) %>%
  ungroup()

# PCA ---------------------------------------------------------------------

for_pca <- count_usable_upis %>% 
  ungroup() %>% 
  arrange(desc(min_cpm)) %>% 
  filter(min_counts >= 100) %>%
  # mutate(sample = paste0(experiment, "_", sample)) %>%
  select(sample, cpm, bc_number) %>% 
  mutate(bc_number = paste0("index_", bc_number)) %>%
  pivot_wider(names_from = bc_number, values_from = "cpm") %>%
  select_if(~ !any(is.na(.))) %>%
  column_to_rownames("sample")

pca_of_samples <- prcomp(for_pca, scale = T) 

prop_variance <- summary(pca_of_samples)$importance[2,]

pca_plot <- pca_of_samples$x %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  separate(sample, into = c("experiment", "fraction", "replicate"), remove = F) %>%
  ggplot(aes(x = -PC1, y = PC2, colour = fraction, shape = experiment)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = paste0("PC1: ", round(prop_variance[1], digits = 3) * 100, "% of variance"),
       y = paste0("PC2: ", round(prop_variance[2], digits = 3) * 100, "% of variance"),
       color = "",
       shape = "") +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))


# Ridgelining -------------------------------------------------------------

mv_scores <- read_tsv(file = "GASR/export_reporter/ppig_fragments/gene_structure_multivalency_123.txt")
mv_scores_masked <- read_tsv(file = "GASR/export_reporter/ppig_fragments/gene_structure_multivalency_and_intronic_masking_123.txt.gz")

upi_mv <- new_barcode_to_struct %>%
  mutate(opt = case_when(exon ~ "Opt", T ~ "DeO")) %>%
  mutate(int_pattern = case_when(intron ~ "Int", T ~ "NoInt")) %>%
  select(bc_number, position, opt, int_pattern) %>%
  group_by(bc_number) %>%
  summarise(structure = paste0(opt, collapse = "-"),
            int_pattern = paste0(int_pattern, collapse = "-")) %>%
  left_join(mv_scores_masked %>% select(-permutation), by = c("structure", "int_pattern"))

ratio_per_replicate <- count_usable_upis %>%
  filter(min_counts >= 200) %>%
  select(experiment, fraction, replicate, cpm, bc_number) %>%
  pivot_wider(names_from = fraction, values_from = cpm) %>%
  mutate(nc_ratio = nucleus/cytoplasm) 

l_plus_ratio <- new_barcode_to_struct %>%
  group_by(bc_number) %>%
  summarise(n_introns = sum(intron),
            n_opt = sum(exon)) %>%
  inner_join(ratio_per_replicate, by = "bc_number") %>%
  left_join(upi_mv) %>%
  drop_na() %>%
  mutate(mv_class = cut(purine_multivalency, 
                        breaks = quantile(purine_multivalency, 
                                          probs = c(0, 0.1, 0.3, 0.6, 0.9, 1)),
                        labels = c("< 10th", "10th - 30th",
                                   "30th - 60th", 
                                   "60th - 90th", 
                                   "90th <"), 
                        include.lowest = T),
         experiment = experiment %>%
           str_replace("h12", "12h dox") %>%
           str_replace("h8", "8h dox") %>%
           str_replace("h4", "4h dox") %>%
           fct_relevel("4h dox", "8h dox")) 

write_tsv(l_plus_ratio, "GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/piggybac_timecourse_gene_ratios.tsv.gz")
l_plus_ratio <- read_tsv("GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/piggybac_timecourse_gene_ratios.tsv.gz") %>%
  mutate(experiment = experiment %>%
           fct_relevel("4h dox", "8h dox"))

ridgeline <- ggplot(l_plus_ratio, 
                    aes(x = log2(nc_ratio),
                        y = mv_class,
                        fill = experiment,
                        # y = factor(n_introns),
                        # fill = experiment,
                    )
) + 
  # facet_wrap(. ~ experiment) +
  geom_density_ridges2(alpha = 0.5, rel_min_height = 0.003, scale = 2) + 
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  theme_classic() +
  coord_cartesian(xlim = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "log2 (nucleus / cytoplasm)", 
       y = "multivalency percentile",
       fill = "") +
  
  theme(axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"))

ridgeline_nolegend <- ridgeline + theme(legend.position = "none")

ridgeline

group_means <- l_plus_ratio %>%
  drop_na() %>%
  group_by(experiment, mv_class, replicate) %>%
  summarise(group_mean_ratio = mean(nc_ratio)) %>%
  ungroup() %>%
  mutate(mv_class = ordered(mv_class)) %>%
  drop_na()

group_means %>%
  group_by(mv_class) %>%
  reframe(model = pairwise.t.test(log2(group_mean_ratio), experiment) %>% tidy()) %>% 
  unnest(model) %>%
  # filter(!((group1 == "12hr dox") & (group2 == "8hr dox"))) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  pairwise_t_tests

write_tsv(pairwise_t_tests, file = "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_timecourse_group_means_multivalency_ttests.tsv")

group_means_plot <- group_means %>%
  ggplot() +
  aes(x = fct_expand(mv_class, ""),
      y = log2(group_mean_ratio),
      fill = experiment,
      group = factor(replicate)) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "multivalency percentile",
       y = "log2 (nucleus / cytoplasm)",
       fill = ""
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black"),) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(drop = F) +
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  coord_flip()

combi_plot <- (ridgeline_nolegend | group_means_plot) + plot_layout(widths = c(4, 2))


# Intron + Multivalency ---------------------------------------------------

l_plus_ratio %>%
  mutate(intron_rich = case_when(n_introns >= 6 ~ "6+ introns",
                                 n_introns >= 4 ~ "4-6 introns",
                                 T ~ "1-3 introns",) %>%
           fct_relevel("1-3 introns", "4-6 introns",)) %>%
  mutate(highly_multivalent = case_when(purine_multivalency >= quantile(purine_multivalency, 0.5) ~ "high multivalency",
                                        # purine_multivalency >= quantile(purine_multivalency, 0.33) ~ "medium GeRM",
                                        T ~ "low multivalency") %>%
           fct_relevel("low multivalency")) %>%
  group_by(intron_rich, highly_multivalent, experiment, replicate) %>%
  summarise(nc_ratio = mean(nc_ratio)) %>%
  ungroup() ->
  intron_multivalency_comparison

intron_multivalency_comparison %>%
  nest(data = -c(experiment, highly_multivalent)) %>%
  mutate(ttest = map(data, ~ { pairwise.t.test(log2(.x$nc_ratio), .x$intron_rich) %>% tidy()})) %>%
  unnest(ttest) %>% 
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  select(-data) ->
  intron_multivalency_ttests

write_tsv(intron_multivalency_ttests, file = "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_timecourse_intron_multivalency_ttests.tsv")

ggplot(intron_multivalency_comparison) +
  aes(x = experiment, y = log2(nc_ratio), fill = intron_rich) +
  geom_point(shape = 21, size = 3, position = position_dodge2(width = 0.9, padding = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = (c("#d1c8e3", "#a4a4ba", "#666d7d"))) +
  # scale_fill_viridis_d(option = "plasma", begin = 0.4, end = 0.8, direction = -1) +
  facet_wrap(. ~ highly_multivalent) +
  theme_classic() +
  # coord_cartesian(ylim = c(-0.3, 0.4)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  labs(y = "log2 (nucleus / cytoplasm)",
       x = "expression time",
       fill = "") ->
  intron_multivalency_plot

intron_multivalency_plot

(ridgeline_nolegend | group_means_plot | intron_multivalency_plot) + plot_layout(widths = c(3, 2, 4)) ->
  big_combi

ggsave("GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_timecourse_main_figure_combiplot.pdf", big_combi,
       width = 9, height = 3.5)

# Example plasmid ---------------------------------------------------------



# Reviewer thing -------------------------------------------------------------------

l_plus_ratio %>%
  separate(int_pattern, into = c(paste0("intpos_", c(1:8)))) %>%
  pivot_longer(cols = contains("intpos"), names_to = "intron_pos", values_to = "intron_presence") %>% 
  group_by(experiment, replicate, intron_pos, intron_presence) %>%
  summarise(mean_nc = mean(nc_ratio)) %>%
  ungroup() %>%
  mutate(intron_pos = intron_pos %>% str_remove("intpos_"),
         intron_presence = case_when(intron_presence == "Int" ~ "Intron",
                                     T ~ "No intron")) ->
  mean_nc_per_intron

l_plus_ratio %>%
  group_by(experiment, replicate) %>%
  summarise(mean_c = mean(cytoplasm),
            mean_n = mean(nucleus),
            mean_nc = mean(nc_ratio),)

ggplot(mean_nc_per_intron %>%
         mutate(intron_presence = intron_presence %>% fct_relevel("No intron")),
       aes(x = intron_pos, y = log2(mean_nc), fill = intron_presence, shape = factor(replicate))) +
  geom_point(position = position_dodge2(width = 0.75), size = 2) +
  facet_wrap(. ~ experiment, ncol = 1) +
  labs(x = "intron position",
       y = "log2 (nucleus / cytoplasm)",
       fill = "") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_shape_manual(values = c(21, 21, 21), guide = "none") +
  theme(legend.position = "top",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")) +
  coord_cartesian(ylim = c(-0.08, 0.16)) +
  scale_fill_manual(values = (c("#d1c8e3", "#666d7d"))) ->
  nc_per_intron_position_plot


# Find an example ---------------------------------------------------------

l_plus_ratio %>%
  select(bc_number, structure) %>%
  distinct() %>%
  group_by(structure) %>%
  summarise(number_architectures = dplyr::n()) %>%
  arrange(desc(number_architectures)) ->
  multiple_architectures

l_plus_ratio %>%
  filter(structure == multiple_architectures$structure[49]) %>% 
  ggplot(aes(x = factor(n_introns + 1), y = log2(nc_ratio), fill = experiment)) +
  geom_point(shape = 21, position = position_dodge2(width = 0.9), size = 3) +
  labs(x = "number of exons",
       fill = "",
       y = "log2 (nucleus / cytoplasm)") +
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) 

l_plus_ratio %>%
  filter(structure == multiple_architectures$structure[49]) %>% 
  ggplot(aes(x = factor(n_introns + 1), y = log2(nc_ratio), fill = experiment)) +
  geom_point(shape = 21, position = position_dodge2(width = 0.9), size = 3) +
  labs(x = "number of exons",
       fill = "",
       y = "log2 (nucleus / cytoplasm)") +
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) ->
  example_reporter_gene

ggsave("GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/example_gene_DeO-Opt-Opt-Opt-Opt-Opt-DeO-DeO.pdf", example_reporter_gene,
       width = 3, height = 3)


aov(log2(nc_ratio) ~ experiment * factor(n_introns), 
    data = l_plus_ratio %>%
      filter(structure == multiple_architectures$structure[49])) %>%
  tidy() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  example_anova

write_tsv(example_anova, "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/example_gene_anova.tsv")

# Intron only analysis ---------------------------------------------------------


int_ridgeline <- ggplot(l_plus_ratio %>%
                          filter(n_introns > 1, n_introns < 8) %>% # Not enough reporters available
                          drop_na(), 
                        aes(x = log2(nc_ratio),
                            y = factor(n_introns),
                            fill = experiment,
                        )
) + 
  geom_density_ridges2(alpha = 0.66, rel_min_height = 0.003, scale = 2) + 
  # scale_fill_viridis_d(option = "cividis") +
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  # scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  coord_cartesian(xlim = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "log2 (nucleus / cytoplasm)", 
       y = "number of introns",
       fill = "") +
  theme(#legend.position = "none",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"))

int_ridgeline_nolegend <- int_ridgeline + theme(legend.position = "none")

int_ridgeline

int_group_means <- l_plus_ratio %>%
  filter(n_introns > 1, n_introns < 8) %>% # Not enough reporters available
  group_by(experiment, n_introns, replicate) %>%
  summarise(group_mean_ratio = mean(nc_ratio)) %>%
  ungroup() %>%
  drop_na() 

int_group_means_plot <- int_group_means %>%
  ggplot() +
  aes(x = fct_expand(factor(n_introns), ""),
      y = log2(group_mean_ratio),
      fill = experiment,
      group = factor(replicate)) +
  geom_point(position = position_dodge(0.5),
             size = 3, shape = 21) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "multivalency percentile",
       y = "log2 (nucleus / cytoplasm)",
       fill = ""
  ) +
  theme(#legend.position = "none",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black"),) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  scale_x_discrete(drop = F) +
  # scale_fill_viridis_d(option = "cividis") +
  scale_fill_manual(values = c("#4E4E4E", "#EE9B21", "#CD2027")) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  coord_flip()

int_combi_plot <- (int_ridgeline_nolegend | int_group_means_plot) + plot_layout(widths = c(4, 2))

ggsave("GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_timecourse_treatment_intron_combiplot.pdf", int_combi_plot,
       width = 5, height = 3.5)

# 
# ggsave("GASR/export_reporter/targeted_sequencing/plots/third_exp_piggybac/introns_ridge_and_group_shortlong.pdf",
#        int_combi_plot, device = "pdf", width = 6, height = 3, units = "in")

lm(log2(group_mean_ratio) ~ scale(n_introns) * experiment_time,
   data = int_group_means %>%
     mutate(experiment_time = experiment %>%
              str_remove("h dox") %>%
              as.numeric())) %>%
  tidy() ->
  intron_number_regression

write_tsv(intron_number_regression, file = "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_intron_number_regression.tsv")

int_group_means %>%
  group_by(n_introns) %>%
  reframe(model = pairwise.t.test(log2(group_mean_ratio), experiment) %>% tidy()) %>% 
  unnest(model) %>%
  # filter(!((group1 == "12hr dox") & (group2 == "8hr dox"))) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  pairwise_t_tests_intron

write_tsv(pairwise_t_tests, file = "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/piggybac_timecourse_intron_number_regression.tsv")


# Per position intron effect ----------------------------------------------

l_plus_ratio %>%
  filter(experiment == "12h dox") %>%
  separate(int_pattern, into = c(paste0("intpos_", c(1:8)))) %>%
  pivot_longer(cols = contains("intpos"), names_to = "intron_pos", values_to = "intron_presence") %>% 
  group_by(replicate, intron_pos, intron_presence) %>%
  summarise(mean_nc = mean(nc_ratio)) %>%
  ungroup() %>%
  mutate(intron_pos = intron_pos %>% str_remove("intpos_"),
         intron_presence = case_when(intron_presence == "Int" ~ "Intron",
                                     T ~ "No intron")) ->
  mean_nc_per_intron

ggplot(mean_nc_per_intron,
       aes(x = intron_pos, y = log2(mean_nc), color = intron_presence, shape = factor(replicate))) +
  geom_point(position = position_dodge2(width = 0.75)) +
  labs(x = "intron position",
       y = "log2 (nucleus / cytoplasm)",
       color = "") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_shape_manual(values = c("circle", "circle", "circle"), guide = "none") +
  theme(#legend.position = "none",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")) +
  scale_color_manual(values= c("red3", "blue4")) +
  coord_cartesian(ylim = c(-0.1, 0.12)) ->
  nc_per_intron_position_plot

ggsave("GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/12hr_individual_intron_positions_nc_plot.pdf",
       nc_per_intron_position_plot, device = "pdf", width = 5, height = 2.2, units = "in")




