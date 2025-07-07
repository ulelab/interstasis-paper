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

segmented_filtered <- read_tsv("GASR/export_reporter/targeted_sequencing/clk1/data/first_experiment_barcodes_upis.txt.gz")
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
  # geom_vline(xintercept = 1600, linetype = "dashed") +
  geom_vline(xintercept = 1800, linetype = "dashed") +
  theme_classic() +
  scale_x_continuous(breaks = breaks_pretty(5),
                     labels = comma) +
  coord_cartesian(xlim = c(0, 10000)) +
  labs(x = "unique barcodes",
       y = "cumulative proportion")


# Annotation --------------------------------------------------------------

barcode_occurence[1:1800,] %>%
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

write_tsv(new_barcode_to_struct, "GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/clk_treatment_barcode_architectures.tsv.gz")
new_barcode_to_struct <- read_tsv("GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/clk_treatment_barcode_architectures.tsv.gz")

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

# for_pca <- count_usable_upis %>% 
#   ungroup() %>% 
#   arrange(desc(min_cpm)) %>% 
#   filter(min_counts >= 100) %>%
#   # mutate(sample = paste0(experiment, "_", sample)) %>%
#   select(sample, cpm, bc_number) %>% 
#   mutate(bc_number = paste0("index_", bc_number)) %>%
#   pivot_wider(names_from = bc_number, values_from = "cpm") %>%
#   select_if(~ !any(is.na(.))) %>%
#   column_to_rownames("sample")
# 
# pca_of_samples <- prcomp(for_pca, scale = T) 
# 
# prop_variance <- summary(pca_of_samples)$importance[2,]
# 
# pca_plot <- pca_of_samples$x %>%
#   as.data.frame() %>%
#   rownames_to_column("sample") %>%
#   separate(sample, into = c("experiment", "fraction", "replicate"), remove = F) %>%
#   ggplot(aes(x = -PC1, y = PC2, colour = fraction, shape = experiment)) +
#   geom_point(size = 3) +
#   scale_color_brewer(palette = "Dark2") +
#   labs(x = paste0("PC1: ", round(prop_variance[1], digits = 3) * 100, "% of variance"),
#        y = paste0("PC2: ", round(prop_variance[2], digits = 3) * 100, "% of variance"),
#        color = "",
#        shape = "") +
#   theme_classic() +
#   theme(axis.text.x = element_text(colour = "black"),
#         axis.text.y = element_text(colour = "black"),
#         axis.ticks = element_line(colour = "black"))

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
  filter(min_counts >= 50) %>%
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
           fct_relevel("dmso")) 

write_tsv(l_plus_ratio, "GASR/export_reporter/targeted_sequencing/minimap_reannotation/data/clk_treatment_gene_ratios.tsv.gz")


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
  scale_fill_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
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

write_tsv(pairwise_t_tests, file = "GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/clk_treatment_group_means_multivalency_ttests.tsv")

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
  scale_fill_manual(values = c("#B4B4B4", "#B9529F")) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  coord_flip()

combi_plot <- (ridgeline_nolegend | group_means_plot) + plot_layout(widths = c(4, 2))

ggsave("GASR/export_reporter/targeted_sequencing/minimap_reannotation/plots/clk_treatment_main_figure_combiplot.pdf", combi_plot,
       width = 5, height = 3.5)


