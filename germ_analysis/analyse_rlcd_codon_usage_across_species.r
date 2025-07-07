rm(list = ls())

library(tidyverse)
library(parallel)
library(corrr)

setwd("~/home/")

# Read data ---------------------------------------------------------------

read_tsv("GASR/germs/data/multiple_species/species_order.tsv") %>%
  arrange(desc(order)) ->
  species_df

species_df %>%
  mutate(species = species %>% fct_relevel(species_df$species),
         species_short = species_short %>% fct_relevel(species_df$species_short)) ->
  species_df

list.files("GASR/germs/data/multiple_species/species_summaries/", full.names = T) %>%
  lapply(read_tsv) %>%
  bind_rows() ->
  species_summaries

species_summaries %>%
  filter(total_genes > 5000,
         propotion_discarded < 0.05) ->
  good_transcriptomes

list.files("GASR/germs/data/multiple_species/species_codons_in_leds/", full.names = T) %>%
  mclapply(read_tsv, mc.cores = 8) %>%
  bind_rows() %>%
  relocate(species) %>% 
  inner_join(species_df) %>%
  filter(species %in% good_transcriptomes$species) ->
  codon_usage_in_leds

list.files("GASR/germs/data/multiple_species/species_codons_outside_leds/", full.names = T) %>%
  mclapply(read_tsv, mc.cores = 8) %>%
  bind_rows() %>%
  relocate(species) %>%
  inner_join(species_df) %>%
  filter(species %in% good_transcriptomes$species) ->
  codon_usage_outside_leds

list.files("~/home/GASR/germs/data/multiple_species/species_codons_per_transcript/", full.names = T) %>%
  mclapply(read_tsv, mc.cores = 8) %>%
  bind_rows() %>%
  relocate(species) %>%
  inner_join(species_df) %>%
  filter(species %in% good_transcriptomes$species) ->
  codon_usage_per_transcript

list.files("GASR/germs/data/multiple_species/species_aa_in_leds/", full.names = T) %>%
  mclapply(read_tsv, mc.cores = 8) %>%
  bind_rows() %>%
  relocate(species) %>%
  inner_join(species_df) %>%
  filter(species %in% good_transcriptomes$species) ->
  aa_usage_in_leds

r_rich_df <- aa_usage_in_leds %>%
  filter(aa == "R") %>% 
  filter(prop > 0.2,
         count >= 10)

# aa_usage_in_leds %>%
#   filter(species %in% good_transcriptomes$species) %>%
#   filter(aa == "R") %>%
#   ggplot(aes(x = count, color = species)) +
#   geom_density() +
#   scale_x_log10()

codons_in_r_rich <- codon_usage_in_leds %>% filter(region %in% r_rich_df$region)
codons_outside_r_rich <- codon_usage_outside_leds %>% filter(region %in% r_rich_df$region)

codons_in_r_rich %>%
  select(species, transcript_id) %>%
  distinct() %>%
  group_by(species) %>%
  summarise(num_tx = dplyr::n()) %>%
  arrange(num_tx)


# Loading entropy ---------------------------------------------------------


list.files("GASR/germs/data/multiple_species/species_entropy/", full.names = T) %>%
  mclapply(read_rds, mc.cores = 8) %>%
  bind_rows() %>%
  relocate(species) ->
  all_entropy

species_summaries %>% mutate(total_genes = total_genes - with_nonstandard_base) %>% dplyr::select(species, total_genes) ->
  gene_count

aa_usage_in_leds %>%
  filter(species %in% good_transcriptomes$species) %>%
  dplyr::select(-c(count, le_region)) %>%
  pivot_wider(names_from = "aa", values_from = "prop") %>%
  replace(is.na(.), 0) %>%
  group_by(species) %>%
  reframe(rs = R + S,
          total_charge = R + E + D + K,
          arg = R,
          net_charge = (R+K) - (E+D)) %>%
  left_join(species_df) ->
  aa_combo_df

ggplot(aa_combo_df, aes(x = species_short, y = total_charge)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

aa_combo_df %>% 
  group_by(species, species_short) %>%
  reframe(total_lcd = dplyr::n(),
            total_charged = sum(total_charge > 0.5)) %>%
  mutate(prop_charged_lcd = total_charged/total_lcd) ->
  prop_lcd_types

ggplot(prop_lcd_types, aes(x = species_short, y = prop_charged_lcd)) +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))

# Numbers of R-rich LCDs --------------------------------------------------

left_join(r_rich_df %>% 
            count(species, name = "r_rich"),
          all_entropy %>%
            filter(species %in% good_transcriptomes$species) %>%
            mutate(n_regions = map(le_region, ~{ length(.x) })) %>%
            filter(n_regions > 0) %>%
            count(species, name = "all_lcds")) %>%
  left_join(all_entropy %>% 
              filter(species %in% good_transcriptomes$species) %>%
              count(species, name = "species_total")) %>%
  mutate(percentage_lcd = all_lcds/species_total,
         percentage_r = r_rich/all_lcds,
         percentage_total_r = r_rich/species_total) %>%
  left_join(species_df) ->
  lcd_counting_df
  
lcd_counting_df %>%
  ggplot() +
  aes(x = species_short, y = percentage_total_r * 100) +
  geom_bar(stat = "identity", fill = "white", color = "black") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "", y = "percentage proteins\ncontaining R-LCD") +
  coord_flip() ->
  percentage_r_lcd_plot

lcd_counting_df %>%
  ggplot() +
  aes(x = species_short, y = percentage_r * 100) +
  geom_bar(stat = "identity", fill = "white", color = "black") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "", y = "percentage LCDs\ncontaining 20% R") +
  coord_flip() ->
  percentage_lcds_r_plot

ggsave("GASR/germs/plots/arginine/multispecies/percentage_lcds_r_plot.pdf",percentage_lcds_r_plot,
       width = 3, height = 3)

ggsave("GASR/germs/plots/arginine/multispecies/percentage_r_lcd_plot.pdf",percentage_r_lcd_plot,
       width = 3, height = 3)



lcd_counting_df %>%
  ggplot() +
  aes(x = species_short, y = percentage_lcd * 100) +
  geom_bar(stat = "identity", fill = "white", color = "black") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "", y = "percentage proteins\ncontaining LCD") +
  coord_flip() 

# Plots within LEDs -------------------------------------------------------------------

codons_in_r_rich %>%
  filter(aa == "R") %>%
  select(species, region, codon, codon_usage) %>%
  # # Scale usage?
  # group_by(species) %>%
  # mutate(codon_usage = scale(codon_usage)) %>%
  # ungroup() %>%
  pivot_wider(names_from = "codon", values_from = "codon_usage") %>%
  nest(data = -species) %>%
  mutate(cor_mat = map(data, ~{
    correlate(.x, method = "spearman") %>%
      # correlate(.x, method = "pearson") %>%
      pivot_longer(cols = -term, names_to = "term_y", values_to = "cor")
  })) %>%
  select(-data) %>%
  unnest(cor_mat) %>%
  mutate(cor = case_when(term == term_y ~ 0,
                         T ~ cor),
         term = fct_relevel(term, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"),
         term_y = fct_relevel(term_y, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC") %>% fct_rev()) %>%
  left_join(species_df %>% select(species, species_short)) ->
  cor_df

# Removing yeast because there are too few domains
ggplot(cor_df %>%
         filter(!(species_short %in% c("S. cerevisiae"))),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 5) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot

# Plotting a smaller number of representative species
ggplot(cor_df %>%
         filter(species_short %in% c("A. thaliana", "C. elegans", "D. rerio", "X. tropicalis", 
                                    "C. porosus", "G. gallus", "B. taurus", "H. sapiens")),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 2) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot_selectd

# ggsave("GASR/germs/plots/arginine/multispecies/correlation_tileplot_selectd.pdf", correlation_tileplot_selectd,
#        width = 4, height = 6)

cor_of_interest <- c("AGA:AGG", "AGA:CGA", "CGG:CGC", "AGA:CGG", "AGA:CGC")

cor_df %>%
  mutate(cor_name = paste0(term, ":", term_y)) %>%
  filter(cor_name %in% cor_of_interest) %>%
  mutate(cor_name = cor_name %>% fct_relevel(cor_of_interest)) ->
  filtered_cor_df

ggplot(filtered_cor_df %>% 
         filter(species != "Saccharomyces_cerevisiae"),
       aes(x = species_short, y = cor, fill = cor_name, color = cor_name)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = cor_name), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  scale_fill_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  scale_color_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  # scale_fill_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  # scale_color_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  labs(x = "", y = "correlation\n(Spearman's rank)", fill = "") ->
  codon_correlation_species_summary

  
# ggsave("GASR/germs/plots/arginine/multispecies/codon_correlation_species_summary.pdf",codon_correlation_species_summary,
#        width = 6, height = 4)

# Looking at codons across all transcripts ------------------------------------------

codon_usage_per_transcript %>%
  filter(aa == "R") %>%
  group_by(species, transcript_id) %>%
  summarise(n_r = sum(codon_count)) %>% 
  ungroup() ->
  number_r_per_transcript

sufficient_r <- number_r_per_transcript %>%
  filter(n_r >= 20, 
         !(transcript_id %in% r_rich_df$transcript_id)) %>%
  select(transcript_id) %>%
  unlist(use.names = F)

codon_usage_per_transcript %>%
  filter(transcript_id %in% sufficient_r, aa == "R") %>%
  select(species, transcript_id, codon, codon_usage) %>%
  # # Scale usage?
  # group_by(species) %>%
  # mutate(codon_usage = scale(codon_usage)) %>%
  # ungroup() %>%
  pivot_wider(names_from = "codon", values_from = "codon_usage") %>%
  nest(data = -species) %>%
  mutate(cor_mat = map(data, ~{
    correlate(.x, method = "spearman") %>%
      # correlate(.x, method = "pearson") %>%
      pivot_longer(cols = -term, names_to = "term_y", values_to = "cor")
  })) %>%
  select(-data) %>%
  unnest(cor_mat) %>%
  mutate(cor = case_when(term == term_y ~ 0,
                         T ~ cor),
         term = fct_relevel(term, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"),
         term_y = fct_relevel(term_y, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC") %>% fct_rev()) %>%
  left_join(species_df %>% select(species, species_short)) ->
  cor_transcript_df

# Removing yeast because there are too few domains
ggplot(cor_transcript_df %>%
         filter(!(species_short %in% c("S. cerevisiae"))),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 5) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot_transcript

# Plotting a smaller number of representative species
ggplot(cor_transcript_df %>%
         filter(species_short %in% c("A. thaliana", "C. elegans", "D. rerio", "X. tropicalis", 
                                     "C. porosus", "G. gallus", "B. taurus", "H. sapiens")),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 2) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot_selectd_transcript

# ggsave("GASR/germs/plots/arginine/multispecies/correlation_tileplot_selectd_outside.pdf", correlation_tileplot_selectd_outside,
#        width = 4, height = 6)


cor_of_interest <- c("AGA:AGG", "AGA:CGA", "CGG:CGC", "AGA:CGG", "AGA:CGC")

cor_transcript_df %>%
  mutate(cor_name = paste0(term, ":", term_y)) %>%
  filter(cor_name %in% cor_of_interest) %>%
  mutate(cor_name = cor_name %>% fct_relevel(cor_of_interest)) ->
  filtered_cor_transcript_df

ggplot(filtered_cor_transcript_df %>% 
         filter(species != "Saccharomyces_cerevisiae"),
       aes(x = species_short, y = cor, fill = cor_name, color = cor_name)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = cor_name), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  scale_color_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  labs(x = "", y = "correlation\n(Spearman's rank)", fill = "") ->
  codon_correlation_species_summary_transcript

library(patchwork)

codon_correlation_species_summary_transcript | codon_correlation_species_summary

# ggsave("GASR/germs/plots/arginine/multispecies/codon_correlation_species_summary.pdf",codon_correlation_species_summary,
#        width = 6, height = 4)

filtered_cor_transcript_df %>%
  dplyr::rename(cor_transcript = cor) %>%
  left_join(filtered_cor_df) %>%
  mutate(correlation_relative = cor - cor_transcript) ->
  relative_cor_df

ggplot(relative_cor_df %>% 
         filter(species != "Saccharomyces_cerevisiae"),
       aes(x = species_short, y = correlation_relative, fill = cor_name, color = cor_name)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = cor_name), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  scale_color_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  labs(x = "", y = "difference in correlation\nbetween R-LCDs and transcripts\n(Spearman's rank)", fill = "") 

# Looking at codons outside LEDs ------------------------------------------

codons_outside_r_rich %>%
  filter(aa == "R") %>%
  select(species, region, codon, codon_usage) %>%
  # # Scale usage?
  # group_by(species) %>%
  # mutate(codon_usage = scale(codon_usage)) %>%
  # ungroup() %>%
  pivot_wider(names_from = "codon", values_from = "codon_usage") %>%
  nest(data = -species) %>%
  mutate(cor_mat = map(data, ~{
    correlate(.x, method = "spearman") %>%
      # correlate(.x, method = "pearson") %>%
      pivot_longer(cols = -term, names_to = "term_y", values_to = "cor")
  })) %>%
  select(-data) %>%
  unnest(cor_mat) %>%
  mutate(cor = case_when(term == term_y ~ 0,
                         T ~ cor),
         term = fct_relevel(term, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC"),
         term_y = fct_relevel(term_y, "AGA", "AGG", "CGA", "CGT", "CGG", "CGC") %>% fct_rev()) %>%
  left_join(species_df %>% select(species, species_short)) ->
  cor_outside_df

# Removing yeast because there are too few domains
ggplot(cor_outside_df %>%
         filter(!(species_short %in% c("S. cerevisiae"))),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 5) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot_outside

# Plotting a smaller number of representative species
ggplot(cor_outside_df %>%
         filter(species_short %in% c("A. thaliana", "C. elegans", "D. rerio", "X. tropicalis", 
                                     "C. porosus", "G. gallus", "B. taurus", "H. sapiens")),
       aes(x = term, y = term_y, fill = cor)) +
  geom_tile(linewidth = 0.5, color = "black") +
  scale_fill_gradient2(
    # low = "blue",
    low = scales::muted("blue"),
    mid = "white",
    high = "red",
    # high = scales::muted("red"),
    midpoint = 0,
    guide = guide_colorbar(frame.colour = "black",
                           ticks.colour = "black", 
                           ticks.linewidth = 0.5, 
                           frame.linewidth = 0.5)) +
  facet_wrap(. ~ species_short, ncol = 2) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(fill = "correlation") ->
  correlation_tileplot_selectd_outside

# ggsave("GASR/germs/plots/arginine/multispecies/correlation_tileplot_selectd_outside.pdf", correlation_tileplot_selectd_outside,
#        width = 4, height = 6)


cor_of_interest <- c("AGA:AGG", "AGA:CGA", "CGG:CGC", "AGA:CGG", "AGA:CGC")

cor_outside_df %>%
  mutate(cor_name = paste0(term, ":", term_y)) %>%
  filter(cor_name %in% cor_of_interest) %>%
  mutate(cor_name = cor_name %>% fct_relevel(cor_of_interest)) ->
  filtered_cor_outside_df

ggplot(filtered_cor_outside_df %>% 
         filter(species != "Saccharomyces_cerevisiae"),
       aes(x = species_short, y = cor, fill = cor_name, color = cor_name)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = cor_name), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "top") +
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  scale_color_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  labs(x = "", y = "correlation\n(Spearman's rank)", fill = "") ->
  codon_correlation_species_summary_outside


codon_correlation_species_summary_outside | codon_correlation_species_summary

# 
# ggsave("GASR/germs/plots/arginine/multispecies/codon_correlation_species_summary.pdf",codon_correlation_species_summary,
#        width = 6, height = 4)


filtered_cor_outside_df %>%
  dplyr::rename(cor_out = cor) %>%
  left_join(filtered_cor_df) %>%
  mutate(correlation_relative = cor - cor_out) ->
  relative_insideout_cor_df

ggplot(relative_insideout_cor_df %>% 
         filter(species != "Saccharomyces_cerevisiae"),
       aes(x = species_short, y = correlation_relative, fill = cor_name, color = cor_name)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = cor_name), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  coord_cartesian(ylim = c(-0.75, 0.5)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  scale_color_manual(values = c("#B9529F", "#777ABA", "#EAAA2A", "#414141", "#7C7B7B")) +
  labs(x = "", y = "difference in correlation\nbetween R-LCDs and transcripts\n(Spearman's rank)", fill = "") 


# Codon usage stuff -------------------------------------------------------


ggplot(transcriptome_codon_usage %>% filter(aa == "R") %>% left_join(species_df),
       aes(x = species_short, y = proportion_usage_transcriptome, fill = codon, color = codon)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = codon), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_fill_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  # scale_color_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "", y = "codon usage", fill = "")

codons_in_r_rich %>% filter(aa == "R") %>% group_by(species, codon) %>%
  summarise(relative_codon_usage_in_lcd = mean(normalised_to_transcriptome)) %>%
  left_join(nice_names) ->
  average_lcd_r_codons

ggplot(average_lcd_r_codons %>% filter(species_short != "S. cerevisiae"),
       aes(x = species_short, y = relative_codon_usage_in_lcd, fill = codon, color = codon)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.3), 
             color = "black") +
  geom_line(show.legend = F,
            aes(group = codon), linewidth = 1,
            position = position_dodge(width = 0.3),
            alpha = 0.25) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  scale_y_log10() +
  # scale_fill_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  # scale_color_manual(values = c("#e3aa87", "#c26648", "#6dad8b", "#404040", "#7a7a7a")) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "", y = "relative codon usage\n(R-LCD / transcriptome)", fill = "")
