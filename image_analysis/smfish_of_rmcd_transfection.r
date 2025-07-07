rm(list = ls())

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(broom)

# setwd("/Users/rupertfaraway/")
setwd("/nemo/lab/ulej/home/users/farawar")

# Load --------------------------------------------------------------------

nuclei <- read_csv("GASR/image_analysis/20240804_Neve_LCD_Transfection/cellprofiler_output/Nuclei.csv")
nucleoplasm <- read_csv("GASR/image_analysis/20240804_Neve_LCD_Transfection/cellprofiler_output/Nucleoplasm.csv")
speckle <- read_csv("GASR/image_analysis/20240804_Neve_LCD_Transfection/cellprofiler_output/TotalSpeckle.csv")
images <- read_csv("GASR/image_analysis/20240804_Neve_LCD_Transfection/cellprofiler_output/Image.csv")

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Image) %>%
  separate(FileName_Image, into = c("lcd", "target", "replicate", "fov", "garbage"), convert = T) %>%
  mutate(target = toupper(target),
         lcd = toupper(lcd),
         lcd = case_when(lcd == "PRPF8" ~ "PRPF38B",
                         lcd == "PRPF38" ~ "PRPF38B",
                         T ~ lcd)) %>%
  filter(lcd != "FILLER",
         lcd != "CDK11") %>%
  select(-c(garbage)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") %>%
  group_by(replicate) %>%
  mutate(scaled_mean_er = Intensity_MeanIntensity_ER / mean(Intensity_MeanIntensity_ER)) %>%
  group_by(target, replicate) %>%
  mutate(scaled_mean_smfish = Intensity_MeanIntensity_smFISH / mean(Intensity_MeanIntensity_smFISH)) %>%
  ungroup() ->
  annotated_nuclei

ggplot(annotated_nuclei) +
  # aes(x = Intensity_IntegratedIntensity_ER) +
  aes(x = Intensity_MeanIntensity_ER, colour = factor(replicate)) +
  facet_wrap(. ~ lcd) +
  geom_density() +
  geom_vline(xintercept = c(0.01, 0.1, .3, 1)) +
  scale_x_log10()

ggplot(annotated_nuclei) +
  aes(x = AreaShape_Area, colour = factor(replicate)) +
  geom_density() +
  geom_vline(xintercept = 6000)

annotated_nuclei %>%
  # Remove nuclei that are suspiciously small, remove nuclei that are on the edges of the image.
  filter(AreaShape_Area > 6000,
         AreaShape_BoundingBoxMinimum_X > 100,
         AreaShape_BoundingBoxMinimum_Y > 100,
         AreaShape_BoundingBoxMaximum_X < 2204,
         AreaShape_BoundingBoxMaximum_Y < 2204,) %>%
  mutate(er_type = case_when(Intensity_MeanIntensity_ER < 0.01 ~ "none",
                             Intensity_MeanIntensity_ER < .1 ~ "low",
                             Intensity_MeanIntensity_ER < .3 ~ "medium",
                             T ~ "high") %>%
           fct_relevel("none", "low", "medium", "high")) %>%
  select(ImageNumber, ObjectNumber, lcd, target, replicate, fov, er_type, nuc_er = Intensity_MeanIntensity_ER, scaled_mean_smfish) ->
  nuclei_type


# Nucleolar loc -----------------------------------------------------------

# Extremely high expression can lead to the R-MCD relocalising to the nuclear speckle

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) ->
  stacked

stacked %>%
  select(lcd, target, replicate, fov, er_type, localisation, ImageNumber, ObjectNumber, nuc_er, intensity = Intensity_MeanIntensity_ER) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  mutate(target = target %>%
           fct_relevel("PSAP", "HNRNPDL", "BRD4", "EIF3A")) ->
  mean_er_ratio

mean_er_ratio[sample(1:nrow(mean_er_ratio), nrow(mean_er_ratio), replace = F),] %>%
  ggplot() +
  aes(x = log2(nuc_er), y = log2(ratio), fill = factor(replicate)) +
  facet_wrap(. ~ lcd) +
  geom_point(stroke = 0.25, aes(fill = factor(replicate)), shape = 21, show.legend = F) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) +
  scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  scale_color_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  theme_classic() +
  coord_cartesian(ylim = c(-2, 4)) +
  theme(axis.text = element_text(colour = "black")) +
  labs(x = expression(log[2] ~ LCD ~ expression),
       y = expression(log[2] ~ LCD ~ "in" ~ speckle ~ vs ~ nucleoplasm)) +
  geom_vline(xintercept = -1.5, linetype = "dashed") ->
  lcd_speckle_enrichment_plot

mean_er_ratio %>% 
  filter(log2(ratio) > -1,
         log2(nuc_er) < -1.5) %>%
  dplyr::select(lcd, target, replicate, fov, ImageNumber, ObjectNumber) ->
  good_lcd_localisation
  
# Ratios in and out of speckle --------------------------------------------

nuclei_type %>%
  inner_join(good_lcd_localisation) %>%
  group_by(lcd, replicate) %>%
  mutate(er_percentile = cut(nuc_er, breaks = quantile(nuc_er, probs = c(0, .5, .9, 1)), labels = c("low", "mid", "high"))) %>%
  dplyr::select(ImageNumber, ObjectNumber, lcd, target, replicate, fov, er_percentile) %>%
  ungroup() %>%
  drop_na() ->
  er_percentile_table

stacked %>%
  select(lcd, target, replicate, fov, er_type, localisation, ImageNumber, ObjectNumber, nuc_er, intensity = Intensity_MeanIntensity_smFISH) %>%
  inner_join(good_lcd_localisation) %>%
  inner_join(er_percentile_table) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  mutate(target = target %>%
           fct_relevel("PSAP", "HNRNPDL", "BRD4", "EIF3A")) ->
  mean_fish_ratio

mean_fish_ratio %>%
  group_by(lcd, target, replicate, er_percentile) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  ungroup() ->
  fish_ratio_group_means


# Plotting ratios for binned groups --------------------------------------------------

ggplot(fish_ratio_group_means) +
  aes(x = target,
      y = log2(mean_ratio),
      fill = er_percentile) +
  facet_wrap(. ~ lcd) +
  geom_point(size = 3, shape = 21, position = position_dodge2(width = 0.5), show.legend = T) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) +
  labs(x = "HCR-FISH target gene",
       y = "log2 mRNA enrichment\n(speckle / nucleoplasm)",
       fill = "") +
  scale_fill_manual(values = c("#b4b4b4", "#aa94b1", "#b9529f")) ->
  binned_groups_plot


ggsave("GASR/image_analysis/20240804_Neve_LCD_Transfection/plots/binned_groups_plot.pdf", binned_groups_plot,
      width = 9, height = 3)

fish_ratio_group_means %>%
  group_by(lcd, target) %>%
  reframe(model = pairwise.t.test(log2(mean_ratio), er_percentile) %>% tidy()) %>% 
  unnest(model) %>%
  group_by(lcd) %>%
  # filter(!((group1 == "12hr dox") & (group2 == "8hr dox"))) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  pairwise_t_tests

write_tsv(pairwise_t_tests, "GASR/image_analysis/20240804_Neve_LCD_Transfection/plots/pairwise_ttests.tsv")

# Plotting example slopes -------------------------------------------------

# psap_brd4_dotplot <-
  # Shuffle the order of the points.
  mean_fish_ratio[sample(1:nrow(mean_fish_ratio), nrow(mean_fish_ratio), replace = F),] %>%
  filter(nuc_er > 0.01) %>%
  ggplot() +
  aes(x = log2(nuc_er), y = log2(ratio), group = factor(replicate)) +
  facet_wrap(lcd ~ target, ncol = 4) +
  geom_point(stroke = 0.25, aes(fill = factor(replicate)), shape = 21, show.legend = F) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  scale_color_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  labs(x = expression(log[2] ~ LCD ~ expression),
       y = expression(log[2] ~ smFISH ~ "in" ~ speckle ~ vs ~ nucleoplasm)) +
  coord_cartesian(ylim = c(-1, 3)) +
  stat_smooth(method = "glm", aes(colour = factor(replicate)),
              geom = "line",
              alpha = 0.75,
              linewidth = 1,
              linetype = "dashed",
              show.legend = F)

# ggsave(psap_brd4_dotplot,
#        filename = "GASR/image_analysis/20240804_Neve_LCD_Transfection/plots/psap_brd4_dotplot.pdf",
#        device = "pdf", units = "in",
#        height = 3, width = 5)

# hnrnpdl_eif3a_dotplot <- 
#   # Shuffle the order of the points.
#   mean_fish_ratio[sample(1:nrow(mean_fish_ratio), nrow(mean_fish_ratio), replace = F),] %>%
#   filter(nuc_er > 0.1,
#          target %in% c("HNRNPDL", "EIF3A")) %>%
#   ggplot() +
#   aes(x = log2(nuc_er), y = log2(ratio), group = factor(replicate)) +
#   facet_wrap(. ~ target) +
#   geom_point(stroke = 0.25, aes(fill = factor(replicate)), shape = 21, show.legend = F) + 
#   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
#   scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
#   scale_color_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
#   theme_classic() +
#   theme(axis.text = element_text(colour = "black")) +
#   labs(x = expression(log[2] ~ PPIG[LCD] ~ expression),
#        y = expression(log[2] ~ smFISH ~ "in" ~ speckle ~ vs ~ nucleoplasm)) +
#   coord_cartesian(ylim = c(-1, 3)) +
#   stat_smooth(method = "glm", aes(colour = factor(replicate)),
#               geom = "line",
#               alpha = 0.75, 
#               linewidth = 1, 
#               linetype = "dashed", 
#               show.legend = F)
# 
# ggsave(hnrnpdl_eif3a_dotplot, 
#        filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/hnrnpdl_eif3a_dotplot.pdf",
#        device = "pdf", units = "in",
#        height = 3, width = 5)
# 
# 
# # Plotting all slopes -----------------------------------------------------
# 
# mean_fish_ratio %>%
#   filter(nuc_er > 0.1, ratio > 0) %>% # Filtering out wildtype cells, as they make the distribution non-normal.
#   mutate(scaled_log_ratio = scale(log2(ratio)),
#          scaled_log_er = scale(log2(nuc_er))) %>%
#   group_by(target, replicate) %>%
#   do(lm(scaled_log_ratio ~ scaled_log_er, data = . ) %>% 
#        tidy) %>%
#   # filter(term == "scaled_log_er") %>%
#   mutate(status = case_when(target %in% c("PSAP", "HNRNPDL") ~ "control",
#                             T ~ "multivalent\npurines")) %>%
#   group_by(term) %>% #A bit janky, but I think the p adjustment should just apply to the ER expression term.
#   mutate(padj = p.adjust(p.value, method = "BH")) %>%
#   ungroup() ->
#   ratio_linear_models
# 
# write_tsv(ratio_linear_models, "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/tables/linear_models.tsv")
# 
# ggplot(ratio_linear_models %>% filter(term == "scaled_log_er")) +
#   aes(x = target,
#       y = estimate,
#       fill = factor(status),
#       group = factor(replicate)) +
#   geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5)) +
#   theme_classic() +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
#   scale_fill_manual(values = c("gray20", "#E86353")) +
#   theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
#         axis.text = element_text(colour = "black")) +
#   labs(x = "",
#        y = "slope of speckle enrichment\ngiven PPIG expression",
#        fill = "") ->
#   all_slopes_plot
# 
# ggsave(all_slopes_plot, 
#        filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/all_slopes_plot.pdf",
#        device = "pdf", units = "in",
#        height = 3, width = 5)
# 
# # Plot R squareds ---------------------------------------------------------
# 
# mean_fish_ratio %>%
#   filter(nuc_er > 0.1, ratio > 0, nuc_er > 0) %>%
#   mutate(scaled_log_ratio = scale(log2(ratio)),
#          scaled_log_er = scale(log2(nuc_er))) %>%
#   group_by(target, replicate) %>%
#   do(lm(scaled_log_ratio ~ scaled_log_er, data = . ) %>% 
#        glance) %>%
#   mutate(status = case_when(target %in% c("PSAP", "HNRNPDL") ~ "control",
#                             T ~ "multivalent\npurines")) ->
#   ratio_linear_glance
# 
# ggplot(ratio_linear_glance) +
#   aes(x = target,
#       y = adj.r.squared,
#       fill = factor(status),
#       group = factor(replicate)) +
#   geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5)) +
#   theme_classic() +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
#   scale_fill_manual(values = c("gray20", "#E86353")) +
#   theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
#         axis.text = element_text(colour = "black")) +
#   labs(x = "",
#        y = "R2 of speckle enrichment\ngiven PPIG expression",
#        fill = "") ->
#   all_r2_plot
# 
# ggsave(all_r2_plot, 
#        filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/all_r2_plot.pdf",
#        device = "pdf", units = "in",
#        height = 3, width = 5)
  