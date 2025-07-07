rm(list = ls())

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(broom)

# setwd("/Users/rupertfaraway/")
setwd("/nemo/lab/ulej/home/users/farawar")

# Load --------------------------------------------------------------------

nuclei <- read_csv("GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/Cellprofiler_Output/Nuclei.csv")
nucleoplasm <- read_csv("GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/Cellprofiler_Output/Nucleoplasm.csv")
speckle <- read_csv("GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/Cellprofiler_Output/TotalSpeckle.csv")
images <- read_csv("GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/Cellprofiler_Output/Image.csv")

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Image) %>%
  separate(FileName_Image, into = c("target", "garbage", "replicate", "fov", "garbage2"), convert = T) %>%
  # Unfortunately, the UBC probes bind all over the genome in strange conserved intronic/intergenic regions.
  # It's actually quite cool, but also very strange and confusing. UBC behaviour in RNAseq is completely different to
  # the behaviour I observe in smFISH, so it seems likely that something weird is happening here with the probes.
  filter(target != "UBC") %>%
  select(-c(garbage, garbage2)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") %>%
  group_by(replicate) %>%
  mutate(scaled_mean_er = Intensity_MeanIntensity_ER / mean(Intensity_MeanIntensity_ER)) %>%
  group_by(target, replicate) %>%
  mutate(scaled_mean_smfish = Intensity_MeanIntensity_smFISH / mean(Intensity_MeanIntensity_smFISH)) %>%
  ungroup() ->
  annotated_nuclei

# ggplot(annotated_nuclei) +
#   # aes(x = Intensity_IntegratedIntensity_ER) +
#   aes(x = scaled_mean_er, colour = factor(replicate)) +
#   facet_wrap(. ~ target) +
#   geom_density() +
#   geom_vline(xintercept = c(0.1, 1, 3, 10)) +
#   scale_x_log10()
# 
# ggplot(annotated_nuclei) +
#   aes(x = AreaShape_Area, colour = factor(replicate)) +
#   geom_density() +
#   geom_vline(xintercept = 4000)

annotated_nuclei %>%
  # Remove nuclei that are suspiciously small, remove nuclei that are on the edges of the image.
  filter(AreaShape_Area > 4000,
         AreaShape_BoundingBoxMinimum_X > 100,
         AreaShape_BoundingBoxMinimum_Y > 100,
         AreaShape_BoundingBoxMaximum_X < 2204,
         AreaShape_BoundingBoxMaximum_Y < 2204,) %>%
  mutate(er_type = case_when(scaled_mean_er < 0.1 ~ "none",
                             scaled_mean_er < 1 ~ "low",
                             scaled_mean_er < 3 ~ "medium",
                             T ~ "high") %>%
           fct_relevel("none", "low", "medium", "high")) %>%
  select(ImageNumber, ObjectNumber, er_type, target, replicate, fov, nuc_er = scaled_mean_er, scaled_mean_smfish) ->
  nuclei_type

# Ratios in and out of speckle --------------------------------------------

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) ->
  stacked

stacked %>%
  select(target, replicate, fov, er_type, localisation, ImageNumber, ObjectNumber, nuc_er, intensity = Intensity_MeanIntensity_smFISH) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  mutate(target = target %>%
           fct_relevel("PSAP", "HNRNPDL", "APP", "SON", "NCL", "BRD4", "EIF3A")) ->
  mean_fish_ratio

# ggplot(mean_fish_ratio) +
#   aes(x = er_type, y = ratio) +
#   facet_wrap(. ~ target) +
#   scale_y_log10() +
#   geom_boxplot() 
  # coord_cartesian(ylim = c(0.5, 5))

mean_fish_ratio %>%
  group_by(target, replicate, er_type) %>%
  summarise(mean_ratio = mean(ratio)) ->
  fish_ratio_group_means


# Plotting ratios for binned groups --------------------------------------------------

ggplot(fish_ratio_group_means) +
  aes(x = er_type,
      y = log2(mean_ratio),
      fill = factor(replicate)) +
  facet_wrap(. ~ target) +
  geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5), show.legend = F) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) +
  labs(x = "PPIG LCD expression",
       y = "log2 mRNA enrichment (speckle / nucleoplasm)") +
  scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) 


# Plotting example slopes -------------------------------------------------

psap_brd4_dotplot <- 
  # Shuffle the order of the points.
  mean_fish_ratio[sample(1:nrow(mean_fish_ratio), nrow(mean_fish_ratio), replace = F),] %>%
  filter(nuc_er > 0.1,
         target %in% c("PSAP", "BRD4")) %>%
  ggplot() +
  aes(x = log2(nuc_er), y = log2(ratio), group = factor(replicate)) +
  facet_wrap(. ~ target) +
  geom_point(stroke = 0.25, aes(fill = factor(replicate)), shape = 21, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  scale_color_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  labs(x = expression(log[2] ~ PPIG[LCD] ~ expression),
       y = expression(log[2] ~ smFISH ~ "in" ~ speckle ~ vs ~ nucleoplasm)) +
  coord_cartesian(ylim = c(-1, 3)) +
  stat_smooth(method = "glm", aes(colour = factor(replicate)),
              geom = "line",
              alpha = 0.75, 
              linewidth = 1, 
              linetype = "dashed", 
              show.legend = F)

ggsave(psap_brd4_dotplot, 
       filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/psap_brd4_dotplot.pdf",
       device = "pdf", units = "in",
       height = 3, width = 5)

hnrnpdl_eif3a_dotplot <- 
  # Shuffle the order of the points.
  mean_fish_ratio[sample(1:nrow(mean_fish_ratio), nrow(mean_fish_ratio), replace = F),] %>%
  filter(nuc_er > 0.1,
         target %in% c("HNRNPDL", "EIF3A")) %>%
  ggplot() +
  aes(x = log2(nuc_er), y = log2(ratio), group = factor(replicate)) +
  facet_wrap(. ~ target) +
  geom_point(stroke = 0.25, aes(fill = factor(replicate)), shape = 21, show.legend = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_fill_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  scale_color_manual(values = c("lightsteelblue3", "indianred3", "gold2")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  labs(x = expression(log[2] ~ PPIG[LCD] ~ expression),
       y = expression(log[2] ~ smFISH ~ "in" ~ speckle ~ vs ~ nucleoplasm)) +
  coord_cartesian(ylim = c(-1, 3)) +
  stat_smooth(method = "glm", aes(colour = factor(replicate)),
              geom = "line",
              alpha = 0.75, 
              linewidth = 1, 
              linetype = "dashed", 
              show.legend = F)

ggsave(hnrnpdl_eif3a_dotplot, 
       filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/hnrnpdl_eif3a_dotplot.pdf",
       device = "pdf", units = "in",
       height = 3, width = 5)


# Plotting all slopes -----------------------------------------------------

mean_fish_ratio %>%
  filter(nuc_er > 0.1, ratio > 0) %>% # Filtering out wildtype cells, as they make the distribution non-normal.
  mutate(scaled_log_ratio = scale(log2(ratio)),
         scaled_log_er = scale(log2(nuc_er))) %>%
  group_by(target, replicate) %>%
  do(lm(scaled_log_ratio ~ scaled_log_er, data = . ) %>% 
       tidy) %>%
  # filter(term == "scaled_log_er") %>%
  mutate(status = case_when(target %in% c("PSAP", "HNRNPDL") ~ "control",
                            T ~ "multivalent\npurines")) %>%
  group_by(term) %>% #A bit janky, but I think the p adjustment should just apply to the ER expression term.
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup() ->
  ratio_linear_models

write_tsv(ratio_linear_models, "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/tables/linear_models.tsv")

ggplot(ratio_linear_models %>% filter(term == "scaled_log_er")) +
  aes(x = target,
      y = estimate,
      fill = factor(status),
      group = factor(replicate)) +
  geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("gray20", "#E86353")) +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text = element_text(colour = "black")) +
  labs(x = "",
       y = "slope of speckle enrichment\ngiven PPIG expression",
       fill = "") ->
  all_slopes_plot

ggsave(all_slopes_plot, 
       filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/all_slopes_plot.pdf",
       device = "pdf", units = "in",
       height = 3, width = 5)

# Plot R squareds ---------------------------------------------------------

mean_fish_ratio %>%
  filter(nuc_er > 0.1, ratio > 0, nuc_er > 0) %>%
  mutate(scaled_log_ratio = scale(log2(ratio)),
         scaled_log_er = scale(log2(nuc_er))) %>%
  group_by(target, replicate) %>%
  do(lm(scaled_log_ratio ~ scaled_log_er, data = . ) %>% 
       glance) %>%
  mutate(status = case_when(target %in% c("PSAP", "HNRNPDL") ~ "control",
                            T ~ "multivalent\npurines")) ->
  ratio_linear_glance

ggplot(ratio_linear_glance) +
  aes(x = target,
      y = adj.r.squared,
      fill = factor(status),
      group = factor(replicate)) +
  geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("gray20", "#E86353")) +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
        axis.text = element_text(colour = "black")) +
  labs(x = "",
       y = "R2 of speckle enrichment\ngiven PPIG expression",
       fill = "") ->
  all_r2_plot

ggsave(all_r2_plot, 
       filename = "GASR/image_analysis/20230302_ER_SC35_smFISH/analysis/plots/all_r2_plot.pdf",
       device = "pdf", units = "in",
       height = 3, width = 5)
  