rm(list = ls())

library(tidyverse)
# library(ggbeeswarm)
library(patchwork)
library(broom)


# Load --------------------------------------------------------------------

nuclei <- read_csv("GASR/image_analysis/20230612_LUC7L3_TRA2B/cellprofiler_output/NucleiMasks.csv")
nucleoplasm <- read_csv("GASR/image_analysis/20230612_LUC7L3_TRA2B/cellprofiler_output/Nucleoplasm.csv")
speckle <- read_csv("GASR/image_analysis/20230612_LUC7L3_TRA2B/cellprofiler_output/TotalSpeckles.csv")
images <- read_csv("GASR/image_analysis/20230612_LUC7L3_TRA2B/cellprofiler_output/Image.csv")

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Images) %>%
  separate(FileName_Images, into = c("stopcodon", "bias", "replicate", "fov", "garbage"), convert = T) %>%
  mutate(construct = paste0(stopcodon, "_", bias)) %>%
  # Unfortunately, the UBC probes bind all over the genome in strange conserved intronic/intergenic regions.
  # It's actually quite cool, but also very strange and confusing. UBC behaviour in RNAseq is completely different to
  # the behaviour I observe in smFISH, so it seems likely that something weird is happening here with the probes.
  select(-c(garbage)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") ->
  annotated_nuclei

ggplot(annotated_nuclei) +
  aes(x = AreaShape_Area, colour = factor(replicate)) +
  facet_wrap(stopcodon ~ bias) +
  geom_density() +
  scale_x_log10() +
  geom_vline(xintercept = 5000)

annotated_nuclei %>%
  filter(AreaShape_Area > 800,
         AreaShape_Area < 4000,
         Intensity_MeanIntensity_TRA2B > 0.003,
         Intensity_MeanIntensity_SC35 > 0.003) %>%
  filter(AreaShape_BoundingBoxMinimum_X > 200, AreaShape_BoundingBoxMaximum_X < 2104,
         AreaShape_BoundingBoxMinimum_Y > 200, AreaShape_BoundingBoxMaximum_Y < 2104,) ->
  filtered_nuclei

filtered_nuclei %>%
  select(ImageNumber, ObjectNumber,
         parent_mGL_i = Intensity_IntegratedIntensity_mGL, parent_mGL_m = Intensity_MeanIntensity_mGL,
         construct, 
         bias, stopcodon, replicate, fov,
         nucleus_size = AreaShape_Area) %>%
  mutate(expressing = parent_mGL_i > 10) ->
  nuclei_type


# Ratios in and out of speckle --------------------------------------------

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) ->
  stacked

stacked %>%
  select(localisation, construct, replicate, fov, ObjectNumber, nucleus_size, Intensity_MeanIntensity_TRA2B, Intensity_IntegratedIntensity_TRA2B) %>%
  filter(localisation == "speckle") %>%
  mutate(speckle_size = Intensity_IntegratedIntensity_TRA2B/Intensity_MeanIntensity_TRA2B,
         prop_speckle = speckle_size / nucleus_size) %>%
  select(-c(localisation, nucleus_size, Intensity_IntegratedIntensity_TRA2B, Intensity_MeanIntensity_TRA2B)) %>%
  filter(prop_speckle > .05, prop_speckle < 0.33) %>%
  drop_na() ->
  speckle_area_filter

stacked %>%
  drop_na() %>%
  select(construct, stopcodon, bias, replicate, fov, localisation, ImageNumber,
         ObjectNumber, intensity = Intensity_MedianIntensity_TRA2B, nucleus_size,
         expressing,
         parent_mGL_i, parent_mGL_m) %>%
  # filter(fov != 1, fov != 10) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  inner_join(speckle_area_filter) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  mutate(construct =  construct %>%
           fct_relevel("NoStop_Poor", "NoStop_Rich", "Stop_Poor", "Stop_Rich"),
         bias = bias %>%
           fct_relevel("Poor"),
         stopcodon = stopcodon %>%
           fct_relevel("NoStop")) ->
  mean_tra2b_ratio

mean_tra2b_ratio %>%
  filter(expressing == T) %>%
  group_by(construct, bias, stopcodon, replicate) %>%
  summarise(mean_ratio = mean(ratio)) ->
  tra2b_ratio_group_means

tra2b_ratio_group_means %>%
  group_by(bias) %>%
  summarise(t.test(mean_ratio ~ stopcodon) %>% tidy()) ->
  tra2b_ratio_ttests

ggplot(mean_tra2b_ratio %>% filter(expressing == T) , 
       aes(x = stopcodon,
           y = ratio,
       )) +
  geom_violin(aes(colour = factor(replicate), fill = factor(replicate)), position = "identity", alpha = 0.2, show.legend = FALSE) +
  # geom_beeswarm(cex = 0.45, alpha = 1, size = .15, aes(colour = factor(replicate),), show.legend = FALSE) +
  geom_point(data = tra2b_ratio_group_means, aes(y = mean_ratio, fill = factor(replicate)), size = 3, shape = 21, show.legend = FALSE) +
  theme_classic() +
  facet_wrap(. ~ bias, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(0.5, 3.5)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3)) +
  labs(fill = "replicate",
       x = "",
       y = "mRNA signal in speckle / nucleoplasm") +
  # scale_colour_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  # scale_fill_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) ->
  tra2b_speckle_violins


ggplot(mean_tra2b_ratio %>% sample_frac(1L) %>% filter(parent_mGL_i > 10) %>%
         mutate(stopcodon = case_when(stopcodon == "NoStop" ~ "LCD Translated",
                                      T ~ "LCD in 3'UTR"),
                bias = case_when(bias == "Poor" ~ "Low Purine",
                                 T ~ "High Purine"))) +
  aes(x = parent_mGL_i, y = ratio, fill = factor(replicate)) +
  facet_wrap(stopcodon ~ bias) +
  geom_point(shape = 21, size = 2, show.legend = F) +
  scale_x_log10() +
  scale_y_log10(breaks = c(1, 2, 3, 4)) +
  coord_cartesian(xlim = c(10, 300), ylim = c(0.9, 4)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  geom_smooth(method = "lm", aes(color = factor(replicate)), se = F, linetype = "dashed", size = 2, show.legend = F) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "replicate",
       x = "mGreenLantern signal (AU)",
       y = "TRA2B signal\n(speckle / nucleoplasm)") ->
  tra2b_mgl_dose_plot

ggsave("GASR/image_analysis/20230612_LUC7L3_TRA2B/plots/tra2b_mgl_dose_plot.pdf", tra2b_mgl_dose_plot,
       width = 5, height = 5)

mean_tra2b_ratio %>% sample_frac(1L) %>% filter(parent_mGL_i > 10) %>%
  mutate(stopcodon = case_when(stopcodon == "NoStop" ~ "LCD Translated",
                               T ~ "LCD in 3'UTR"),
         bias = case_when(bias == "Poor" ~ "Low Purine",
                          T ~ "High Purine")) %>%
  mutate(construct = paste0(stopcodon, ", ", bias)) %>%
  group_by(construct, replicate) %>%
  summarise(model = lm(scale(log2(ratio)) ~ scale(log2(parent_mGL_i))) %>% tidy()) %>%
  unnest(model) %>%
  filter(term != "(Intercept)") ->
  models

write_tsv(models, "GASR/image_analysis/20230612_LUC7L3_TRA2B/plots/tra2b_mgl_model_summaries.tsv")

ggplot(models %>% mutate(construct = construct %>% str_replace(", ", "\n")),
       aes(x = construct, y = estimate, fill = factor(replicate))) +
  geom_point(shape = 21, size = 3, show.legend = F, position = position_dodge2(width = 0.6)) + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "", y = "slope", fill = "replicate") ->
  tra2b_mgl_slope_plot

ggsave("GASR/image_analysis/20230612_LUC7L3_TRA2B/plots/tra2b_mgl_slope_plot.pdf", tra2b_mgl_slope_plot,
       width = 4.5, height = 3.25)

models %>%
  ungroup() %>%
  reframe(ttest = pairwise.t.test(estimate, construct) %>% tidy()) %>%
  unnest(ttest) ->
  slope_comparisons

write_tsv(slope_comparisons, "GASR/image_analysis/20230612_LUC7L3_TRA2B/plots/tra2b_mgl_slopecomparisons.tsv")


