rm(list = ls())

library(tidyverse)
# library(ggbeeswarm)
library(patchwork)
library(broom)

setwd("/Users/rupertfaraway/")

# Load --------------------------------------------------------------------

nuclei <- read_csv("WorkImages/20230201_TRA2B_PPIG_SC35/cellprofiler_output/Speckle_vs_Nucleoplasm/NucleiMasks.csv")
nucleoplasm <- read_csv("WorkImages/20230201_TRA2B_PPIG_SC35/cellprofiler_output/Speckle_vs_Nucleoplasm/Nucleoplasm.csv")
speckle <- read_csv("WorkImages/20230201_TRA2B_PPIG_SC35/cellprofiler_output/Speckle_vs_Nucleoplasm/TotalSpeckles.csv")
images <-  read_csv("WorkImages/20230201_TRA2B_PPIG_SC35/cellprofiler_output/Speckle_vs_Nucleoplasm/Image.csv")

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Images) %>%
  separate(FileName_Images, into = c("target", "condition", "replicate", "fov", "garbage"), convert = T) %>%
  mutate(condition = str_replace(condition, pattern = "control", replacement = "DMSO") %>% fct_relevel("DMSO")) %>%
  select(-c(garbage)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") ->
  annotated_nuclei

ggplot(annotated_nuclei) +
  aes(x = AreaShape_Area, colour = factor(replicate)) +
  facet_wrap(. ~ condition) +
  geom_density() +
  geom_vline(xintercept = 1000)

ggplot(annotated_nuclei) +
  aes(x = Intensity_MeanIntensity_PPIG, colour = factor(replicate)) +
  facet_wrap(. ~ condition) +
  geom_density() +
  scale_x_log10() +
  geom_vline(xintercept = .01)

annotated_nuclei %>%
  filter(AreaShape_Area > 1000,) %>%
  filter(AreaShape_BoundingBoxMinimum_X > 100, AreaShape_BoundingBoxMaximum_X < 2204,
         AreaShape_BoundingBoxMinimum_Y > 100, AreaShape_BoundingBoxMaximum_Y < 2204,) ->
  filtered_nuclei

filtered_nuclei %>%
  mutate(nonexpressing = (condition == "dox") & (Intensity_MeanIntensity_PPIG < 0.01)) %>%
  select(ImageNumber, ObjectNumber,
         target, condition, replicate, fov, nonexpressing,
         parent_er = Intensity_MeanIntensity_PPIG,
         nucleus_size = AreaShape_Area) ->
  nuclei_type

# Ratios in and out of speckle --------------------------------------------

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) ->
  stacked

stacked %>%
  filter(!nonexpressing,
         Intensity_MedianIntensity_TRA2B > 0) %>%
  select(target, condition, replicate, fov, localisation, ImageNumber, ObjectNumber, intensity = Intensity_MedianIntensity_TRA2B, nucleus_size) %>%
  # filter(fov != 1, fov != 10) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(big_nuc = nucleus_size > median(nucleus_size)) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  drop_na() %>%
  mutate(condition = condition %>%
           fct_relevel("DMSO")) ->
  mean_ratio

mean_ratio %>%
  group_by(target, condition, replicate) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  ungroup() ->
  ratio_group_means

ratio_group_means %>%
  summarise(t.test(mean_ratio ~ condition) %>% tidy()) ->
  ttest

ggplot(mean_ratio, 
       aes(x = condition,
           y = ratio,
       )) +
  geom_violin(aes(colour = factor(replicate), fill = factor(replicate)), position = "identity", alpha = 0.2, show.legend = FALSE) +
  # geom_beeswarm(cex = 0.45, alpha = 1, size = .15, aes(colour = factor(replicate),), show.legend = FALSE) +
  geom_point(data = fish_ratio_group_means, aes(y = mean_ratio, fill = factor(replicate)), size = 3, shape = 21, show.legend = FALSE) +
  theme_classic() +
  # facet_wrap(. ~ target) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(1, 5)) +
  scale_y_log10(breaks = c(1, 2, 3, 4)) +
  labs(fill = "replicate",
       x = "",
       y = "TRA2B signal\n(speckle / nucleoplasm)") +
  scale_colour_manual(values = c("#a1b0c5", "#d9cd8d")) +
  scale_fill_manual(values = c("#a1b0c5",  "#d9cd8d")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) +
  annotate("text", x = 1.5, y = 4, label = "*", size = 16) ->
  tra2b_dox_vs_dmso_plot

ggsave("WorkImages/20230201_TRA2B_PPIG_SC35/plots/tra2b_dox_vs_dmso_ratio_plot.pdf",
       tra2b_dox_vs_dmso_plot, width = 3, height = 3)


# Dose --------------------------------------------------------------------

stacked %>%
  filter(condition == "dox",
         # Filter weird shit
         parent_er > 0, Intensity_MedianIntensity_TRA2B > 0) %>%
  select(replicate, fov, localisation, ImageNumber, ObjectNumber, intensity = Intensity_MedianIntensity_TRA2B, parent_er, nucleus_size) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  drop_na() %>%
  mutate(big_nuc = nucleus_size > median(nucleus_size)) %>%
  mutate(ratio = speckle/nucleoplasm) ->
  dox_mean_ratio

ggplot(dox_mean_ratio %>% sample_frac(1L) %>% filter(parent_er > 5e-3)) +
  aes(x = parent_er, y = ratio, fill = replicate) +
  geom_point(shape = 21, size = 2, show.legend = F) +
  scale_x_log10() +
  scale_y_log10(breaks = c(1, 2, 3, 4)) +
  coord_cartesian(xlim = c(5e-3, 2e-01), ylim = c(1, 5)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  geom_smooth(method = "lm", aes(color = replicate), se = F, linetype = "dashed", size = 2, show.legend = F) +
  scale_colour_manual(values = c("#a1b0c5", "#d9cd8d")) +
  scale_fill_manual(values = c("#a1b0c5",  "#d9cd8d")) +
  labs(fill = "replicate",
        x = "PPIG LCD signal (AU)",
       y = "TRA2B signal\n(speckle / nucleoplasm)") ->
  tra2b_dox_dose_plot

ggsave("WorkImages/20230201_TRA2B_PPIG_SC35/plots/tra2b_dox_dose_plot.pdf",
       tra2b_dox_dose_plot, width = 3, height = 3)


