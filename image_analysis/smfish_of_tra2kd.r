rm(list = ls())

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(broom)

# setwd("/Users/rupertfaraway/")
setwd("/nemo/lab/ulej/home/users/farawar")

# Load --------------------------------------------------------------------

nuclei <- read_csv("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/cellprofile_output/Nuclei.csv")
cytoplasm <- read_csv("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/cellprofile_output/FakeCyto.csv")
nucleoplasm <- read_csv("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/cellprofile_output/Nucleoplasm.csv")
speckle <- read_csv("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/cellprofile_output/TotalSpeckle.csv")
images <- read_csv("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/cellprofile_output/Image.csv")

# Blah --------------------------------------------------------------------

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Image) %>%
  separate(FileName_Image, into = c("condition", "target", "replicate", "fov", "garbage"), convert = T) %>%
  select(-c(garbage)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") %>%
  group_by(condition, replicate) %>%
  mutate(scaled_er = Intensity_MeanIntensity_ER / mean(Intensity_MeanIntensity_ER)) %>%
  ungroup() ->
  annotated_nuclei

ggplot(annotated_nuclei) +
  aes(x = AreaShape_Area, colour = factor(condition)) +
  geom_density() +
  geom_vline(xintercept = 8000)

ggplot(annotated_nuclei) +
  aes(x = Intensity_MeanIntensity_TRA2B, color = factor(condition), linetype = factor(replicate)) +
  # aes(x = scaled_mean_er, colour = factor(replicate)) +
  # facet_wrap(. ~ target) +
  geom_density() +
  geom_vline(xintercept = c(0.006, 0.01)) +
  scale_x_log10()

annotated_nuclei %>%
  # Remove nuclei that are suspiciously small, remove nuclei that are on the edges of the image.
  filter(AreaShape_Area > 8000,
         AreaShape_BoundingBoxMinimum_X > 50,
         AreaShape_BoundingBoxMinimum_Y > 50,
         AreaShape_BoundingBoxMaximum_X < 2254,
         AreaShape_BoundingBoxMaximum_Y < 2254,) %>%
  select(ImageNumber, ObjectNumber, condition, replicate, fov, nuc_er = scaled_er, nuc_tra2b = Intensity_MeanIntensity_TRA2B) ->
  full_nuclei

# Filter nuclei based on speckle segmentation -----------------------------

left_join(
  speckle %>%
    mutate(speckle_area = Intensity_IntegratedIntensity_TRA2B / Intensity_MeanIntensity_TRA2B) %>%
    select(ImageNumber, nucleus_number = Parent_Nuclei, speckle_area),
  nuclei %>%
    select(ImageNumber, nucleus_number = ObjectNumber, nucleus_area = AreaShape_Area,)) %>%
  left_join(image_key) %>%
  mutate(proportion_speckle = speckle_area / nucleus_area) ->
  proportion_of_nucleus_speckle

ggplot(proportion_of_nucleus_speckle, aes(x = proportion_speckle, color = condition)) +
  geom_density() +
  geom_vline(xintercept = c(0.04, 0.2))

proportion_of_nucleus_speckle %>%
  filter(proportion_speckle >= 0.04, proportion_speckle <= 0.2) %>%
  select(ImageNumber, condition, replicate, fov, ObjectNumber = nucleus_number) ->
  nuclei_with_normal_speckles

ggplot(annotated_nuclei %>% inner_join(nuclei_with_normal_speckles) %>% inner_join(nuclei_type)) +
  aes(x = scaled_er, color = factor(condition)) +
  # aes(x = scaled_mean_er, colour = factor(replicate)) +
  # facet_wrap(. ~ target) +
  geom_density() +
  geom_vline(xintercept = c(.75, 1.5)) +
  scale_x_log10()

full_nuclei %>%
  inner_join(nuclei_with_normal_speckles) %>%
  mutate(er_type = case_when(nuc_er <= quantile(nuc_er, 0.33) ~ "low",
                             nuc_er <= quantile(nuc_er, 0.66) ~ "medium",
                             T ~ "high") %>%
           fct_relevel("low", "medium", "high")) ->
  nuclei_type
  
# Ratios in and out of speckle --------------------------------------------

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) %>%
  inner_join(nuclei_with_normal_speckles) ->
  stacked

stacked %>%
  select(condition, replicate, fov, localisation, ImageNumber, ObjectNumber, er_type, nuc_er, nuc_tra2b, intensity = Intensity_MeanIntensity_smFISH) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(ratio = speckle/nucleoplasm) ->
  mean_fish_ratio

mean_fish_ratio %>%
  # filter( ( (tra2b_type == "low") & (condition == "TRA2KD") ) | ( (tra2b_type != "low") & (condition == "Control") ) ) %>%
  # filter(er_type == "high") %>%
  group_by(condition, replicate, er_type) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  ungroup() ->
  fish_ratio_group_means

ggplot(fish_ratio_group_means) +
  aes(x = er_type, y = mean_ratio, fill = condition) +
  geom_point(position = position_dodge2(width = 0.5), size = 3, shape = 21) +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "PPIG LCD expression",
       y = "log2 EIF3A mRNA enrichment\n(speckle / nucleoplasm)",
       fill = "") +
  scale_fill_manual(values = c("#B6549A", "#b4b4b4")) ->
  speckle_retention_plot

fish_ratio_group_means %>%
  group_by(er_type) %>%
  summarise(ttest = t.test(mean_ratio ~ condition) %>% tidy()) %>%
  unnest(ttest) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) ->
  speckle_retention_ttests

ggsave("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/speckle_retention_plot.pdf", speckle_retention_plot,
       height = 3, width = 3.5)

write_tsv(speckle_retention_ttests, "GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/speckle_retention_ttests.tsv")
write_tsv(fish_ratio_group_means, "GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/speckle_retention_data.tsv")

# TRA2 KD plot ------------------------------------------------------------

nuclei_type %>%
  group_by(condition, replicate) %>%
  summarise(mean_tra2b = mean(nuc_tra2b)) %>%
  ungroup() ->
  nuclei_tra2b_means

ggplot(nuclei_tra2b_means) +
  aes(x = condition, y = mean_tra2b, fill = condition) +
  geom_point(position = position_dodge2(width = 0.5), size = 3, shape = 21) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "nuclear TRA2B signal (A.U.)",
       fill = "") +
  scale_fill_manual(values = c("#B6549A", "#b4b4b4")) +
  coord_cartesian(ylim = c(0, 0.016)) ->
  tra2_kd_plot

nuclei_tra2b_means %>%
  summarise(ttest = t.test(mean_tra2b ~ condition) %>% tidy()) %>%
  unnest(ttest) ->
  tra2_kd_ttest

ggsave("GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/tra2_kd_plot.pdf", tra2_kd_plot,
       height = 3, width = 3)

write_tsv(tra2_kd_ttest, "GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/tra2_kd_ttest.tsv")
write_tsv(nuclei_tra2b_means, "GASR/image_analysis/20250502_Neve_TRA2_KD_FISH/plots/tra2_kd_data.tsv")

