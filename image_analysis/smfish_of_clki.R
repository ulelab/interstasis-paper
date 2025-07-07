rm(list = ls())

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(broom)

setwd("/Users/rupertfaraway/")

# Load --------------------------------------------------------------------

nuclei <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/Nuclei.csv")
nucleoplasm <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/Nucleoplasm.csv")
exp_nuc <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/ExpSpeckNucleoplasm.csv")
speckle <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/TotalSpeckle.csv")
exp_speckle <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/ExpTotalSpeckle.csv")
images <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/Image.csv")

# Wrangle -----------------------------------------------------------------

images %>%
  select(ImageNumber, FileName_Image) %>%
  separate(FileName_Image, into = c("target", "condition", "replicate", "fov", "garbage"), convert = T) %>%
  select(-c(garbage)) ->
  image_key

nuclei %>%
  inner_join(image_key, by = "ImageNumber") ->
  annotated_nuclei

ggplot(annotated_nuclei) +
  aes(x = AreaShape_Area, colour = factor(replicate)) +
  facet_wrap(. ~ condition) +
  geom_density() +
  geom_vline(xintercept = 5000)

annotated_nuclei %>%
  filter(AreaShape_Area > 5000) %>%
  filter(Location_Center_X > 200, Location_Center_X < 2104,
         Location_Center_Y > 200, Location_Center_Y < 2104,) ->
  filtered_nuclei
  
filtered_nuclei %>%
  select(ImageNumber, ObjectNumber,
         target, condition, replicate, fov,
         nucleus_size = AreaShape_Area) ->
  nuclei_type

# Ratios in and out of speckle --------------------------------------------

list("nucleoplasm" = nucleoplasm, "speckle" = speckle) %>%
  bind_rows(.id = "localisation") %>%
  inner_join(nuclei_type, by = c("ImageNumber", "ObjectNumber")) %>%
  mutate(bg_norm_smfish = Intensity_MeanIntensity_smFISH - Intensity_MinIntensity_smFISH,
         bg_norm_dapi = Intensity_MeanIntensity_DAPI - Intensity_MinIntensity_DAPI,
         bg_norm_son = Intensity_MeanIntensity_Son - Intensity_MinIntensity_Son) ->
  stacked

stacked %>%
  select(target, condition, replicate, fov, localisation, ImageNumber, ObjectNumber, intensity = bg_norm_smfish, nucleus_size) %>%
  # filter(fov != 1, fov != 10) %>%
  pivot_wider(names_from = localisation, values_from = intensity) %>%
  mutate(big_nuc = nucleus_size > median(nucleus_size)) %>%
  mutate(ratio = speckle/nucleoplasm) %>%
  mutate(target = target %>%
           fct_relevel("PSAP", "HNRNPDL", "BRD4", "EIF3A"),
         condition = condition %>%
           fct_relevel("DMSO")) ->
  mean_fish_ratio

mean_fish_ratio %>%
  group_by(target, condition, replicate) %>%
  summarise(mean_ratio = mean(ratio)) ->
  fish_ratio_group_means

fish_ratio_group_means %>%
  group_by(target) %>%
  summarise(t.test(mean_ratio ~ condition) %>% tidy()) ->
  mrna_speckle_ttests

ggplot(mean_fish_ratio, 
                   aes(x = condition,
                       y = ratio,
                   )) +
  geom_violin(aes(colour = factor(replicate), fill = factor(replicate)), position = "identity", alpha = 0.2, show.legend = FALSE) +
  # geom_beeswarm(cex = 0.45, alpha = 1, size = .15, aes(colour = factor(replicate),), show.legend = FALSE) +
  geom_point(data = fish_ratio_group_means, aes(y = mean_ratio, fill = factor(replicate)), size = 3, shape = 21, show.legend = FALSE) +
  theme_classic() +
  facet_wrap(. ~ target, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(0.5, 3.5)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3)) +
  labs(fill = "replicate",
       x = "",
       y = "mRNA signal in speckle / nucleoplasm") +
  scale_colour_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  scale_fill_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black")) ->
  mrna_speckle_violins
  
ggsave("WorkImages/20230518_CLKi_smFISH/plots/mrna_speckle_violins.pdf", mrna_speckle_violins,
       width = 6, height = 3)

write_tsv(mrna_speckle_ttests, "WorkImages/20230518_CLKi_smFISH/plots/mrna_speckle_ttests.tsv")


# Speckles per nucleus ----------------------------------------------------

spn <- read_csv("WorkImages/20230518_CLKi_smFISH/cellprofile_output/MaskedSpeckles.csv") %>%
  left_join(image_key) %>%
  drop_na()

names(image_key)

spn %>%
  ggplot() +
  aes(x = Intensity_MeanIntensity_DAPI, fill = condition) +
  geom_density(alpha = 0.5)

spn %>%
  inner_join(filtered_nuclei %>% 
               select(ImageNumber, target, condition,
                      replicate, fov, parent_size = AreaShape_Area,
                      parent_meandapi = Intensity_MeanIntensity_DAPI,
                      parent_meanfish = Intensity_MeanIntensity_smFISH,
                      Parent_Nuclei = ObjectNumber),
             by = c("ImageNumber", "target", "condition",
                    "replicate", "fov", "Parent_Nuclei")) %>%
  mutate(condition = condition %>% fct_relevel("DMSO")) %>%
  group_by(target, condition, replicate, fov, Parent_Nuclei, parent_size) %>%
  summarise(sp_son_mean_int = mean(Intensity_MeanIntensity_Son),
            sp_dapi_mean_int = mean(Intensity_MeanIntensity_DAPI),
            sp_dapi_mean_normalised = mean(Intensity_MeanIntensity_DAPI/parent_meandapi),
            sp_fish_mean_normalised = mean(Intensity_MeanIntensity_smFISH/parent_meanfish),
            total_speckle_size = sum(AreaShape_Area),
            average_speckle_size = mean(AreaShape_Area),
            proportion_nucleus_speckle = sum(AreaShape_Area/parent_size),
            eccentricity = mean(AreaShape_Eccentricity),
            stdev_son = mean(Intensity_StdIntensity_Son),
            stdev_smfish = mean(Intensity_StdIntensity_smFISH),
            number_of_speckles = dplyr::n(),) %>%
  mutate(number_of_speckles_per_pixel = number_of_speckles/parent_size) %>%
  ungroup() ->
  spn_summaries

spn_summaries %>%
  group_by(target, condition, replicate, fov) %>%
  summarise(number_of_nuclei = dplyr::n(),
            mean_ratio = mean(average_speckle_size)) %>%
  ungroup() %>%
  ggplot(aes(x = number_of_nuclei, y = mean_ratio, color = condition)) +
  geom_point() +
  # facet_wrap(. ~ target) +
  scale_y_log10()



spn_summaries %>%
  group_by(target, condition, replicate,) %>%
  summarise(sp_son_mean_int = mean(sp_son_mean_int),
            sp_dapi_mean_int = mean(sp_dapi_mean_int),
            sp_dapi_mean_normalised = mean(sp_dapi_mean_normalised),
            sp_fish_mean_normalised = mean(sp_fish_mean_normalised),
            total_speckle_size = mean(total_speckle_size),
            average_speckle_size = mean(average_speckle_size),
            proportion_nucleus_speckle = mean(proportion_nucleus_speckle),
            number_of_speckles = mean(number_of_speckles),
            eccentricity = mean(eccentricity),
            stdev_son = mean(stdev_son),
            stdev_smfish = mean(stdev_smfish),
            number_of_speckles_per_pixel = mean(number_of_speckles_per_pixel)) %>%
  ungroup() %>%
  sample_frac(1L) ->
  spn_per_replicate_means
  

ggplot(spn_summaries) +
  aes(x = number_of_speckles_per_pixel, fill = condition) +
  geom_density(alpha = 0.5)


ggplot(spn_per_replicate_means) +
  aes(x = condition, y = number_of_speckles_per_pixel) +
  geom_bar(stat = "summary", fun = "mean", colour = "black", fill = "white", width = 0.6) +
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar",
               width = 0.3) + 
  # geom_point(position = position_dodge2(width = .5), shape = 21, size = 3) +
  # coord_cartesian(ylim = c(0, 0.0035)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10)) +
  labs(x = "", y = "number of speckles\nper nuclear pixel") ->
  speckles_per_nuclear_area_plot

ggplot(spn_per_replicate_means) +
  aes(x = condition, y = average_speckle_size) +
  geom_bar(stat = "summary", fun = "mean", colour = "black", fill = "white", width = 0.6) +
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar",
               width = 0.3) + 
  # geom_point(position = position_dodge2(width = .5), shape = 21, size = 3) +
  coord_cartesian(ylim = c(0, 60)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10)) +
  labs(x = "", y = "average speckle size (pixels)") ->
  speckle_size_plot

ggplot(spn_per_replicate_means) +
  aes(x = condition, y = eccentricity) +
  geom_bar(stat = "summary", fun = "mean", colour = "black", fill = "white", width = 0.6) +
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar",
               width = 0.3) + 
  # geom_point(position = position_dodge2(width = .5), shape = 21, size = 3) +
  coord_cartesian(ylim = c(0.5, 0.7)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10)) +
  labs(x = "", y = "eccentricity of speckles") ->
  speckle_eccentricity_plot

t.test(number_of_speckles_per_pixel ~ condition, data = spn_per_replicate_means) %>% tidy() ->
  speckle_number_ttest

t.test(eccentricity ~ condition, data = spn_per_replicate_means) %>% tidy() ->
  eccentricity_ttest

t.test(average_speckle_size ~ condition, data = spn_per_replicate_means) %>% tidy() ->
  speckle_size_ttest

ggsave("WorkImages/20230518_CLKi_smFISH/plots/speckle_eccentricity_plot.pdf", speckle_eccentricity_plot,
       width = 3, height = 3)

ggsave("WorkImages/20230518_CLKi_smFISH/plots/speckle_size_plot.pdf", speckle_size_plot,
       width = 3, height = 3)

ggsave("WorkImages/20230518_CLKi_smFISH/plots/speckles_per_nuclear_area_plot.pdf", speckles_per_nuclear_area_plot,
       width = 3, height = 3)

write_tsv(speckle_number_ttest, "WorkImages/20230518_CLKi_smFISH/plots/speckle_eccentricity_plot_ttest.tsv")
write_tsv(eccentricity_ttest, "WorkImages/20230518_CLKi_smFISH/plots/eccentricity_ttest.tsv")
write_tsv(speckle_size_ttest, "WorkImages/20230518_CLKi_smFISH/plots/speckle_size_ttest.tsv")


