rm(list = ls())

# Collection of useful packages for dealing with data and plotting.
library(tidyverse)

# Useful for cleaning up the format of statistical tests
library(broom)

library(scico)

# Load --------------------------------------------------------------------

# Loading the names of the images and their 'ImageNumber'.
images <- read_csv("Documents/cellprofile_analysis/20230627_Neve_ER/NEW_ePB_output/Image.csv") %>%
  # Select only these columns:
  select(FileName_image, ImageNumber) %>%
  # Split the filename up into new columns (keep the original filename for now with remove = F).
  separate(FileName_image, into = c("epb", "polya", "time", "replicate", "tif"), remove = F) %>%
  # Get rid of the useless new columns.
  select(-c(epb, polya, tif)) %>%
  # Assign an order to the times
  mutate(time = time %>% fct_relevel("0hr", "4hr", "8hr"))

# Reading in the cell, nuclear and cytoplasmic measurements.
cells <- read_csv("Documents/cellprofile_analysis/20230627_Neve_ER/NEW_ePB_output/cell_objects.csv")
nuclei <- read_csv("Documents/cellprofile_analysis/20230627_Neve_ER/NEW_ePB_output/nuclear_objects.csv")
cyto <- read_csv("Documents/cellprofile_analysis/20230627_Neve_ER/NEW_ePB_output/cyto_objects.csv")


# Quality control ---------------------------------------------------------

# We need to remove objects that don't really make any sense:
# Objects should have a single nucleus and a single cytoplasm.
# The cells shouldn't be right at the edge of the image.
# The sizes of these things should be normal.

# We also want to add in the information about the drug dosage.

# The image is 2304x2304. Let's take stuff that doesn't go within 100 pixels of the edge (a bit arbitrary).
cells %>%
  # Join together the cell table and the image info table.
  left_join(images) %>%
  filter(AreaShape_BoundingBoxMinimum_X > 100,
    AreaShape_BoundingBoxMaximum_X < 2204,
    AreaShape_BoundingBoxMinimum_Y > 100,
    AreaShape_BoundingBoxMaximum_Y < 2204) %>%
  # Let's also make sure that the number of nuclei and cytoplasms is 1 each. (== means equals, = is used for setting values of variables and options)
  filter(Children_cyto_objects_Count == 1,
         Children_nuclear_objects_Count == 1) ->
  cells_filtered

# We seem to lose a shitload of objects when we filter out things that are close to the edge.

# Here I am just plotting the mean intensity of the poly_A signal in the cells.
# I want to check that there is a reasonably smooth distribution. 
# I split up the cells based on their dose, and plot them in different colors.
# It seems like the intensity of poly-A mRNA in the cell does go down.
ggplot(cells_filtered, aes(x = Intensity_MeanIntensity_polyA, color = time)) +
  geom_density()

# As does the total poly-A signal.
ggplot(cells_filtered, aes(x = Intensity_IntegratedIntensity_polyA, color = time)) +
  geom_density()


# Building a big table of everything --------------------------------------

# Ultimately, we want to have a big table where we know the various nuclear, cytoplasmic and whole cell measurements for each cell.
# We will need to rename all of the columns to reflect whether they are cell/cyto/nuc measurements.
# The things that will be used to relate all the measurements are the identifying columns:
# ImageNumber, ObjectNumber, dose, replicate

# Looking at the column names:
names(cells_filtered)

# Let's also get rid of the useless non-measurement columns like 'Children_cyto_objects_Count' (served its purpose)
cells_filtered %>%
  select(ImageNumber, ObjectNumber, replicate, time,
         # You can rename columns while selecting them.
         cell_area = AreaShape_Area,
         # Using tidyselect commands to save me time in selecting measurement columns
         contains("Intensity")) %>%
  # Rename everything that has 'Intensity' using this command I found on the internet
  # and that I don't fully understand. Basically combines 'cell_' with all column names, 
  # provided that column name contains "Intensity".
  rename_with(~str_c("cell_", .), .cols = contains("Intensity")) ->
  cells_filtered_named

# Now let's do this renaming thing to the cytoplasmic and nuclear tables, so that we can combine them all.
# The info about dose and replicate number are already in the cell table, so we don't need to add it here.
# We also already filtered the cells at the edge, so when we join the tables together, the nuclei and cytoplasm
# that are too close to the edge will also get filtered.

names(nuclei)

nuclei %>%
  select(ImageNumber, ObjectNumber, nucleus_area = AreaShape_Area,
         contains("Intensity")) %>%
  rename_with(~str_c("nucleus_", .), .cols = contains("Intensity")) ->
  nuclei_named

cyto %>%
  select(ImageNumber, ObjectNumber, cyto_area = AreaShape_Area,
         contains("Intensity")) %>%
  rename_with(~str_c("cyto_", .), .cols = contains("Intensity")) ->
  cyto_named

# Now to join the tables. We can use an 'inner join', where we match the rows of two tables
# based on a set of columns they share. In this case, we are matching based on the image number
# and the object number. If one of the tables has an image/object pair that isn't in the other table,
# then this gets removed. Because we filtered a lot of cells earlier, this inner join means we will
# only take nuclei/cytoplasms that match up to the remaining cells.
# We will do this in two steps - first join the nuclei, then the cytoplasms.

cells_filtered_named %>%
  inner_join(nuclei_named, by = c("ImageNumber", "ObjectNumber")) %>%
  inner_join(cyto_named, by = c("ImageNumber", "ObjectNumber")) ->
  big_table



# One last filtering ------------------------------------------------------

# One last thing to check - how big are the nuclei and cytoplasms? 
# We might have some weird tiny nuclei coming from mitotic cells, or
# cells where the segmentation failed for some reason.

names(big_table)

ggplot(big_table, aes(x = nucleus_area)) +
  geom_density() +
  geom_vline(xintercept = mean(big_table$nucleus_area))

ggplot(big_table, aes(x = cyto_area)) +
  geom_density() +
  geom_vline(xintercept = mean(big_table$nucleus_area))

# This looks pretty good, although there are a small number of very tiny nuclei.
# Let's filter out any cells where the nucleus is super tiny (we decide arbitrarily)
# Why not an area of 10,000?

# Let's also look at the cytoplasms and cell areas.

ggplot(big_table, aes(x = cyto_area)) +
  geom_density()

ggplot(big_table, aes(x = cell_area)) +
  geom_density() +
  # I'm adding a vertical line to the plot at around the average size of a nucleus, for comparison.
  geom_vline(xintercept = mean(big_table$nucleus_area))

# This all seems quite reasonable. Let's just filter out the tiny nuclei, even though there aren't many.

big_table %>%
  filter(nucleus_area > 10000,
         cyto_area > 10000,) ->
  big_table_filtered

# Actual analysis ---------------------------------------------------------

# Two main questions: 
# what is the nuclear/cytoplasmic mRNA ratio in response to dose?
# does mRNA get more 'speckly' in response to dose?

# For the first question, lets calculate some ratios. Let's start with the total intensity ratio.
big_table_filtered$integrated_nc_ratio <- 
  big_table_filtered$nucleus_Intensity_IntegratedIntensity_polyA / big_table_filtered$cell_Intensity_IntegratedIntensity_polyA

big_table_filtered$mean_nc_ratio <- 
  big_table_filtered$nucleus_Intensity_MeanIntensity_polyA / big_table_filtered$cell_Intensity_MeanIntensity_polyA


big_table_filtered$nuclear_intensity_ratio <- 
  big_table_filtered$nucleus_Intensity_UpperQuartileIntensity_polyA / big_table_filtered$nucleus_Intensity_MedianIntensity_polyA


big_table_filtered %>%
  ggplot() +
  aes(x = nucleus_Intensity_IntegratedIntensity_ER,
      y = integrated_nc_ratio,
      fill = time) +
  geom_point(shape = 21, size = 3) +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "lm", se = F, aes(color = time),
              linetype = 'dashed', linewidth = 1, show.legend = F) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  geom_hline(yintercept = big_table_filtered %>% 
               filter(time == "0hr") %>%
               summarise(median(integrated_nc_ratio)) %>% 
               unlist(use.names = F),
             linetype = "dashed",
             alpha = 0.5) +
  scale_fill_scico_d(direction = -1, begin = 0.15, end = 0.8) +
  scale_color_scico_d(direction = -1, begin = 0.15, end = 0.8) +
  labs(x = "PPIG LCD signal (AU)",
       y = "poly-A FISH signal (nucleus / cytoplasm)",
       fill = "dox induction") ->
  fish_ratio_plot
fish_ratio_plot

big_table_filtered %>%
  ggplot() +
  aes(x = nucleus_Intensity_IntegratedIntensity_ER,
      y = nuclear_intensity_ratio,
      fill = time) +
  geom_point(shape = 21, size = 3) +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "lm", se = F, aes(color = time),
              linetype = 'dashed', linewidth = 1, show.legend = F) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  scale_fill_scico_d(direction = -1, begin = 0.15, end = 0.8) +
  scale_color_scico_d(direction = -1, begin = 0.15, end = 0.8) +
  labs(x = "PPIG LCD signal (AU)",
       y = "poly-A FISH nuclear intensity (upper quartile / median)",
       fill = "dox induction") ->
  fish_speckly_plot
  









# First with a density plot (I just like it).
ggplot(big_table_filtered,
       aes(x = integrated_nc_ratio, color = time)) +
  geom_density()

# Looks quite good. Let's plot it as a boxplot (easy to read, and compact).

ggplot(big_table_filtered,
       aes(x = time, y = integrated_nc_ratio)) +
  geom_boxplot() +
  # axis labels
  labs(x = "expression time",
       # \n is interpretted as a line break.
       y = "nuclear / cytoplasmic\noligo-dT signal") +
  # Messing with the appearance.
  theme_classic() +
  scale_y_log10() +
  theme(axis.text = element_text(color = "black")) ->
  nc_ratio_plot

nc_ratio_plot

# Let's do some pairwise significance tests.
big_table_filtered %>%
  summarise(pairwise = pairwise.t.test(integrated_nc_ratio, time) %>% tidy()) %>%
  unnest(cols = everything()) %>%
  # I'm only interested in comparisons to 0. I limit the comparisons before correcting for multiple testing,.
  # filter(group2 == "0uM") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  nc_stats

# As a quick way to address 'speckliness' lets look at how much brighter the upper quartile of polyA signal is than the median.
big_table_filtered$nuclear_intensity_ratio <- big_table_filtered$nucleus_Intensity_UpperQuartileIntensity_polyA / big_table_filtered$nucleus_Intensity_MedianIntensity_polyA

ggplot(big_table_filtered, aes(x = nuclear_intensity_ratio, color = dose)) +
  geom_density()

# That looks clean to me - let's use it and make another boxplot
ggplot(big_table_filtered,
       aes(x = dose, y = nuclear_intensity_ratio)) +
  geom_boxplot() +
  # axis labels
  labs(x = "CLK-IN-T3 dose",
       y = "upper quartile / median\noligo-dT signal") +
  # Messing with the appearance.
  theme_classic() +
  theme(axis.text = element_text(color = "black")) ->
  nucint_ratio_plot

# Let's do some pairwise significance tests.
big_table_filtered %>%
  summarise(pairwise = pairwise.t.test(nuclear_intensity_ratio, dose) %>% tidy()) %>%
  unnest(cols = everything()) %>%
  # I'm only interested in comparisons to 0. I limit the comparisons before correcting for multiple testing,.
  filter(group2 == "0uM") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  nucint_stats

ggsave("Documents/cellprofile_analysis/20230627_Neve_CLKi/plots/nucint_ratio_plot.pdf", nucint_ratio_plot,
       width = 3, height = 3)

ggsave("Documents/cellprofile_analysis/20230627_Neve_CLKi/plots/nc_ratio_plot.pdf", nc_ratio_plot,
       width = 3, height = 3)

write_tsv(nucint_stats, "Documents/cellprofile_analysis/20230627_Neve_CLKi/plots/nucint_ratio_stats.tsv")
write_tsv(nc_stats, "Documents/cellprofile_analysis/20230627_Neve_CLKi/plots/nc_ratio_stats.tsv")

