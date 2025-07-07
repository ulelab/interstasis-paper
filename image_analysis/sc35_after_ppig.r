# This line just removes any old data that R has loaded, giving us a fresh workspace.
rm(list = ls())

# The tidyverse is a set of packages for R that are used for things like manipulating tables of data and making plots.
# Two of the most important packages in the tidyverse are dplyr (for tables) and ggplot2 (for plots).
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(broom)

# setwd("/Users/rupertfaraway/")

# Some basics -------------------------------------------------------------

# In R, you use <- to assign some value to a new variable. 
# Telling R x <- 6 means 'make a variable called x that has the value 6.
x <- 6
x
x - 2

# A variable could be more complicated: it could be a word, list of things, or a table, for example.
y <- "aaaaaaaaaaaaaaaaaa"
y

z <- c(x, y)
z

# A little note here - when we look at x, we see 6, but when we look at z, we now see "6" in quotes.
# R has converted the number/value 6 to the character "6" because it doesn't want us to have a vector with both numbers and string.
# So while we can get "6" out of z like so:
z[1]

# We can't subtract from it anymore. It's now a string, not a numeric value. This may seem very arcane, but it is an important distinction.
# z[1] - 2

# Loading in data ---------------------------------------------------------

# First, we read in the .csv files using the read_csv() function (from readr in the tidyverse) and assign the data to new variables.
# There are three files we'll need:
#   The information about each nuclear speckle - this is our main dataset
#   The information about each nucleus - we can use this to categorise cells by export reporter expression.
#   The information about each image - we need this to determine the replicate number and condition that things belong to.

speckles <- read_csv("GASR/image_analysis/20221130_NeveSRRM2ExpRep_2/cellprofiler_output/ShapeIntensity_MaskedSpeckles.csv")
nuclei <- read_csv("GASR/image_analysis/20221130_NeveSRRM2ExpRep_2/cellprofiler_output/ShapeIntensity_NucObjects.csv")
images <- read_csv("GASR/image_analysis/20221130_NeveSRRM2ExpRep_2/cellprofiler_output/ShapeIntensity_Image.csv")

# Sanitise the data about images ------------------------------------------------------------------

# "Wrangling" the data into the shape that you want to do your analysis or make your figures is 90%+ of the work in any project.

# The head command prints out the first few rows of a table. 
head(images)

# We can see that the table is too long to see all of the column names, so we can print the column names using names().
names(images)

# There is a lot of useless stuff - we really just want to know what the file name was for each image.
# We can see that there are columns called 'FileName_Images' and 'ImageNumber'.
# We just want to have a table that we can use to convert the image numbers into information about condition and replicate.
# First, we can use select() to take the columns that we want.

images_clean <- select(images, FileName_Images, ImageNumber)

images_clean

# Now we want to add columns that tell us whether the image is from is a control or dox sample and which replicate it is.
# For these images, we know that Neve took 10 images per well. 
# We also know the first 10 images (dox samples) had the wrong settings, so we'll want to remove them.

# I'll use the word() command to split up the file name and take one part.

# word() works like this:
word("my_dog_is_named_andrew", 5, sep = "_")
# Here we broke up the string based on where we could see "_", then we took the 5th element ("andrew").
# We can also give word() a vector of strings.
word(
  c("my_dog_is_named_andrew",
  "my_dog_is_named_sabrina",
  "my_dog_is_named_alfonso"),
  5,
  sep = "_")

# We can look at the FileName_Images column from images_clean table we made:
images_clean["FileName_Images"]
# or like this:
images_clean$FileName_Images
# we can see that if we split the names up using "_", we would find the dox/control information in the third column.

# We can add a new column to the dataframe by using the $ symbol.
images_clean$condition <- word(images_clean$FileName_Images, 3, sep = "_")

head(images_clean)

# There are many ways we could convert the ImageNumbers to the replicate numbers now.
# We want a vector that goes 1, 1, 1, 1, ... 2, 2, 2, .. with each number being repeated 10 times.
# We can first make a vector using the function c(), then we can repeat each element of the vector using rep().
# We know that the first set of 10 dox images need to be thrown out, so let's give them a totally different number.

# Make the vector.
replicates_1_value <- c(1, 2, 3, 200, 1, 2, 3)
# Repeat each element of this vector 10 times:
replicates_vector <- rep(replicates_1_value, each = 10)

replicates_vector

# And make this vector a new column in our data

images_clean$replicate <- replicates_vector

# Now we can get rid of any rows that have the replicate number 200.
# The way to do this using base R would be to get information about which rows have replicate equal
# to 200 using the operator "==", then take the other rows (taking the opposite using the operator !)
# to a new dataframe.

bad_rows <- images_clean$replicate == 200
bad_rows
!bad_rows

# We can take the rows using these square brackets []. We put the comma to show that these are rows, not columns.
# Note: we could simulaneously subset rows and columns like this: dataframe[1:10, 4] (take rows 1:10 of column 4).
images_clean_filtered <- images_clean[!bad_rows,]

# I've intentionally coded this in a very step-by-step way, where we make variables for each step of the process, then
# pass those variables into new functions to make new variables, etc. This is useful to understand what is 
# happening at each step, and if anything goes wrong at any step, you can trace it back.

# If I was actually writing this code, I would probably write something much quicker like this:

# Using the %>% operator (from magrittr in tidyverse) passes the images variable as the first 'argument' to
# the next function (in this case, select()).
images_table <- images %>%
  select(FileName_Images, ImageNumber) %>%
  # The mutate function is from dplyr in tidyverse, and it makes new columns in a dataframe.
  mutate(condition = word(FileName_Images, 3, sep = "_"), # Mutate knows that FileName_Images is a column in my dataframe.
         replicate = rep(c(1,2,3,200,1,2,3), each = 10)) %>%
  # The filter() command (also from dplyr) knows that replicate is a column in my dataframe.
  # Filter takes forward all of the data where the column 'replicate' is not equal to 200.
  filter(replicate != 200)

# We can see that this gives the same result:
identical(images_table, images_clean_filtered)


# Sanitise the data about nuclei ------------------------------------------

# One thing we can see in the images is that some cells in the dox condition don't have ER expression.
# I want to be able to take this into account in my analysis, so I want to make a table with information about
# each nucleus's export reporter expression.

# I also might want to filter out nuclei that are on the edges of the image.

names(nuclei)
  
nuclei_table <- nuclei %>%
  # The select command allows us to choose columns but also to give them different names.
  # I'm renaming the columns that have the same name in the 'speckles' table so that when we 
  # put everything together, we won't have two columns with the same name.
  # I've also named renamed 'ObjectNumber' to have the same name as it does in the speckle table.
  select(ImageNumber,
         "Nuc_Max_X" = AreaShape_BoundingBoxMaximum_X,
         "Nuc_Min_X" = AreaShape_BoundingBoxMinimum_X,
         "Nuc_Max_Y" = AreaShape_BoundingBoxMaximum_Y,
         "Nuc_Min_Y" = AreaShape_BoundingBoxMinimum_Y,
         "Parent_NucObjects" = ObjectNumber, 
         "Nuc_Intensity_MeanIntensity_ExpRep" = Intensity_MeanIntensity_ExpRep,
         "Nuc_Intensity_MeanIntensity_Srrm2" = Intensity_MeanIntensity_Srrm2,
         "Nuc_AreaShape_Area" = AreaShape_Area,
         "Nuc_Intensity_UpperQuartileIntensity_ExpRep" = Intensity_UpperQuartileIntensity_ExpRep,
         "Nuc_Intensity_LowerQuartileIntensity_ExpRep" = Intensity_LowerQuartileIntensity_ExpRep) %>%
  # The total amount of export reporter or SRRM2 should be the mean intensity multiplied by the total size of the nucleus.
  mutate(Nuc_Total_ExpRep = Nuc_Intensity_MeanIntensity_ExpRep * Nuc_AreaShape_Area,
         Nuc_Total_Srrm2 = Nuc_Intensity_MeanIntensity_Srrm2 * Nuc_AreaShape_Area)

# The dimensions of the images are 2304 by 2304. Let's remove nuclei that fall within 200 pixels of the edge.
edge_distance <- 200

nuclei_table_trimmed <- nuclei_table %>%
  filter(Nuc_Max_X < (2304 - edge_distance),
         Nuc_Min_X > edge_distance,
         Nuc_Max_Y < (2304 - edge_distance),
         Nuc_Min_Y > edge_distance,)

# Now, let's combine this information with the information we got about the image numbers earlier.
# We do this by performing a 'join' of the data. We choose a column that's shared between two images,
# which is the 'ImageNumber' column. The join will add the columns from one table to another table, 
# matching them up based on this joining column.
# There is a super useful cheatsheat for dplyr that shows you things like this:
# https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
# It's on page 2 on the right.

# We want to do an 'inner join' - this will discard any rows that don't have a match between the two tables.
# We're doing this because we want to get rid of all the nuclei that are from those first 10 dox images.

nuclei_annotated <- inner_join(nuclei_table_trimmed, images_table) # dplyr figures out that we have a shared column 'ImageNumber' to join by.

nuclei_annotated

# Let's use ggplot to make a quick density plot that shows the different levels of ER signal in the dox and control conditions.
# ggplot works by adding different layers of the plot together.
# First we start the ggplot by giving it our data.
ggplot(nuclei_annotated) +
  # Next we tell ggplot what the different columns to use are, and which parts of the plot they will make.
  # For this plot, we can show the ER signal on the x axis, and to group their fill colour by the condition (dox vs control).
  aes(x = Nuc_Total_ExpRep, fill = condition) +
  # Now we tell ggplot to make a density plot using this information.
  # I use the option 'alpha' to change the opacity of the fill colour, so you can see both of the densities.
  geom_density(alpha = 0.5)

# You can tweak the ggplot endlessly.
# I don't like the grey background, so I can change the theme to a simpler one.
# I can also set the axis limits with coord_cartesian(), and change the axis titles with labs().

ggplot(nuclei_annotated) +
  aes(x = Nuc_Total_ExpRep, fill = condition) +
  geom_density(alpha = 0.5) +
  labs(x = "total nuclear mScarlet", fill = "") +
  coord_cartesian(xlim = c(0, 2000)) +
  theme_minimal()

# Or I can make the x-axis logarithmic.
ggplot(nuclei_annotated) +
  aes(x = Nuc_Total_ExpRep, fill = condition) +
  geom_density(alpha = 0.5) +
  labs(x = "total nuclear mScarlet", fill = "") +
  # coord_cartesian(xlim = c(0, 2000)) +
  scale_x_log10() +
  theme_minimal()

# Or I can make the x-axis logarithmic and cut off the very very tiny values.
ggplot(nuclei_annotated) +
  aes(x = Nuc_Total_ExpRep, fill = condition) +
  geom_density(alpha = 0.5) +
  labs(x = "total nuclear mScarlet", fill = "") +
  coord_cartesian(xlim = c(3, 4000)) +
  scale_x_log10() +
  theme_minimal()

# Another thing we can look at - does dox induction change the size of the nuclei?
ggplot(nuclei_annotated) +
  aes(x = Nuc_AreaShape_Area, fill = condition) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 5000) +
  labs(x = "total nuclear area", fill = "") +
  theme_minimal()

# Just for fun, let's look at how the nuclear size and mean intensity of the mScarlet look together.
# I'm not using total mScarlet because this contains nuclear size as part of its calculation, so it would definitely be correlated.
ggplot(nuclei_annotated) +
  # I use 'colour' instead of 'fill' for dots, as they don't have a fill by default.
  aes(x = Nuc_AreaShape_Area, y = Nuc_Intensity_MeanIntensity_ExpRep, colour = condition) +
  # geom_point() tells ggplot to make a scatterplot.
  geom_point() +
  # We saw earlier than intensities are easier to look at on a log scale. This time, intensities are on the y-axis
  scale_y_log10() +
  theme_minimal()

# The effect isn't super clear, but we can see that it's probably not an effect caused by dox alone.
ggplot(nuclei_annotated) +
  aes(x = Nuc_Total_Srrm2, fill = condition) +
  geom_density(alpha = 0.5) +
  labs(x = "total nuclear SRRM2", fill = "") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 2000))

# # One more thing: how does dox induction change the total amount of SRRM2 in the nucleus?
nuclei_annotated %>%
  group_by(condition, replicate) %>%
  summarise(Nuc_Total_Srrm2 = mean(Nuc_Total_Srrm2)) %>%
  ungroup() ->
  nuclei_group_means

t.test(Nuc_Total_Srrm2 ~ condition, data = nuclei_group_means) %>% tidy() -> sc35_ttest

ggplot(nuclei_annotated, 
       aes(x = condition,
           y = Nuc_Total_Srrm2,
           # colour = factor(replicate),
           fill = factor(replicate)
       )
) +
  geom_violin(position = "identity", aes(colour = factor(replicate)), alpha = 0.2, show.legend = F) +
  geom_point(data = nuclei_group_means, size = 3, shape = 21, show.legend = F) +
  scale_colour_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  scale_fill_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  theme_classic() +
  labs(fill = "replicate",
       x = "",
       y = "average SC35 intensity per cell (A.U.)") +
  # coord_cartesian(ylim = c(0, 300)) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(color = "black")) +
  annotate("text", x = 1.5, y = 2500, label = paste0("p = ", round(sc35_ttest$p.value, digits = 3))) ->
  srrm2_intensity_plot

# nuclei_srrm2_plot <- ggplot(nuclei_annotated) +
#   aes(x = Nuc_Total_Srrm2, fill = condition) +
#   geom_density(alpha = 0.5) +
#   labs(x = "total nuclear SRRM2", fill = "") +
#   theme_minimal() +
#   coord_cartesian(xlim = c(0, 2000))
# 
# ggsave("WorkImages/20221130_NeveSRRM2ExpRep_2/plots/srrm2_per_speckle.pdf", nuclei_srrm2_plot, 
#        device = "pdf", units = "in",
#        width = 5, height = 4)


# This is useful to know when we interpret changes to the number, size, etc. of speckles

# Analysing things about speckles -----------------------------------------------------------------

# Lets combine our table of nuclear information with our table of speckle information.
# We'll use an inner join, so any speckles with a 'Parent_NucObjects' that belongs to a nucleus we filtered out.
merged_table <- inner_join(speckles, nuclei_annotated)

# Now we have info for each speckle about what condition and replicate it comes from, as well as
# the size of the nucleus, the intensity of the ER signal, etc.
names(merged_table)

# Remove nuclei that show leaky expression in control, or fail to express with dox.
merged_table %>%
  filter(AreaShape_Area > 1,
         Nuc_AreaShape_Area > 5000,
         ((condition == "control") & (Nuc_Total_ExpRep <= 100)) | 
         ((condition == "dox") & (Nuc_Total_ExpRep >= 100))) ->
  merged_table

# There are thousands of possible plots to make, but the most obvious things we can look at are:
# How do speckle size and intensity change with dox induction?

# We want to know this information on a nucleus-by-nucleus basis.
# To do this, we want to look at each nucleus one by one, then summarise the information.
# We can do this very quickly by grouping rows in the dataframe together based on their nucleus number.
# Again, the dplyr cheaetsheet is very helpful to understand this more intuitively.
# https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
# Look at the grouping section in the bottom of the second page.


merged_table %>% names()

speckle_average_table <- merged_table %>%
  # IMPORTANT - when we summarise the grouped data, we will lose any variables 
  # that we are not using to group the data, even if those variables don't break the groups up at all.
  # So here, we can keep the condition and replicate columns, even though they have totally different
  # groups of nuclei.
  group_by(condition, replicate, Parent_NucObjects, Nuc_AreaShape_Area, Nuc_Total_ExpRep) %>%
  # Now we can summarise different information on a nucleus-by-nucleus basis.
  # Let's take the mean signal inside each speckle and the mean size of each speckle.
  # Let's also take the number of speckles, which is equal to the number of rows per nucleus, using the function n().
  # Let's also keep the nuclear area, so we can compare the area of speckles to the area of the nucleus.
  summarise(average_speckle_intensity = mean(Intensity_MeanIntensity_Srrm2),
            average_speckle_area = mean(AreaShape_Area),
            total_speckle_area = sum(AreaShape_Area),
            number_of_speckles = n(),) %>%
  # It's good practice to remove the grouping variables, because this can unexpectedly interfere with
  # stuff later on.
  ungroup() %>%
  # The proportion of the 
  mutate(proportion_speckle_area = total_speckle_area/Nuc_AreaShape_Area,
         speckles_per_pixel = number_of_speckles/Nuc_AreaShape_Area)

# Now let's see how these things differ between the conditions.
ggplot(speckle_average_table) +
  aes(x = average_speckle_intensity, fill = condition) +
  geom_density(alpha = 0.5) +
  theme_minimal()

# We could also plot the different replicates seperately by giving them different linetypes.
# I had to convert the replicate variable into a 'factor' (remember the differene between 6 and "6").
# ggplot doesn't like to group things based on a numeric value.
ggplot(speckle_average_table) +
  aes(x = average_speckle_intensity, fill = condition, linetype = factor(replicate)) +
  geom_density(alpha = 0.3) +
  theme_minimal()

# There is some variability, but all the replicates from dox are lower than all the replicates from control on average.
ggplot(speckle_average_table) +
  aes(x = average_speckle_area, fill = condition, linetype = factor(replicate)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(0,200))

ggplot(speckle_average_table) +
  aes(x = number_of_speckles/Nuc_AreaShape_Area, fill = condition, linetype = factor(replicate)) +
  geom_density(alpha = 0.3) +
  theme_minimal()

ggplot(speckle_average_table) +
  aes(x = proportion_speckle_area, fill = condition, linetype = factor(replicate)) +
  geom_density(alpha = 0.3) +
  theme_minimal()

# We know the SRRM2 levels don't change between these conditions, but we can see that it's distributed very differently.
# dox+ cells have:
#   Fewer speckles
#   Dimmer speckles
#   Bigger speckles
#   The speckles take up less of the nucleus??

# Summarising replicates and doing statistics -----------------------------

speckle_group_means <- speckle_average_table %>%
  group_by(condition, replicate) %>%
  summarise(average_speckle_area = mean(average_speckle_area),
            average_speckle_intensity = mean(average_speckle_intensity),
            speckles_per_pixel = mean(speckles_per_pixel)
            ) %>%
  ungroup()

int_ttest <- t.test(average_speckle_intensity ~ condition, data = speckle_group_means) %>% tidy()
area_ttest <- t.test(average_speckle_area ~ condition, data = speckle_group_means) %>% tidy()
speckperpix_ttest <- t.test(speckles_per_pixel ~ condition, data = speckle_group_means) %>% tidy()

# speck_int_plot <- ggplot(speckle_average_table, 
#        aes(x = condition,
#            y = average_speckle_intensity,
#            fill = factor(replicate)
#            )
#        ) +
#   # geom_jitter()
#   # geom_beeswarm(cex = 0.75, alpha = 0.5, size = .25, aes(colour = factor(replicate),), show.legend = FALSE) +
#   geom_violin(position = "identity", aes(colour = factor(replicate)), alpha = 0.2, show.legend = F) +
#   geom_point(data = speckle_group_means, size = 3, shape = 21) +
#   scale_colour_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
#   scale_fill_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
#   theme_classic() +
#   labs(fill = "replicate",
#        x = "",
#        y = "average speckle intensity (au)") +
#   theme(axis.text.x = element_text(size = 10),
#         axis.text = element_text(color = "black")) +
#   annotate("text", x = 1.5, y = 0.2, label = paste0("p = ", round(int_ttest$p.value, digits = 3)))


speck_size_plot <- ggplot(speckle_average_table, 
       aes(x = condition,
           y = average_speckle_area,
           # colour = factor(replicate),
           fill = factor(replicate)
       )
  ) +
  geom_violin(position = "identity", aes(colour = factor(replicate)), alpha = 0.2, show.legend = F) +
  geom_point(data = speckle_group_means, size = 3, shape = 21, show.legend = F) +
  scale_colour_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  scale_fill_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  theme_classic() +
  labs(fill = "replicate",
       x = "",
       y = "average speckle size (pixels)") +
  coord_cartesian(ylim = c(0, 250)) +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(color = "black")) +
  annotate("text", x = 1.5, y = 200, label = paste0("p = ", round(area_ttest$p.value, digits = 3)))

speck_num_plot <- ggplot(speckle_average_table, 
       aes(x = condition,
           y = speckles_per_pixel,
           # colour = factor(replicate),
           fill = factor(replicate)
       )
) +
  geom_violin(position = "identity", aes(colour = factor(replicate)), alpha = 0.2, show.legend = F) +
  geom_point(data = speckle_group_means, size = 3, shape = 21, show.legend = F) +
  scale_colour_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  scale_fill_manual(values = c("#878787", "#b2b2b2", "#575756")) +
  theme_classic() +
  labs(fill = "replicate",
       x = "",
       y = "number of speckles per pixel") +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(color = "black")) +
  annotate("text", x = 1.5, y = 0.0055, label = paste0("p = ", round(speckperpix_ttest$p.value, digits = 7)))
  # coord_cartesian(ylim = c(0, 300))

merged_plots <- (srrm2_intensity_plot | speck_size_plot | speck_num_plot) + plot_layout(guides = "collect")

ggsave("GASR/image_analysis/20221130_NeveSRRM2ExpRep_2/plots/merged_plots_speckle_params.pdf", merged_plots,
       device = "pdf", units = "in",
       width = 9, height = 3)

# Dose dependence ---------------------------------------------------------


ggplot(speckle_average_table %>% 
         filter(condition == "dox",) %>%
         sample_frac(1L)) +
  aes(x = Nuc_Total_ExpRep, y = average_speckle_area, fill = factor(replicate)) +
  geom_point(shape = 21, size = 2, show.legend = F) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(100, 3000), ylim = c(10, 300)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  geom_smooth(method = "lm", aes(color = factor(replicate)), se = F, linetype = "dashed", size = 2, show.legend = F) +
  scale_colour_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  scale_fill_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  labs(fill = "replicate",
       x = "PPIG LCD signal (AU)",
       y = "average speckle size (pixels)") ->
  speckle_size_dose_plot

ggplot(speckle_average_table %>% 
         filter(condition == "dox") %>%
         sample_frac(1L)) +
  aes(x = Nuc_Total_ExpRep, y = speckles_per_pixel, fill = factor(replicate)) +
  geom_point(shape = 21, size = 2, show.legend = F) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(100, 3000), ylim = c(1e-04, 3.5e-03)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  geom_smooth(method = "lm", aes(color = factor(replicate)), se = F, linetype = "dashed", size = 2, show.legend = F) +
  scale_colour_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  scale_fill_manual(values = c("#a1b0c5", "#bd4f4d", "#e0bf07")) +
  labs(fill = "replicate",
       x = "PPIG LCD signal (AU)",
       y = "average speckle size (pixels)") ->
  speckle_per_pix_dose_plot

merged_dose <- (speckle_size_dose_plot | speckle_per_pix_dose_plot)

ggsave("WorkImages/20221130_NeveSRRM2ExpRep_2/plots/merged_dose_plot.pdf", merged_dose,
       device = "pdf", units = "in",
       width = 6.5, height = 3)

