setwd("/camp/lab/ulej/home/users/farawar/")

rm(list = ls())

library(tidyverse)
library(cowplot)
library(patchwork)
library(lemon)

# mESC -----------------------------------------------------------------

#Loading and formatting.

sample_names <- list.files("GASR/tmap/mesc_mitosis_1/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

file_paths <- list.files("GASR/tmap/mesc_mitosis_1/results/regions/", pattern = ".txt", full.names = T)

tables <- lapply(file_paths, function(x) read_tsv(x)) %>% 
  set_names(sample_names) %>% 
  bind_rows(.id = "sample") %>%
  separate(sample, into = c("target", "phase", "fraction", "replicate"), remove = F) %>%
  filter(fraction != "Wildtype") %>%
  mutate(subtype_2 = recode(subtype,
                            `ncRNA lncRNA` = "other ncRNA",
                            `intron lncRNA` = "other ncRNA",
                            `ncRNA sRNA` = "other ncRNA",
                            `ncRNA miRNA` = "other ncRNA",
                            `ncRNA mt_tRNA` = "other ncRNA",
                            `ncRNA mt_rRNA` = "other ncRNA",
                            `ncRNA Mt_rRNA` = "other ncRNA",
                            `ncRNA tRNA` = "other ncRNA",
                            `ncRNA rRNA_pseudogene` = "other ncRNA")) %>%
  group_by(sample, target, phase, fraction, replicate, subtype_2) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(sample, target, phase, fraction, replicate) %>%
  mutate(percentage = 100*count/sum(count)) %>%
  ungroup() %>%
  mutate(target = fct_relevel(target, "Prpf8", "SmB", "Son",  "Tra2b", "Srrm2")) %>%
  mutate(subtype_2 = fct_relevel(subtype_2, "UTR5 mRNA", "CDS mRNA", "UTR3 mRNA", "intron pre-mRNA",
                                 "ncRNA rRNA", "ncRNA snRNA",  "ncRNA snoRNA", "other ncRNA"))

total_counts <- tables %>%
  group_by(sample, target, phase, fraction, replicate) %>%
  summarise(count = sum(count))

#Total counts plot
ggplot(total_counts, aes(x = sample, y = count, fill = phase)) +
  facet_wrap(target ~ fraction, scales = "free_y") +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  scale_fill_viridis_d(begin = 0.3, end = 0.7) +
  labs(y = "Number of mapped reads", 
       x = "", 
       fill = "") +
  theme_minimal_grid(line_size = 0.1, font_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        legend.position = "top",
        legend.justification=0.4,
        legend.text = element_text(size = 10, hjust = 0)) +
  coord_flip()

#Subtype percentages plot
ggplot(tables %>%
         filter(target %in% c("Tra2b"),
                phase == "Asynchronous") %>%
         arrange(target), 
       aes(x = subtype_2, y = percentage, fill = fraction)) +
  facet_wrap(. ~ target) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", width = 0.5, fun.data = mean_se, position = position_dodge(width = 0.9), size = 1) +
  scale_fill_manual(values = c("grey70", "#6aa358", "#78689e")) +
  scale_y_continuous(limits=c(0,60), expand = c(0,0)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.9), show.legend = F) + 
  scale_shape(solid = FALSE) + 
  labs(y = "Percentage of reads", 
       x = "", 
       fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 10, hjust = 0)) ->
  subtypes_tra2b

ggsave("GASR/paper_figures/clip_subtypes_tra2b.pdf", subtypes_tra2b
       width = 4, height = 3)
