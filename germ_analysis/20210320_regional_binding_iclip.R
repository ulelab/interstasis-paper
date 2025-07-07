setwd("/camp/lab/ulej/home/users/farawar/")

rm(list = ls())

library(tidyverse)
library(cowplot)
library(patchwork)
library(lemon)

# mESC Mitosis -----------------------------------------------------------------

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
         filter(target %in% c("Tra2b", "Son", "Srrm2"),
                phase == "Asynchronous") %>%
         mutate(target = target %>% fct_relevel("Tra2b", "Srrm2")) %>%
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
  subtypes_tra2b_srrm2_son

ggsave("GASR/paper_figures/clip_subtypes_tra2b_srrm2_son.pdf", subtypes_tra2b_srrm2_son,
       width = 6, height = 3)

# mESC nuc/cyto -----------------------------------------------------------------

#Loading and formatting.

sample_names <- list.files("GASR/tmap/PRPF8_mESC_1/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

file_paths <- list.files("GASR/tmap/PRPF8_mESC_1/results/regions/", pattern = ".txt", full.names = T)


list.files("GASR/tmap/PRPF8_mESC_1/results/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

tables <- lapply(file_paths, function(x) read_tsv(x)) %>% 
  set_names(sample_names) %>% 
  bind_rows(.id = "sample") %>%
  separate(sample, into = c("target", "nothing", "location", "fraction", "replicate", "nothing2", "nothing3"), remove = F) %>%
  filter(location %in% c("nucleus", "cytoplasm")) %>%
  select(-c(nothing, nothing2, nothing3)) %>%
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
  group_by(sample, target, location, fraction, replicate, subtype_2) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(sample, target, location, fraction, replicate) %>%
  mutate(percentage = 100*count/sum(count)) %>%
  ungroup() %>%
  mutate(target = case_when(target == "prpf8" ~ "Prpf8", T ~ "Mock") %>% fct_relevel("Prpf8")) %>%
  mutate(location = case_when(location == "nucleus" ~ "Nucleus", T ~ "Cytoplasm")) %>%
  mutate(fraction = case_when(fraction == "highmw" ~ "High MW", fraction == "lowmw" ~ "Low MW", T ~ "All")) %>%
  mutate(sample = paste(target, location, fraction, replicate, sep = " ")) %>%
  mutate(subtype_2 = fct_relevel(subtype_2, "UTR5 mRNA", "CDS mRNA", "UTR3 mRNA", "intron pre-mRNA",
                                 "ncRNA rRNA", "ncRNA snRNA",  "ncRNA snoRNA", "other ncRNA"))

total_counts <- tables %>%
  group_by(sample, target, location, fraction, replicate) %>%
  summarise(count = sum(count))

#Total counts plot
ggplot(total_counts, aes(x = sample, y = count, fill = location)) +
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
ggplot(tables, aes(x = subtype_2, y = percentage, fill = location)) +
  facet_wrap(target ~ fraction) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", width = 0.5, fun.data = mean_se, position = position_dodge(width = 0.9), size = 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7) +
  scale_y_continuous(limits=c(0,90), expand = c(0,0)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.9), show.legend = F) + 
  scale_shape(solid = FALSE) + 
  labs(y = "Percentage of reads", 
       x = "", 
       fill = "") +
  theme_minimal_grid(line_size = 0.1, font_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        legend.position = "top",
        legend.justification=0.5,
        legend.text = element_text(size = 10, hjust = 0))

# mESC nuc/cyto -----------------------------------------------------------------

#Loading and formatting.

sample_names <- list.files("GASR/tmap/PRPF8_mESC_1/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

file_paths <- list.files("GASR/tmap/PRPF8_mESC_1/results/regions/", pattern = ".txt", full.names = T)

tables <- lapply(file_paths, function(x) read_tsv(x)) %>% 
  set_names(sample_names) %>% 
  bind_rows(.id = "sample") %>%
  separate(sample, into = c("target", "nothing", "location", "fraction", "replicate", "nothing2", "nothing3"), remove = F) %>%
  filter(location %in% c("nucleus", "cytoplasm")) %>%
  select(-c(nothing, nothing2, nothing3)) %>%
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
  group_by(sample, target, location, fraction, replicate, subtype_2) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(sample, target, location, fraction, replicate) %>%
  mutate(percentage = 100*count/sum(count)) %>%
  ungroup() %>%
  mutate(target = case_when(target == "prpf8" ~ "Prpf8", T ~ "Mock") %>% fct_relevel("Prpf8")) %>%
  mutate(location = case_when(location == "nucleus" ~ "Nucleus", T ~ "Cytoplasm")) %>%
  mutate(fraction = case_when(fraction == "highmw" ~ "High MW", fraction == "lowmw" ~ "Low MW", T ~ "All")) %>%
  mutate(sample = paste(target, location, fraction, replicate, sep = " ")) %>%
  mutate(subtype_2 = fct_relevel(subtype_2, "UTR5 mRNA", "CDS mRNA", "UTR3 mRNA", "intron pre-mRNA",
                                 "ncRNA rRNA", "ncRNA snRNA",  "ncRNA snoRNA", "other ncRNA"))

total_counts <- tables %>%
  group_by(sample, target, location, fraction, replicate) %>%
  summarise(count = sum(count))

#Total counts plot
ggplot(total_counts, aes(x = sample, y = count, fill = location)) +
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
ggplot(tables, aes(x = subtype_2, y = percentage, fill = location)) +
  facet_wrap(target ~ fraction) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", width = 0.5, fun.data = mean_se, position = position_dodge(width = 0.9), size = 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7) +
  scale_y_continuous(limits=c(0,90), expand = c(0,0)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.9), show.legend = F) + 
  scale_shape(solid = FALSE) + 
  labs(y = "Percentage of reads", 
       x = "", 
       fill = "") +
  theme_minimal_grid(line_size = 0.1, font_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        legend.position = "top",
        legend.justification=0.5,
        legend.text = element_text(size = 10, hjust = 0))

# PRPF8 and SMB in ENCODE lines-----------------------------------------------------------------

#Loading and formatting.

prpf8_sample_names <- list.files("GASR/tmap/PRPF8_ENCODElines/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

prpf8_file_paths <- list.files("GASR/tmap/PRPF8_ENCODElines/results/regions/", pattern = ".txt", full.names = T)

smb_sample_names <- list.files("GASR/tmap/SmB_ENCODElines/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

smb_file_paths <- list.files("GASR/tmap/SmB_ENCODElines/results/regions/", pattern = ".txt", full.names = T)

tables <- list("SmB" = lapply(smb_file_paths, function(x) read_tsv(x)) %>% 
                 set_names(smb_sample_names) %>% 
                 bind_rows(.id = "sample"),
               "PRPF8" = lapply(prpf8_file_paths, function(x) read_tsv(x)) %>% 
                 set_names(prpf8_sample_names) %>% 
                 bind_rows(.id = "sample")) %>%
  bind_rows(.id = "target") %>%
  separate(sample, into = c("cell_line", "condition", "replicate"), remove = F) %>%
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
  group_by(sample, target, cell_line, condition, replicate, subtype_2) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(sample, target, cell_line, condition, replicate) %>%
  mutate(percentage = 100*count/sum(count)) %>%
  ungroup() %>%
  mutate(cell_line = case_when(cell_line == "hepg2" ~ "HepG2", T ~ "K562")) %>%
  mutate(condition = case_when(condition == "med" ~ "Med. Lysis", condition == "mild" ~ "Mild Lysis",
                               condition == "high" ~ "High MW", T ~ "Low MW")) %>%
  mutate(sample = paste(target, cell_line, condition, replicate, sep = " ")) %>%
  mutate(subtype_2 = fct_relevel(subtype_2, "UTR5 mRNA", "CDS mRNA", "UTR3 mRNA", "intron pre-mRNA",
                                 "ncRNA rRNA", "ncRNA snRNA",  "ncRNA snoRNA", "other ncRNA"))

total_counts <- tables %>%
  group_by(sample, target, cell_line, condition, replicate) %>%
  summarise(count = sum(count))

#Total counts plot
ggplot(total_counts, aes(x = sample, y = count, fill = condition)) +
  facet_wrap(cell_line ~ target, scales = "free_y") +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
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
ggplot(tables, aes(x = subtype_2, y = percentage, fill = condition)) +
  facet_wrap(cell_line ~ target) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", width = 0.5, fun.data = mean_se, position = position_dodge(width = 0.9), size = 1) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_y_continuous(limits=c(0,75), expand = c(0,0)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.9), show.legend = F) + 
  scale_shape(solid = FALSE) + 
  labs(y = "Percentage of reads", 
       x = "", 
       fill = "") +
  theme_minimal_grid(line_size = 0.1, font_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        legend.position = "top",
        legend.justification=0.5,
        legend.text = element_text(size = 10, hjust = 0))

# SmB Mitosis thing -------------------------------------------------------


#Loading and formatting.

sample_names <- list.files("GASR/tmap/SmB_cellcycle/results/regions/", pattern = ".txt") %>%
  sapply(function(x) str_remove(x, ".subtypes.txt")) %>%
  unname()

file_paths <- list.files("GASR/tmap/SmB_cellcycle/results/regions/", pattern = ".txt", full.names = T)

tables <- lapply(file_paths, function(x) read_tsv(x)) %>% 
  set_names(sample_names) %>% 
  bind_rows(.id = "sample") %>%
  separate(sample, into = c("phase", "replicate"), remove = F) %>%
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
  group_by(sample, phase, replicate, subtype_2) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(sample, phase, replicate) %>%
  mutate(percentage = 100*count/sum(count)) %>%
  ungroup() %>%
  mutate(phase = phase %>% fct_relevel("G1", "S", "G2", "M")) %>%
  mutate(subtype_2 = fct_relevel(subtype_2, "UTR5 mRNA", "CDS mRNA", "UTR3 mRNA", "intron pre-mRNA",
                                 "ncRNA rRNA", "ncRNA snRNA",  "ncRNA snoRNA", "other ncRNA"))

total_counts <- tables %>%
  group_by(sample, phase, replicate) %>%
  summarise(count = sum(count))

#Total counts plot
ggplot(total_counts, aes(x = sample, y = count, fill = phase)) +
  facet_wrap( ~ phase, scales = "free_y") +
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
ggplot(tables, aes(x = subtype_2, y = percentage, fill = phase)) +
  # facet_wrap(target ~ fraction) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", width = 0.5, fun.data = mean_se, position = position_dodge(width = 0.9), size = 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7) +
  scale_y_continuous(limits=c(0,62), expand = c(0,0)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.9), show.legend = F) + 
  scale_shape(solid = FALSE) + 
  labs(y = "Percentage of reads", 
       x = "", 
       fill = "") +
  theme_minimal_grid(line_size = 0.1, font_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        legend.position = "top",
        legend.justification=0.5,
        legend.text = element_text(size = 10, hjust = 0))

