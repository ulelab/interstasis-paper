rm(list = ls())

list.of.packages <- c("stringi", "tidyverse", "plyranges",
                      "parallel", "tictoc", "patchwork",
                      "viridis", "broom")

for(i in list.of.packages) {
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}


# Genome transcriptome map ------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", 
                               col_types = cols())

germ_clusters <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "cds"), sep = ":|-|@", remove = F) %>%
  dplyr::select(-cds) %>%
  mutate(strand = "+")

germ_clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) ->
  cluster_rename_df

# Set factor order
cluster_rename_df$cluster_name_short <- fct_relevel(cluster_rename_df$cluster_name_short, cluster_rename_df$cluster_name_short)

germ_gr <- makeGRangesFromDataFrame(germ_clusters, keep.extra.columns = T, seqnames.field = "transcript_id")

genome_transcriptome_map <- read_tsv("GASR/all_clip/longtx_exons.gff3",
                                     col_names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
                                     col_types = "cccddcccc",
                                     comment = "#") %>%
  mutate(transcript_id = word(attributes, 2, sep = "transcript_id=|;gene_type"),
         exon_number = as.numeric(word(attributes, 2, sep = "exon_number=|;exon_id"))) %>%
  select(-c(attributes, phase, source, type)) %>%
  mutate(width = end - start) %>%
  arrange(transcript_id, exon_number) %>%
  group_by(transcript_id) %>%
  mutate(tx_start = c(1, (1 + cumsum(width)[-length(transcript_id)])),
         tx_end = tx_start + width,
         exon_end = end,
         exon_start = start,
         tx_strand = "+") %>%
  ungroup() %>%
  dplyr::rename(genome_strand = strand,
                genome_start = start,
                genome_end = end)

tx_loc_genome_info <- makeGRangesFromDataFrame(genome_transcriptome_map,
                                               keep.extra.columns = T, 
                                               seqnames.field = "transcript_id",
                                               start.field = "tx_start",
                                               end.field = "tx_end",
                                               strand = "tx_strand")

genome_loc_tx_info <- makeGRangesFromDataFrame(genome_transcriptome_map,
                                               keep.extra.columns = T, 
                                               seqnames.field = "chr",
                                               start.field = "genome_start",
                                               end.field = "genome_end",
                                               strand = "genome_strand")
                                               

# GeRM exons --------------------------------------------------------------

germ_exons_tx <- find_overlaps(tx_loc_genome_info, germ_gr, minoverlap = 50) %>% unique()

germ_exons_df <- germ_exons_tx %>% as.data.frame()

germ_exons_df %>%
  dplyr::rename(transcript_id = seqnames,
                tx_start = start,
                tx_end = end,
                tx_strand = strand) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T,
                           seqnames.field = "chr",
                           start.field = "genome_start",
                           end.field = "genome_end",
                           strand.field = "genome_strand") ->
  germ_exons_genome

list(germ_exons_df %>% 
  dplyr::select(cluster, cluster_name, width) %>%
    left_join(cluster_rename_df),
genome_transcriptome_map %>% 
  mutate(cluster = 100,
         cluster_name_short = "non-multivalent",
         cluster_name = "non-multivalent") %>%
  dplyr::select(cluster, cluster_name, cluster_name_short, width)) %>%
  bind_rows() %>%
  mutate(cluster_name_short = cluster_name_short %>%
           fct_relevel(as.character(cluster_rename_df$cluster_name_short))) %>%
  filter(cluster_name_short != "None",
         width >= 50)->
  exon_width_df
       

ggplot(exon_width_df) +
  aes(y = width, x = cluster_name_short, fill = cluster_name_short) +
  geom_violin(show.legend = F, alpha = 0.5) +
  geom_boxplot(width = 0.3, show.legend = F, outlier.size = 0.25, alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_log10() +
  geom_hline(yintercept = median(genome_transcriptome_map$width), linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c(cluster_rename_df$cluster_colour[-1], "grey80") ) +
  labs(x = "", y = "exon length") +
  coord_cartesian(ylim = c(50, 20000)) ->
  germ_exon_lengths

ggsave("GASR/germs/plots/germ_exons/germ_exon_lengths.pdf", germ_exon_lengths, height = 3, width = 4)



pairwise.t.test(x = log(exon_width_df$width), g = exon_width_df$cluster_name_short) %>% tidy() %>%
 filter((group1 == "non-multivalent") | (group2 == "non-multivalent")) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  germ_exon_lengths_ttests

write_tsv(germ_exon_lengths_ttests, "GASR/germs/plots/germ_exons/germ_exon_lengths_ttests.tsv")

# TRA2 KD -----------------------------------------------------------------

se_rmats <- read_tsv("GASR/rnaseq/Best_et_al_2014_Tra2_kd/rmats/SE.MATS.JCEC.txt") %>% .[,-1] #duplicate ID column must be removed
mxe_rmats <- read_tsv("GASR/rnaseq/Best_et_al_2014_Tra2_kd/rmats/MXE.MATS.JCEC.txt") %>% .[,-1]
ri_rmats <- read_tsv("GASR/rnaseq/Best_et_al_2014_Tra2_kd/rmats/RI.MATS.JCEC.txt") %>% .[,-1] 

se_grange <- makeGRangesFromDataFrame(se_rmats, keep.extra.columns = T, seqnames.field = "chr", 
                         start.field = "exonStart_0base", end.fiel = "exonEnd", strand.field = "strand")

se_rmats$germ_exon <- case_when(
  count_overlaps(se_grange, germ_exons_genome %>% 
                   filter(cluster %in% c(1,2,7))) > 0 ~ "GA multivalency",
          T ~ "non-multivalent")

se_rmats$good_exon <- case_when(
  count_overlaps(se_grange, genome_loc_tx_info) > 0 ~ "good",
  T ~ "bad")

ggplot(se_rmats) +
  aes(x = IncLevelDifference, y = -log10(PValue), color = germ_exon) +
  geom_point()

ggplot(se_rmats %>% filter(good_exon == "good")) +
  aes(x = -IncLevelDifference, color = factor(germ_exon)) +
  geom_density() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("red3", "black")) +
  labs(x = "dPSI after TRA2 double KD", color = "") ->
  exon_skipping_t2kd_plot

ggsave("GASR/germs/plots/germ_exons/exon_skipping_t2kd_plot.pdf", exon_skipping_t2kd_plot, height = 3, width = 5)  

ggplot(se_rmats %>% 
        filter(good_exon == "good") %>%
         mutate(mean_psi = lapply(IncLevel1, function(x) x %>% strsplit(",") %>% unlist() %>% as.numeric() %>% mean()) %>%
                  unlist()) %>%
         drop_na()) +
  aes(x = mean_psi, color = factor(germ_exon)) +
  geom_density() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual(values = c("red3", "black")) +
  labs(x = "exon PSI in control", color = "") ->
  exon_inclusion_plot

ggsave("GASR/germs/plots/germ_exons/exon_inclusion_plot.pdf", exon_inclusion_plot, height = 3, width = 5)  


# Across tissues ----------------------------------------------------------

vast_table <- read_tsv("GASR/published_summary_data/vast_db/PSI_TABLE-hg38.tab.gz")

vast_table_gr <- vast_table %>%
  filter(grepl(pattern = "EX", EVENT)) %>% select(COORD) %>%
  separate(COORD, into = c("seqname", "start", "end")) %>%
  makeGRangesFromDataFrame()

vast_table_gr$one_of_my_exons <- count_overlaps(vast_table_gr, genome_loc_tx_info) > 0
vast_table_gr$germ_exon <- case_when(
  count_overlaps(vast_table_gr, germ_exons_genome %>% 
                   filter(cluster %in% c(1,2,7))) > 0 ~ "GA multivalency",
  T ~ "non-multivalent")

vast_table_gr %>% as.data.frame() %>%
  mutate(COORD = paste0(seqnames, ":", start, "-", end)) %>%
  dplyr::select(COORD, one_of_my_exons, germ_exon) ->
  vast_exon_info

vast_table %>%
  left_join(vast_exon_info, by= "COORD") %>%
  filter(one_of_my_exons == TRUE) %>%
  select(-contains("-Q")) %>%
  select(-c(FullCO)) %>%
  pivot_longer(cols = c(6:150), names_to = "tissue", values_to = "psi") ->
  vast_table_with_germ

vast_table_with_germ %>%
  filter(COMPLEX %in% c("S", "C1", "C2", "C3", "MIC", "ANN")) %>%
  drop_na() %>%
  group_by(germ_exon, tissue) %>%
  summarise(mean_psi = mean(psi),
            median_psi = median(psi),
            percentage_95 = mean(psi > 0.95)) ->
  mean_psi_tissues

ggplot(mean_psi_tissues) +
  aes(x = tissue, y = mean_psi, fill = germ_exon) +
  geom_point(shape = 21) +
  coord_cartesian(ylim = c(75, 100))
  
ggplot(mean_psi_tissues) +
  aes(x = germ_exon, y = percentage_95 * 100, fill = germ_exon) +
  geom_violin(show.legend = F, alpha = 0.5) +
  geom_boxplot(width = 0.3, outlier.size = 0.25, show.legend = F, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red3", "black")) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "", y = "percentage of exons\nconstitutive in tissue\n(PS > .95)") ->
  percentage_exons_constitutive_tissue

ggsave("GASR/germs/plots/germ_exons/percentage_exons_constitutive_tissue.pdf", percentage_exons_constitutive_tissue,
       width = 3, height = 3)  

t.test(mean_psi_tissues$percentage_95 ~ mean_psi_tissues$germ_exon) %>% tidy() ->
  percentage_exons_constitutive_tissue_ttest

write_tsv(percentage_exons_constitutive_tissue_ttest, "GASR/germs/plots/germ_exons/percentage_exons_constitutive_tissue_ttest.tsv")
  
  
