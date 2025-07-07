rm(list = ls())
setwd("/camp/lab/ulej/home/users/farawar/")

# Packages ----------------------------------------------------------------

list.of.packages <- c("stringi", "tidyverse", "parallel", "ggrepel", "tictoc", "patchwork", "viridis", "BiocManager", "factoextra", "broom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/")

bioconductor.packages <- c("plyranges")
new.bioconductor <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()["Package"])]
for(i in new.bioconductor){
  BiocManager::install(i, lib = .libPaths()[1])
}

for(i in c(list.of.packages, bioconductor.packages)){
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}

# Functions ---------------------------------------------------------------

read_crosslinks <- function(xl_path){
  # xl_path <- "/camp/lab/ulej/home/users/kuretk/data/merge_eclip_replicates_20210906/merged_files/HepG2-TRA2A-merged.xl.bed.gz"
  read_tsv(xl_path, 
           col_names = c("chr", "start", "end", "nothing", "count", "strand"), 
           col_types = "cddcdc") %>%
    select(-nothing) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    anchor_5p() %>%
    mutate(width = 1)
}

convert_to_transcriptome <- function(xl_path, gtm){
  # xl_path <- eclip_files[which(eclip_names == eclip_name)]
  
  # gtm <- genome_transcriptome_map
  
  find_overlaps_directed(read_crosslinks(xl_path), gtm) %>%
    as.data.frame() %>%
    mutate(exon_position = case_when(strand == "+" ~ start - exon_start,
                                     strand == "-" ~ exon_end - start)) %>%
    mutate(transcript_position = tx_pos + exon_position) %>%
    select(transcript_id, transcript_position, count) %>%
    makeGRangesFromDataFrame(seqnames.field = "transcript_id",
                             start.field = "transcript_position",
                             end.field = "transcript_position",
                             keep.extra.columns = T)
}

percentile_per_kmer_eclip <- function(eclip_name, region_filter) {
  
  # eclip_name <- "K562-TRA2A"
  # region_filter <- "CDS"
  print(eclip_name)
  
  output_name <- paste0("GASR/germs/data/clip/individual_kmer_multivalency_percentiles/", 
                        region_filter, "_",
                        eclip_name, ".txt.gz")
  
  # Reload already output data (if loop failed midway through when run before) (it did)
  if(file.exists(output_name)){
    return(read_tsv(output_name))
  }
  
  clip <- convert_to_transcriptome(eclip_files[which(eclip_names == eclip_name)],
                                   genome_transcriptome_map)
  
  
  germs_ranges %>%
    filter(region == region_filter) %>%
    mutate(xlink_counts = count_overlaps(., clip)) %>%
    as.data.frame() %>%
    group_by(kmer) %>%
    summarise(sum_xl = sum(xlink_counts),
              xl_per_kmer = sum_xl/dplyr::n(),
              percentile = sum(kmer_percentile * xlink_counts)/sum(xlink_counts),
              average_scaled_mv = sum(scaled_mv * xlink_counts)/sum(xlink_counts)) ->
    percentile_table
  
  write_tsv(percentile_table, output_name) 
  
  return(percentile_table)
  
}

percentile_per_kmer_iclip <- function(iclip_name, region_filter) {
  
  # iclip_name <- "HEK293-HNRNPH1"
  # region_filter <- "UTR3"
  print(iclip_name)
  
  output_name <- paste0("GASR/germs/data/clip/individual_kmer_multivalency_percentiles/", 
                        region_filter, "_",
                        iclip_name, ".txt.gz")
  
  # Reload already output data (if loop failed midway through when run before) (it did)
  if(file.exists(output_name)){
    return(read_tsv(output_name))
  }
  
  clip <- convert_to_transcriptome(iclip_files[which(iclip_names == iclip_name)],
                                   genome_transcriptome_map)
  
  
  germs_ranges %>%
    filter(region == region_filter) %>%
    mutate(xlink_counts = count_overlaps(., clip)) %>%
    as.data.frame() %>%
    group_by(kmer) %>%
    summarise(sum_xl = sum(xlink_counts),
              xl_per_kmer = sum_xl/dplyr::n(),
              percentile = sum(kmer_percentile * xlink_counts)/sum(xlink_counts),
              average_scaled_mv = sum(scaled_mv * xlink_counts)/sum(xlink_counts)) ->
    percentile_table
  
  write_tsv(percentile_table, output_name) 
  
  return(percentile_table)
  
}

# Parameters --------------------------------------------------------------

k_len <- 5

# Get multivalent regions -------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", 
                               col_types = cols())


germs <- read_tsv(file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")

# For normalising the multivalency
minkm <- min(germs$kmer_multivalency)
medkm <- median(germs$kmer_multivalency)

# Probably won't use smoothed values here, but good to keep it consistent. Ignore the 0s from the edges of transcripts.
minskm <- min(germs$smoothed_kmer_multivalency[germs$smoothed_kmer_multivalency > 0])
medskm <- median(germs$smoothed_kmer_multivalency[germs$smoothed_kmer_multivalency > 0])

# Normalise the multivalency
germs <- germs %>%
  mutate(med_norm = 
           (kmer_multivalency - minkm)/
           (medkm - minkm),
         med_norm_smooth = 
           (smoothed_kmer_multivalency - minskm)/
           (medskm - minskm))

germs_ranges <- left_join(germs, transcript_details, by = "transcript_id") %>%
  dplyr::select(-c(peak_start, peak_end, peak_location)) %>%
  group_by(transcript_id) %>%
  mutate(start = row_number()) %>%
  ungroup() %>%
  mutate(end = start + 4,
         region = case_when(start < cds_start ~ "UTR5", #Yes, I am aware that I am not taking k-length into account here.
                            start > cds_end ~ "UTR3", #How would you prefer I classify boundary cases?
                            T ~ "CDS")) %>%
  group_by(region, kmer) %>%
  mutate(kmer_percentile = rank(med_norm) / dplyr::n(),
         scaled_mv = scale(med_norm)) %>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = T,
                           seqnames.field = "transcript_id")

# Get crosslink locations -------------------------------------------------

eclip_folder <- "/camp/lab/ulej/home/users/kuretk/data/merge_eclip_replicates_20210906/merged_files/"

# For some reason, the TDP data is a file of size 0? It binds introns/3UTRs anyway, so fuck it.
eclip_files <- list.files(eclip_folder, pattern = ".gz", full.names = T) %>% str_subset("K562-TARDBP", negate = T)
eclip_names <- list.files(eclip_folder, pattern = ".gz", full.names = F) %>% word(1, sep = "-merged") %>% str_subset("K562-TARDBP", negate = T)

sus_eclip <- read_lines("GASR/lists/eclip_datasets_clara_is_suspicious_of.txt")

iclip_files <- c("GASR/tmap/EJC_Hauer/results/xlinks_merged_genome/casc3.xlinks.bed",
                 "GASR/tmap/EJC_Hauer/results/xlinks_merged_genome/eif4a3.xlinks.bed",
                 "GASR/tmap/EJC_Hauer/results/xlinks_merged_genome/rnps1.xlinks.bed",
                 "GASR/tmap/HNRNPH1_Braun/results/xlinks_merged_genome/hnrnph1_wt.xlinks.bed",
                 "GASR/tmap/TRA2B_Best_et_al/results/xlinks_merged_genome/TRA2B.xlinks.bed",
                 "GASR/tmap/SR_Krchnakova_Hela/results/xlinks_merged_genome/srsf2.xlinks.bed",
                 "GASR/tmap/SR_Krchnakova_Hela/results/xlinks_merged_genome/srsf5.xlinks.bed",
                 "GASR/tmap/SR_Krchnakova_Hela/results/xlinks_merged_genome/srsf6.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l_clip.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l2_clip.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l3_clip.xlinks.bed")

iclip_names <- c("Hela-CASC3",
                 "Hela-EIF4A3",
                 "Hela-RNPS1",
                 "HEK293-HNRNPH1",
                 "MDAMB231-TRA2B",
                 "Hela-SRSF2",
                 "Hela-SRSF5",
                 "Hela-SRSF6",
                 "K562-LUC7L",
                 "K562-LUC7L2",
                 "K562-LUC7L3")

# Make a map to translate genomic coords to transcriptomic coords ------------------------------

genome_transcriptome_map <- read_tsv("GASR/all_clip/longtx_exons.gff3",
                                     col_names = c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
                                     col_types = "cccddcccc",
                                     comment = "#") %>%
  mutate(transcript_id = word(attributes, 2, sep = "transcript_id=|;gene_type"),
         exon_number = as.numeric(word(attributes, 2, sep = "exon_number=|;exon_id"))) %>%
  select(-c(attributes, phase, source, type)) %>%
  mutate(width = end - start) %>%
  arrange(transcript_id, exon_number) %>%
  group_by(transcript_id) %>%
  mutate(tx_pos = c(1, (1 + cumsum(width)[-length(transcript_id)])),
         exon_end = end,
         exon_start = start) %>% # We need these columns as metadata
  makeGRangesFromDataFrame(keep.extra.columns = T)

# Calculate eclip preferences ------------------------------------------------------------

tic("Big calculation")
list_eclip_kmer_percentiles_cds <- lapply(eclip_names, function(x) percentile_per_kmer_eclip(x, "CDS"))
toc()

big_eclip_kmer_percentiles_cds <- list_eclip_kmer_percentiles_cds %>%
  set_names(eclip_names) %>%
  bind_rows(.id = "dataset")

write_tsv(big_eclip_kmer_percentiles_cds,
          "GASR/germs/data/clip/crosslinks_and_kmer_multivalencies_cds.txt.gz")


tic("UTR3 Calculation")
list_eclip_kmer_percentiles_utr3 <- mclapply(eclip_names, function(x) percentile_per_kmer_eclip(x, "UTR3"), mc.cores = 8)
toc()

big_eclip_kmer_percentiles_utr3 <- list_eclip_kmer_percentiles_utr3 %>%
  set_names(eclip_names) %>%
  bind_rows(.id = "dataset")

write_tsv(big_eclip_kmer_percentiles_utr3,
          "GASR/germs/data/clip/crosslinks_and_kmer_multivalencies_utr3.txt.gz")

# Calculate iclip preferences ------------------------------------------------------------

tic("Big calculation")
list_iclip_kmer_percentiles_cds <- lapply(iclip_names, function(x) percentile_per_kmer_iclip(x, "CDS"))
toc()

big_iclip_kmer_percentiles_cds <- list_iclip_kmer_percentiles_cds %>%
  set_names(iclip_names) %>%
  bind_rows(.id = "dataset") %>%
  mutate(gene_name = word(dataset, 2, sep = "-"))

write_tsv(big_iclip_kmer_percentiles_cds,
          "GASR/germs/data/clip/iclip_crosslinks_and_kmer_multivalencies_cds.txt.gz")

tic("UTR3 Calculation")
list_iclip_kmer_percentiles_utr3 <- mclapply(iclip_names, function(x) percentile_per_kmer_iclip(x, "UTR3"), mc.cores = 8)
toc()

big_iclip_kmer_percentiles_utr3 <- list_iclip_kmer_percentiles_utr3 %>%
  set_names(iclip_names) %>%
  bind_rows(.id = "dataset") %>%
  mutate(gene_name = word(dataset, 2, sep = "-"))

write_tsv(big_iclip_kmer_percentiles_utr3,
          "GASR/germs/data/clip/iclip_crosslinks_and_kmer_multivalencies_utr3.txt.gz")

# Reload ------------------------------------------------------------------

big_eclip_kmer_percentiles_cds <- read_tsv("GASR/germs/data/clip/crosslinks_and_kmer_multivalencies_cds.txt.gz") %>%
  separate(dataset, into = c("cell_line", "gene_name"), remove = F)

big_eclip_kmer_percentiles_utr3 <- read_tsv("GASR/germs/data/clip/crosslinks_and_kmer_multivalencies_utr3.txt.gz") %>%
  separate(dataset, into = c("cell_line", "gene_name"), remove = F)

big_iclip_kmer_percentiles_cds <- read_tsv("GASR/germs/data/clip/iclip_crosslinks_and_kmer_multivalencies_cds.txt.gz") %>%
  separate(dataset, into = c("cell_line", "gene_name"), remove = F)

big_iclip_kmer_percentiles_utr3 <- read_tsv("GASR/germs/data/clip/iclip_crosslinks_and_kmer_multivalencies_utr3.txt.gz") %>%
  separate(dataset, into = c("cell_line", "gene_name"), remove = F)

# CDS analysis ------------------------------------------------------------

# eclip
big_eclip_kmer_percentiles_cds %>%
  group_by(dataset) %>%
  summarise(crosslinks_in_cds = sum(sum_xl)) -> 
  crosslinks_per_dataset_cds

big_eclip_kmer_percentiles_cds %>%
  group_by(dataset) %>%
  mutate(strong_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 T ~ "weak"),
         graded_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 (rank(-xl_per_kmer) <= 200) & (rank(-xl_per_kmer) > 50) ~ "medium",
                                 T ~ "weak"),) %>%
  # group_by(dataset, strong_kmer) %>%
  group_by(dataset, gene_name, graded_kmer) %>%
  summarise(mean_perc = mean(percentile)) %>%
  pivot_wider(names_from = "graded_kmer", values_from = "mean_perc") %>%
  mutate(ratio = strong/weak) %>%
  left_join(crosslinks_per_dataset_cds) ->
  eclip_kmer_percentile_ratios_cds

crosslinks_per_dataset_cds %>%
  filter(crosslinks_in_cds > quantile(crosslinks_in_cds, 0.33)) %>%
  select(dataset) %>%
  unlist(use.names = F) ->
  datasets_with_enough_counts_to_use_for_further_data_analysis_cds

# Public clip
big_iclip_kmer_percentiles_cds %>%
  group_by(dataset) %>%
  summarise(crosslinks_in_cds = sum(sum_xl)) -> 
  crosslinks_per_dataset_cds

big_iclip_kmer_percentiles_cds %>%
  group_by(dataset) %>%
  mutate(strong_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 T ~ "weak"),
         graded_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 (rank(-xl_per_kmer) <= 200) & (rank(-xl_per_kmer) > 50) ~ "medium",
                                 T ~ "weak"),) %>%
  # group_by(dataset, strong_kmer) %>%
  group_by(dataset, gene_name, graded_kmer) %>%
  summarise(mean_perc = mean(percentile)) %>%
  pivot_wider(names_from = "graded_kmer", values_from = "mean_perc") %>%
  mutate(ratio = strong/weak) %>%
  left_join(crosslinks_per_dataset_cds) ->
  iclip_kmer_percentile_ratios_cds

crosslinks_per_dataset_cds %>%
  filter(crosslinks_in_cds > quantile(crosslinks_in_cds, 0.33)) %>%
  select(dataset) %>%
  unlist(use.names = F) ->
  datasets_with_enough_counts_to_use_for_further_data_analysis_cds

# UTR3 analysis ------------------------------------------------------------

# ECNODE
big_eclip_kmer_percentiles_utr3 %>%
  group_by(dataset) %>%
  summarise(crosslinks_in_utr3 = sum(sum_xl)) -> 
  crosslinks_per_dataset_utr3

big_eclip_kmer_percentiles_utr3 %>%
  group_by(dataset) %>%
  mutate(strong_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 T ~ "weak"),
         graded_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 (rank(-xl_per_kmer) <= 200) & (rank(-xl_per_kmer) > 50) ~ "medium",
                                 T ~ "weak"),) %>%
  # group_by(dataset, strong_kmer) %>%
  group_by(dataset, gene_name, graded_kmer) %>%
  summarise(mean_perc = mean(percentile)) %>%
  pivot_wider(names_from = "graded_kmer", values_from = "mean_perc") %>%
  mutate(ratio = strong/weak) %>%
  left_join(crosslinks_per_dataset_utr3) ->
  eclip_kmer_percentile_ratios_utr3

crosslinks_per_dataset_utr3 %>%
  filter(crosslinks_in_utr3 > quantile(crosslinks_in_utr3, 0.33)) %>%
  select(dataset) %>%
  unlist(use.names = F) ->
  datasets_with_enough_counts_to_use_for_further_data_analysis_utr3

# Public clip
big_iclip_kmer_percentiles_utr3 %>%
  group_by(dataset) %>%
  summarise(crosslinks_in_utr3 = sum(sum_xl)) -> 
  crosslinks_per_dataset_utr3

big_iclip_kmer_percentiles_utr3 %>%
  group_by(dataset) %>%
  mutate(strong_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 T ~ "weak"),
         graded_kmer = case_when(rank(-xl_per_kmer) <= 50 ~ "strong",
                                 (rank(-xl_per_kmer) <= 200) & (rank(-xl_per_kmer) > 50) ~ "medium",
                                 T ~ "weak"),) %>%
  # group_by(dataset, strong_kmer) %>%
  group_by(dataset, gene_name, graded_kmer) %>%
  summarise(mean_perc = mean(percentile)) %>%
  pivot_wider(names_from = "graded_kmer", values_from = "mean_perc") %>%
  mutate(ratio = strong/weak) %>%
  left_join(crosslinks_per_dataset_utr3) ->
  iclip_kmer_percentile_ratios_utr3

crosslinks_per_dataset_utr3 %>%
  filter(crosslinks_in_utr3 > quantile(crosslinks_in_utr3, 0.33)) %>%
  select(dataset) %>%
  unlist(use.names = F) ->
  datasets_with_enough_counts_to_use_for_further_data_analysis_utr3

# Plots? ------------------------------------------------------------------

eclip_kmer_percentile_ratios_cds %>%
  ggplot() +
  aes(x = crosslinks_in_cds, y = ratio) +
  geom_point() +
  geom_vline(xintercept = median(eclip_kmer_percentile_ratios_cds$crosslinks_in_cds))


eclip_kmer_percentile_ratios_utr3 %>%
  # filter(dataset %in% datasets_with_enough_counts_to_use_for_further_data_analysis) %>%
  ggplot() +
  aes(x = strong) +
  geom_density() +
  scale_x_log10()

big_eclip_kmer_percentiles_cds %>%
  filter(dataset == "K562-UPF1") %>%
  # filter(gene_name == "PPIG") %>%
  ggplot() +
  aes(x = xl_per_kmer, y = percentile * 100) +
  geom_point(fill = "white", shape = 21, size = 1.5) +
  scale_x_log10(breaks = c(0.07, 0.1, 0.15)) +
  coord_cartesian(ylim = c(25, 75)) +
  theme_classic() +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  labs(x = "crosslinks per kmer",
       y = "mean GeRM percentile",
       title = "UPF1") +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5)) ->
  UPF1_mv_pref_plot

ggsave("GASR/germs/plots/UPF1_eCLIP_multivalency_preference_dotplot.pdf", UPF1_mv_pref_plot,
       width = 2.5, height = 2.5)


big_iclip_kmer_percentiles_cds %>%
  filter(dataset == "MDAMB231-TRA2B") ->
  t2b_dataset


t2b_dataset %>%
  ggplot() +
  aes(x = xl_per_kmer, y = percentile * 100) +
  geom_point(fill = "white", shape = 21, size = 1.5) +
  geom_label_repel(data = t2b_dataset %>%
                    arrange(desc(xl_per_kmer)) %>%
                    slice(1:5),
                  aes(label = kmer), max.iter = Inf, size = 3, force = 10,
                  alpha = 0.75) +
  scale_x_log10() +
  theme_classic() +
  coord_cartesian(ylim = c(20, 80)) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  labs(x = "crosslinks per kmer",
       y = "mean GeRM percentile",
       title = "TRA2B") +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5)) ->
  TRA2B_mv_pref_plot

ggsave("GASR/germs/plots/TRA2B_iCLIP_multivalency_preference_dotplot.pdf", TRA2B_mv_pref_plot,
       width = 2.5, height = 2.5)


# Save lists of multivalent binders ---------------------------------------

eclip_kmer_percentile_ratios_cds %>%
  ungroup() %>%
  filter(crosslinks_in_cds > quantile(crosslinks_in_cds, 0.33),
         ratio > 1.1) %>%
  select(dataset) %>%
  unlist(use.names = F) %>%
  write_lines("GASR/germs/data/clip/multivalent_eclip_datasets_cds.txt")


eclip_kmer_percentile_ratios_utr3 %>%
  ungroup() %>%
  filter(crosslinks_in_utr3 > quantile(crosslinks_in_utr3, 0.33),
         ratio > 1.1) %>%
  select(dataset) %>%
  unlist(use.names = F) %>%
  write_lines("GASR/germs/data/clip/multivalent_eclip_datasets_utr3.txt")



# This section didn't go in the paper but I think it was kinda a cool idea. I don't think the data was good enough at this point to be super conclusive though.

# Loading alphafold data ------------------------------------------------------------

gencode_to_swissprot <- read_tsv("GASR/lists/gencode.v29.metadata.SwissProt.txt",
                                 col_names = c("transcript_id", "uniprot_id", "id2"),
                                 col_type = "c")

af2_disorder <- read_tsv("General/alpha_fold_disorder/af2_human_disorder_pred.tsv.gz",
                         skip = 1,
                         col_names = c("name", "pos", "aa", "lddt", "disorder", "rsa", "ss", "disorder_25", "binding_25_0581"),
                         col_types = "cdcdddcdd")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

transcript_details_uniprot <- inner_join(gencode_to_swissprot, transcript_details)

transcript_af2_disorder <- af2_disorder %>%
  group_by(name) %>%
  nest(data = c(pos, aa, lddt, disorder, rsa, ss, disorder_25, binding_25_0581)) %>%
  ungroup() %>%
  mutate(uniprot_id = word(name, 2, sep = fixed("-"))) %>%
  dplyr::select(-name) %>%
  inner_join(transcript_details_uniprot)

#Does the length of the alpha fold vector correctly match the length of the CDS? If not, discard.
transcript_af2_disorder <- transcript_af2_disorder %>%
  ungroup() %>%
  mutate(disorder_nrow = map(data, ~{nrow(.x)}) %>% unlist()) %>%
  filter((disorder_nrow + 1) == (cds_length/3)) %>%
  dplyr::select(-disorder_nrow)
                                     
transcript_af2_disorder %>%
  mutate(median_lddt = map(data, ~{median(.x$lddt)}) %>% unlist(),
         mean_lddt = map(data, ~{mean(.x$lddt)}) %>% unlist(),
         p20_lddt = map(data, ~{quantile(.x$lddt, 0.2)}) %>% unlist(),
         perc_dis = map(data, ~{sum(.x$lddt < 0.6)/length(.x$lddt)}) %>% unlist(),) %>%
  select(gene_name, median_lddt, p20_lddt, mean_lddt, perc_dis) ->
  median_lddt_per_gene


# AlphaFold and binding plots ---------------------------------------------------

eclip_kmer_percentile_ratios_cds %>%
  filter(!(dataset %in% sus_eclip)) %>%
  group_by(gene_name) %>%
  summarise(strong = mean(strong), ratio = mean(ratio)) %>%
  ungroup() %>%
  inner_join(median_lddt_per_gene) %>% 
  mutate(ratio_class = cut(strong, breaks = c(0, quantile(strong, c(.9, 1))), include_lowest = T,
                           labels = c("bottom 90%", "top 10%"))) ->
  eclip_multivalency_and_protein_order

wilcox.test(perc_dis ~ ratio_clas, eclip_multivalency_and_protein_order)


lm(perc_dis ~ strong, eclip_multivalency_and_protein_order) %>%
  summary()

eclip_multivalency_and_protein_order %>% 
  ggplot(aes(x = ratio_class, y = perc_dis)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "tendency to bind\nto multivalent kmers",
       y = "proportion protein disordered (alphafold)")





