rm(list = ls())
setwd("/nemo/lab/ulej/home/users/farawar/")

# Packages ----------------------------------------------------------------

list.of.packages <- c("stringi", "tidyverse", "parallel", "tictoc", "patchwork", "viridis", "BiocManager", "factoextra", "broom")
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

generate_ratio_tables <- function(xl_path, gtm){
  # xl_path <- clip_files[18]
  # gtm <- genome_transcriptome_map
  
  tic(xl_path)
  
  # Get the transcriptome crosslinks.
  txl <- convert_to_transcriptome(xl_path, gtm)
  
  allseqlevels <- unique(c(seqlevels(txl), seqlevels(mv_region_gr)))
  
  seqlevels(txl) <- allseqlevels
  seqlevels(mv_region_gr) <- allseqlevels
  
  # Get the density of crosslinks in each multivalent region.
  mv_xl_counts <- find_overlaps(mv_region_gr, txl) %>%
    as.data.frame() %>%
    full_join(mv_region_gr %>%
                as.data.frame(), 
              by = c("seqnames", "start", "end", "region_type", "peak_identifier", "width", "strand",
                     "cluster", "cluster_name", "cluster_colour", "gene_name", 
                     "gene_id", "cds_start", "cds_length", "tx_length", "cds_end")) %>%
    dplyr::rename(transcript_id = seqnames) %>%
    replace_na(list(count = 0)) %>%
    group_by(peak_identifier, gene_name, transcript_id, region_type, cluster, cluster_name, cluster_colour, width) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(mv_xl_per_nt = count/width)
  
  # This makes everything much simpler to count.
  uncounted_txl <- txl %>% as.data.frame() %>% uncount(count) %>% makeGRangesFromDataFrame()
  
  # Get the density of crosslinks 
  tx_xl_counts <- transcript_details %>%
    mutate(utr3_length = tx_length - cds_end) %>%
    select(transcript_id, gene_name, utr5_length = cds_start, cds_length, utr3_length, total_length = tx_length) %>%
    mutate(utr5_counts = count_overlaps(all_utr5_gr, uncounted_txl),
           cds_counts = count_overlaps(all_cds_gr, uncounted_txl),
           utr3_counts = count_overlaps(all_utr3_gr, uncounted_txl),
           total_counts = count_overlaps(all_tx_gr, uncounted_txl)) %>%
    mutate(utr5_xl_per_nt = utr5_counts / utr5_length,
           cds_xl_per_nt = cds_counts / cds_length,
           utr3_xl_per_nt = utr3_counts / utr3_length,
           total_xl_per_nt = total_counts / total_length)
  
  protein_region_counts <- tibble(region_type = c("UTR5", "CDS", "UTR3", "Total"),
                                  count = c(sum(tx_xl_counts$utr5_counts),
                                            sum(tx_xl_counts$cds_counts),
                                            sum(tx_xl_counts$utr3_counts),
                                            sum(tx_xl_counts$total_counts)))
  
  # Merge on the details about tx counts
  # Should I actually be subtracting the mv counts?
  # Problem - often UTR annotation is wrong. I am taking the longest possible UTR annotation. So 3UTR variance will be especially high.
  
  # How best to actually address this? Negative binomial regression of xl densities? How to weight by sample size (number of xlinks) and not bias towards long tx?
  ratio_table <- mv_xl_counts %>%
    left_join(tx_xl_counts, by = c("gene_name", "transcript_id")) %>%
    # filter(total_xl_per_nt > quantile(tx_xl_counts$total_xl_per_nt, 0.5)) %>% # Only look at transcripts with high count density
    mutate(ratio = case_when(region_type == "UTR3" ~ mv_xl_per_nt / utr3_xl_per_nt,
                             region_type == "UTR5" ~ mv_xl_per_nt / utr5_xl_per_nt,
                             region_type == "CDS" ~ mv_xl_per_nt / cds_xl_per_nt),
           boring_ratio = mv_xl_per_nt / total_xl_per_nt)
  
  # group_by(region_type, cluster, cluster_colour, cluster_name) %>%
  # summarise(median_ratio = median(ratio), 
  #           number_in_cluster = dplyr::n(),
  #           .groups = "drop")
  
  toc()
  
  return(list("ratios" = ratio_table,
              "counts" = protein_region_counts))
}

# Parameters --------------------------------------------------------------

k_len <- 5




# Get multivalent regions -------------------------------------------------

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", 
                               col_types = cols())

cds <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")

mv_region_df <- cds %>%
  bind_rows(.id = "region_type") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "region"), sep = ":|-|@", remove = F) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  dplyr::select(-c(V1, V2)) %>%
  mutate(end = as.numeric(end) + (k_len - 1),
         start = as.numeric(start))

mv_region_gr <- mv_region_df %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id",
                           keep.extra.columns = T)

mv_region_df %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) ->
  cluster_rename_df

# Get transcript regions --------------------------------------------------

all_utr5_gr <- transcript_details %>%
  mutate(start = 1) %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id", 
                           end.field = "cds_start")

all_cds_gr <- transcript_details %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id", 
                           start.field = "cds_start", 
                           end.field = "cds_end")

all_utr3_gr <- transcript_details %>%
  mutate(cds_end = cds_end + 1) %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id", 
                           start.field = "cds_end", 
                           end.field = "tx_length")

all_tx_gr <- transcript_details %>%
  mutate(start = 1) %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id", 
                           end.field = "tx_length")


# Get crosslink locations -------------------------------------------------

eclip_folder <- "/camp/lab/ulej/home/users/kuretk/data/merge_eclip_replicates_20210906/merged_files/"
good_eclip <- read_lines("GASR/germs/data/clip/multivalent_eclip_datasets_cds.txt")



eclip_files <- list.files(eclip_folder, pattern = ".gz", full.names = T) %>% str_subset(stri_paste(good_eclip, collapse = "|"))
eclip_names <- list.files(eclip_folder, pattern = ".gz", full.names = F) %>% word(1, sep = "-merged") %>% str_subset(stri_paste(good_eclip, collapse = "|"))


iclip_files <- c("GASR/tmap/HNRNPH1_Braun/results/xlinks_merged_genome/hnrnph1_wt.xlinks.bed",
                 "GASR/tmap/TRA2B_Best_et_al/results/xlinks_merged_genome/TRA2B.xlinks.bed",
                 "GASR/tmap/SR_Krchnakova_Hela/results/xlinks_merged_genome/srsf2.xlinks.bed",
                 "GASR/tmap/SR_Krchnakova_Hela/results/xlinks_merged_genome/srsf6.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l_clip.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l2_clip.xlinks.bed",
                 "GASR/tmap/LUC7L123_Daniels/results/xlinks_merged_genome/luc7l3_clip.xlinks.bed")

iclip_names <- c("HEK293-HNRNPH1",
                 "MDAMB231-TRA2B",
                 "Hela-SRSF2",
                 "Hela-SRSF6",
                 "K562-LUC7L",
                 "K562-LUC7L2",
                 "K562-LUC7L3")

clip_files <- c(eclip_files, iclip_files)
clip_names <- c(eclip_names, iclip_names)

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

# Calculate ratios ------------------------------------------------------------

tic("Total elapsed time")
ratios_and_counts <- lapply(clip_files, function(x) generate_ratio_tables(x, genome_transcriptome_map))
toc()

all_ratios <- map(ratios_and_counts, ~{ .x %>% pluck(1) }) %>%
  set_names(clip_names) %>%
  bind_rows(.id = "sample")

write_tsv(all_ratios, file = "GASR/germs/data/clip/crosslink_ratios_for_clip_in_cds_regions.tsv.gz")

all_counts <- map(ratios_and_counts, ~{ .x %>% pluck(2) }) %>%
  set_names(clip_names) %>%
  bind_rows(.id = "sample")

write_tsv(all_counts, file = "GASR/germs/data/clip/all_clip_sample_counts.txt.gz")

# Load ratios -------------------------------------------------------------
 
all_ratios <- read_tsv(file = "GASR/germs/data/clip/crosslink_ratios_for_clip_in_cds_regions.tsv.gz")
all_counts <- read_tsv(file = "GASR/germs/data/clip/all_clip_sample_counts.txt.gz")

# Cluster samples and clusters for CDS ---------------------------------------------------------

# Manual filter (samples that don't have much preference, or don't make sense based on literature about their binding motifs)
manual_filter <- c("HepG2-SND1", "K562-LIN28B", "K562-SAFB", "HepG2-IGF2BP3",
  "HepG2-UPF1", "K562-EIF4G2", "HepG2-XRN2", "K562-HNRNPC",
  "HepG2-NOLC1", "K562-PUM1", "K562-SLTM", "HepG2-SLTM",
  "HepG2-RBFOX2", "HepG2-NKRF")

# Only top 1/3 of samples by CDS counts
whitelist_samples <- all_counts %>% 
  filter(region_type == "CDS") %>%
  filter(count > quantile(count, 0)) %>%
  dplyr::select(sample) %>%
  unlist(use.names = F)

transcript_whitelist <- all_ratios %>% dplyr::select(transcript_id, sample, total_counts) %>%
  group_by(transcript_id) %>%
  summarise(all_counts = sum(total_counts)) %>%
  filter(all_counts > quantile(all_counts, 0.5)) %>%
  dplyr::select(transcript_id) %>%
  unlist(use.names = F)

all_ratios_cds <- all_ratios %>%
  mutate(boring_ratio = mv_xl_per_nt / cds_xl_per_nt) %>%
  dplyr::select(-ratio) %>%
  filter(cluster > 0,
         !(sample %in% manual_filter),
         transcript_id %in% transcript_whitelist,
         is.finite(boring_ratio)) %>% 
  drop_na() %>%
  select(sample, cluster, cluster_name, cluster_colour, boring_ratio) %>%
  group_by(sample, cluster, cluster_name, cluster_colour) %>%
  filter(dplyr::n() > 20) %>%
  summarise(log2_ratio = log2(mean(boring_ratio))) %>%
  ungroup()

hclust_samples <- all_ratios_cds %>% 
  select(sample, cluster_name, log2_ratio) %>%
  pivot_wider(names_from = cluster_name, values_from = log2_ratio) %>%
  column_to_rownames("sample") %>%
  replace(is.na(.), 0) %>%
  # drop_na() %>%
  eclust("hclust")

fviz_dend(hclust_samples)

hclust_clusters <- all_ratios_cds %>% 
  select(sample, cluster_name, log2_ratio) %>%
  pivot_wider(names_from = sample, values_from = log2_ratio) %>%
  column_to_rownames("cluster_name") %>%
  replace(is.na(.), 0) %>%
  # drop_na() %>%
  eclust("hclust", k.max = 4)

fviz_dend(hclust_clusters)

lapply(hclust_clusters$labels[hclust_clusters$order], function(x) {
  cluster_rename_df$cluster_name_short[which(cluster_rename_df$cluster_name == x)]
}) %>% unlist() ->
  cluster_name_short_order

# hclust of the clusters kinda sucks, prefer to just use the order of the multivalency snail from 
# my UMAP

cluster_rename_df[c(3,2,8,4,5,7,9),] ->
  cluster_rename_df_ordered

all_ratios_cds_ordered <- all_ratios_cds %>%
  mutate(sample = sample %>% fct_relevel(hclust_samples$labels[hclust_samples$order])) %>%
  left_join(cluster_rename_df) %>%
  mutate(cluster_name_short = cluster_name_short %>%
           fct_relevel(cluster_rename_df_ordered$cluster_name_short)) %>%
  arrange(sample, cluster_name_short)

first_half_samples <- all_ratios_cds_ordered %>%
  arrange(sample, cluster_name_short) %>%
  dplyr::select(sample) %>%
  distinct() %>%
  dplyr::slice(1:round(dplyr::n()/2)) %>%
  unlist(use.names = F)

heatmap_plot_wrapped <- ggplot(all_ratios_cds_ordered %>%
         arrange(sample, cluster_name) %>%
         mutate(first_half = (sample %in% first_half_samples)),
       aes(x = cluster_name_short, y = sample, fill = log2_ratio)) + 
  geom_tile(colour = "black", size = 0.3) +
  facet_wrap(. ~ first_half, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  # scale_fill_viridis() +
  theme_classic() +
  labs(fill = "log2 fold enrichment") +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
        legend.title.align = 0.5,
        legend.title = element_text(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(filename = "GASR/germs/plots/clip_crosslinks_in_germs_clusters_heatmap_wrapped.pdf", heatmap_plot_wrapped,
       device = "pdf", units = "in", width = 8, height = 5)

heatmap_plot_long <- ggplot(all_ratios_cds_ordered %>%
                                 arrange(sample, cluster_name),
                               aes(x = cluster_name_short, y = sample, fill = log2_ratio)) + 
  geom_tile(colour = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", 
                                              ticks.linewidth = 0.5, 
                                              frame.linewidth = 0.5)) +
  theme_classic() +
  labs(fill = "log2 fold enrichment") +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
        legend.title.align = 1,
        legend.title = element_text(),
        legend.position = "top",
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(filename = "GASR/germs/plots/clip_crosslinks_in_germs_clusters_heatmap_long.pdf", heatmap_plot_long,
       device = "pdf", units = "in", width = 3.5, height = 8)
