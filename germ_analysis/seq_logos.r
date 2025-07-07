rm(list = ls())

library(tidyverse)
library(Biostrings)
# library(seqLogo)
# library(universalmotif)
library(ggseqlogo)
library(patchwork)

# load relevant data ------------------------------------------------------

germ_scores <- read_tsv("GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz") %>%
  group_by(transcript_id) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  filter(smoothed_kmer_multivalency > 0)

# For normalising the multivalency
minkm <- min(germ_scores$kmer_multivalency)
medkm <- median(germ_scores$kmer_multivalency)

# Normalise the multivalency
germ_scores <- germ_scores %>%
  mutate(kmer_multivalency = 
           (kmer_multivalency - minkm)/
           (medskm - minkm)) %>%
  mutate(kmer_percentile = rank(kmer_multivalency)/length(kmer_multivalency))

tx_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt")
nts <- Biostrings::readDNAStringSet("GASR/lists/longest_gencode29.fa")


# Codon usage -------------------------------------------------------------


aa_df <- tibble(aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

base_nts <- c("A", "C", "G", "T")

codon_matrix <- do.call(expand.grid, rep(list(base_nts), 3))

codon_table <- tibble(codon = do.call(paste0, c(codon_matrix))) %>%
  mutate(aa = as.character(translate(DNAStringSet(codon), no.init.codon = T))) %>%
  group_by(aa) %>%
  mutate(codon_number = c(1:length(aa))) %>%
  ungroup() %>%
  arrange(aa)

tx_details %>% 
  mutate(seq = nts[match(tx_details$transcript_id, names(nts))] %>% as.character()) %>%
  mutate(cds_seq = pmap(list(cds_start, cds_end, seq), 
                        function(cds_start, cds_end, transcript_seq){
                          str_sub(transcript_seq, cds_start, cds_end)
                        }) %>% 
           unlist()) %>%
  mutate(codons = purrr::map(cds_seq, function(cds_seq){
    substring(cds_seq, seq(1, nchar(cds_seq), 3), seq(3, nchar(cds_seq), 3))
  })) ->
  seq_df

transcriptome_codon_usage <- seq_df$codons %>% 
  unlist() %>%
  table() %>%
  as_tibble() %>%
  set_names(c("codon", "count")) %>%
  left_join(codon_table, by = "codon") %>%
  arrange(aa, codon) %>%
  group_by(aa) %>%
  mutate(proportion_usage_transcriptome = count / sum(count)) %>%
  dplyr::select(-c(count, codon_number))

# Load mutation shit ------------------------------------------------------

mutated_multivalency <- read_tsv("GASR/germs/data/mutated_multivalency_correct.tsv.gz")

mutated_multivalency %>%
  # There are a few NA values around because they are too close to the start/end of the transcript - not sure how these got through.
  drop_na() %>%
  mutate(native = (native - minkm)/
           (medkm - minkm),
         mutated = (mutated - minkm)/
           (medkm - minkm),
         diff = native - mutated,
         ratio = native / mutated) ->
  mutated_multivalency

bullshit <- mutated_multivalency %>% group_by(transcript_id) %>% summarise(bullshit = "*" %in% aa) %>% filter(bullshit)

transcript_conservation <- read_tsv("genomes/conservation/conservation_in_transcriptome_coordinates.tsv.gz")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_length - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

mutated_multivalency_conservation <- mutated_multivalency %>%
  dplyr::rename("position" = "transcript_position") %>%
  mutate(codon2_position = case_when(position_in_codon == 1 ~ position + 1,
                                     position_in_codon == 3 ~ position - 1)) %>%
  left_join(transcript_conservation, by = c("transcript_id", "position")) %>%
  left_join(transcript_conservation %>% dplyr::rename(codon2_position = position, 
                                               codon2_phylop_17_primates = phylop_17_primates,
                                               codon2_phylop_470_mammals = phylop_470_mammals,
                                               codon2_phylop_100_vertebrates = phylop_100_vertebrates), 
            by = c("transcript_id", "codon2_position")) %>%
  inner_join(tx_details) %>%
  # Add information about local multivalency
  filter(!(transcript_id %in% bullshit$transcript_id)) %>%
  mutate(codon2_diff_17p = phylop_17_primates - codon2_phylop_17_primates,
         codon2_diff_100v = phylop_100_vertebrates - codon2_phylop_100_vertebrates,
         codon2_diff_470m = phylop_470_mammals - codon2_phylop_470_mammals,) %>%
  group_by(codon, position_in_codon) %>%
  # Scale by subtraction because phylop is log scaled.
  mutate(scaled_17p = phylop_17_primates - mean(phylop_17_primates),
         scaled_100v = phylop_100_vertebrates - mean(phylop_100_vertebrates),
         scaled_470m = phylop_470_mammals - mean(phylop_470_mammals),
         codon2_scaled_17p = codon2_diff_17p - mean(codon2_diff_17p),
         codon2_scaled_100v = codon2_diff_100v - mean(codon2_diff_100v),
         codon2_scaled_470m = codon2_diff_470m - mean(codon2_diff_470m),
  ) %>%
  ungroup() 


mutated_multivalency_conservation %>%
  group_by(transcript_id, gene_name) %>%
  summarise(mean_ratio = mean(ratio),
            mean_s100v = mean(scaled_100v),
            mean_c2s100v = mean(codon2_scaled_100v)) %>% view()


# conservation plotter function -------------------------------------------------------------------


conservation_plotter <- function(gene_name, begin, end, only_mutable){
  
  # gene_name = "RBM39"
  # begin = 150
  # end = 350
  
  fav_txid <- tx_details$transcript_id[tx_details$gene_name == gene_name]
  fav_roi_start <- begin
  fav_roi_end <- end
  
  
  # position >= 800, position <= 1000
  
  fav_seq <- nts[[fav_txid]]
  fav_cons <- transcript_conservation %>% 
    filter(transcript_id == fav_txid)
  fav_mutable_pos = mutated_multivalency_conservation %>% 
    filter(transcript_id == fav_txid) %>%
    dplyr::select(position)
  
  
  tibble(position = c(1:nchar(fav_seq)),
         nt = str_split(fav_seq %>% as.character(), "") %>% unlist(),
         conservation = fav_cons$phylop_100_vertebrates) %>%
    mutate(mutable = position %in% fav_mutable_pos$position) %>%
    filter(position >= fav_roi_start, position <= fav_roi_end) ->
    fav_cons_per_pos
  
  if(only_mutable) {
    fav_cons_per_pos %>%
      mutate(conservation = case_when(mutable == T ~ conservation, mutable == F ~ 0)) ->
      fav_cons_per_pos
  }
  
  
  fav_cons_per_pos %>%
    dplyr::select(-mutable) %>%
    pivot_wider(names_from = "position",
                values_from = "conservation") %>%
    replace(is.na(.), 0) %>%
    arrange(nt) %>%
    column_to_rownames("nt") %>%
    as.matrix() ->
    fav_matrix
  
  fav_cons_per_pos %>%
    mutate(conservation = 1) %>%
    dplyr::select(-mutable) %>%
    pivot_wider(names_from = "position",
                values_from = "conservation") %>%
    replace(is.na(.), 0) %>%
    arrange(nt) %>%
    column_to_rownames("nt") %>%
    as.matrix() ->
    fav_matrix_all_1
  
  
  ggseqlogo(fav_matrix, method='custom', seq_type='dna') +
    theme_classic() +
    coord_cartesian(xlim = c(1, (fav_roi_end - fav_roi_start - 1))) +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          axis.title.x = element_blank()) +
    labs(y = "conservation\n100 vert.") ->
    fav_motif
  
  ggseqlogo(fav_matrix_all_1, method='custom', seq_type='dna') +
    theme_classic() +
    coord_cartesian(xlim = c(1, (fav_roi_end - fav_roi_start - 1))) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) ->
    fav_motif_all_1
  
  
  # Random codons -----------------------------------------------------------
  
  fav_start <- tx_details$cds_start[tx_details$transcript_id == fav_txid]
  fav_end <- tx_details$cds_end[tx_details$transcript_id == fav_txid]
  
  fav_cds <- fav_seq[fav_start:fav_end]
  fav_prot <- translate(fav_cds)
  
  set.seed(12345)
  
  tibble(aa = fav_prot %>% as.character() %>% str_split("") %>% unlist(),
         codon = substring(fav_cds %>% as.character(), seq(1, nchar(fav_cds), 3), seq(3, nchar(fav_cds), 3))) %>%
    mutate(new_codon = map(codon, ~{
      # .x = "TGG"
      my_aa = transcriptome_codon_usage$aa[transcriptome_codon_usage$codon == .x]
      
      if(sum(transcriptome_codon_usage$aa == my_aa) < 2) {
        return(.x)
      }
      
      relevant_aa <- transcriptome_codon_usage %>% filter(aa == my_aa, codon != .x)
      sample(relevant_aa$codon, 1, prob = relevant_aa$proportion_usage_transcriptome)
      
    }) %>% unlist()) %>%
    dplyr::select(new_codon) %>%
    unlist() %>%
    paste(collapse = "") ->
    fav_cds_newcodons
  
  tibble(nt = str_split(fav_cds_newcodons %>% as.character(), "") %>% unlist(),
         conservation = 1) %>%
    mutate(position = row_number() + (fav_start - 1)) %>%
    filter(position >= fav_roi_start, position <= fav_roi_end) ->
    fav_scramble_df
  
  fav_scramble_df %>%
    pivot_wider(names_from = "position",
                values_from = "conservation") %>%
    replace(is.na(.), 0) %>%
    arrange(nt) %>%
    column_to_rownames("nt") %>%
    as.matrix() ->
    fav_matrix_newcodons
  
  ggseqlogo(fav_matrix_newcodons, method='custom', seq_type='dna') +
    theme(axis.text.x = element_blank()) +
    theme_classic() +
    coord_cartesian(xlim = c(1, (fav_roi_end - fav_roi_start - 1))) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) ->
    fav_newcodons_motif
  
  mutated_multivalency %>%
    filter(transcript_id == fav_txid, transcript_position >= fav_roi_start, transcript_position <= fav_roi_end - 2) %>%
    dplyr::rename("position" = "transcript_position") %>% 
    ggplot(aes(x = position, y = ratio)) +
    scale_y_log10() +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    
    geom_segment(aes(xend = position, y = 1, yend = ratio)) +
    geom_point(shape = 21, size = 1, fill = "white") +
    coord_cartesian(xlim = c(fav_roi_start,fav_roi_end - 2)) +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          legend.position = "none") ->
    fav_ratio_plot
  
  (fav_motif_all_1 / fav_motif / fav_ratio_plot / fav_newcodons_motif) + plot_layout(heights = c(1, 4, 2, 1)) ->
    fav_combined_conservation_logos
  
  return(fav_combined_conservation_logos)
}

conservation_plotter("LUC7L3", 895, 975, T) ->
  LUC7L3_mut

conservation_plotter("LUC7L3", 895, 975, F) ->
  LUC7L3_all

conservation_plotter("UPF3B", 700, 780, T) ->
  UPF3B_mut

conservation_plotter("UPF3B", 700, 780, F) ->
  UPF3B_all

conservation_plotter("CCDC61", 650, 750, T) ->
  CCDC61_mut

conservation_plotter("CCDC61", 650, 750, F) ->
  CCDC61_all

conservation_plotter("TMEM200B", 255, 335, T) ->
  TMEM200B_mut

conservation_plotter("TMEM200B", 255, 335, F) ->
  TMEM200B_all

ggsave("GASR/germs/plots/example_genes/conservation/LUC7L3_mut.pdf", LUC7L3_mut,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/LUC7L3_all.pdf", LUC7L3_all,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/UPF3B_mut.pdf", UPF3B_mut,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/UPF3B_all.pdf", UPF3B_all,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/CCDC61_mut.pdf", CCDC61_mut,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/CCDC61_all.pdf", CCDC61_all,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/TMEM200B_mut.pdf", TMEM200B_mut,
       width = 8, height = 3)

ggsave("GASR/germs/plots/example_genes/conservation/TMEM200B_all.pdf", TMEM200B_all,
       width = 8, height = 3)