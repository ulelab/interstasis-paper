rm(list = ls())

library(tidyverse)
library(germs)
library(parallel)
library(Biostrings)
library(broom)

# Parameters --------------------------------------------------------------

k_length <- 5
window_size <- 123
smoothing_size <- 123
lambda <- 1
scaling_function <- function(x) exp(lambda * x)
threshold_percentile <- 0.98

# Fetch transcripts -------------------------------------------------------

sequences <- readDNAStringSet("genomes/hs/fasta/longest_gencode29.fa") %>% as.character()

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_length - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

sequences <- sequences[match(transcript_details$transcript_id, names(sequences))] #Reorder sequences to match transcript detail table

# Create codon and aa tables ----------------------------------------------

aa_df <- tibble(aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

nts <- c("A", "C", "G", "T")

codon_matrix <- do.call(expand.grid, rep(list(nts), 3))

codon_table <- tibble(codon = do.call(paste0, c(codon_matrix))) %>%
  mutate(aa = as.character(translate(DNAStringSet(codon), no.init.codon = T))) %>%
  group_by(aa) %>%
  mutate(codon_number = c(1:length(aa))) %>%
  ungroup() %>%
  arrange(aa)

# Make HDM and PDV --------------------------------------------------------

hdm <- create_hamming_distance_matrix(k_length, lambda = lambda, scale_fun = scaling_function)
pdv <- create_positional_distance_vector(window_size, k_length)

# Load previous GeRMS --------------------------------------------------------------------

germs <- read_tsv(file = "GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz")

germs <- germs %>%
  group_by(transcript_id) %>%
  mutate(position = c(1:dplyr::n())) %>%
  ungroup()

germs_with_cds_clusters <- germs %>%
  group_by(transcript_id) %>%
  mutate(has_cluster = sum(peak_location == "CDS") > 0) %>%
  ungroup() %>%
  filter(has_cluster) %>%
  dplyr::select(-has_cluster)

cds_clusters <- read_tsv(file = "GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")

median_smoothed_germs <- germs %>%
  filter(smoothed_kmer_multivalency > 0) %>%
  left_join(transcript_details) %>%
  group_by(transcript_id) %>%
  filter(position >= cds_start,
         position <= cds_end) %>%
  ungroup() %>%
  summarise(median_mv = median(smoothed_kmer_multivalency)) %>%
  unlist(use.names = F)

# For normalising the multivalency
minkm <- min(germs$kmer_multivalency)
medkm <- median(germs$kmer_multivalency)

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


# Reshuffle and calculate -------------------------------------------------

codon_reshuffler <- function(input_txid) {
  #Reshuffle codons while retaining codon usage. Requires protein sequence (as string)
  #and codon sequence (as vector of strings).
  
  #A dataframe 'codon_table' containing all codons numbered within each amino acid must be defined globally.
  
  set.seed(sample(1:100, 1))
  
  # input_gene_name <- "BOD1L1"
  # input_txid <- transcript_details$transcript_id[transcript_details$gene_name == input_gene_name]
  
  input_cds_start <- transcript_details$cds_start[transcript_details$transcript_id == input_txid]
  input_cds_end <- transcript_details$cds_end[transcript_details$transcript_id == input_txid]
  input_tx_len <- transcript_details$tx_length[transcript_details$transcript_id == input_txid]
  
  input_cds <- sequences[input_txid] %>% substr(input_cds_start, input_cds_end)
  input_5utr <- sequences[input_txid] %>% substr(1, input_cds_start - 1)
  input_3utr <- sequences[input_txid] %>% substr(input_cds_end + 1, input_tx_len)
  
  input_codons <- substring(input_cds, seq(1, nchar(input_cds), 3), seq(3, nchar(input_cds), 3))
  input_aa <- input_cds %>% DNAStringSet() %>% translate() %>% as.character() %>% unname()
  
  return(
    c(input_5utr,
      tibble(aa = str_split(input_aa, "") %>% unlist(),
             codon = input_codons) %>%
        left_join(codon_table, by = c("aa", "codon")) %>%
        group_by(aa) %>%
        mutate(codon_number = codon_number[sample(length(codon_number))]) %>% #Randomly reorder the codon numbers.
        ungroup() %>%
        dplyr::select(-codon) %>%
        left_join(codon_table, by = c("aa", "codon_number")) %>%
        dplyr::select(codon) %>%
        unlist(use.names = F),
      input_3utr)  %>%
      paste0(collapse = "")
  )
}

codon_shuffled_local_multivalency_potential <- function(input_txid, number_of_shuffles){
  # input_txid = "ENST00000334701.11"
  # number_of_shuffles = 5
  
  lapply(c(1:number_of_shuffles), function(x) {
    # x=1
    codon_reshuffler(input_txid) %>% 
      calculate_kmer_multivalencies_df(.,
                                       input_txid,
                                       k_length,
                                       window_size,
                                       smoothing_size,
                                       hdm,
                                       pdv) %>%
      mutate(position = c(1:dplyr::n()))
    }) %>%
    set_names(c(1:number_of_shuffles)) %>%
    bind_rows(.id = "iteration") %>%
    group_by(sequence_name, position) %>%
    summarise(mean_shuffled_kmer_multivalency = mean(kmer_multivalency),
      mean_shuffled_smoothed_kmer_multivalency = mean(smoothed_kmer_multivalency)) %>%
    ungroup() %>%
    arrange(position)
}

library(tictoc)
tic("shuffling and calculating")
germs_shuffled <- mclapply(unique(germs$transcript_id),
                                             function(x) codon_shuffled_local_multivalency_potential(x, 10), mc.cores = 16) %>%
  bind_rows() %>%
  dplyr::rename(transcript_id = sequence_name) %>%
  left_join(germs, by = c("transcript_id", "position")) %>%
  mutate(med_norm =
           (kmer_multivalency - minkm)/
           (medkm - minkm),
         med_norm_smooth =
           (smoothed_kmer_multivalency - minskm)/
           (medskm - minskm),
         med_norm_shuff =
           (mean_shuffled_kmer_multivalency - minkm)/
           (medkm - minkm),
         med_norm_shuff_smooth =
           (mean_shuffled_smoothed_kmer_multivalency - minskm)/
           (medskm - minskm))
toc()

germs_shuffled %>%
mutate(med_norm =
         (kmer_multivalency - minkm)/
         (medkm - minkm),
       med_norm_smooth =
         (smoothed_kmer_multivalency - minskm)/
         (medskm - minskm),
       med_norm_shuff =
         (mean_shuffled_kmer_multivalency - minkm)/
         (medkm - minkm),
       med_norm_shuff_smooth =
         (mean_shuffled_smoothed_kmer_multivalency - minskm)/
         (medskm - minskm)) ->
  germs_shuffled

# write_tsv(germs_shuffled,
#                     file = "GASR/germs/data/smoothed_kmer_multivalency_native_and_shuffled.tsv.gz")

# germs_shuffled <- read_tsv("GASR/germs/data/smoothed_kmer_multivalency_native_and_shuffled.tsv.gz")

# Effect on low complexity domains ----------------------------------------

lcds <- read_tsv("GASR/germs/data/low_entropy_clusters.tsv.gz",
               col_types = cols()) %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end"), sep = ":|-", remove = F, convert = T) %>%
  left_join(transcript_details, by = "transcript_id") %>%
  mutate(start = ((start * 3) - 3) + cds_start, # A region that starts at position 1 should end up with its start being equal to the cds_start.
         end = (end * 3) + cds_start + 2 + (k_length - 1)) # A region that ends on the last amino acid should end up with its end being equal to the cds_end - 3 (no stop codon in entropy) 

left_join(germs_shuffled, 
          lcds[c("transcript_id", "start", "end")],
          relationship = "many-to-many") %>%
  # the join makes multiple copies of genes when they encode 2+ LCDs, but we will get rid of them later with the distinct() function call.
  mutate(is_lcd = case_when((position >= start) & (position <= end) ~ "LCD",
                           T ~ "non-LCD") %>%
           fct_relevel("non-LCD")) %>%
  mutate(id = case_when(is_lcd == "LCD" ~ paste0(transcript_id, ":", start, "-", end),
                        T ~ "no")) %>%
  select(-c(start, end)) %>%
  distinct() ->
  germs_shuffled_lcd

# Need to filer out UTRs which are not of interest

germs_shuffled_lcd %>%
  left_join(transcript_details) %>%
  filter(position > cds_start,
         position < (cds_end - (k_length - 1))) ->
  germs_shuffled_lcd_cds

germs_shuffled_lcd_cds %>%
  group_by(transcript_id, is_lcd, id) %>%
  summarise(native = mean(med_norm_smooth),
            shuffled = mean(med_norm_shuff_smooth)) %>%
  select(-c(transcript_id, id)) %>%
  pivot_longer(cols = c(native, shuffled), names_to = "type", values_to = "germ") ->
  long_lcd_comparison

write_tsv(long_lcd_comparison, "GASR/germs/data/local_shuffling_summary.tsv.gz")

long_lcd_comparison <- read_tsv("GASR/germs/data/local_shuffling_summary.tsv.gz")

ggplot(long_lcd_comparison,
       aes(x = germ, color = is_lcd, linetype = type)) +
  stat_ecdf(n = 1000) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.8) +
  coord_cartesian(xlim = c(0,6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic() +
  labs(x = "mean GeRM in region",
       y = "cumulative density",
       color = "", linetype = "") +
  scale_color_manual(values = c("grey70", "black")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) ->
  lcd_shuffling_germ_plot

long_lcd_comparison %>%
  group_by(is_lcd) %>%
  summarise(ttest = t.test(germ ~ type) %>% tidy()) %>%
  unnest(ttest) ->
  lcd_shuffling_germ_ttests

ggsave("GASR/germs/plots/lcd_shuffling_germ_plot.pdf", lcd_shuffling_germ_plot,
       device = "pdf", units = "in",
       height = 2, width = 4)

write_tsv(lcd_shuffling_germ_ttests, "GASR/germs/plots/lcd_shuffling_germ_ttests.tsv")
write_tsv(long_lcd_comparison, "GASR/germs/plots/lcd_shuffling_germ_data.tsv")

germs_shuffled_lcd_cds %>%
  mutate(ratio = med_norm_smooth / med_norm_shuff_smooth) %>%
  group_by(transcript_id, is_lcd, id) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  ungroup() %>%
  left_join(transcript_details) ->
  lcd_ratio_comparison

# This kinda sucks because most non-LCD regions will have a more stable variance because they can be whole proteins, so the result isn't easily interpretable
# 
# ggplot(lcd_ratio_comparison,
#        aes(x = mean_ratio, color = is_lcd)) +
#   stat_ecdf() +
#   # geom_vline(xintercept = median(germs_shuffled_lcd_cds$med_norm_smooth), linetype = "dashed", alpha = 0.8) +
#   # geom_vline(xintercept = p10_entropy, linetype = "dashed", alpha = 0.8, color = "red3") +
#   # coord_cartesian(xlim = c(0,6)) +
#   # scale_x_continuous(breaks = scales::pretty_breaks(3)) +
#   scale_x_log10() +
#   theme_classic() +
#   labs(x = "mean GeRM in region",
#        y = "cumulative density",
#        color = "", linetype = "") +
#   scale_color_manual(values = c("grey70", "black")) +
#   theme(axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black"),
#         axis.ticks = element_line(color = "black"),
#   ) 

ggplot(lcd_ratio_comparison,
       aes(x = mean_ratio)) +
  stat_ecdf() +
  facet_wrap(. ~ is_lcd) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.8) +
  # geom_vline(xintercept = p10_entropy, linetype = "dashed", alpha = 0.8, color = "red3") +
  # coord_cartesian(xlim = c(0,6)) +
  # scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  scale_x_log10() +
  theme_classic() +
  labs(x = "ratio of native:shuffled GeRM in region",
       y = "cumulative density",
       color = "", linetype = "") +
  scale_color_manual(values = c("grey70", "black")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  )


# Jure's shuffling --------------------------------------------------------

jure_germ_in_lcd <- read_tsv('GASR/germs/plots/jure_shuffle/from_jure/meanGerm_LCD.tsv')
jure_germ_out_lcd <- read_tsv('GASR/germs/plots/jure_shuffle/from_jure/meanGerm_notLCD.tsv')

list(jure_germ_in_lcd, jure_germ_out_lcd) %>%
  bind_rows() %>%
  select(transcript_id = sequence_name, peak_id, is_lcd = data, mean_unshuffled_smooth, mean_shuffled_smooth) %>%
  pivot_longer(cols = c(mean_unshuffled_smooth, mean_shuffled_smooth), names_to = "type", values_to = "germ") %>%
  mutate(germ = 
           (germ - minskm)/
           (medskm - minskm),
         type = case_when(type == "mean_unshuffled_smooth" ~ "native",
                          T ~ "global shuffle") %>%
           fct_relevel("native")) ->
  jure_df

list(long_lcd_comparison,
     jure_df %>% select(-peak_id) %>% filter(type == "global shuffle")) %>%
  bind_rows() %>%
  ungroup() %>%
  mutate(type = type %>% fct_relevel("native", "global shuffle"))->
  merged_df



median(jure_df$germ[jure_df$is_lcd == "non-LCD"])

ggplot(jure_df,
       aes(x = germ, color = is_lcd, linetype = type)) +
  stat_ecdf(n = 1000) +
  geom_vline(xintercept = jure_df %>% 
               filter(is_lcd == "non-LCD", type == "native") %>% 
               summarise(median(germ)) %>%
               unlist(use.names = F), 
             linetype = "dashed", alpha = 0.8) +
  coord_cartesian(xlim = c(0,6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic() +
  labs(x = "mean GeRM in region",
       y = "cumulative density",
       color = "", linetype = "") +
  scale_color_manual(values = c("#479d99", "#b4b4b4")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) ->
  jure_lcd_shuffle_plot

jure_df %>%
  group_by(is_lcd) %>%
  summarise(ttest = t.test(germ ~ type) %>% tidy()) %>%
  unnest(ttest) ->
  jure_lcd_shuffling_germ_ttests

ggsave("GASR/germs/plots/jure_shuffle/jure_lcd_shuffle_plot.pdf", jure_lcd_shuffle_plot,
       device = "pdf", units = "in",
       height = 2, width = 4)

write_tsv(jure_lcd_shuffling_germ_ttests, "GASR/germs/plots/jure_shuffle/jure_lcd_shuffling_germ_ttests.tsv")

write_tsv(jure_df, "GASR/germs/plots/jure_shuffle/jure_lcd_shuffling_germ_table.tsv")

# Older shite -------------------------------------------------------------



cds_clusters_shuffled_vs_native <- germs_with_cds_clusters_shuffled %>%
  filter(peak_location == "CDS",
         smoothed_kmer_multivalency > 0) %>% #Incase the peak is near the start or end of the transcript, where we have 0s for smoothed scores
  mutate(peak_identifier = paste0(transcript_id, ":", peak_start, "-", peak_end, "@", peak_location)) %>%
  group_by(peak_identifier) %>%
  summarise(mean_native_smoothed_multivalency = mean(med_norm_smooth),
            mean_shuffled_smoothed_multivalency = mean(med_norm_shuff_smooth)) %>%
  left_join(cds_clusters) %>%
  mutate(ratio_native_shuffled = mean_shuffled_smoothed_multivalency / mean_native_smoothed_multivalency,
         diff_native_shuffled = mean_shuffled_smoothed_multivalency - mean_native_smoothed_multivalency) %>%
  drop_na()

write_tsv(cds_clusters_shuffled_vs_native, file = "GASR/germs/data/germs_CDS_umap_with_native_and_shuffled_mv.tsv")
cds_clusters_shuffled_vs_native <- read_tsv("GASR/germs/data/germs_CDS_umap_with_native_and_shuffled_mv.tsv")


cds_clusters_shuffled_vs_native %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  arrange(cluster) %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df


cds_clusters_shuffled_vs_native %>%
  left_join(cluster_rename_df) %>%
  arrange(cluster_name_short) %>%
  filter(new_order > 0) ->
  cds_clusters_shuffled_vs_native_named

ggplot(cds_clusters_shuffled_vs_native_named) +
  aes(x = log2(ratio_native_shuffled), colour = cluster_name_short) +
  stat_ecdf() +
  scale_color_manual(values = cluster_rename_df$cluster_colour[-1]) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black")) +
  labs(colour = "",
       x = "log2 (native / shuffled) GeRM score in region",
       y = "cumulative density") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) ->
  shuffled_ratio_ecdf_plot

ggsave("GASR/germs/plots/ecdf_shuffled_ratio_germs.pdf", shuffled_ratio_ecdf_plot,
       device = "pdf", units = "in",
       height = 3, width = 4.5)


# Without clusters --------------------------------------------------------

cds_clusters_shuffled_vs_native %>%
  select(mean_native_smoothed_multivalency, mean_shuffled_smoothed_multivalency) %>%
  pivot_longer(cols = everything(), names_to = "type", values_to = "multivalency") %>%
  mutate(type = type %>%
           str_replace("mean_native_smoothed_multivalency", "native") %>% 
           str_replace("mean_shuffled_smoothed_multivalency", "shuffled")) ->
  simple_comparison

ggplot(simple_comparison,
       aes(x = multivalency, color = type)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic() +
  labs(x = "GeRM score",
       y = "cumulative density",
       color = "") +
  scale_color_manual(values = c("black", "grey70")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) 

ggplot(cds_clusters_shuffled_vs_native,
       aes(x = ratio_native_shuffled)) +
  stat_ecdf() +
  # coord_cartesian(xlim = c(0, 6)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic() +
  labs(x = "ratio of GeRM score\n(shuffled / native)",
       y = "cumulative density",
       color = "") +
  # scale_color_manual(values = c("black", "grey70")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
  ) 


# Stats -------------------------------------------------------------------

cluster_rename_df %>% 
  filter(new_order != 0) %>% 
  select(cluster_name_short) %>% 
  unlist(use.names = F) ->
  testing_clusters

testing_clusters %>%
  lapply(., function(x) {
    germ_shuff <- cds_clusters_shuffled_vs_native_named$ratio_native_shuffled[cds_clusters_shuffled_vs_native_named$cluster_name_short == x]
    
    return(t.test(log2(germ_shuff)) %>% tidy())
  }) %>%
  set_names(testing_clusters) %>%
  bind_rows(.id = "cluster_name") %>%
  mutate(padj = p.adjust(p.value, method = "BH")) ->
  shuff_ttests

write_tsv(shuff_ttests, "GASR/germs/data/stats/shuffling_in_germ_regions.tsv")
