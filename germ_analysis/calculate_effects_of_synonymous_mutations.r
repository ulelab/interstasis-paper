list.of.packages <- c("germs", "Biostrings", "tidyverse", "parallel", "tictoc",
                      "patchwork", "dbscan", "umap", "viridis") 

for(i in list.of.packages){
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}

rm(list = ls())
setwd("/camp/lab/ulej/home/users/farawar/")

# Functions ---------------------------------------------------------------

calculate_mutated_multivalency <- function(index) {
  tic(paste0("Index: ", index))
  # index <- which(transcript_details$gene_name == "SON")
  # print(index)
  
  sequence <- sequences[index]
  split_seq <- sequence %>% str_split("") %>% unlist()
  cds <- str_sub(sequence, transcript_details$cds_start[index], (transcript_details$cds_end[index]-3))
  start_boundary = window_size / 2
  end_boundary = transcript_details$tx_length[index] - (window_size / 2)
  
  aa_sequence <- cds %>% 
    DNAStringSet() %>%
    translate() %>%
    as.character()
  
  tibble(codon = strsplit(cds, "(?<=.{3})", perl = TRUE)[[1]]) %>%
    mutate(codon_position = row_number()) %>%
    inner_join(codon_table, relationship = "many-to-many", by = "codon") %>%
    mutate(transcript_position = (3 *(codon_position - 1)) + position_in_codon + (transcript_details$cds_start[index] - 1)) %>%
    filter(transcript_position > start_boundary,
           transcript_position < end_boundary) ->
    substitution_positions
  
  # Add in the native codons.
  bind_rows(list(substitution_positions %>% mutate(alternative_codon = codon) %>%
    mutate(substitute_nt = pmap(list(position_in_codon, codon), function(a, b) { str_split(b, "") %>% unlist() %>% pluck(a) }) %>%
             unlist()) %>%
    distinct(),
    substitution_positions)) %>%
      arrange(codon_position) ->
    substitution_positions_native
  
  substitution_positions_native %>%
    mutate(flanking_sequence = paste0(str_sub(sequence, transcript_position - flanking_sequence_length, transcript_position - 1), 
                                      substitute_nt,
                                      str_sub(sequence, transcript_position + 1 , transcript_position + flanking_sequence_length))) %>%
    mutate(germ_vector = germs::list_kmer_multivalencies(flanking_sequence, k_len = k_length, window_size = window_size, hamming_distances = hdm, positional_distances = pdv)) %>%
    mutate(max_multivalency = map(germ_vector, ~{ .x[(flanking_sequence_length - 3) : (flanking_sequence_length + 1)] %>% max() }) %>%
             unlist()) %>% 
    select(transcript_position, codon_position, aa, codon, alternative_codon, position_in_codon, substitute_nt, max_multivalency) ->
    substitutions_mv
  
  substitutions_mv %>%
    mutate(native = case_when(codon == alternative_codon ~ "native",
                              T ~ "mutated")) %>%
    group_by(transcript_position, codon_position, aa, codon, position_in_codon, native) %>%
    summarise(average_max_mv = mean(max_multivalency)) %>%
    ungroup() %>%
    pivot_wider(names_from = "native", values_from = "average_max_mv") %>%
    mutate(diff = native - mutated,
           ratio = native / mutated) ->
    diffy
  toc()
  
  return(diffy)
}

# Parameters --------------------------------------------------------------

k_length <- 5
window_size <- 123
smoothing_size <- 123
lambda <- 1
scaling_function <- function(x) exp(lambda * x)
threshold_percentile <- 0.98

flanking_sequence_length = ((window_size - 1)/2) + (k_length - 1)

# Fetch transcripts -------------------------------------------------------

sequences <- readDNAStringSet("genomes/hs/fasta/longest_gencode29.fa") %>% as.character()

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_length - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

sequences <- sequences[match(transcript_details$transcript_id, names(sequences))] #Reorder sequences to match transcript detail table

# Make HDM and PDV --------------------------------------------------------

hdm <- create_hamming_distance_matrix(k_length, lambda = lambda, scale_fun = scaling_function)
pdv <- create_positional_distance_vector(window_size, k_length)

# Create codon and aa tables ----------------------------------------------

aa_df <- tibble(aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

nts <- c("A", "C", "G", "T")

codon_matrix <- do.call(expand.grid, rep(list(nts), 3))

codon_table_basic <- tibble(codon = do.call(paste0, c(codon_matrix))) %>%
  mutate(aa = as.character(translate(DNAStringSet(codon), no.init.codon = T))) %>%
  group_by(aa) %>%
  # mutate(codon_number = c(1:length(aa))) %>%
  ungroup() %>%
  arrange(aa) 
   
codon_table <- codon_table_basic %>%
  mutate(options = lapply(codon_table_basic$codon, function(x) {
    # x = "ATG"
    # print(x)
    aa = codon_table_basic$aa[codon_table_basic$codon == x]
    split_codon = str_split(x, "") %>% unlist()
    codons = codon_table_basic$codon[codon_table_basic$aa == aa]
    
    if(length(codons) == 1) {
      return(NULL)
    }
    
    codon_distances = lapply(codons, function(y) {
      # y = "AGG"
      y = str_split(y, "") %>% unlist()
      match_vector = c(y[1] != split_codon[1],
                       y[2] != split_codon[2],
                       y[3] != split_codon[3])
      as.numeric(match_vector)
    })
    
    tibble(alternative_codon = codons,
           match_pos = codon_distances) %>%
      mutate(codon_distance = map(match_pos, ~{ sum(.x) }) %>%
               unlist()) %>%
      filter(codon_distance == 1) %>%
      mutate(position_in_codon = map(match_pos, ~{ which(.x == 1) }) %>% unlist()) %>%
      mutate(substitute_nt = pmap(list(position_in_codon, alternative_codon), function(a, b) { str_split(b, "") %>% unlist() %>% pluck(a) }) %>%
               unlist()) %>%
      select(-c(codon_distance, match_pos)) %>%
    return()
  })) %>%
  unnest(options)


# Calculate effets of mutations ----------------------------------------------------------------

library(tictoc)
options(dplyr.summarise.inform = FALSE)

tic("calcluate mv")
mutated_multivalency <- lapply(c(1:length(sequences)),
                               function(x) calculate_mutated_multivalency(x)) %>%
  set_names(names(sequences)) %>%
  bind_rows(.id = "transcript_id") %>%
  left_join(transcript_details[c("transcript_id", "gene_name")])
toc()

write_tsv(mutated_multivalency, "GASR/germs/data/mutated_multivalency_correct.tsv.gz")
