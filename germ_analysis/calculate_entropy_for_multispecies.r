#!/usr/bin/env Rscript
rm(list = ls())

# Pass the filenames as arguments to this script
input_file = commandArgs(trailingOnly=TRUE)

# setwd("~/home/genomes/CDS/")
# input_file <- "Danio_rerio.GRCz11.cds.all.fa"

# Init --------------------------------------------------------------------

library(tidyverse)
library(Biostrings)
library(parallel)

# Params ------------------------------------------------------------------

win_size <- 41
human_98 <- 0.6746539

# Functions ---------------------------------------------------------------

aa_entropy <- function(aa, win_size){
  #Calculate local entropy in a sliding window of {win_size} along an amino acid sequence {aa}.
  
  # aa <- protein_entropy$aa_sequence[4885]
  # win_size <- 41
  
  #Turn the string in to all 20mers.
  short_seqs <- substring(aa, 
                          seq(1, nchar(aa) - win_size, 1), #Exclude stop codon.
                          seq(win_size, nchar(aa) - 1 , 1)) #Exclude stop codon.
  
  #Split the 20mers up, then transpose them so that entropy can be calculated. 
  ent <- as.matrix(short_seqs) %>%
    str_split("", simplify = T) %>%
    t() %>% 
    MolecularEntropy("AA")
  
  return(ent$H) #Return the vector of entropy values
  
}

#Stupid shit for the AA entropy 
AminoAcids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AAbyGroup = c("aliphatic", "cysteine", "acidic", "acidic", "aromatic", "aliphatic", "basic", "aliphatic", "basic", "aliphatic", "aliphatic", "aminic", "proline", "aminic", "basic", "hydroxylated", "hydroxylated", "aliphatic", "aromatic", "aromatic")
AAGroups = c("acidic","aliphatic", "aminic",  "aromatic", "basic","cysteine","hydroxylated", "proline")
small = c(1,2,3,6,12,13,16,17,18)
polar = c(2,3,4,7,9,12,14,15,16,17,19,20)
hydrophobic = c(1, 2,5,6,7,8,9,10,11,17,18,19,20)

#This function is lifted from the old HDMD package, which is no longer available on CRAN.
MolecularEntropy = function(x, type){
  
  # x <- blah[,1]
  type <- "AA"
  
  H = NULL
  
  if(!is.matrix(x)){
    if(is.vector(x)){ 
      x = strsplit(x, "")
      x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
    }
    if(is.list(x)){
      x = sapply(x, strsplit, "")
      x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
    }
  }
  
  if (type == "DNA") {
    DNA = c("A", "C", "G", "T")
    if(!all(x %in% DNA)){print("Warning: Data set contains non-nucleotide elements")}
    
    counts = apply(x, 2, function(z){table(factor(z, levels = DNA))})
    
    freq = apply(counts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
    
    H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0)) })
    H = H/log(4)
    
  }
  
  if (type == "AA") {
    if(!all(x %in% AminoAcids)){print("Warning: Data set contains non-Amino Acid elements")}
    
    counts = apply(x, 2, function(z){table(factor(z, levels = AminoAcids))})
    
    freq = apply(counts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
    
    H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0))} )
    H = H/log(20)
    
  }
  
  if (type == "GroupAA") {
    if(!all(x %in% AminoAcids)){print("Warning: Data set contains non-Amino Acid elements")}
    counts = apply(x, 2, function(z){table(factor(z, levels = AminoAcids))})
    
    grpcounts = apply(counts, 2, function(site){c(site[3]+site[4], site[1]+site[6]+site[8]+site[10]+ site[11]+site[18], site[12]+site[14], site[5]+site[19]+site[20], site[7]+site[9]+site[15], site[2], site[16]+site[17],site[13])})
    rownames(grpcounts)=AAGroups
    
    freq = apply(grpcounts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
    counts = grpcounts
    
    H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0)) })
    H = H/log(8)
    
  }
  
  result <- list(counts = counts, freq=freq, H=H)
  return(result)
}


call_le_regions <- function(ent, threshold){
  
  # ent = multivalency_df$ent[multivalency_df$transcript_id == "ENST00000322723.8"][[1]]
  # theshold <- ent_threshold
  
  ent_pos <- ent <= threshold
  
  #If there are no multivalent positions, return an empty list.
  if(sum(ent_pos) == 0){
    return(list())
  }
  
  ent_rle <- rle(ent_pos)
  
  ent_end <- cumsum(ent_rle$lengths) #End positions for all runs.
  ent_start <- c(1, ent_end[-length(ent_end)] + 1)[ent_rle$values] #Start positions for runs above threshold.
  ent_end <- ent_end[ent_rle$values] #End positions for runs above threshold.
  
  # Often, we end up finding multiple multivalent regions which, when considering the window size, are made up of
  # overlapping regions of the transcript. If these regions have an overlap of at least 2/3 the window size, they
  # are merged.
  
  ent_dist <- (ent_end[-length(ent_end)] + win_size) - ent_start[-1] #Find the distance between clusters. Many
  
  olap <- ent_dist >= win_size/3
  
  ent_start <- ent_start[!c(F, olap)] #Keep the first start of the overlaps
  ent_end <- ent_end[!c(olap, F)] # Keep the last end of the overlaps.
  
  #Make a list of mv starts and ends.
  return(
    lapply(c(1:length(ent_end)), function(x){
      c(ent_start[x], ent_end[x])
    })
  )
}

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

# load --------------------------------------------------------------------

species_name <- input_file %>% word(1, sep = fixed("."))


if(file.exists(paste0("~/home/GASR/germs/data/multiple_species/species_codons_per_transcript/", 
                      species_name, ".codons_per_transcript.tsv.gz"))) {
  quit(save = "no")
}


fa <- readDNAStringSet(input_file)

cds_df <- tibble(fa_head = names(fa),
       cds_len = width(fa),
       cds_seq = as.character(fa)) %>%
  separate(fa_head, 
           into = c("transcript_id", "is_cds","genomic_loc", "gene_id", "gene_biotype", "transcript_biotype", "gene_name"),
           sep = " ",
           extra = "drop") %>%
  dplyr::select(-is_cds) %>%
  filter(gene_biotype == "gene_biotype:protein_coding" &
           transcript_biotype == "transcript_biotype:protein_coding") %>%
  mutate(start_codon = substr(cds_seq, 1, 3),
         stop_codon = substr(cds_seq, nchar(cds_seq) - 2, nchar(cds_seq)),
         frame = (nchar(cds_seq) %% 3) == 0) %>%
  mutate(canonical_start = start_codon == "ATG",
         canonical_stop = stop_codon %in% c("TGA", "TAA", "TAG")) %>%
  mutate(species = species_name)

non_dna <- letters %>% toupper() %>% str_subset(pattern = "A|C|T|G", negate = T) %>% paste(collapse = "|")

longest_cds <- cds_df %>%
  filter(frame == T,
         canonical_start == T,
         canonical_stop == T) %>%
  group_by(gene_id) %>%
  arrange(desc(cds_len), .by_group = T) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(frame, canonical_start, canonical_stop, start_codon, stop_codon)) %>% #Columns contain no information at this point.
  filter(cds_len >= 200) %>% #Filter very short transcripts.
  filter(!grepl("Mito|MT|mitochonrion", x = genomic_loc)) %>% #Mitochondrial translation has strange codons.
  filter(!grepl("dubious", x = gene_name)) %>% # Some genes are just called 'dubious'. We don't want them.
  mutate(has_nondna_bases = str_detect(cds_seq, non_dna)) 

tibble(total_genes = nrow(longest_cds),
       with_nonstandard_base = sum(longest_cds$has_nondna_bases)) %>%
  mutate(propotion_discarded = with_nonstandard_base / total_genes,
         species = species_name) ->
  summary_stats

write_tsv(summary_stats, paste0("~/home/GASR/germs/data/multiple_species/species_summaries/", species_name, "summary.tsv"))

longest_cds %>% #Proteins with weird ass characters in sequence can't be translated. There are not many, luckily.
  filter(!has_nondna_bases) %>%
  mutate(protein_seq = cds_seq %>% DNAStringSet() %>% translate() %>% as.character() %>% unname()) ->
  longest_cds

if(file.exists(paste0("~/home/GASR/germs/data/multiple_species/species_entropy/", species_name, ".entropy.rds.gz"))) {
  
  all_entropy <- read_rds(paste0("~/home/GASR/germs/data/multiple_species/species_entropy/", species_name, ".entropy.rds.gz"))
  
} else {
  
  longest_cds %>%
    mutate(entropy = mclapply(protein_seq, function(x) { aa_entropy(x, win_size) }, mc.cores = 8)) ->
    all_entropy
  
  all_entropy %>%
    # sample_n(1000) %>%
    mutate(le_region = mclapply(entropy, function(x) call_le_regions(x, human_98), mc.cores = 16)) %>%
    mutate(codons = map(cds_seq, ~{ strsplit(.x, "(?<=.{3})", perl = TRUE)[[1]] })) ->
    all_entropy
  
  write_rds(all_entropy, paste0("~/home/GASR/germs/data/multiple_species/species_entropy/", species_name, ".entropy.rds.gz"))
  
}

# Codons in leds ----------------------------------------------------------

# Make a table containing each low entropy domain, and the proportion of the total made up by each amino acid.
aa_in_le_regions <- all_entropy %>%
  unnest(le_region) %>%
  mutate(region = map2(transcript_id, # Make a name for the multivalent region (tx_id:startkmer-endkmer)
                       le_region,
                       ~{ paste0(.x, ":", .y[1], "-", .y[2] + (win_size - 1)) }) %>%
           unlist(),
         aa = map2(protein_seq, le_region,
                   ~{ substring(.x, .y[1], .y[2] + (win_size - 1)) %>%
                       strsplit("") })) %>%
  select(species, transcript_id, le_region, region, aa) %>%
  unnest(aa) %>% unnest(aa) %>% # Double unnest fuckery.
  group_by(species, transcript_id, le_region, region, aa) %>%
  summarise(count = dplyr::n()) %>%
  mutate(prop = count/sum(count)) %>%
  ungroup()

write_tsv(aa_in_le_regions, 
          paste0("~/home/GASR/germs/data/multiple_species/species_aa_in_leds/", species_name, "aa_in_leds.tsv.gz"))

transcriptome_codon_usage <- all_entropy %>% 
  select(species, codon = codons) %>%
  unnest(codon) %>%
  count(species, codon, name = "count") %>%
  left_join(codon_table, by = "codon") %>%
  arrange(species, aa, codon) %>%
  group_by(species, aa) %>%
  mutate(proportion_usage_transcriptome = count / sum(count)) %>%
  dplyr::select(-c(count, codon_number)) %>%
  ungroup()

# Codon usages ------------------------------------------------------------------

codons_in_leds <- all_entropy %>%
  select(-c(genomic_loc, gene_id, transcript_biotype, gene_biotype, cds_len)) %>%
  unnest(le_region) %>%
  mutate(region = map(le_region, ~{ paste0(.x[1], "-", (.x[2] + (win_size - 1) ) ) }) %>% unlist() %>% paste0(transcript_id, ":", .)) %>% 
  mutate(codons_in_domain = pmap(list(le_region, codons),
                                 function(le_region, codons){
                                   codons[le_region[1]:(le_region[2] + (win_size - 1))]
                                 })) %>%
  mutate(codons_frequency = map(codons_in_domain, ~{
    tibble(codon = .x) %>%
      group_by(codon) %>%
      count() %>%
      ungroup() %>%
      right_join(codon_table, by = "codon")
  })) %>% 
  dplyr::select(-c(le_region, codons, codons_in_domain)) %>%
  unnest(codons_frequency) %>%
  replace_na(list(n = 0))

codon_usage_in_leds <- codons_in_leds %>%
  group_by(species, region, aa) %>%
  mutate(codon_usage = n/sum(n)) %>%
  ungroup() %>%
  left_join(transcriptome_codon_usage) %>%
  mutate(normalised_to_transcriptome = codon_usage/proportion_usage_transcriptome)

write_tsv(codon_usage_in_leds %>% select(-c(cds_seq, protein_seq, entropy, n)), 
          paste0("~/home/GASR/germs/data/multiple_species/species_codons_in_leds/", species_name, ".codons_in_leds.tsv.gz"))


# Codons per transcript ---------------------------------------------------

# species_codons_per_transcript

all_entropy %>%
  select(-c(genomic_loc, gene_id, transcript_biotype, gene_biotype, cds_len, le_region)) %>%
  mutate(codons_frequency = map(codons, ~{
    tibble(codon = .x) %>%
      group_by(codon) %>%
      count() %>%
      ungroup() %>%
      right_join(codon_table, by = "codon")
  })) %>% 
  dplyr::select(-c(codons)) %>%
  unnest(codons_frequency) %>%
  replace_na(list(n = 0)) ->
  codons_in_transcript

codon_usage_per_transcript <- codons_in_transcript %>%
  group_by(species, transcript_id, aa) %>%
  mutate(codon_usage = n/sum(n)) %>%
  ungroup() %>%
  left_join(transcriptome_codon_usage) %>%
  mutate(normalised_to_transcriptome = codon_usage/proportion_usage_transcriptome,
         protein_length = nchar(protein_seq)) %>%
  dplyr::rename(codon_count = n)

write_tsv(codon_usage_per_transcript %>% select(-c(cds_seq, protein_seq, entropy)), 
          paste0("~/home/GASR/germs/data/multiple_species/species_codons_per_transcript/", species_name, ".codons_per_transcript.tsv.gz"))

# Codons outside of LCDs --------------------------------------------------

codons_outside_leds <- all_entropy %>%
  select(-c(genomic_loc, gene_id, transcript_biotype, gene_biotype, cds_len)) %>%
  unnest(le_region) %>%
  mutate(region = map(le_region, ~{ paste0(.x[1], "-", (.x[2] + (win_size - 1) ) ) }) %>% unlist() %>% paste0(transcript_id, ":", .)) %>% 
  mutate(codons_in_domain = pmap(list(le_region, codons),
                                 function(le_region, codons){
                                   c(codons[1:(le_region[1] - 1)], codons[le_region[2] + win_size:length(codons)]) 
                                 })) %>%
  mutate(codons_frequency = map(codons_in_domain, ~{
    tibble(codon = .x) %>%
      group_by(codon) %>%
      count() %>%
      ungroup() %>%
      right_join(codon_table, by = "codon")
  })) %>% 
  dplyr::select(-c(le_region, codons, codons_in_domain)) %>%
  unnest(codons_frequency) %>%
  replace_na(list(n = 0))

codon_usage_outside_leds <- codons_outside_leds %>%
  group_by(species, region, aa) %>%
  mutate(codon_usage = n/sum(n)) %>%
  ungroup() %>%
  left_join(transcriptome_codon_usage) %>%
  mutate(normalised_to_transcriptome = codon_usage/proportion_usage_transcriptome)

write_tsv(codon_usage_outside_leds %>% select(-c(cds_seq, protein_seq, entropy, n)), 
          paste0("~/home/GASR/germs/data/multiple_species/species_codons_outside_leds/", species_name, ".codons_outside_leds.tsv.gz"))
