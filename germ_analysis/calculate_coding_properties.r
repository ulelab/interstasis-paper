rm(list = ls())

library(tidyverse)
library(parallel)
library(Biostrings)
library(ggsci)

# Params -----------------------------------------------------------------

k_length <- 5
window_size <- 123
smoothing_size <- 123

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

# Fetch transcripts -------------------------------------------------------

sequences <- readDNAStringSet("genomes/hs/fasta/longest_gencode29.fa") %>% as.character()

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

sequences <- sequences[match(transcript_details$transcript_id, names(sequences))] #Reorder sequences to match transcript detail table

# Load clusters --------------------------------------------------------------------

cds_clusters <- read_tsv("GASR/germs/data/germs_CDS_umap_clusters.tsv.gz")

annotated_clusters <- cds_clusters %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "cds"), sep = ":|-|@", remove = F, convert = T) %>%
  dplyr::select(-cds) %>%
  mutate(end = end + 4) %>% #Adjust for the k length.
  left_join(transcript_details)

# Fetch coding information ------------------------------------------------

aa_usage <- transcript_details %>%
  mutate(protein_sequence = sequences[transcript_id] %>% 
           str_sub(cds_start, cds_end) %>% 
           DNAStringSet() %>%
           translate() %>%
           as.character()) %>%
  select(protein_sequence) %>%
  unlist(use.names = F) %>%
  paste0() %>%
  str_split("") %>%
  unlist(use.names = F) %>%
  table() %>%
  as_tibble() %>%
  set_names(c("aa", "count")) %>%
  filter(aa != "*") %>%
  mutate(proteome_proportion = count/sum(count)) %>%
  select(-count)

aa_clusters <- annotated_clusters %>%
  mutate(start = case_when(start < cds_start ~ as.integer(cds_start),
                           T ~ start),
         end = case_when(end > cds_end ~ cds_end,
                         T ~ end)) %>%
  mutate(aa_start = round((start - cds_start) / 3) + 1) %>%
  mutate(aa_end = aa_start + round((end - start) / 3)) %>%
  mutate() %>%
  mutate(protein_sequence = sequences[transcript_id] %>% 
           str_sub(cds_start, cds_end) %>% 
           DNAStringSet() %>%
           translate() %>%
           as.character() %>%
           str_sub(aa_start, aa_end))

fav_aa <- aa_clusters %>%
  mutate(most_common_aa = map(protein_sequence, ~{ .x %>% str_split("") %>% 
      table() %>%
      as_tibble %>%
      set_names("aa", "occurence") %>% 
      arrange(desc(occurence)) %>%
      dplyr::slice(1L) %>%
      dplyr::select(aa) %>%
      unlist(use.names = F)}) %>% unlist())

kept_aa <- fav_aa %>% group_by(most_common_aa) %>% count() %>% filter(n >= 25) %>% dplyr::select(most_common_aa) %>% unlist(use.names = F)

fav_aa %>% 
  group_by(cluster_name, most_common_aa) %>%
  count() %>% 
  group_by(cluster_name) %>% mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  mutate(most_common_aa = most_common_aa %>% fct_relevel("R", "K", "E", "D",
                                                         "G", "A", "L", "P",
                                                         "S", "T", "Q")) %>%
  arrange(most_common_aa)

# plot_colors <- c(pal_simpsons()(10), "#cfcfcf")

plot_colors <- c("#577399",
                 "#BDD5EA",
                 "#E86252",
                 "#ad4c7e",
                 "#5D576B",
                 "#DFA06E",
                 "#D8D4D5",
                 "#297373",
                 "#E9D758",
                 "#B592A0",
                 "#A1C181")


fav_aa_umap_plot <-
  ggplot(fav_aa %>%
         filter(cluster > 0,
                most_common_aa %in% kept_aa) %>%
           mutate(most_common_aa = most_common_aa %>% fct_relevel("R", "K", "E", "D",
                                                                  "G", "A", "L", "P",
                                                                  "S", "T", "Q")) %>%
           sample_frac(1L) %>%
           arrange(most_common_aa == "R"),
       aes(x = V1, y = V2, fill = most_common_aa)) +  
  scale_fill_manual(values = plot_colors) +
  geom_point(size = 3, pch = 21) +
  guides(fill = guide_legend(ncol = 2, 
                             title = element_blank(),
                             override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
  )

# Don't emphasise arginines
fav_aa_umap_plot_not_r_ordered <-
  ggplot(fav_aa %>%
           filter(cluster > 0,
                  most_common_aa %in% kept_aa) %>%
           mutate(most_common_aa = most_common_aa %>% fct_relevel("R", "K", "E", "D",
                                                                  "G", "A", "L", "P",
                                                                  "S", "T", "Q")) %>%
           sample_frac(1L),
         aes(x = V1, y = V2, fill = most_common_aa)) +  
  scale_fill_manual(values = plot_colors) +
  geom_point(size = 3, pch = 21) +
  guides(fill = guide_legend(ncol = 2, 
                             title = element_blank(),
                             override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
  )

ggsave(filename = "GASR/germs/plots/germs_CDS_fav_aa.pdf",
       fav_aa_umap_plot,
       device = "pdf", units = "in",
       height = 4.5, width = 5.9) #Yeah that's right, 5.9.

ggsave(filename = "GASR/germs/plots/germs_CDS_fav_aa_not_r_ordered.pdf",
       fav_aa_umap_plot_not_r_ordered,
       device = "pdf", units = "in",
       height = 4.5, width = 5.9) #Yeah that's right, 5.9.


# Amino acids in clusters -------------------------------------------------

aa_clusters %>%
  dplyr::select(cluster, cluster_name, cluster_colour) %>%
  distinct() %>%
  mutate(cluster_name_short = c("None", "G-rich Pur.", "A-rich Pur.",
                                "G-rich G/C", "C-rich", "CUG-rep", "CU-rich",
                                "Pur. + C", "CAG-rep")) %>%
  mutate(new_order = c(0, 2, 1, 4, 6, 5, 7, 3, 8)) %>%
  arrange(new_order) %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel(cluster_name_short)) ->
  cluster_rename_df

# cluster_order <- aa_clusters %>%
#   filter(cluster > 0) %>%
#   select(cluster_name, cluster) %>%
#   distinct()
# 
# cluster_order <- cluster_order[c(2,7,1,4,6,8,3,5),]

aa_biases <- aa_clusters %>%
  filter(cluster > 0) %>%
  group_by(cluster, cluster_name, cluster_colour) %>%
  # summarise(aas = paste(unlist(protein_sequence), collapse = "")) %>%
  summarise(aa_distribution = map(aas,
                                  ~{ .x %>% 
                                      str_split("") %>%
                                      unlist(use.names = F) %>%
                                      table() %>%
                                      as_tibble() %>%
                                      set_names("aa", "occurence") %>%
                                      filter(aa != "*") %>%
                                      mutate(proportion = occurence/sum(occurence)) %>%
                                      right_join(aa_usage)})) %>%
  unnest(aa_distribution) %>%
  ungroup() %>%
  mutate(relative_abundance = proportion/proteome_proportion) %>%
  left_join(cluster_rename_df)


# SLOW AS FUCK
# aa_biases_per_region <- aa_clusters %>%
#   filter(cluster > 0) %>%
#   group_by(cluster, cluster_name, cluster_colour) %>%
#   summarise(aa_distribution = map(protein_sequence,
#                                   ~{ .x %>% 
#                                       str_split("") %>%
#                                       unlist(use.names = F) %>%
#                                       table() %>%
#                                       as_tibble() %>%
#                                       set_names("aa", "occurence") %>%
#                                       filter(aa != "*") %>%
#                                       mutate(proportion = occurence/sum(occurence)) %>%
#                                       right_join(aa_usage, by = "aa")})) %>%
#   unnest(aa_distribution) %>%
#   mutate(relative_abundance = proportion/proteome_proportion)
# 
# aa_biases %>%
#   filter(aa %in% c("E", "R", "K", "D")) %>%
#   ggplot() +
#   aes(x = cluster_name, y = relative_abundance) +
#   geom_bar(stat = "identity", width = 0.5) +
#   facet_wrap(. ~ aa) +
#   scale_y_log10() +
#   theme_classic() +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   labs(x = "", y = "fold enrichment \n over proteome") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"),
#         axis.text.y = element_text(colour = "black"))


aa_biases %>%
  mutate(aa = aa %>% fct_relevel("E", "D", "K", "R", "P", "A", "G", "S", "Q")) %>%
  filter(aa %in% c("E", "D", "K", "R", "P", "A", "G", "S", "Q")) %>%
  arrange(new_order) %>%
  ggplot() +
  aes(x = cluster_name_short, y = aa, fill = log2(relative_abundance)) +
  geom_tile(colour = "black", linewidth = 0.75) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", ticks.linewidth = 0.5, frame.linewidth = 0.5)) +
  labs(fill = "log2\nenrichment\nover\nproteome") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(vjust = 0.5, colour = "black"),
        legend.title.align =  0.5) ->
  aa_bias_per_region_heatmap


ggsave(filename = "GASR/germs/plots/aa_bias_per_region_heatmap.pdf",
       aa_bias_per_region_heatmap,
       device = "pdf", units = "in",
       height = 3.5, width = 4)


# Arginine amounts --------------------------------------------------------

aa_clusters %>%
  filter(cluster > 0) %>%
  left_join(cluster_rename_df) %>%
  group_by(cluster, cluster_name, cluster_colour) %>%
  # summarise(aas = paste(unlist(protein_sequence), collapse = "")) %>%
  mutate(aa_distribution = map(protein_sequence,
                                  ~{ .x %>% 
                                      str_split("") %>%
                                      unlist(use.names = F) %>%
                                      table() %>%
                                      as_tibble() %>%
                                      set_names("aa", "occurence") %>%
                                      filter(aa != "*") %>%
                                      mutate(proportion = occurence/sum(occurence))})) %>%
  ungroup() %>%
  select(peak_identifier, new_order, cluster_name_short, cluster_colour, aa_distribution) %>%
  unnest(aa_distribution) %>%
  select(-occurence) %>%
  pivot_wider(names_from = aa, values_from = proportion) %>%
  replace(is.na(.), 0) ->
  aa_proportion_per_region

ggplot(aa_proportion_per_region, aes(R, color = cluster_name_short)) + 
  geom_density() +
  facet_wrap(. ~ cluster_name_short)

aa_proportion_per_region %>%
  mutate(is_r = R >= (2 * aa_usage$proteome_proportion[aa_usage$aa == "R"]),
         r_over_k = R > K) %>%
  group_by(cluster_name_short) %>%
  summarise(n_r = sum(is_r),
            n_other = dplyr::n() - sum(is_r)) %>%
  pivot_longer(cols = c(n_r, n_other)) %>%
  mutate(name = name %>% str_replace("n_r", "R-rich") %>% str_replace("n_other", "Other")) ->
  prop_argy

argy_annots <- prop_argy %>%
  mutate(name = name %>% str_replace("R-rich", "R")) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  mutate(position = Other + R + 100,
         name = "R-rich")

prop_argy %>%
  ggplot(aes(x = cluster_name_short, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  geom_text(data = argy_annots, aes(x = cluster_name_short, y = position, label = R)) +
  scale_fill_manual(values = c("white", "black")) +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 400, 800, 1200)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black")) +
  labs(x = "",
       y = "number of regions",
       fill = "") ->
  r_rich_prop_plot

ggsave("GASR/germs/plots/arginine/r_double_proteome_per_germ_cluster.pdf", r_rich_prop_plot,
       height = 4, width = 4)

