# Setup -------------------------------------------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(pbapply)
library(patchwork)
library(plyranges)
library(roll)
library(Biostrings)

# entropy stuff -----------------------------------------------------------

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

# Functions ---------------------------------------------------------------

icount_bedloader <- function(directory, filter_term){
  #Take a directory containing .bed files from icount-style xlinks sites.
  #Samples can be filtered using filter_term. If no filtering is desired, set filter_term to "".
  
  # directory <- "test/"
  
  bedfiles <- list.files(directory, pattern = ".bed", full.names = T)
  bednames <- list.files(directory, pattern = ".bed", full.names = F)
  
  #Filter the samples based on the filtering
  bedfiles <- bedfiles[grepl(filter_term, bednames)]
  bednames <- bednames[grepl(filter_term, bednames)]
  
  bed <- pblapply(bedfiles, function(x){
    fread(x)
  }, cl = 16)
  
  bedlist <- pblapply(bed, function(x){
    #For negative stranded reads, add the width of the read (sum the cigar number) to the start position
    #to get the actual 5' position of the read.
    
    x$start[x$V6 == "+"] <- x$V2[x$V6 == "+"] + 1 #BED starts are 0-based.
    x$start[x$V6 == "-"] <- x$V3[x$V6 == "-"] #But BED end positions are 1-based.
    
    return(
      GRanges(seqnames = x$V1,
              ranges = IRanges(start = x$start, width = 1),
              strand = x$V6,
              count = x$V5))
  }, cl = 16)
  
  names(bedlist) <- lapply(bednames, function(x){
    strsplit(x, ".", fixed = T)[[1]][1]
  })
  
  bedgr <- unlist(GRangesList(bedlist))
  
  bedgr$sample <- names(bedgr)
  
  return(bedgr)
}

xlink_profile <- function(iclip, gene_name, smooth_width, plot_title){
  
  # smooth_width = 50
  # gene_name = "RBM39"
  # iclip <- all_clip
  # plot_title <- "NCL"
  
  lib_sizes <- iclip %>% group_by(sample) %>%
    summarise(xlinks = sum(count)/1000000)
  
  if(!smooth_width %% 2){
    smooth_width <- smooth_width + 1
  }
  
  smooth_weights <- dnorm(1:smooth_width, mean = ceiling(smooth_width/2), sd = smooth_width/5)
  smooth_weights <- smooth_weights/sum(smooth_weights) #Ensure the kernel adds to 1.
  
  iclip <- select(iclip, seqnames, start, end, width, strand, count, sample, target)
  
  sample_key <- iclip %>% 
    select(sample, target) %>%
    distinct()
  
  map <- data.frame(start = rep(c(1:(tx_details$tx_length[tx_details$gene_name == gene_name] + (smooth_width - 1))), 
                                each = length(iclip$sample %>% unique)),
                    sample = rep(iclip$sample %>% unique, 
                                 (tx_details$tx_length[tx_details$gene_name == gene_name]) + (smooth_width - 1))) %>%
    left_join(sample_key)
  
  transcript_id <- tx_details$transcript_id[tx_details$gene_name == gene_name]
  
  region_xlinks <- iclip %>%
    filter(seqnames == transcript_id) %>%
    left_join(lib_sizes, by = "sample") %>%
    mutate(count = count/xlinks) %>%
    dplyr::select(-c(end, seqnames, strand, width, xlinks)) %>%
    full_join(map, by = c("start", "sample", "target")) %>%
    replace_na(list(count = 0)) %>%
    arrange(start) %>%
    group_by(sample, target) %>%
    mutate(smooth_counts = c(rep(0, (smooth_width - 1)/2), 
                             count %>% roll_mean(width = smooth_width, weights = smooth_weights, online = F) %>% .[-c(1:smooth_width-1)],
                             rep(0, (smooth_width - 1)/2))) %>%
    group_by(start, target) %>%
    summarise(mean_count = mean(smooth_counts),
              se_count = sd(smooth_counts)/sqrt(dplyr::n())) %>%
    mutate(se_top = mean_count + se_count,
           se_bottom = mean_count - se_count) %>%
    #Trim off the edges introduced for the purposes of smoothing.
    mutate(start = start - (smooth_width - 1)/2) %>%
    filter(start >= 1,
           start <= tx_details$tx_length[tx_details$gene_name == gene_name])
  
  ggplot(region_xlinks,
         aes(x = start)) +
    geom_line(aes(y = mean_count, 
                  colour = target), 
              linewidth = 0.5, 
              alpha = 0.7, 
              show.legend = F) +
    geom_ribbon(aes(ymin = se_top,
                    ymax = se_bottom,
                    fill = target),
                alpha = 0.3,
                show.legend = F) +
    theme_classic() +
    facet_grid(target ~ ., scales = "free_y") +
    labs(y = "CLIP CPM") +
    scale_colour_manual(values = c("red4", "orangered3", "blue4", "purple4", "green4", "gray50")) +
    scale_fill_manual(values = c("red4", "orangered3", "blue4", "purple4", "green4", "gray50")) +
    theme(axis.line.x = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "plain"),
          legend.title = element_blank()) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2, min.n = 2))
  
}


germ_profile <- function(gene_name) {
  
  germ_scores_shuff %>% 
    filter(transcript_id == tx_details$transcript_id[tx_details$gene_name == gene_name]) %>%
    mutate(pos = position) %>%
    ggplot() +
    aes(x = pos) +
    geom_line(aes(y = med_norm_smooth)) +
    geom_line(aes(y = med_norm_shuff), linetype = "dashed") +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black")) +
    labs(x = "", y = "smoothed norm.\nGeRM score") +
    # coord_cartesian(xlim = c(1, tx_details_seq$tx_length[tx_details_seq$gene_name == gene_name])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3, min.n = 2))
  
}

germ_profile_no_shuffle <- function(gene_name) {
  
  germ_scores %>% 
    filter(transcript_id == tx_details$transcript_id[tx_details$gene_name == gene_name]) %>%
    mutate(pos = position) %>%
    ggplot() +
    aes(x = pos) +
    geom_line(aes(y = med_norm_smooth)) +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black")) +
    labs(x = "", y = "smoothed norm.\nGeRM score") +
    # coord_cartesian(xlim = c(1, tx_details_seq$tx_length[tx_details_seq$gene_name == gene_name])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3, min.n = 2))
  
}

disorder_profile <- function(gene_name2){
  
  transcript_af2_disorder %>%
    filter(gene_name == gene_name2) %>%
    unnest(data) %>%
    ungroup() ->
    gene_disorder_df
  
  data.frame(disorder = rep(gene_disorder_df$disorder_25, each = 3),
             plddt = rep(gene_disorder_df$lddt, each = 3),
             rsa = rep(gene_disorder_df$rsa, each = 3),
             pos = gene_disorder_df$cds_start[1]:(gene_disorder_df$cds_end[1] - 3)) %>%
    ggplot() +
    aes(x = pos, y = plddt) +
    geom_line() +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black")) +
    labs(x = "", y = "pLDDT") +
    # coord_cartesian(xlim = c(1, gene_disorder_df$tx_length[1])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3, min.n = 2))
}

plot_annotation <- function(gene_name) {
  
  data.frame(xmi = c(0, 
                     tx_details$cds_start[tx_details$gene_name == gene_name],
                     tx_details$cds_end[tx_details$gene_name == gene_name]),
             xma = c(tx_details$cds_start[tx_details$gene_name == gene_name], 
                     tx_details$cds_end[tx_details$gene_name == gene_name],
                     tx_details$tx_length[tx_details$gene_name == gene_name]),
             ymi = c(0.4, 0.2, 0.4),
             yma = c(0.6, 0.8, 0.6)) %>%
    ggplot(aes(xmin = xmi,
               xmax = xma,
               ymin = ymi,
               ymax = yma)) +
    geom_rect(fill = "grey30", colour = "grey30", size = 1) +
    theme_classic() +
    coord_cartesian(ylim = c(0,1)) +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(colour = "black")) +
    annotate("text", x = (tx_details$cds_start[tx_details$gene_name == gene_name] +
                            tx_details$cds_end[tx_details$gene_name == gene_name])/2, 
             y = 0.5, label = gene_name, colour = "white", fontface = 2, size = 3) +
    labs(x = "Position along transcript",
         y = "")
  
}

aa_profile <- function(gene_name, window){
  
  # gene_name = "RBM27"
  # window = 20
  
  aas <- tx_details_seq$aa[tx_details_seq$gene_name == gene_name] %>% strsplit("") %>% unlist()
  
  charged_aa <- aas %in% c("D", "E", "K", "R") %>% as.numeric()
  
  rolled_aa <- roll_mean(c(rep(0, window/2 - 1), charged_aa, rep(0, window/2)),
                         width = window)[-c(1:(window-1))]
  
  charged_df <- data.frame(charged = rep(rolled_aa, each = 3), 
                           nothing = "nothing") %>%
    mutate(position = 
             tx_details_seq$cds_start[tx_details_seq$gene_name == gene_name]:
             tx_details_seq$cds_end[tx_details_seq$gene_name == gene_name])
  
  ggplot(charged_df, aes(position, nothing, fill = charged)) +
    geom_tile(show.legend = F) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.direction = "horizontal",
          legend.title = element_text(hjust = 0.5)) +
    labs(fill = "Proportion charged\namino acids") +
    coord_cartesian(xlim = c(1, tx_details_seq$tx_length[tx_details_seq$gene_name == gene_name])) +
    # scale_fill_viridis_c()
    scale_fill_gradient2(low = "white",
                         mid = "#ffce47",
                         high = "red", 
                         midpoint = median(charged_df$charged),
                         limits = c(0, max(charged_df$charged)))
  
}

entropy_profile <- function(gene_name2, window){
  
  # gene_name2 = "APP"
  # window = 41
  
  charge_colour = "darkorange2"
  
  # Get the net charge and smooth
  aas_stop <- tx_details_seq$aa[tx_details_seq$gene_name == gene_name2] %>% strsplit("") %>% unlist()
  aas <- aas_stop[-length(aas_stop)] #Drop the stop codon
  
  charged_aa <- aas %in% c("D", "E", "K", "R") %>% as.numeric()
  
  rolled_aa <- roll_mean(charged_aa, window)[-c(1:(window-1))] #NAs are at the start from roll_mean
  
  # Get entropy
  entropy_values <- aa_entropy(tx_details_seq$aa[tx_details_seq$gene_name == gene_name2], window)
  
  values_df <- data.frame(position = c(1:(length(entropy_values) * 3)) +
               tx_details_seq$cds_start[tx_details_seq$gene_name == gene_name2] + # shift positions to start codon
               ((window - 1) / 2) * 3, # shift positions to account for the lost values due to the smoothing window
             entropy = rep(entropy_values, each = 3),
             charged_aa = rep(rolled_aa, each = 3))
  
  ent_range = max(values_df$entropy) - min(values_df$entropy)
  charge_range = max(values_df$charged_aa) - min(values_df$charged_aa)
  
  multiplier = charge_range / ent_range
  subtraction = min(values_df$entropy * multiplier) - min(values_df$charged_aa * multiplier)
  
  ggplot(values_df) +
    aes(x = position) +
    geom_line(aes(y = entropy)) +
    geom_line(aes(y = (charged_aa + subtraction) / multiplier), colour = charge_colour) +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          # axis.title.y.right = element_text(color = charge_colour),
          # axis.text.y.right =  element_text(color = charge_colour),
          axis.line.y.right = element_line(color = charge_colour),
          axis.ticks.y.right = element_line(color = charge_colour),
          ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3, min.n = 2),
                       name = "amino acid\nentropy",
                       sec.axis = sec_axis(~ (. * multiplier) - subtraction,
                                           name = "prop. charged\namino acids",
                                           breaks = scales::pretty_breaks(n = 4, min.n = 2)))
  
  
}

assemble_plot <- function(gene_name) {
  
  # gene_name <- "HSP90AA1"
  cds_length <- tx_details$cds_length[tx_details$gene_name == gene_name]
  start_pos <- tx_details$cds_start[tx_details$gene_name == gene_name] - (cds_length / 20)
  end_pos <- tx_details$cds_end[tx_details$gene_name == gene_name] + (cds_length / 20)
  
  profile <- xlink_profile(all_clip, gene_name, 50, "")
  
  annotation_plot <- plot_annotation(gene_name)
  
  entropy_prof <- entropy_profile(gene_name, 41)
  
  germ_prof <- germ_profile_no_shuffle(gene_name)
  
  # disorder_prof <- disorder_profile(gene_name)
  
  ((profile +
      coord_cartesian(xlim = c(start_pos, end_pos))) /
      (germ_prof +
         coord_cartesian(xlim = c(start_pos, end_pos))) /
      # (disorder_prof +
      #    coord_cartesian(xlim = c(start_pos, end_pos))) /
      (entropy_prof +
         coord_cartesian(xlim = c(start_pos, end_pos))) / 
      (annotation_plot +
         coord_cartesian(xlim = c(start_pos, end_pos),
                         ylim = c(0, 1)))) +
    plot_layout(heights = c(10, 2, 2, 1))
  
}

assemble_plot_noclip <- function(gene_name) {
  
  # gene_name <- "NCL"
  cds_length <- tx_details$cds_length[tx_details$gene_name == gene_name]
  start_pos <- tx_details$cds_start[tx_details$gene_name == gene_name] - (cds_length / 20)
  end_pos <- tx_details$cds_end[tx_details$gene_name == gene_name] + (cds_length / 20)
  
  annotation_plot <- plot_annotation(gene_name)
  
  # aa_prof <- aa_profile(gene_name, 20)
  
  germ_prof = germ_profile(gene_name)
  
  disorder_prof <- disorder_profile(gene_name)
  
  entropy_profile <- entropy_profile(gene_name, 41)
  
  # wobble_cons_prof <- wobble_cons_profile(gene_name, 20)
  
  (
    (germ_prof +
       coord_cartesian(xlim = c(start_pos, end_pos))) /
      (entropy_profile +
         coord_cartesian(xlim = c(start_pos, end_pos))) /
      (disorder_prof +
         coord_cartesian(xlim = c(start_pos, end_pos))) /
      # (aa_prof +
      #    coord_cartesian(xlim = c(start_pos, end_pos))) / 
      (annotation_plot +
         coord_cartesian(xlim = c(start_pos, end_pos),
                         ylim = c(0, 1)))) +
    plot_layout(heights = c(3, 3, 3, 1))
  
}

assemble_plot_nodisorder <- function(gene_name) {
  
  # gene_name <- "NCL"
  cds_length <- tx_details$cds_length[tx_details$gene_name == gene_name]
  start_pos <- tx_details$cds_start[tx_details$gene_name == gene_name] - (cds_length / 20)
  end_pos <- tx_details$cds_end[tx_details$gene_name == gene_name] + (cds_length / 20)
  
  annotation_plot <- plot_annotation(gene_name)
  
  # aa_prof <- aa_profile(gene_name, 20)
  
  germ_prof = germ_profile(gene_name)
  
  # disorder_prof <- disorder_profile(gene_name)
  
  entropy_profile <- entropy_profile(gene_name, 41)
  
  # wobble_cons_prof <- wobble_cons_profile(gene_name, 20)
  
  (
    (germ_prof +
       coord_cartesian(xlim = c(start_pos, end_pos))) /
      (entropy_profile +
         coord_cartesian(xlim = c(start_pos, end_pos))) /
      # (disorder_prof +
      #    coord_cartesian(xlim = c(start_pos, end_pos))) /
      # (aa_prof +
      #    coord_cartesian(xlim = c(start_pos, end_pos))) / 
      (annotation_plot +
         coord_cartesian(xlim = c(start_pos, end_pos),
                         ylim = c(0, 1)))) +
    plot_layout(heights = c(3, 3, 1))
  
}

assemble_plot_noclip_nonothin <- function(gene_name) {
  
  # gene_name <- "NCL"
  cds_length <- tx_details$cds_length[tx_details$gene_name == gene_name]
  start_pos <- tx_details$cds_start[tx_details$gene_name == gene_name] - (cds_length / 20)
  end_pos <- tx_details$cds_end[tx_details$gene_name == gene_name] + (cds_length / 20)
  
  annotation_plot <- plot_annotation(gene_name)
  
  # aa_prof <- aa_profile(gene_name, 20)
  
  germ_prof = germ_profile(gene_name)
  
  # disorder_prof <- disorder_profile(gene_name)
  
  # entropy_profile <- entropy_profile(gene_name, 41)
  
  # wobble_cons_prof <- wobble_cons_profile(gene_name, 20)
  
  (
    (germ_prof +
       coord_cartesian(xlim = c(start_pos, end_pos))) /
      # (aa_prof +
      #    coord_cartesian(xlim = c(start_pos, end_pos))) / 
      (annotation_plot +
         coord_cartesian(xlim = c(start_pos, end_pos),
                         ylim = c(0, 1)))) +
    plot_layout(heights = c(3, 1))
  
}

# Load data ---------------------------------------------------------------

germ_scores <- read_tsv("GASR/germs/data/kmer_multivalency_within_annotated_germs_peaks.tsv.gz") %>%
  group_by(transcript_id) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  filter(smoothed_kmer_multivalency > 0) 

# For normalising the multivalency
minskm <- min(germ_scores$smoothed_kmer_multivalency)
medskm <- median(germ_scores$smoothed_kmer_multivalency)

# Normalise the multivalency
germ_scores <- germ_scores %>%
  mutate(med_norm_smooth = 
           (smoothed_kmer_multivalency - minskm)/
           (medskm - minskm))

germ_scores_shuff <- read_tsv("GASR/germs/data/smoothed_kmer_multivalency_native_and_shuffled_for_cds_cluster_transcripts.tsv.gz") %>%
  filter(smoothed_kmer_multivalency > 0) %>%
  mutate(med_norm_smooth = 
           (smoothed_kmer_multivalency - minskm)/
           (medskm - minskm),
         med_norm_shuff = 
           (mean_shuffled_smoothed_kmer_multivalency - minskm)/
           (medskm - minskm))
         
tx_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt")
nts <- Biostrings::readDNAStringSet("GASR/lists/longest_gencode29.fa")

germ_scores_shuff_overview <- read_tsv("GASR/germs/data/germs_CDS_umap_with_native_and_shuffled_mv.tsv") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "cds"), sep = "-|:|@") %>%
  left_join(tx_details, by = "transcript_id") %>%
  mutate(multi_diff = mean_native_smoothed_multivalency - mean_shuffled_smoothed_multivalency)

tx_details_seq <- data.frame(seq = as.character(nts),
                             transcript_id = names(nts)) %>%
  right_join(tx_details) %>%
  mutate(cds = str_sub(seq, cds_start, cds_end)) %>%
  mutate(aa = cds %>% DNAStringSet() %>% translate() %>% as.character())

# Loading disorder ------------------------------------------------------------

gencode_to_swissprot <- read_tsv("GASR/lists/gencode.v29.metadata.SwissProt.txt",
                                 col_names = c("transcript_id", "uniprot_id", "id2"),
                                 col_type = "c")

af2_disorder <- read_tsv("General/alpha_fold_disorder/af2_human_disorder_pred.tsv.gz",
                         skip = 1,
                         col_names = c("name", "pos", "aa", "lddt", "disorder", "rsa", "ss", "disorder_25", "binding_25_0581"),
                         col_types = "cdcdddcdd")

transcript_details_uniprot <- inner_join(gencode_to_swissprot, tx_details)

transcript_af2_disorder <- af2_disorder %>%
  group_by(name) %>%
  nest(data = c(pos, aa, lddt, disorder, rsa, ss, disorder_25, binding_25_0581)) %>%
  ungroup() %>%
  mutate(uniprot_id = word(name, 2, sep = fixed("-"))) %>%
  dplyr::select(-name) %>%
  inner_join(transcript_details_uniprot, multiple = "all")

#Does the length of the alpha fold vector correctly match the length of the CDS? If not, discard.
transcript_af2_disorder <- transcript_af2_disorder %>%
  ungroup() %>%
  mutate(disorder_nrow = map(data, ~{nrow(.x)}) %>% unlist()) %>%
  filter((disorder_nrow + 1) == (cds_length/3)) %>%
  dplyr::select(-disorder_nrow)

# Plots without CLIP ----------------------------------------------------------

# R-rich, likes A
LUC7L3 <- assemble_plot_nodisorder("LUC7L3")
ggsave("GASR/germs/plots/example_genes/conservation/LUC7L3_profile.pdf", LUC7L3,
       width = 8, height = 2.5)

UPF3B <- assemble_plot_nodisorder("UPF3B")
ggsave("GASR/germs/plots/example_genes/conservation/UPF3B_profile.pdf", UPF3B,
       width = 8, height = 2.5)


# R-rich, but likes GC
CCDC61 <- assemble_plot_nodisorder("CCDC61")
ggsave("GASR/germs/plots/example_genes/conservation/CCDC61_profile.pdf", CCDC61,
       width = 8, height = 2.5)

TMEM200B <- assemble_plot_nodisorder("TMEM200B")
ggsave("GASR/germs/plots/example_genes/conservation/TMEM200B_profile.pdf", TMEM200B,
       width = 8, height = 2.5)


# Loading CLIP data -------------------------------------------------------

sr_hela <- icount_bedloader("GASR/tmap/SR_Krchnakova_Hela/results/xlinks_tome/", "srsf2|srsf6") %>%
  unname() %>%
  as.data.frame() %>%
  separate(sample, into = c("target", "replicate"), remove = F) %>%
  mutate(target = toupper(target))

sreclip <- icount_bedloader("GASR/tmap/SR_eCLIP/results/xlinks_tome/", "") %>%
  unname() %>%
  as.data.frame() %>%
  separate(sample, into = c("target", "celltype", "replicate"), remove = F) %>%
  mutate(target = toupper(target)) %>%
  filter(replicate != "mock")

tra2b <- icount_bedloader("GASR/tmap/TRA2B_Best_et_al/results/xlinks_tome/", "") %>%
  unname() %>%
  as.data.frame() %>%
  separate(sample, into = c("target", "replicate"), remove = F) %>%
  mutate(target = toupper(target))

all_clip <- list(sreclip, tra2b) %>% bind_rows()

# Plots -------------------------------------------------------------------

app <- assemble_plot("APP")

hsp90aa1 <- assemble_plot("HSP90AA1")

ncl <- assemble_plot("NCL")

eif3a <- assemble_plot("EIF3A")

psap <- assemble_plot("PSAP")

hnrnpdl <- assemble_plot("HNRNPDL")

ggsave("GASR/germs/plots/example_genes/clip/HSP90AA1_sr_prots.pdf",
       hsp90aa1,
       device = "pdf", units = "in",
       width = 6, height = 6)

ggsave("GASR/germs/plots/example_genes/clip/APP_sr_prots.pdf",
       app,
       device = "pdf", units = "in",
       width = 6, height = 6)

ggsave("GASR/germs/plots/example_genes/clip/NCL_sr_prots.pdf",
       ncl,
       device = "pdf", units = "in",
       width = 6, height = 6)

ggsave("GASR/germs/plots/example_genes/clip/EIF3A_sr_prots.pdf",
       eif3a,
       device = "pdf", units = "in",
       width = 6, height = 6)

ggsave("GASR/germs/plots/example_genes/clip/HNRNPDL_sr_prots.pdf",
       hnrnpdl,
       device = "pdf", units = "in",
       width = 6, height = 6)

ggsave("GASR/germs/plots/example_genes/clip/PSAP_sr_prots.pdf",
       psap,
       device = "pdf", units = "in",
       width = 6, height = 6)
