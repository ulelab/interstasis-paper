# Setup -------------------------------------------------------------------

rm(list = ls())
setwd("/camp/lab/ulej/home/users/farawar")

library(data.table)
library(tidyverse)
library(pbapply)
library(cowplot)
library(patchwork)
library(plyranges)
library(rtracklayer)
library(roll)

library(ShortRead)
library(ggpubr)

# Functions ---------------------------------------------------------------

icount_bedloader <- function(directory){
  #Take a directory containing .bed files from iMaps's xlinks sites.
  
  # directory <- "test/"
  
  bedfiles <- list.files(directory, pattern = ".bed", full.names = T)
  bednames <- list.files(directory, pattern = ".bed", full.names = F)
  
  bed <- pblapply(bedfiles, function(x){
    fread(x)
  })
  
  bedlist <- pblapply(bed, function(x){
    
    x$start[x$V6 == "+"] <- x$V2[x$V6 == "+"] + 1 #BED starts are 0-based.
    x$start[x$V6 == "-"] <- x$V3[x$V6 == "-"] #But BED end positions are 1-based.
    
    return(
      GRanges(seqnames = x$V1,
              ranges = IRanges(start = x$start, width = 1),
              strand = x$V6,
              count = x$V5))
  })
  
  names(bedlist) <- bednames %>% str_remove(".bed")
  
  bedgr <- unlist(GRangesList(bedlist))
  
  bedgr <- keepStandardChromosomes(bedgr, pruning.mode = "coarse")
  
  seqlevelsStyle(bedgr) <- "UCSC"
  
  bedgr$sample <- names(bedgr)
  
  return(bedgr)
}

lib_size <- function(bf){
  #Takes a GRanges from bamloader/bedloader and returns a df of library sizes for each sample (divided by 1 million).
  bf %>% unname() %>% as.data.frame() %>% group_by(sample) %>% summarise(xlinks = sum(count) / 1000000) %>% return()
}

exonmap <- function(reads, model){
  #Takes output of bedloader and maps x-link positions to exons.
  
  #Where do the xlinks land in the exons?
  key <- findOverlaps(reads, model, ignore.strand = F)
  
  #Take exons that have hits - repeated for each crosslink in that exon.
  locs <- model[subjectHits(key)]
  
  #Assign crosslink positions to exons.
  locs$xlink <- start(ranges(reads))[queryHits(key)]
  
  #Assign metacolumns from iCLIP to location list.
  for(i in names(mcols(reads))){
    mcols(locs)[i] <- mcols(reads)[queryHits(key), i]
  }
  # 
  # locs$sample <- reads$sample[queryHits(key)]
  # 
  # locs$junction <- reads$junction[queryHits(key)]
  # 
  # locs$span <- reads$span[queryHits(key)]
  
  return(locs)
}

metaprofile <- function(xlinks, exons, ss, map_range, smooth_width, grouping, grouping2, normalise_grouping2){
  
  #xlinks - granges from bedloader of xlink positions and counts, as well as annotated with sample grouping variables.
  #smooth_width - width of the smoothing kernel. Must be odd, or everything will break.
  
  #normalise grouping2 only works if grouping2 is a column in the exons input, not crosslinks.
  
  # xlinks = asynch_clip_tg %>%
  #   filter(target == "Srrm2") %>%
  #   filter(fraction == "Low")
  # exons = mm_internal_ga_ex
  # ss = "5ss"
  # map_range = c(50,0)
  # smooth_width = 11
  # grouping = "target"
  # grouping2 = "motif_class"
  # normalise_grouping2 = T
  
  #Prevent fuckery.
  if(!smooth_width %% 2){
    smooth_width <- smooth_width + 1
  }
  
  if(normalise_grouping2 & grouping2 != "none") {
    motif_class_norm <- exons %>% as.data.frame() %>% group_by(across(grouping2)) %>% summarise(class_size = dplyr::n()) %>%
      mutate(class_size = class_size/sum(class_size))
  } else if(grouping2 != "none") {
    #A tragic solution to my bloated code.
    motif_class_norm <- tibble(tempname = mcols(xlinks)[,grouping2] %>% unique(),
                               class_size = 1) %>% set_names(c(grouping2, "class_size"))
  } else {
    motif_class_norm <- tibble()
  }
  
  #Extend map range to account for smoothing
  extended_map_range <- map_range + ((smooth_width - 1)/2)
  
  library_size_df <- lib_size(xlinks)
  
  if(ss == "5ss"){
    
    new_exons <- exons %>%
      anchor_3p() %>%
      mutate(width = 1) %>%
      anchor_3p() %>%
      stretch(extended_map_range[1]) %>%
      anchor_5p() %>%
      stretch(extended_map_range[2])
    
  } else if(ss == "3ss") {
    
    new_exons <- exons %>%
      anchor_5p() %>%
      mutate(width = 1) %>%
      anchor_5p() %>%
      stretch(extended_map_range[2]) %>%
      anchor_3p() %>%
      stretch(extended_map_range[1])
    
  }
  
  clip_df <- exonmap(xlinks, new_exons) %>%
    as.data.frame() %>%
    mutate(pos = case_when(strand == "+" ~ (xlink - start) - extended_map_range[1],
                           strand == "-" ~ (end - xlink) - extended_map_range[1]))
  
  
  #Determine all possible permutations of grouping variables.
  if(grouping2 == "none") {
    
    expanded_factors <- clip_df %>%
      group_by(across(c("sample", grouping))) %>%
      summarise()
    
  } else {
    
    expanded_factors <- clip_df %>%
      group_by(across(c("sample", grouping, grouping2))) %>%
      summarise()
    
  }
  
  map <- data.frame(pos = -extended_map_range[1]:extended_map_range[2] %>% rep(each = nrow(expanded_factors)),
                    sample = expanded_factors$sample) %>%
    left_join(expanded_factors, relationship = "many-to-many")
  
  
  #Set up smoothing parameters.
  
  smooth_weights <- dnorm(1:smooth_width, mean = ceiling(smooth_width/2), sd = smooth_width/5)
  smooth_weights <- smooth_weights/sum(smooth_weights) #Ensure the kernel adds to 1.
  
  #Group by position, count reads for each group/sample, then normalise by grouping factors.
  if(grouping2 == "none") {
    positions <- clip_df %>%
      ungroup() %>%
      group_by(across(c("pos", "sample", grouping))) %>%
      summarise(count = sum(count)) %>%
      ungroup() %>%
      full_join(map, by = c("pos", "sample", grouping)) %>%
      replace_na(list(count = 0)) %>%
      left_join(library_size_df, by = "sample") %>%
      mutate(count = count/xlinks) %>%
      select(-xlinks) %>%
      arrange(pos) %>%
      group_by(sample) %>%
      mutate(smooth_counts = c(rep(0, (smooth_width - 1)/2), 
                               count %>% roll_mean(width = smooth_width, weights = smooth_weights, online = F) %>% .[-c(1:smooth_width-1)],
                               rep(0, (smooth_width - 1)/2))) %>%
      group_by(across(c("pos", grouping))) %>%
      summarise(mean_count = mean(smooth_counts),
                se_count = sd(smooth_counts)/sqrt(dplyr::n())) %>%
      mutate(se_top = mean_count + se_count,
             se_bottom = mean_count - se_count)
    
    
  } else {
    positions <- clip_df %>%
      ungroup() %>%
      group_by(across(c("pos", "sample", grouping, grouping2))) %>%
      summarise(count = sum(count)) %>%
      ungroup() %>%
      full_join(map, by = c("pos", "sample", grouping, grouping2)) %>%
      replace_na(list(count = 0)) %>%
      left_join(library_size_df, by = "sample") %>%
      mutate(count = count/xlinks) %>%
      left_join(motif_class_norm, by = grouping2) %>%
      mutate(count = count/class_size) %>%
      select(-c(xlinks, class_size)) %>%
      arrange(pos) %>%
      group_by(across(c("sample", grouping2))) %>%
      mutate(smooth_counts = c(rep(0, (smooth_width - 1)/2), 
                               count %>% roll_mean(width = smooth_width, weights = smooth_weights, online = F) %>% .[-c(1:smooth_width-1)],
                               rep(0, (smooth_width - 1)/2))) %>%
      group_by(across(c("pos", grouping, grouping2))) %>%
      summarise(mean_count = mean(smooth_counts),
                se_count = sd(smooth_counts)/sqrt(dplyr::n())) %>%
      mutate(se_top = mean_count + se_count,
             se_bottom = mean_count - se_count)
    
  }
  
  positions %>%
    filter(pos >= -map_range[1], 
           pos <= map_range[2]) %>% 
    return()
  
}

kmer_generator <- function(nucleotides, k, edge_trim1, edge_trim2){
  # For dinucleotides supplied in the format c("N", "N"), generate a matrix of all
  # combinations of these nucleotides with length k, then remove the edge_trim1 first 
  # and edge_trim2 last combinations (these will be biased towards a single nucleotide).
  # Supply the final combinations as a '|' separated string for regex.
  
  # nucleotides <- c("G", "A")
  # k = 5
  # edge_trim = 1
  
  kmer_matrix <- expand.grid(c(rep(list(nucleotides), k))) #Make a matrix of all combinations of the nucleotides.
  kmers <- apply(format(kmer_matrix), 1, paste, collapse = "") #Paste the rows of the matrix together.
  
  return(
    paste(
      kmers[(1+edge_trim1):(length(kmers)-edge_trim2)], collapse = "|")) #Paste the kmers together, seperate by | for regex purposes and trim edge cases if desired.
}

exon_motif_classifier <- function(exons, motifs, class_breaks) {
  
  # exons <- mm_internal_ex
  
  class_breaks <- c(0, class_breaks, 1)
  
  class_labels <- lapply(1:(length(class_breaks)-1), function(x) {
    paste0(100 * class_breaks[x], "% - ", 100 * class_breaks[x+1], "%")
  }) %>% unlist() %>% as.factor()
  
  classified_exons <- exons %>%
    as.data.frame() %>%
    mutate(seq = getSeq(gen_fasta, exons) %>% as.character()) %>%
    mutate(motif_count = str_count(seq, motifs)) %>%
    mutate(motif_density = motif_count/width) %>%
    mutate(motif_class = cut(motif_density, 
                             breaks = quantile(motif_density, 
                                               probs = class_breaks), 
                             include.lowest = T, 
                             labels = class_labels)) %>%
    group_by(transcript_id) %>%
    mutate(transcript_motif_class = rev(class_labels)[1] %in% motif_class) %>%
    ungroup() %>%
    select(-c(seq)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    return()
}

# Loading exons for mouse -------------------------------------------------

gcm22 <- read_delim("genomes/mm/annotation/gencode.vM22.annotation.gff3",
                    delim = "\t",
                    col_names = c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
                    comment = "#")

mouse_tx <- read_lines("GASR/lists/longest_proteincoding_transcript_mm.txt")

gen_fasta <- FaFile("genomes/mm/fasta/GRCm38.p6.genome.fa")

#This is really really slow and I'm sorry.
#All exons from the longest protein coding isoform of each coding gene (excluding mitochondrial transcripts).

mm_proco_ex <- gcm22 %>%
  filter(type == "exon") %>%
  mutate(transcript_id = sapply(attributes,
                                function(x){
                                  str_split(string = x,
                                            pattern = "transcript_id=|;gene_type")[[1]][2]})) %>%
  filter(transcript_id %in% mouse_tx) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  select(-c(score, type, source, phase, attributes))

mm_internal_ex <- mm_proco_ex %>%
  as.data.frame() %>%
  group_by(transcript_id) %>%
  mutate(exon_number = row_number()) %>%
  filter(exon_number > 1,
         exon_number < max(exon_number)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

motifs <- kmer_generator(c("G", "A"), 5, 11, 1)

mm_internal_ga_ex <- exon_motif_classifier(exons = mm_internal_ex,
                                           motifs = motifs,
                                           class_breaks = c(.5,.9, .97)) 



# SON SRRM2 TRA2B ---------------------------------------------------------


# Load crosslinks ---------------------------------------------------------

mesc_mito <- icount_bedloader("GASR/tmap/mesc_mitosis_1/results/xlinks/")
sams <- mesc_mito$sample %>% Rle()

mesc_mito$target <- sams@values %>% str_split("_", simplify = T) %>% .[,1] %>% rep(sams@lengths)
mesc_mito$phase <- sams@values %>% str_split("_", simplify = T) %>% .[,2] %>% rep(sams@lengths)
mesc_mito$fraction <- sams@values %>% str_split("_", simplify = T) %>% .[,3] %>% rep(sams@lengths)

mesc_mito %>% filter(phase == "Asynchronous") %>% filter(!(target %in% c("SmB", "Prpf8"))) ->
  asynch_clip


mesc_mito_tg <- icount_bedloader("GASR/tmap/mesc_mitosis_1/results/xlinks_tome_genome/")
sams_tg <- mesc_mito_tg$sample %>% Rle()

mesc_mito_tg$target <- sams_tg@values %>% str_split("_", simplify = T) %>% .[,1] %>% rep(sams_tg@lengths)
mesc_mito_tg$phase <- sams_tg@values %>% str_split("_", simplify = T) %>% .[,2] %>% rep(sams_tg@lengths)
mesc_mito_tg$fraction <- sams_tg@values %>% str_split("_", simplify = T) %>% .[,3] %>% rep(sams_tg@lengths)

mesc_mito_tg %>% filter(phase == "Asynchronous") %>% filter(!(target %in% c("SmB", "Prpf8"))) ->
  asynch_clip_tg


# Percentage crosslinks spliced -------------------------------------------

mm_internal_ga_ex %>%
  anchor_3p() %>%
  mutate(width = 1) %>%
  anchor_3p() %>%
  mutate(width = 40) ->
  mm_internal_ga_ex_5ss_40

calculate_proportion_spliced <- function(clip_target, clip_fraction, exon_regions, minimum_xlinks){
  
  # exon_regions <- mm_internal_ga_ex_5ss_50
  # clip_target <- "Tra2b"
  # clip_fraction <- "All"
  # minimum_xlinks <- 10
  
  comparison_all <- asynch_clip %>% filter(target == clip_target,
                                           fraction == clip_fraction)
  
  comparison_tg <- asynch_clip_tg %>% filter(target == clip_target,
                                             fraction == clip_fraction)
  
  lapply(unique(comparison_all$sample), function(x){
    # x = unique(comparison_all$sample)[2]
    
    comparison_all_sample <- comparison_all %>% filter(sample == x)
    comparison_all_tg <- comparison_tg %>% filter(sample == x)
    
    
    tibble(type = exon_regions$motif_class,
           motif_density = exon_regions$motif_density,
           all_counts = count_overlaps(exon_regions, comparison_all_sample),
           spliced_counts = count_overlaps(exon_regions, comparison_all_tg)) %>%
      filter(all_counts > minimum_xlinks) %>%
      mutate(spliced_proportion = spliced_counts/all_counts) %>%
      mutate(sample = x)
  }) %>%
    bind_rows() %>%
    mutate(fraction = clip_fraction,
           target = clip_target)
  
}

min_xl = 10

list(
  calculate_proportion_spliced("Tra2b", "All", mm_internal_ga_ex_5ss_40, min_xl),
  # calculate_proportion_spliced("Srrm2", "High", mm_internal_ga_ex_5ss_40, min_xl),
  # calculate_proportion_spliced("Srrm2", "Low", mm_internal_ga_ex_5ss_40, min_xl),
  # calculate_proportion_spliced("Son", "High", mm_internal_ga_ex_5ss_40, min_xl),
  # calculate_proportion_spliced("Son", "Low", mm_internal_ga_ex_5ss_40, min_xl)) 
  ) %>%
  bind_rows() ->
  all_spliced_percentages

all_spliced_percentages %>%
  group_by(target, fraction, type, sample) %>%
  summarise(mean_spliced_proportion = mean(spliced_proportion)) %>%
  ungroup() %>%
  ggplot(aes(x = target, y = mean_spliced_proportion, fill = type)) +
  geom_point(shape = 21, size = 3, position = position_dodge2(width = 0.9)) +
  facet_wrap(target ~ fraction, scales = "free_x", ncol = 5) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_fill_manual(values = c("#FABA6B", "#E67C43", "#FA4C2D", "#E60E33")) +
  labs(x = "",
       y = "mean proportion of spliced reads",
       fill = "exon GA-richness percentile") ->
  ga_richness_plot

ggsave("GASR/germs/plots/mouse/percentage_spliced_reads_at_5ss.pdf", ga_richness_plot,
       height = 3, width = 5)


