# Setup -------------------------------------------------------------------

rm(list=ls())

setwd("/camp/lab/ulej/home/users/farawar")
library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)
library(pbapply)
library(plyranges)
library(roll)

# Functions  --------------------------------------------------------------

old_metaprofile <- function(xlinks, peaks, map_range, smooth_width, grouping, grouping2){
  
  # xlinks <- prpf8_mesc %>% filter(target == "Mock", location == "Nucleus")
  # peaks <- mesc_prpf8_low_peaks
  # map_range = 50
  # smooth_width = 5
  # grouping = "location"
  # grouping2 = "none"
  
  # Prevent fuckery.
  
  if(!smooth_width %% 2){
    smooth_width <- smooth_width + 1
  }
  xlinks <- xlinks %>% filter(strand == "+")
  
  #Extend map range to account for smoothing
  extended_map_range <- map_range + ((smooth_width - 1)/2)
  
  library_size_df <- lib_size(xlinks)
  
  new_peaks <- peaks %>%
    mutate(start = start + round(width/2)) %>%
    mutate(width = 1) %>%
    stretch(extended_map_range * 2)
  
  clip_df <- region_map(xlinks, new_peaks) %>%
    as.data.frame() %>%
    mutate(pos = (xlink - start) - extended_map_range)
  
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
  
  map <- data.frame(pos = -extended_map_range:extended_map_range %>% rep(each = nrow(expanded_factors)),
                    sample = expanded_factors$sample) %>%
    left_join(expanded_factors)
  
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
      dplyr::select(-xlinks) %>%
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
      dplyr::select(-xlinks) %>%
      arrange(pos) %>%
      group_by(sample) %>%
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
    filter(pos >= -map_range,
           pos <= map_range) %>%
    return()
}

metaprofile <- function(xlinks, peaks, map_range, smooth_width, grouping_var){
  
  # This metaprofiler takes CLIP replicates, and plots their crosslinking distribution as mean + SE
  # around the midpoint of different GeRM regions (grouping). It normalises by the number of GeRM regions used.
  
  # xlinks <- asynch_clip %>% filter(target == "Tra2b")
  # peaks <- germs_gr
  # map_range = 400
  # smooth_width = 51
  # grouping_var = "cluster_name"
  
  # Prevent fuckery.
  
  if(!smooth_width %% 2){
    smooth_width <- smooth_width + 1
  }
  
  xlinks <- xlinks %>% filter(strand == "+")
  
  #Extend map range to account for smoothing
  extended_map_range <- map_range + ((smooth_width - 1)/2)
  
  library_size_df <- lib_size(xlinks)
  grouping_size_df <- as.data.frame(peaks) %>% count(across(all_of(grouping_var)), name = "group_size")
  
  new_peaks <- peaks %>%
    mutate(start = start + round(width/2)) %>%
    mutate(width = 1) %>%
    stretch(extended_map_range * 2)
  
  clip_df <- region_map(xlinks, new_peaks) %>%
    as.data.frame() %>%
    mutate(pos = (xlink - start) - extended_map_range)
  
  #Determine all possible permutations of grouping variables.
  expanded_factors <- clip_df %>%
    group_by(across(all_of(c("sample", grouping_var)))) %>%
    summarise()
  
  library_size_df$sample
  
  map <- data.frame(pos = -extended_map_range:extended_map_range %>% rep(each = length(library_size_df$sample)),
                    sample = library_size_df$sample) %>%
    left_join(expanded_factors, by = "sample", relationship = "many-to-many")
  
  #Set up smoothing parameters.
  smooth_weights <- dnorm(1:smooth_width, mean = ceiling(smooth_width/2), sd = smooth_width/5)
  smooth_weights <- smooth_weights/sum(smooth_weights) #Ensure the kernel adds to 1.
  
  #Group by position, count reads for each group/sample, then normalise by grouping factors.
  positions <- clip_df %>%
    ungroup() %>%
    group_by(across(all_of(c("pos", "sample", grouping_var)))) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    full_join(map, by = c("pos", "sample", grouping_var)) %>%
    replace_na(list(count = 0)) %>%
    left_join(library_size_df, by = "sample") %>%
    mutate(count = count/xlinks) %>%
    # dplyr::select(-xlinks) %>%
    left_join(grouping_size_df, by = grouping_var) %>%
    mutate(count = count/group_size) %>%
    # dplyr::select(-group_size) %>%
    arrange(pos) %>%
    ungroup() %>%
    group_by(across(all_of(c("sample", grouping_var)))) %>%
    mutate(smooth_counts = c(rep(0, (smooth_width - 1)/2),
                             count %>% roll_mean(width = smooth_width, weights = smooth_weights, online = F) %>% .[-c(1:smooth_width-1)],
                             rep(0, (smooth_width - 1)/2))) %>%
    group_by(across(all_of(c("pos", grouping_var)))) %>%
    summarise(mean_count = mean(smooth_counts),
              se_count = sd(smooth_counts)/sqrt(dplyr::n())) %>%
    mutate(se_top = mean_count + se_count,
           se_bottom = mean_count - se_count)
  
  positions %>%
    filter(pos >= -map_range,
           pos <= map_range) %>%
    return()
}


region_map <- function(reads, model){
  #Where do the xlinks land in the regions?
  key <- findOverlaps(reads, model, ignore.strand = F)
  #Take regions that have hits.
  locs <- model[subjectHits(key)]
  #Assign crosslink positions to regions
  locs$xlink <- start(ranges(reads))[queryHits(key)]
  #Assign metacolumns from iCLIP to location list.
  for(i in names(mcols(reads))){
    mcols(locs)[i] <- mcols(reads)[queryHits(key), i]
  }
  return(locs)
}

lib_size <- function(bf){
  #Takes a GRanges from bamloader/bedloader and returns a df of library sizes for each sample.
  outdf <- data.frame(sample = names(table(bf$sample)),
                      xlinks = as.vector(unlist(unname(table(bf$sample)))))
  outdf$xlinks <- outdf$xlinks/1000000
  return(outdf)
}

icount_bedloader_tome <- function(directory){
  #Take a directory containing .bed files from iMaps's xlinks sites.
  bedfiles <- list.files(directory, pattern = ".bed", full.names = T)
  bednames <- list.files(directory, pattern = ".bed", full.names = F)
  bed <- pblapply(bedfiles, cl = 6, function(x){
    fread(x)
  })
  bedlist <- pblapply(bed, cl = 6, function(x){
    #For negative stranded reads, add the width of the read (sum the cigar number) to the start position
    #to get the actual 5' position of the read.
    x$start[x$V6 == "+"] <- x$V2[x$V6 == "+"] + 1 #BED starts are 0-based.
    x$start[x$V6 == "-"] <- x$V3[x$V6 == "-"] #But BED end positions are 1-based.
    return(
      GRanges(seqnames =  x$V1, #gsub("\\..*","", x$V1)
              ranges = IRanges(start = x$start, width = 1),
              strand = x$V6,
              count = x$V5))
  })
  names(bedlist) <- lapply(bednames, function(x){
    strsplit(x, ".", fixed = T)[[1]][1]
  })
  bedgr <- unlist(GRangesList(bedlist))
  bedgr$sample <- names(bedgr)
  return(bedgr)
}



# Loading mouse annotations ------------------------------------------------------------

mm_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_mm_details.txt")

mm_cds_gr <- tibble(seqnames = mm_details$transcript_id,
                    start = mm_details$cds_start,
                    end = mm_details$cds_end,
                    strand = "+") %>%
  makeGRangesFromDataFrame()


read_tsv("GASR/germs/data/mouse/germs_CDS_umap_clusters.tsv.gz") %>%
  separate(peak_identifier, into = c("transcript_id", "start", "end", "region"), sep = ":|-|@") ->
  germs_clusters_mouse

germs_clusters_mouse %>%
  mutate(cluster_name_short = case_when(cluster == 0 ~ "None",
                                        cluster == 1 ~ "C-rich",
                                        cluster == 2 ~ "GC-rich",
                                        cluster == 3 ~ "CUG-rep",
                                        cluster == 4 ~ "CAG-rep",
                                        cluster == 5 ~ "Purine-rich")) %>%
  makeGRangesFromDataFrame(seqnames.field = "transcript_id", keep.extra.columns = T) ->
  germs_gr

# Loading my iCLIP data ---------------------------------------------------

mesc_mito <- icount_bedloader_tome("GASR/tmap/mesc_mitosis_1/results/xlinks_tome/")
sams <- mesc_mito$sample %>% Rle()

mesc_mito$target <- sams@values %>% str_split("_", simplify = T) %>% .[,1] %>% rep(sams@lengths)
mesc_mito$phase <- sams@values %>% str_split("_", simplify = T) %>% .[,2] %>% rep(sams@lengths)
mesc_mito$fraction <- sams@values %>% str_split("_", simplify = T) %>% .[,3] %>% rep(sams@lengths)

mesc_mito %>% filter(phase == "Asynchronous") %>% filter(!(target %in% c("SmB", "Prpf8"))) ->
  asynch_clip

# Plot --------------------------------------------------------------------

# Tra2b
metaprofile(asynch_clip %>% filter(target == "Tra2b"),
            germs_gr %>% filter(cluster > 0),
            400, 51, "cluster_name_short") %>%
  mutate(cluster_name_short = cluster_name_short %>% fct_relevel("Purine-rich")) ->
  tra2b_maps

ggplot(tra2b_maps,
       aes(x = pos)) +
  geom_line(aes(y = mean_count, 
                colour = cluster_name_short)) +
  geom_ribbon(aes(ymin = se_bottom,
                  ymax = se_top,
                  fill = cluster_name_short),
              alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(y = "CPM per region",
       x = "distance from region midpoint",
       title = "Tra2b") +
  scale_fill_manual(values = c("red3", "blue4", "green4", "mediumpurple3", "orange1")) +
  scale_color_manual(values = c("red3", "blue4", "green4", "mediumpurple3", "orange1")) +
  theme(legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.8)) ->
  tra2b_plot


ggsave("GASR/germs/plots/mouse/tra2b_in_germs_clusters.pdf",
       tra2b_plot, width = 3, height = 2.5)
