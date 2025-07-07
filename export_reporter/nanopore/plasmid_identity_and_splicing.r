rm(list = ls())

library(tidyverse)
library(Biostrings)
library(plyranges)
library(stringdist)
library(patchwork)
library(furrr)

# Load genome--------------------------------------------------------------------

reference_genome <- readDNAStringSet("GASR/export_reporter/nanopore/minimap2/all_reporter_sequences_with_introns.fa")

intron_locations <- read_csv(file = "GASR/export_reporter/best_introns/output/reporter.high_purine.output.csv")[4:19,] %>%
  separate(key, into = c("intron", "intron_number", "type"), convert = T) %>%
  # For some reason I missed a couple nucleotides at the start of this file (also 0 indexing convert to 1)
  mutate(value = as.numeric(value) + 3) %>%
  # mutate(value = as.numeric(value) + 21) %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::select(-intron)

intron_locations %>%
  mutate(seqnames = "reporter") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) ->
  introns_gr

exons_gr <- setdiff_ranges(GRanges(seqnames = "reporter", ranges = IRanges(start = 2, end = 3141)), introns_gr) %>%
    # This seems to make the ranges correct - something to do with how setdiff is working.
    anchor_3p() %>% stretch(1)

introns_gr <- introns_gr %>%
  # This seems to make the ranges correct - something to do with how setdiff is working.
  anchor_5p() %>% stretch(1)

# Load genomic reads ------------------------------------------------------
# 
# dna_reads <- read_tsv("GASR/export_reporter/nanopore/minimap2/plasmid/plasmid_alignment.filtered.table.tsv",
#                       col_names = c("read_id", "flag", "opt_pattern", "start", "cigar", "sequence")) %>%
#   filter(flag <= 16)
# 
# # dna_reads_s <- dna_reads %>% sample_frac(0.01)
# 
# # dna_reads %>% filter(read_id == "fd2e2dd4-7704-4dfc-8869-846af733ced7") %>% .["sequence"] %>% unlist()
# # dna_reads %>% filter(read_id == "fd2e2dd4-7704-4dfc-8869-846af733ced7") %>% .["cigar"] %>% unlist()
# 
# # Barcode extraction ------------------------------------------------------
# 
# # Cigar elements that consume reference:
# # M, D, N, =, X
# 
# us_seq = toupper("attacgtctcgatcgacaac")
# ds_seq = "GGATCCGCAGGCCTCTGCTA"
# 
# dna_reads %>%
#   mutate(cigar_parsed = map(cigar, ~{
#     tibble(cig_length = .x %>%
#              str_split("[:letter:]") %>% 
#              unlist() %>% .[c(1:(length(.) - 1))] %>%
#              as.numeric(),
#            cig_type = .x %>% str_extract_all("[:letter:]") %>% unlist()) ->
#       cig_table
#   }),
#   alignment_length = map(cigar_parsed, ~{
#     .x %>%
#       filter(cig_type %in% c("M", "D", "N", "=", "X")) %>%
#       summarise(l = sum(cig_length)) %>%
#       unlist(use.names = F)
#   }) %>%
#     unlist(),
#   end = start + alignment_length - 1) ->
#   dna_reads_parsed
# 
# dna_reads_parsed %>%
#   ggplot(aes(x = end)) +
#   geom_density()
# 
# dna_reads_parsed %>% 
#   # (almost) full length reads only
#   filter(start <= 175, end >= 3111) %>%
#   mutate(last_200 = str_sub(sequence, nchar(sequence) - 200, nchar(sequence))) %>%
#   mutate(us_find = map(last_200, ~{
#     afind(.x, us_seq) %>% as.data.frame() %>% set_names("us_loc", "us_dist", "us_match")
#   }),
#   ds_find = map(last_200, ~{
#     afind(.x, ds_seq) %>% as.data.frame() %>% set_names("ds_loc", "ds_dist", "ds_match")
#   })) %>%
#   unnest(c(us_find, ds_find)) %>%
#   filter(us_dist <= 4, ds_dist <= 4) %>%
#   mutate(barcode = str_sub(last_200, us_loc + 20, ds_loc - 1),
#          barcode_length = nchar(barcode)) ->
#   dna_barcode_str_matches
# 
# # DNA Intron detection --------------------------------------------------------
# 
# dna_barcode_str_matches %>%
#   mutate(junction_positions = pmap(.l = list(cigar_parsed = cigar_parsed,
#                                              read_start = start),
#                                    .f = function(cigar_parsed, read_start) {
#                                      cigar_parsed %>%
#                                        # Take the cigar types that consume the reference
#                                      filter(cig_type %in% c("M", "D", "N", "=", "X")) %>%
#                                        mutate(is_n = cig_type == "N") %>%
#                                        # Stretches that run between junctions are 'exons' in the read.
#                                        mutate(region = rep(c(1:(is_n %>% rle() %>% .["values"] %>% unlist() %>% length())), 
#                                                            times = is_n %>% rle() %>% .["lengths"] %>% unlist())) %>%
#                                        group_by(region, is_n) %>%
#                                        reframe(region_length = sum(cig_length)) %>%
#                                        mutate(start = c(1,(cumsum(region_length)[-dplyr::n()] + 1)),
#                                               end = cumsum(region_length)) %>%
#                                        mutate(start = start + read_start,
#                                               end = end + read_start)
#                                    })) %>%
#   mutate(intron_pattern = map(junction_positions, ~{
#     
#     # If there are no splice junctions, every intron is present
#     if(sum(.x$is_n) == 0) {
#       return(paste(rep(1, 8), collapse = "-"))
#     }
#     
#     .x %>%
#       filter(is_n) %>%
#       mutate(seqnames = "reporter") %>%
#       makeGRangesFromDataFrame() ->
#       read_junctions
#     
#     # Check if there are junctions that span over an exon for some reason
#     # This would suggest something weird happened during assembly.
#     if(sum(count_overlaps(exons_gr,
#                           read_junctions,
#                           minoverlap = 50)) > 0) {
#       return("misassembled")
#     }
#     
#     # Which reporter junctions overlap with the introns in the read
#     count_overlaps(introns_gr,
#                    read_junctions,
#                    minoverlap = 50) ->
#       junction_overlaps
#     
#     # If the intron is spliced out of the DNA read, then it isn't present in the reporter.
#     paste((-junction_overlaps + 1), collapse = "-") %>%
#       return()
#     
#   }) %>% unlist()) ->
#   dna_barcode_junctions
# 
# dna_barcode_junctions %>%
#   filter(intron_pattern != "misassembled",
#          # BC length 17-23
#          barcode_length > 16, barcode_length < 24) %>%
#   select(read_id, flag, opt_pattern, intron_pattern, barcode) ->
#   dna_barcode_junctions_short
#   
# write_tsv(dna_barcode_junctions_short, "GASR/export_reporter/nanopore/minimap2/processed_data/plasmid_barcodes_and_architecture.tsv.gz")

dna_barcode_junctions_short <- read_tsv("GASR/export_reporter/nanopore/minimap2/processed_data/plasmid_barcodes_and_architecture.tsv.gz")

dna_barcode_junctions_short %>%
  group_by(barcode, opt_pattern, intron_pattern) %>%
  summarise(number_reads = dplyr::n()) %>%
  group_by(barcode) %>%
  mutate(barcode_count = sum(number_reads),
         barcode_types = dplyr::n()) %>%
  ungroup() %>%
  # filter(barcode_count > 2) %>%
  arrange(barcode, number_reads) %>%
  group_by(barcode) %>%
  mutate(cumsum_reads = cumsum(number_reads)) %>%
  dplyr::slice(1L) %>%
  filter(cumsum_reads >= 0.75 * barcode_count) %>%
  ungroup() ->
  consensus_barcodes

write_tsv(consensus_barcodes, "GASR/export_reporter/nanopore/minimap2/processed_data/plasmid_consensus_barcodes.tsv.gz")

# Plotting representation  ---------------------------------------------------------------

dna_barcode_junctions_short %>%
  mutate(arch_df = future_pmap(.l = list(opt_pattern = opt_pattern, intron_pattern = intron_pattern),
                              .f = function(opt_pattern, intron_pattern){

                                # opt_pattern = dna_barcode_junctions_short$opt_pattern[1]
                                # intron_pattern = dna_barcode_junctions_short$intron_pattern[1]
                                
                                opt_pattern %>%
                                  str_split("-") %>%
                                  unlist() %>%
                                  str_detect("Opt") ->
                                  opt_presence
                                
                                intron_pattern %>%
                                  str_split("-") %>%
                                  unlist() %>% 
                                  as.numeric() %>% 
                                  as.logical() ->
                                  intron_presence
                                
                                tibble(chunk = c(1:8),
                                       opt = opt_presence,
                                       int = intron_presence) %>%
                                  return()
                                
                              })) %>%
  dplyr::select(read_id, arch_df) %>%
  unnest(arch_df) %>%
  group_by(chunk) %>%
  count(opt, int) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(opt = case_when(opt ~ "GA-rich",
                         T ~ "GA-poor") %>%
           fct_relevel("GA-poor"),
         int = case_when(int ~ "intron",
                         T ~ "intronless") %>%
           fct_relevel("intron")) ->
  chunk_proportion
  

ggplot(chunk_proportion, aes(x = int, y = prop, fill = opt)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(. ~ chunk, nrow = 1) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black",)) +
  scale_fill_manual(values = c("#B4B4B4", "#B9529F")) +
  labs(x = "", y = "proportion of plasmids", fill = "")
  
chunk_proportion %>%
  mutate(opt_int = paste0(opt, ", ", int) %>%
           fct_relevel("GA-poor, intron", "GA-poor, intronless", 
                       "GA-rich, intron", "GA-rich, intronless")) ->
  stacked_chunk_proportion

ggplot(stacked_chunk_proportion, aes(x = "hi", y = prop, fill = opt_int)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(. ~ chunk, nrow = 1) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black",),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#B4B4B4", "#878787", "#B9529F", "#7d316a")) +
  labs(x = "", y = "proportion of plasmids", fill = "") ->
  stacked_chunk_prop_plot

# RNA barcodes----------------------------------------------------------


# Load RNA reads --------------------------------------------------------------

rna_reads <- read_tsv("GASR/export_reporter/nanopore/minimap2/alignment.sorted.filtered.table.tsv",
                      col_names = c("read_id", "flag", "opt_pattern", "start", "cigar", "sequence")) 

# rna_reads_s <- rna_reads %>% filter(opt_pattern == "Opt-Opt-Opt-Opt-Opt-Opt-Opt-Opt")

# Barcode extraction ------------------------------------------------------

# Cigar elements that consume reference:
# M, D, N, =, X

us_seq = toupper("attacgtctcgatcgacaac")
ds_seq = "GGATCCGCAGGCCTCTGCTA"

library(tictoc)

plan(multisession, workers = 4)

tic("Parsing RNA cigars")
rna_reads %>%
  # slice(1:20000) %>%
  mutate(cigar_parsed = future_map(cigar, ~{
    tibble(cig_length = .x %>%
             str_split("[:letter:]") %>% 
             unlist() %>% .[c(1:(length(.) - 1))] %>%
             as.numeric(),
           cig_type = .x %>% str_extract_all("[:letter:]") %>% unlist()) ->
      cig_table
  }),
  alignment_length = future_map(cigar_parsed, ~{
    .x %>%
      filter(cig_type %in% c("M", "D", "N", "=", "X")) %>%
      summarise(l = sum(cig_length)) %>%
      unlist(use.names = F)
  }) %>%
    unlist(),
  end = start + alignment_length - 1) ->
  rna_reads_parsed
toc()

# rna_reads_parsed %>%
#   ggplot(aes(x = start)) +
#   geom_density()

tic("Parsing RNA barcodes")

rna_reads_parsed %>% 
  # (almost) full length reads only
  filter(start <= 175, end >= 3111) %>%
  mutate(last_200 = str_sub(sequence, nchar(sequence) - 200, nchar(sequence))) %>%
  mutate(us_find = future_map(last_200, ~{
    afind(.x, us_seq) %>% as.data.frame() %>% set_names("us_loc", "us_dist", "us_match")
  }),
  ds_find = future_map(last_200, ~{
    afind(.x, ds_seq) %>% as.data.frame() %>% set_names("ds_loc", "ds_dist", "ds_match")
  })) %>%
  unnest(c(us_find, ds_find)) %>%
  filter(us_dist <= 4, ds_dist <= 4) %>%
  mutate(barcode_rna = str_sub(last_200, us_loc + 20, ds_loc - 1),
         barcode_rna_length = nchar(barcode_rna)) %>%
  dplyr::rename(opt_pattern_rna = opt_pattern) ->
  rna_barcode_str_matches

toc()

# rna_barcode_str_matches %>%
#   ggplot(aes(x = barcode_rna_length)) +
#   geom_density()

tic("Matching RNA barcodes")
rna_barcode_str_matches %>%
  dplyr::count(barcode_rna) %>%
  arrange(desc(n)) %>%
  mutate(match_index = future_map(barcode_rna, ~{
    
    # .x <- rna_barcode_str_matches$barcode[40]
    # message(paste0(barcode, "\n"))
    
    distances <- stringdist(.x, consensus_barcodes$barcode)
    
    matches <- which((distances == min(distances)) & (distances < 4))
    
    if(length(matches) == 0) {
      return(NA)
    }
    
    if(length(matches) == 1) {
      return(matches)
    }
    
    # If there are multiple matches, I check to see if there is a consensus between those matches.
    # I only take the barcode matches where there is a consensus for gene structure.
    
    architectures <- tibble(opt_pattern = consensus_barcodes$opt_pattern[matches], 
                            intron_pattern = consensus_barcodes$intron_pattern[matches],
                            match_index = matches)
    
    architectures %>%
      count(opt_pattern, intron_pattern) %>%
      mutate(proportion = n/sum(n)) %>%
      filter(proportion >= 0.66) ->
      architecture_proportions
    
    # No consensus, no service
    if(nrow(architecture_proportions) == 0) {
      return(NA)
    }
    
    architectures %>%
      filter(opt_pattern == architecture_proportions$opt_pattern,
             intron_pattern == architecture_proportions$intron_pattern) %>%
      dplyr::slice(1) %>%
      dplyr::select(match_index) %>%
      unlist(use.names = F) %>%
      return()
    
  }) %>%
    unlist()) ->
  rna_barcode_matches
toc()

# rna_barcode_matches<- rna_barcode_matches %>% dplyr::rename(barcode_rna = barcode)

rna_barcode_matches %>%
  drop_na() %>%
  mutate(opt_pattern_dna = consensus_barcodes$opt_pattern[match_index],
         int_pattern_dna = consensus_barcodes$intron_pattern[match_index]) ->
  rna_barcode_matches_arch

rna_barcode_str_matches %>%
  inner_join(rna_barcode_matches_arch, by = "barcode_rna") %>%
  # Make sure the RNA mapping matches the expected DNA mapping.
  # This is not true for ~5% of reads.
  filter(opt_pattern_rna == opt_pattern_dna) ->
  rna_reads_arch_annotated

write_rds(rna_reads_arch_annotated, "GASR/export_reporter/nanopore/minimap2/processed_data/rna_reads_architecture_annotated.rds")

# stringdist(rna_barcode_str_matches$barcode[1], dna_barcode_junctions_short$barcode)


# Is the whole transcript correctly spliced? ------------------------------

introns_gr %>%
  anchor_3p() %>%
  mutate(width = 0) %>%
  anchor_center() %>%
  mutate(width = 20) ->
  ss3_introns_gr

introns_gr %>%
  anchor_5p() %>%
  mutate(width = 0) %>%
  anchor_center() %>%
  mutate(width = 20) ->
  ss5_introns_gr

tic("Getting read junction positions")
rna_reads_arch_annotated %>%
  # slice(1:1000) %>%
  dplyr::select(read_id, cigar_parsed, start, opt_pattern_dna, int_pattern_dna) %>%
  mutate(junction_positions = future_pmap(.l = list(cigar_parsed = cigar_parsed,
                                             read_start = start),
                                   .f = function(cigar_parsed, read_start) {
                                     cigar_parsed %>%
                                       # Take the cigar types that consume the reference
                                       filter(cig_type %in% c("M", "D", "N", "=", "X")) %>%
                                       mutate(is_n = cig_type == "N") %>%
                                       # Stretches that run between junctions are 'exons' in the read.
                                       mutate(region = rep(c(1:(is_n %>% rle() %>% .["values"] %>% unlist() %>% length())), 
                                                           times = is_n %>% rle() %>% .["lengths"] %>% unlist())) %>%
                                       group_by(region, is_n) %>%
                                       reframe(region_length = sum(cig_length)) %>%
                                       mutate(start = c(1,(cumsum(region_length)[-dplyr::n()] + 1)),
                                              end = cumsum(region_length)) %>%
                                       mutate(start = start + (read_start - 1),
                                              end = end + (read_start - 1))
                                   })) -> 
  rna_reads_junction_positions
toc()

write_rds(rna_reads_junction_positions, "GASR/export_reporter/nanopore/minimap2/processed_data/rna_reads_junction_positions.rds")

tic("Checking read splicing")
rna_reads_junction_positions %>%
  # slice_sample(n = 10000) %>% 
  mutate(is_it_well_spliced = future_pmap(.l = list(junction_positions = junction_positions,
                                             opt_pattern_dna = opt_pattern_dna,
                                             int_pattern_dna = int_pattern_dna),
                                   .f = function(junction_positions, opt_pattern_dna, int_pattern_dna) {
                                     
                                     # st <- Sys.time()
                                     
                                     # junction_positions <- rna_reads_junction_positions$junction_positions[[500]]
                                     # opt_pattern_dna <- rna_reads_junction_positions$opt_pattern_dna[[500]]
                                     # int_pattern_dna <- rna_reads_junction_positions$int_pattern_dna[[500]]
                                     
                                     opt_pattern_dna %>%
                                       str_split("-") %>%
                                       unlist() %>%
                                       str_detect("Opt") ->
                                       opt_presence
                                     
                                     int_pattern_dna %>%
                                       str_split("-") %>%
                                       unlist() %>% 
                                       as.numeric() %>% 
                                       as.logical() ->
                                       intron_presence
                                     
                                     if(nrow(junction_positions) == 1) {
                                       # To deal with very rare cases where the are no introns,
                                       # I am making a fake intron in the middle of nowhere so I don't have to change my code
                                       read_introns <- GRanges(seqnames = "reporter",
                                                               ranges = IRanges(start = 50000, end = 50100))
                                     } else {
                                       junction_positions %>%
                                         filter(is_n) %>%
                                         mutate(seqnames = "reporter") %>%
                                         makeGRangesFromDataFrame() ->
                                         read_introns
                                     }
                                     
                                     junction_positions %>%
                                       filter(!is_n) %>%
                                       mutate(seqnames = "reporter") %>%
                                       makeGRangesFromDataFrame() ->
                                       read_exons
                                     
                                     # reg_def_time = as.numeric(Sys.time() - st)
                                     # 
                                     # st <- Sys.time()
                                     
                                     # Intron is spliced in at all?
                                     count_overlaps(introns_gr, read_introns, minoverlap = 20) ->
                                       introns_spliced
                                     
                                     # Exon is spliced in? (at least a bit)
                                     count_overlaps(exons_gr, read_exons, minoverlap = 50) ->
                                       exons_spliced
                                     
                                     read_5ss = read_introns
                                     read_3ss = read_introns
                                     end(read_5ss) = start(read_5ss)
                                     start(read_3ss) = end(read_3ss)
                                     
                                     # Is the 5'SS correct
                                     count_overlaps(ss5_introns_gr, read_5ss) ->
                                       ss5_correct
                                     
                                     # Is the 3'SS correct
                                     count_overlaps(ss3_introns_gr, read_3ss) ->
                                       ss3_correct
                                     
                                     # 
                                     # good_splice_time = as.numeric(Sys.time() - st)
                                     # 
                                     # st <- Sys.time()
                                     
                                     tibble(chunk = c(1:8),
                                            ir = introns_spliced == 0,
                                            a5ss = (introns_spliced == 1) & (ss5_correct == 0),
                                            a3ss = (introns_spliced == 1) & (ss3_correct == 0),
                                            correct = (ir | a5ss | a3ss) == FALSE,
                                            genomic_intron = intron_presence,
                                            genomic_opt = opt_presence) %>% 
                                       nest(intron_splicing = c(chunk, ir, a5ss, a3ss, correct, genomic_intron, genomic_opt)) %>%
                                       mutate(correctly_spliced = sum(introns_spliced, exons_spliced, ss3_correct, ss5_correct) == 33) ->
                                       read_splicing_output
                                     
                                     # make_table_time = as.numeric(Sys.time() - st)
                                     # 
                                     # tibble(reg_def_time = reg_def_time,
                                     #        good_splice_time = good_splice_time,
                                     #        make_table_time = make_table_time
                                     #        ) %>%
                                     #   return()
                                     
                                     return(read_splicing_output)
                                   })) %>%
  unnest(is_it_well_spliced) ->
  is_splicing_correct
toc()

write_rds(is_splicing_correct, "GASR/export_reporter/nanopore/minimap2/processed_data/rna_reads_splicing_correct.rds")

is_splicing_correct %>%
  dplyr::select(read_id, intron_splicing) %>%
  unnest(intron_splicing) %>% 
  filter(genomic_intron) %>%
  mutate(correct_5ss = !ir & !a5ss,
         correct_3ss = !ir & !a3ss) %>%
  dplyr::select(read_id, chunk, genomic_opt, correct_5ss, a5ss, ir, a3ss, correct_3ss) %>%
  pivot_longer(cols = c("correct_5ss", "a5ss", "ir", "a3ss", "correct_3ss"), names_to = "splice_type", values_to = "splice_presence") %>%
  mutate(splice_type = case_when(splice_type == "correct_5ss" ~ "5'SS",
                                 splice_type == "correct_3ss" ~ "3'SS",
                                 splice_type == "a5ss" ~ "A5'SS",
                                 splice_type == "a3ss" ~ "A3'SS",
                                 splice_type == "ir" ~ "IR") %>%
           fct_relevel("5'SS",
                       "A5'SS",
                       "IR",
                       "A3'SS",
                       "3'SS")) %>%
  group_by(chunk, genomic_opt) %>%
  mutate(number_reads = length(unique(read_id))) %>%
  group_by(chunk, genomic_opt, number_reads, splice_type) %>%
  reframe(count = sum(splice_presence)) %>%
  mutate(proportion_event = count / number_reads,
         genomic_opt = case_when(genomic_opt ~ "GA-rich",
                                 T ~ "GA-poor") %>%
           fct_relevel("GA-poor")) ->
  splicing_type_by_chunk

ggplot(splicing_type_by_chunk, 
       aes(x = splice_type, y = proportion_event, fill = genomic_opt)) +
  geom_point(shape = 21, size = 3, position = position_dodge2(width = 0.5)) +
  facet_wrap(. ~ chunk, nrow = 1) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", angle = 90, hjust = 1, vjust = 0)) +
  scale_fill_manual(values = c("#B4B4B4", "#B9529F")) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", alpha = 0.2) +
  labs(x = "", y = "proportion of reads\nwith splicing event", fill = "chunk\ntype") ->
  junction_usage_plot
  

is_splicing_correct %>%
  group_by(opt_pattern_dna) %>%
  summarise(percentage_correctly_spliced = sum(correctly_spliced) / dplyr::n()) ->
  perc_spliced_per_genetype

ggplot(perc_spliced_per_genetype, aes(x = percentage_correctly_spliced)) +
  geom_histogram(fill = "grey30") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black",)) +
  labs(x = "proportion of reads correctly spliced") ->
  proportion_reads_spliced_plot 

(proportion_reads_spliced_plot | junction_usage_plot) + plot_layout(widths = c(1, 5)) ->
  combined_splicing_plot

ggsave("GASR/export_reporter/nanopore/minimap2/plots/combined_splicing_from_longread_plot.pdf", combined_splicing_plot,
       width = 10, height = 3)

ggsave("GASR/export_reporter/nanopore/minimap2/plots/proportion_transcripts_correctly_spliced.pdf", proportion_reads_spliced_plot,
       width = 3, height = 3)

ggsave("GASR/export_reporter/nanopore/minimap2/plots/junction_usage_plot.pdf", junction_usage_plot,
       width = 8, height = 3)

(stacked_chunk_prop_plot / junction_usage_plot) ->
  chunk_types_and_splicing_plot

ggsave("GASR/export_reporter/nanopore/minimap2/plots/chunk_types_and_junction_usage_plot.pdf", chunk_types_and_splicing_plot,
       width = 10, height = 5)
