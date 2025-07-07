rm(list = ls())

library(tidyverse)
library(broom)
library(furrr)
plan(multisession, workers = 2)
library(roll)
library(patchwork)

# Nuclear retention -------------------------------------------------------

er_int <- read_tsv("GASR/export_reporter/quantseq/neve_timecourse_1/tables/publication_summary_data_interaction.tsv")

er_int %>%
  filter(padj < 0.05, log2FoldChange > 0) ->
  retained

# Disorder ----------------------------------------------------------------

gencode_to_swissprot <- read_tsv("GASR/lists/gencode.v29.metadata.SwissProt.txt",
                                 col_names = c("transcript_id", "uniprot_id", "id2"),
                                 col_type = "c")

af2_disorder <- read_tsv("General/alpha_fold_disorder/af2_human_disorder_pred.tsv.gz",
                         skip = 1,
                         col_names = c("name", "pos", "aa", "lddt", "disorder", "rsa", "ss", "disorder_25", "binding_25_0581"),
                         col_types = "cdcdddcdd")

transcript_details <- read_tsv("GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols())

transcript_details_uniprot <- inner_join(gencode_to_swissprot, transcript_details)

transcript_af2_disorder <- af2_disorder %>%
  group_by(name) %>%
  nest(data = c(pos, aa, lddt, disorder, rsa, ss, disorder_25, binding_25_0581)) %>%
  ungroup() %>%
  mutate(uniprot_id = word(name, 2, sep = fixed("-"))) %>%
  dplyr::select(-name) %>%
  inner_join(transcript_details_uniprot)

#Does the length of the alpha fold vector correctly match the length of the CDS? If not, discard.
transcript_af2_disorder <- transcript_af2_disorder %>%
  ungroup() %>%
  mutate(disorder_nrow = map(data, ~{nrow(.x)}) %>% unlist()) %>%
  filter((disorder_nrow + 1) == (cds_length/3)) %>%
  dplyr::select(-disorder_nrow)

# Choose threshold --------------------------------------------------------

transcript_af2_disorder %>%
  filter(transcript_id %in% er_int$transcript_id) %>%
  mutate(gene_class = case_when(transcript_id %in% retained$transcript_id ~ "retained",
                                T ~ "control") %>%
           fct_relevel("control")) ->
  af2_disorder_erintgenes
  
# 0.5 seems like a nice cutoff for plDDT
af2_disorder_erintgenes %>%
  unnest(data) %>%
  ggplot(aes(x = lddt)) +
  geom_density()


# IDRs --------------------------------------------------------------------

idr_smooth_width = 40

af2_disorder_erintgenes %>%
  # .[1:20,] %>%
  mutate(idr_info = future_map(data, ~{
    # af2_disorder_erintgenes$data[[3]] -> .x
    
    roll_mean(.x$lddt, width = idr_smooth_width) -> rolled_lddt
    
    # Set NA to 1 (excluded from being disordered)
    rolled_lddt[is.na(rolled_lddt)] <- 1
    
    rle(rolled_lddt <= 0.5) -> idrs
    
    cumsum(idrs$lengths)
    
    c()
    
    if(sum(idrs$values) == 0) {
      tibble(max_idr_length = 0,
             total_idr_length = 0,
             percentage_protein_disordered = 0) %>%
        return()
    } else {
      # A single positive window means idr_smooth_width window had mean of 0.5 pLDDT - to account I am adding the window back
      idrs$lengths[idrs$values] + (idr_smooth_width/2) -> idr_lengths
      
      tibble(max_idr_length = max(idr_lengths),
             total_idr_length = sum(idr_lengths),
             percentage_protein_disordered = total_idr_length / length(rolled_lddt)) %>%
        return()
    }
    
    
  })) %>%
  unnest(idr_info) ->
  protein_disorder_info

transcript_af2_disorder %>%
  mutate(idr_info = map(data, ~{
    
    roll_mean(.x$lddt, width = idr_smooth_width) -> rolled_lddt
    
    # Remove the NAs created by the smoothing window
    rolled_lddt[!is.na(rolled_lddt)] -> rolled_lddt
    
    rle(rolled_lddt <= 0.5) -> idrs
    
    cumsum(idrs$lengths) + (idr_smooth_width - 1) -> region_ends
    
    c(1, 1 + region_ends[-length(region_ends)]) -> region_starts
    
    
    
    if(sum(idrs$values) == 0) {
      tibble(idr_start = NA,
             idr_end = NA) %>%
        return()
    } else {
      tibble(idr_start = region_starts[idrs$values],
             idr_end = region_ends[idrs$values]) %>%
        return()
    }
    
    
  })) %>%
  unnest(idr_info) ->
  protein_disordered_regions

write_tsv(protein_disordered_regions %>% select(-c(data, cds_start, cds_length, tx_length, cds_end)), "GASR/germs/data/alphafold_disorder/table_of_locations_for_jure.txt.gz")

af2_disorder_erintgenes %>%
  mutate(percentage_disordered = future_map(data, ~{
    sum(.x$lddt <= 0.5) / length(.x$lddt)
  }) %>%
    unlist()) ->
  percentage_protein_disordered

percentage_protein_disordered %>%
  ggplot(aes(x = gene_class, y = percentage_disordered, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "proportion protein disordered\n(alphafold pLDDT < 0.5)",
       color = "") ->
  proportion_protein_disordered_plot

protein_disorder_info %>%
  ggplot(aes(x = gene_class, y = percentage_protein_disordered, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "proportion protein disordered",
       color = "") ->
  proportion_protein_disordered_plot

protein_disorder_info %>%
  ggplot(aes(x = gene_class, y = max_idr_length, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  coord_cartesian(ylim = c(0,900)) +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "longest disordered region",
       color = "") ->
  longest_disordered_plot
  
t.test(percentage_protein_disordered ~ gene_class, data = protein_disorder_info) %>% tidy() -> prop_disordered_ttest
t.test(max_idr_length ~ gene_class, data = protein_disorder_info) %>% tidy() -> idr_length_ttest

write_tsv(prop_disordered_ttest, "GASR/germs/plots/protein_features/prop_disordered_ttest.tsv")

percentage_protein_disordered %>%
  select(-data) ->
  percentage_protein_disordered

write_tsv(percentage_protein_disordered, "GASR/germs/plots/protein_features/percentage_protein_disordered_sourcedata.tsv")

# Proportion charged ------------------------------------------------------

library(Biostrings)
library(roll)

prots <- readAAStringSet("GASR/lists/longest_gencode29_prots.fa")

aa_window_size = 40
charge_threshold = 0.4

tibble(transcript_id = names(prots),
       protein = as.character(prots)) %>%
  filter(transcript_id %in% er_int$transcript_id) %>%
  # sample_frac(0.02) %>%
  mutate(prop_charged = future_map(protein, 
                            ~{
                              # .x = as.character(prots$ENST00000568838.1)
                              zeros <- rep(0, nchar(.x))
                              
                              one_locs <- .x %>% str_locate_all("K|E|D|R") %>% as.data.frame() %>% select(start) %>% unlist(use.names = F)
                              
                              zeros[one_locs] <- 1
                              
                              prop_charged <- roll_mean(zeros, aa_window_size)[-c(1:(aa_window_size-1))]
                              
                              return(sum(prop_charged > charge_threshold) / length(prop_charged))
                            }) %>% unlist(),
         
         gene_class = case_when(transcript_id %in% retained$transcript_id ~ "retained",
                                       T ~ "control") %>%
                  fct_relevel("control")) ->
  prop_charged

prop_charged %>%
  ggplot(aes(x = gene_class, y = prop_charged, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "proportion protein highly charged",
       color = "") ->
  proportion_charged_plot

t.test(prop_charged ~ gene_class, data = prop_charged) %>% tidy() -> prop_charged_ttest

write_tsv(prop_charged_ttest, "GASR/germs/plots/protein_features/prop_charged_ttest.tsv")

prop_charged %>%
  select(-protein) ->
  prop_charged

write_tsv(prop_charged, "GASR/germs/plots/protein_features/prop_charged_sourcedata.tsv")

aa_window_size = 40
charge_threshold = 0.4

tibble(transcript_id = names(prots),
       protein = as.character(prots)) %>%
  filter(transcript_id %in% er_int$transcript_id) %>%
  # sample_frac(0.02) %>%
  mutate(prop_polarnoncharged = future_map(protein, 
                                   ~{
                                     # .x = as.character(prots$ENST00000568838.1)
                                     zeros <- rep(0, nchar(.x))
                                     
                                     one_locs <- .x %>% str_locate_all("S|T|N|Q") %>% as.data.frame() %>% select(start) %>% unlist(use.names = F)
                                     
                                     zeros[one_locs] <- 1
                                     
                                     prop_polarnoncharged <- roll_mean(zeros, aa_window_size)[-c(1:(aa_window_size-1))]
                                     
                                     return(sum(prop_polarnoncharged > charge_threshold) / length(prop_polarnoncharged))
                                   }) %>% unlist(),
         
         gene_class = case_when(transcript_id %in% retained$transcript_id ~ "retained",
                                T ~ "control") %>%
           fct_relevel("control")) ->
  prop_polarnoncharged

prop_polarnoncharged %>%
  ggplot(aes(x = gene_class, y = prop_polarnoncharged, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "proportion protein highly polar (non-charged)",
       color = "") ->
  prop_polarnoncharged_plot

t.test(prop_polarnoncharged ~ gene_class, data = prop_polarnoncharged) %>% tidy() -> prop_polarnoncharged_ttest

write_tsv(prop_polarnoncharged_ttest, "GASR/germs/plots/protein_features/prop_polarnoncharged_ttest.tsv")

prop_polarnoncharged %>%
  select(-protein) ->
  prop_polarnoncharged

write_tsv(prop_polarnoncharged, "GASR/germs/plots/protein_features/prop_polarnoncharged_sourcedata.tsv")

(proportion_charged_plot | proportion_protein_disordered_plot | longest_disordered_plot) + plot_layout(guides = "collect") ->
  combined_charge_disorder_plot

ggsave("GASR/germs/plots/protein_features/combined_charge_disorder_plot.pdf", combined_charge_disorder_plot,
       width = 5, height = 3)


ggsave("GASR/germs/plots/protein_features/prop_polarnoncharged_plot.pdf", prop_polarnoncharged_plot,
       width = 1.66, height = 3)


# Jure's stuff ------------------------------------------------------------

jure_table <- read_tsv("GASR/germs/plots/protein_features/P03_H01_A02_mastertable_fil.tsv") %>%
  mutate(gene_class = case_when(Groups_analysis == "interstatic" ~ "retained",
                                T ~ "control"))

jure_table %>% names()


jure_table %>%
  ggplot(aes(x = gene_class, y = FuzDrop, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "phase separation propensity\n(FuzDrop score)",
       color = "") ->
  fuzdrop_plot

t.test(FuzDrop ~ gene_class, data = jure_table, na.action = "na.omit") %>% tidy() -> fuzdrop_ttest

write_tsv(fuzdrop_ttest, "GASR/germs/plots/protein_features/fuzdrop_ttest.tsv")

jure_table %>%
  ggplot(aes(x = gene_class, y = DosPS, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "phase separation propensity\n(DosPS score)",
       color = "") ->
  dosps_plot

t.test(DosPS ~ gene_class, data = jure_table, na.action = "na.omit") %>% tidy() -> dosps_ttest

write_tsv(dosps_ttest, "GASR/germs/plots/protein_features/dosps_ttest.tsv")

fuzdrop_plot | dosps_plot ->
  phase_separation_combiplot

ggsave("GASR/germs/plots/protein_features/phase_separation_combiplot.pdf", phase_separation_combiplot,
       width = 3.75, height = 3)

jure_table %>%
  ggplot(aes(x = gene_class, y = pHaplo, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "dosage sensitivity in deletion\n(haploinsufficiency)",
       color = "") ->
  phaplo_plot

t.test(pHaplo ~ gene_class, data = jure_table, na.action = "na.omit") %>% tidy() -> phaplo_ttest

write_tsv(phaplo_ttest, "GASR/germs/plots/protein_features/phaplo_ttest.tsv")

jure_table %>%
  ggplot(aes(x = gene_class, y = pTriplo, color = gene_class)) +
  stat_boxplot(geom ='errorbar', width = 0.2, show.legend = F) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "dosage sensitivity in duplication\n(triplosensitivity)",
       color = "") ->
  ptriplo_plot

t.test(pTriplo ~ gene_class, data = jure_table, na.action = "na.omit") %>% tidy() -> ptriplo_ttest

write_tsv(ptriplo_ttest, "GASR/germs/plots/protein_features/ptriplo_ttest.tsv")


phaplo_plot | ptriplo_plot ->
  dosage_sensitivity_plot

ggsave("GASR/germs/plots/protein_features/dosage_sensitivity_plot.pdf", dosage_sensitivity_plot,
       width = 3.75, height = 3)

jure_table %>%
  ggplot(aes(x = gnomAD_pLI, color = gene_class)) +
  geom_density() +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(color = "",
       # y = "dosage sensitivity in duplication\n(triplosensitivity)",
       x = "probability of loss-intolerance\n(gnomAD)") ->
  gnomad_density_plot

jure_table %>%
  dplyr::select(gene_class, gnomAD_pLI) %>%
  drop_na() %>%
  group_by(gene_class) %>%
  summarise(proportion_high_gnomaD = sum(gnomAD_pLI > 0.5)/dplyr::n()) %>%
  ggplot(aes(x = gene_class, y = proportion_high_gnomaD, color = gene_class)) +
  geom_bar(stat = "identity", fill = NA) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("#B4B4B4", "#B9529F")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs(x = "",
       y = "proportion high loss-intolerance\n(gnomAD)",
       color = "") ->
  high_gnomad_plot

(gnomad_density_plot | high_gnomad_plot) + plot_layout(guides = "collect", widths = c(2,1)) ->
  gnomad_combi_plot
  
ggsave("GASR/germs/plots/protein_features/gnomad_combi_plot.pdf", gnomad_combi_plot,
       width = 6, height = 3)
