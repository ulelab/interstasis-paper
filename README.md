# interstasis-paper
A repository containing the code for Faraway et al. (2025).

The GeRM package, used in this work to score Generalised RNA Multivalency, is available at: [https://github.com/ulelab/germ](https://github.com/ulelab/germ).


## GeRM analysis
- create_germ_table.r
  - Analyses GeRM profiles across longest protein coding transcript per gene, calls peaks, clusters peaks using UMAP/HDBSCAN and plots clusters.
- example_profile_plots.r
  - Plots GeRM profiles and amino acid entropy as seen for example genes like LUC7L3, as well as examples of CLIP crosslinks found later in the paper.
- seq_logos.r
  - Plots the seq logos and conservation information for the example gene plots.
- alphafold_disorder.r
  - Analysis of the Alphafold predictions of protein regions encoded by CDS GeRM regions.
- entropy_in_germ_regions.r
  - Analysis of the amino acid entropy of protein regions encoded by CDS GeRM regions.
- shuffled_cds_multivalency.r
  - The influence of codon shuffling on GeRM scores in regions encoding low complexity domains.
- calculate_effects_of_synonymous_mutations.r
  - Calculating the degree to which synonymous mutations influence local GeRM score across the transcriptome (used for downstream analysis of conservation of GeRM).
- cluster_low_entropy.r
  - Defining low complexity domains (as well as clustering them based on their amino acid identity, which isn’t really included in the paper, even though it was kind of cool).
- conservation_of_multivalency.r
  - Relationship between conservation of mutable positions and the degree to which they contribute to multivalency.
- calculate_coding_properties.r
  - Calculate the proportion of each amino acid encoded within each CDS GeRM region and replot UMAP, as well as supplementary heatmap.
- germs_cluster_ontologies.r
  - Gene ontology analysis of the genes with different types of GeRM region. Plotting of ontology heatmaps.
- arginine_codons.r
  - Analysis of the arginine codons within different R-LCDs and the gene ontologies of R-LCDs with different codon biases.
- analyse_rlcd_codon_usage_across_species.r
  - Analyse the correlation of arginine codon usage within R-LCDs for many different species.
- calculate_entropy_for_multispecies.r
  - Calculate the amino acid entropy and call LCDs for a given species.
- launch_jobs_for_multispecies.r
  - Launch entropy calculation jobs for all species.
- proportion_of_mv_per_kmer_per_cluster.r
  - For each CDS GeRM cluster identified, what are the k-mers that contribute the most to the total GeRM score in the clusters.
- er_interaction_retained_protein_properties.r
  - Analysis of the properties of the proteins encoded by transcripts that are retained in the nucleus following PPIG expression.
- are_germ_exons_spliced.r
  - Analysis of where exons that host GeRM sequences are subject to alternative splicing after TRA2 knockdown.
- multivalency_preference_of_clip_xlinks.r
  - Assessing whether CLIP crosslinks are more likely to fall in k-mers when they have relatively high or low GeRM scores. Used to identify RBPs that prefer to bind to motifs in multivalent contexts (these datasets are subsequently used to look at binding in GeRM - subtypes).
- crosslinks_within_germs_clusters.r
  - Looking at enrichments of CLIP crosslinks within different GeRM regions relative to the rest of the CDS, then plotting heatmaps.
- WHERE IS THE TRA2B iCLIP ANALYSIS?
- WHERE IS THE VASTDB CODE?
## Transcriptomics
- deseq_ppig_3prime_end.r
  - Perform DESeq2 analysis on the 3’ end data from PPIG overexpression timecourse. Write out summary tables.
- deseq_clk_inhibition.r
  - Perform DESeq2 analysis on the 3’ end data from CLK inhibition with CLK-IN-T3. Write out summary tables. Also, perform gene ontology analysis.
- germs_metrics_for_classifying_response.r
  - Analyse DESeq2 output from PPIG overexpression and CLK inhibition with respect to various GeRM properties of transcripts. Perform gene ontology analysis on the transcripts differentially localised upon dox induction. Analyse the coding content of those - differentially localised transcripts.
## Image analysis
- smfish_of_ppig_overexpression.cpproj
  - CellProfiler pipeline to analyse smFISH signal in nuclear speckles after PPIG overexpression.
- smfish_of_ppig_overexpression.r
  - Analysis of smFISH signal enrichment in nuclear speckles versus the nucleoplasm after PPIG overexpression.
- smfish_of_clki.cpproj
  - CellProfiler pipeline to analyse smFISH signal in nuclear speckles after CLK inhibition.
- smfish_of_clki.r
  - Analysis of smFISH signal enrichment in nuclear speckles versus the nucleoplasm after CLK inhibition. Includes analysis of speckle morphology.
- smfish_of_rmcd_transfection.cpproj
  - CellProfiler pipeline to analyse smFISH signal in nuclear speckles after transfection of various R-MCDs.
- smfish_of_rmcd_transfection.r
  - Analysis of smFISH signal enrichment in nuclear speckles versus the nucleoplasm after transfection of various R-MCDs.
- smfish_of_rmcd_tra2kd.cpproj
  - CellProfiler pipeline to analyse smFISH signal in nuclear speckles after knockdown of TRA2A and B.
- smfish_of_rmcd_transfection.r
  - Analysis of smFISH signal enrichment in nuclear speckles versus the nucleoplasm after knockdown of TRA2A and B.
- sc35_after_ppig.cpproj
  - Segmentation of nuclear speckles by SC35 immunofluorescence with CellProfiler.
- sc35_after_ppig.r
  - A very weird script that started as a demonstration of how to use R to analyse CellProfiler data, but ended up generating the plots that we used in the final paper.
- tra2b_after_luc7l3.cpproj
  - CellProfiler segmentation of nuclear speckles and quantification of TRA2B immunofluorescence signal.
- tra2b_after_luc7l3.r
  - Analysis of TRA2B colocalisation with nuclear speckles after transfection of LUC7L3’s R-MCD.
- polyA_after_ppig.R
  - Analysis of poly-A FISH signal in nuclear speckles and nucleoplasm after PPIG overexpression.
- tra2b_after_ppig.cpproj
  - CellProfiler segmentation of nuclear speckles and quantification of TRA2B immunofluorescence signal.
- tra2b_after_ppig.r
  - Analysis of TRA2B colocalisation with nuclear speckles after PPIG overexpression.
## Export reporter
### Sequence generation
- codon_optimisation_with_spliceai.ipynb
  - Takes the PPIG-mScarlet sequence, inserts the intron sequences and performs codon optimisation (and intronic mutagenesis) to maximise or minimise multivalency while preserving the efficiency of the splice sites.
### Nanopore
- plasmid_identity_and_splicing.r
  - Analyses nanopore data from the reporter plasmid pool to relate plasmid barcodes to their gene architectures, then analyses nanopore data from RNA extracted from cells expressing the pool to analyse the efficiency of the splicing of the reporter transcripts.
### Targeted sequencing
- clk_treatment_analysis.r
  - Analysis of the barcodes in the sequencing experiment by relating them to the gene structures in the plasmid pool. In this experiment, the reporter expression was induced from a polyclonal pool with PiggyBac integration of the plasmid pool, either in the - presence or absence of CLK-IN-T3.
- clk_treatment_barcode_parsing.r
  - Parsing the raw sequencing reads to extract the unique barcodes. Deduplicating reads based on UMIs. In this experiment, the reporter expression was induced from a polyclonal pool with PiggyBac integration of the plasmid pool, either in the presence or absence - of CLK-IN-T3.
- first_transfection_analysis.r
  - Analysis of the barcodes in the sequencing experiment by relating them to the gene structures in the plasmid pool. In this experiment, the plasmid pool was transiently transfected, and RNA was harvested after 16 hours.
- first_transfection_barcode_parsing.r
  - Parsing the raw sequencing reads to extract the unique barcodes. Deduplicating reads based on UMIs. In this experiment, the plasmid pool was transiently transfected, and RNA was harvested after 16 hours.
- piggybac_timecourse_analysis.r
  - Analysis of the barcodes in the sequencing experiment by relating them to the gene structures in the plasmid pool. In this experiment, the reporter expression was induced for 4, 8 or 12 hours from a polyclonal pool with PiggyBac integration of the plasmid pool.
- piggybac_timecourse_barcode_parsing.r
  - Parsing the raw sequencing reads to extract the unique barcodes. Deduplicating reads based on UMIs. In this experiment, the reporter expression was induced for 4, 8 or 12 hours from a polyclonal pool with PiggyBac integration of the plasmid pool.
- second_transfection_analysis.R
  - Analysis of the barcodes in the sequencing experiment by relating them to the gene structures in the plasmid pool. In this experiment, the plasmid pool was transiently transfected, and RNA was harvested after 8 or 24 hours.
 
## Bind-N-Seq Analysis 
The code for computing the binding potential of RBPs based on *in vitro* motifs is available in **Bind-N-Seq_BindingPotentialAnalysis/**.  
Input data and output files are deposited on [Zenodo](https://doi.org/10.5281/zenodo.15682421) due to their large size.

## Random Forest Modelling
The code for training a random forest model to predict nuclear retention of endogenous transcripts upon induction of R-MCD is available in **RandomForestModel/**.  
Input data and output files are deposited on [Zenodo](https://doi.org/10.5281/zenodo.15682421) due to their large size.
- second_transfection_barcode_parsing.r
  - Parsing the raw sequencing reads to extract the unique barcodes. Deduplicating reads based on UMIs. In this experiment, the plasmid pool was transiently transfected, and RNA was harvested after 8 or 24 hours.- 

