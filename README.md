# fungicides_metabarcoding

A repository containing the data and R code for the associated paper titled "Effects of fungicides on aquatic fungi and bacteria: a comparison of morphological and molecular approaches from a microcosm experiment", currently in submission.

Created by Romana K. Salis


## Content Overview


### Rcode:

1. DADA2_analysis_WP1_ITS.R - R code used for sequence processing of the ITS data using the DADA2 pipeline, subsetting the data and normalising read counts
2. DADA2_analysis_WP1_16S.R - R code used for sequence processing of the 16S data


3. data_analysis.R - R code used the data analyses of the 5 datasets, including: species richness calculations, RDA analyses, comparisons of individual species, differential abundance analysis of ASVs, analysis of decomposed leaf mass and bacterial density, code for creating the figures and performing the statistical analyses.

### Files for the sequence processing:

4. sh_general_release_dynamic_04.02.2020_LR.fasta - customised database with long read sequences
5. unite.taxa1_LR_rescued.csv - ITS taxonomic assignment including rescued taxonomy
6. ITS.Fto_CT1.rds - R data file: phyloseq object of entire ITS dataset (DATASET 3)
7. ITS.Fto_CT1_AH.rds  - R data file: phyloseq object of aquatic hyphomycetes ASVs only (DATASET 4)
8. 1ITS.Fto_CT1_AHsp.rds - R data file: phyloseq object with aquatic hyphomycetes summarised to species level (DATASET 2)
9. B16S.Fto_CT1.rds - R data file: phyloseq object of 16S dataset (DATASET 5)

### Files for data analysis:

10. Raw_Conidia_DE.csv - taxon abundance table for morphology identified aquatic hyphomycetes (DATASET 1)
11. richness_conidia_DE.csv - morphology aquatic hyphomycete richness (DATASET 1)
12. richness_conidia.csv - morphology aquatic hyphomycete richness - treatment means + CI (DATASET 1)


13. asv_table_ITS.Fto_CT1_AHsp.csv - taxon abundance table for the ITS-identified aquatic hyphomycete species (DATASET 2)
14. asv_table_ITS.Fto_CT1_AHsp.ra_withsampledata_IDs.csv - relative Aquatic Hyphomycete species abundances with sample data, including species IDs (DATASET 2)
15. ITSAH_sp_richness.csv - aquatic hyphomycete species richness for the ITS dataset (DATASET 2)


16. asv_table_ITS.Fto_CT1.csv - taxon abundance table for the entire ITS dataset - ASVs (DATASET 3)
17. ITSall_ASV_richness.csv - ASV richness for the ITS dataset (DATASET 3)
18. sample_data_ITS.Fto_CT1.csv - sample data for the ITS dataset (DATASET 3)
19. asv_table_ITS.Fto_CT1.ra_withsampledata.csv - relative ITS ASV abundances with sample data (DATASET 3)
20. ITS_cycle1_DESeq_sigtab_ID.csv, ITS_cycle2_DESeq_sigtab_ID.csv, ITS_cycle3_DESeq_sigtab_ID.csv - files for DESeq analysis (DATASET 3)
21. asvtax_table_ITS.Fto_CT1_sp.guilds.prop.csv - ASVs assigned to functional groups fungal trophic modes (DATASET 3)


22. asv_table_ITS.Fto_CT1_AH.csv - taxon abundance table for the ITS-identified aquatic hyphomycete ASVs (DATASET 4)
23. asv_table_ITS.Fto_CT1_AH.ra_withsampledata.csv - relative Aquatic Hyphomycete ASV abundances with sample data (DATASET 4)
24. ITSAH_ASV_richness.csv - aquatic hyphomycete ASV richness for the ITS dataset (DATASET 4)


25. asv_table_B16S.Fto_CT1.csv - taxon abundance table for the 16S bacterial ASVs (DATASET 5)
26. B16S_ASV_richness.csv - ASV richness for the 16S dataset (DATASET 5)
27. sample_data_B16S.Fto_CT1.csv - sample data for the 16S dataset (DATASET 5)
28. asv_table_B16S.Fto_CT1.ra_withsampledata.csv - relative 16S ASV abundances with sample data (DATASET 5)
29. 16S_cycle1_DESeq_sigtab_ID.csv, 16S_cycle2_DESeq_sigtab_ID.csv, 16S_cycle3_DESeq_sigtab_ID.csv  - files for DESeq analysis (DATASET 5)
30. 16S_functions.csv - functions assigned to the bacterial ASVs using the Functional Annotation of Prokaryotic Taxa database version 1.2.6 (FAPROTAX)


31. rawdata_DE.csv - raw data including bacterial density and leaf mass
