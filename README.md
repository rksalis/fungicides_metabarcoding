# fungicides_metabarcoding

A repository containing the data and R code for the associated paper titled "Effects of fungicides on aquatic fungi and bacteria: a comparison of morphological and molecular approaches from a microcosm experiment", currently in submission.

Created by Romana K. Salis


Content Overview

Rcode:
1. DADA2_analysis_WP1_ITS.R - R code used for sequence processing of the ITS data using the DADA2 pipeline, subsetting the data and normalising read counts
2. DADA2_analysis_WP1_16S.R - R code used for sequence processing of the 16S data

3. data_analysis.R - R code used the data analyses of the 5 datasets, including: species richness calculations, RDA analyses, comparisons of individual species, differential abundance analysis of ASVs, analysis of decomposed leaf mass and bacterial density, code for creating the figures and performing the statistical analyses.

Files for the sequence processing:
4. sh_general_release_dynamic_04.02.2020_LR.fasta - customised database with long read sequences
5. unite.taxa1_LR_rescued.csv - ITS taxonomic assignment including rescued taxonomy
6. ITS.Fto_CT1.rds - R data file: phyloseq object of entire ITS dataset (DATASET 3)
7. ITS.Fto_CT1_AH.rds  - R data file: phyloseq object of aquatic hyphomycetes ASVs only (DATASET 4)
8. 1ITS.Fto_CT1_AHsp.rds - R data file: phyloseq object with aquatic hyphomycetes summarised to species level (DATASET 2)
9. B16S.Fto_CT1.rds - R data file: phyloseq object of 16S dataset (DATASET 5)

Files for data analysis:
10. Raw_Conidia_DE.csv - taxon abundance table for morphology identified aquatic hyphomycetes (DATASET 1)
11. richness_conidia.csv - morphology aquatic hyphomycete richness - treatment means + CI (DATASET 1)

12. asv_table_ITS.Fto_CT1.csv - taxon abundance table for the entire ITS dataset - ASVs (DATASET 3)
13. ITSall_ASV_richness.csv - ASV richness for the ITS dataset (DATASET 3)
14. sample_data_ITS.Fto_CT1.csv - sample data for the ITS dataset (DATASET 3)
15. asv_table_ITS.Fto_CT1.ra_withsampledata.csv - relative ITS ASV abundances with sample data (DATASET 3)

16. asv_table_ITS.Fto_CT1_AH.csv - taxon abundance table for the ITS-identified aquatic hyphomycete ASVs (DATASET 4)
17. asv_table_ITS.Fto_CT1_AH.ra_withsampledata.csv - relative Aquatic Hyphomycete ASV abundances with sample data (DATASET 4)

18. asv_table_ITS.Fto_CT1_AHsp.csv - taxon abundance table for the ITS-identified aquatic hyphomycete species (DATASET 2)
19. asv_table_ITS.Fto_CT1_AHsp.ra_withsampledata_IDs.csv - relative Aquatic Hyphomycete species abundances with sample data, including species IDs (DATASET 2)

20. asv_table_B16S.Fto_CT1.csv - taxon abundance table for the 16S bacterial ASVs (DATASET 5)
21. sample_data_B16S.Fto_CT1.csv - sample data for the 16S dataset (DATASET 5)
22. asv_table_16S.Fto_CT1.ra_withsampledata.csv - relative 16S ASV abundances with sample data (DATASET 5)

23. rawdata_DE.csv - raw data including bacterial density and leaf mass
