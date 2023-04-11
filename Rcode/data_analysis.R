## data analysis of Conidia, ITS and 16S datasets 
#### subset for paper 1: TU-1 vs control, without shredders
## R.K.Salis and V.C. Schreiner, April 2023

#1. Species richness calculations and plots
  #1.1 - Species richness conidia morphology
  #1.2 - ITS ASV richness all
  #1.3 - ITS ASV richness Aquatic Hyphomycetes
  #1.4 - ITS SPECIES richness Aquatic Hyphomycetes
  #1.5 - ASV richness 16S

#2. Statistical analysis of richness glms with poisson distribution
  #2.1 - analysis of Conidia species richness (morphology)
  #2.2 - analysis of ITS ASV richness - ALL
  #2.3 - analysis of ITS ASV richness - Aquatic Hyphomycetes
  #2.4 - analysis of ITS species richness - Aquatic Hyphomycetes
  #2.5 - analysis of 16S ASV richness

# 3. RDA analyses #####
  #3.1 - RDA for ITS all fungi (ASV level)
  #3.2 - RDA - ITS Aquatic Hyphomycetes (ASV level) cycles 1,2,3
  #3.3 - RDA - ITS Aquatic Hyphomycetes (species level) cycles 1,2,3
  #3.4 - RDA - Conidia Aquatic Hyphomycetes (morphology) cycles 2 and 3
  #3.5 - RDA - 16S Bacteria (ASV level) cycles 1,2,3

# 4. Individual species comparison #####
#Looking closer at following species:
  #4.1 Alatospora acuminata
  #4.2 Articulospora tetracladia
  #4.3 Flagellospora curvula
  #4.4 Lemonniera aquatica
  #4.5 Tetrachaetum elegans
  #4.6 Tetracladium marchalianum

# 5. DESeq Analyses - ITS and 16S #####
  #5.1 - ITS cycle 1 
  #5.2 - ITS cycle 2
  #5.3 - ITS cycle 3
  #5.4 - 16S cycle 1
  #5.5 - 16S cycle 2
  #5.6 - 16S cycle 3

# 6. Leaf decomposition + bacterial density - data VS ####
  #6.1 - decomposed leaf mass
  #6.2 - bacterial density

# 7. Plots
  #7.1 - Figure 2: Species richness for fungi and bacteria, decomposed leaf mass and bacterial density across cycles and treatments.
  #7.2 - Figure 3: RDA plots
    #7.2.1 - ITS all (ASV level) RDA plot
    #7.2.2 - ITS AH (ASV level) RDA plot
    #7.2.3 - ITS AH (species level) RDA plot
    #7.2.4 - Conidia RDA plot
    #7.2.5 - bacteria 16S (ASV level) RDA plot
  #7.3 - Figure S7: Aquatic Hyphomycetes taxa plots

#install/load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("S4Vectors")
install.packages("ggplot2")
install.packages("vegan")
install.packages("lsmeans")
install.packages("effects")
install.packages("cowplot")
install.packages("Rmisc")
install.packages("reshape2")
library(phyloseq); packageVersion("phyloseq")
library(S4Vectors); packageVersion("S4Vectors")
library(ggplot2); packageVersion("ggplot2")
library(lsmeans); packageVersion("lsmeans")
library(vegan); packageVersion("vegan")
library(effects)
library(cowplot)
library(Rmisc)
library(lme4)
require(reshape2)


theme_set(theme_bw())

# two different work directives (Apple, Windows)
os_sys <- switch(Sys.info()[['sysname']],
                 Windows= {print("Windows")},
                 Darwin  = {print("Darwin")}) 

prj <- if (os_sys == "Windows") {
  setwd("")
} else {
  setwd("") 
}

prj <- getwd()

##### 1. Species richness calculations and plots #####
#1.1 - Species richness conidia morphology
richness_conidia <- read.table(file.path(prj, "datafiles/richness_conidia.csv"), sep = ",",  header = TRUE, na.strings = c("NR", "NA")) 
richness_conidia
p1_AH_Morphspecies <- ggplot(richness_conidia, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = richness.mean, 
                                                                                           ymin = richness.lower, ymax = richness.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ 
  ylab("Species richness")+ xlab("")+ ylim(0,20)+ ggtitle("a) Aquatic Hyphomycetes (morph. species)")+ theme_classic()+ 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), 
        legend.position = c(0.85, 0.8), strip.text.x = element_text(face = "bold" ), text = element_text(size=14), 
        legend.title = element_blank())
p1_AH_Morphspecies

#1.2 - ITS ASV richness all
ITSall_ASV1 <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
# VS: I adjustes the read in of files that they can stay in the separate folder
ITSall_ASV <- ITSall_ASV1[, -1]
ITSall_sampledata <- read.table(file.path(prj, "datafiles/sample_data_ITS.Fto_CT1.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
richness <- rowSums(ITSall_ASV > 0) 
S<-specnumber(ITSall_ASV) 
H <- diversity(ITSall_ASV) #shannon diversity
simpson<-diversity(ITSall_ASV,index="simpson") #simpsons 1 - D
inversimp<-diversity(ITSall_ASV,index="invsimpson") #inverse simpsons 1/D
J<-H/log(S) #Pielou's evenness
evenness<-inversimp/S #simpson's evenness (1/D)/S
diversity <- data.frame(richness, S, H, simpson, inversimp, J, evenness)
diversity
ITSall_ASV_richness <- data.frame(ITSall_sampledata,diversity)
ITSall_ASV_richness
write.csv(ITSall_ASV_richness, file = "datafiles/ITSall_ASV_richness.csv")
ITSall_richness <- group.CI(richness ~ fungicide_treatment + cycle, data=ITSall_ASV_richness, ci = 0.95)
ITSall_richness

p2_all_asv <- ggplot(ITSall_richness, aes(x = fungicide_treatment)) + 
  geom_pointrange(aes(y = richness.mean, ymin = richness.lower, ymax = richness.upper, 
                      shape= fungicide_treatment))+ facet_wrap(~cycle)+  ylab('ASV richness')+ 
  xlab("")+ ylim(0,250)+ ggtitle("d) Entire ITS data set (ASVs)")+ theme_classic()+ theme(axis.text.x = element_blank(), 
                                                                                          axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), 
                                                                                          legend.position="none", strip.text.x = element_text(face = "bold" ), text = element_text(size=14))
p2_all_asv

#1.3 - ITS ASV richness Aquatic Hyphomycetes
ITSAH_ASV1 <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1_AH.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
ITSAH_ASV <- ITSAH_ASV1[, -1]
richness <- rowSums(ITSAH_ASV > 0) 
S<-specnumber(ITSAH_ASV) 
H <- diversity(ITSAH_ASV) #shannon diversity
simpson<-diversity(ITSAH_ASV,index="simpson") #simpsons 1 - D
inversimp<-diversity(ITSAH_ASV,index="invsimpson") #inverse simpsons 1/D
J<-H/log(S) #Pielou's evenness
evenness<-inversimp/S #simpson's evenness (1/D)/S
diversity <- data.frame(richness, S, H, simpson, inversimp, J, evenness)
diversity
ITSAH_ASV_richness <- data.frame(ITSall_sampledata,diversity)
ITSAH_ASV_richness
write.csv(ITSAH_ASV_richness, file = "datafiles/ITSAH_ASV_richness.csv")
ITSAH_richness <- group.CI(richness ~ fungicide_treatment + cycle, data=ITSAH_ASV_richness, ci = 0.95)
ITSAH_richness

p3_AH_asv <- ggplot(ITSAH_richness, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = richness.mean, 
                                                                                        ymin = richness.lower, ymax = richness.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('ASV richness')+ 
  xlab("")+ ylim(0,80)+ ggtitle("c) Aquatic Hyphomycetes (ITS ASVs)")+ theme_classic()+ theme(axis.text.x = element_blank(), 
                                                                                              axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
                                                                                              strip.text.x = element_text(face = "bold" ), text = element_text(size=14))
p3_AH_asv

#1.4 - ITS SPECIES richness Aquatic Hyphomycetes
ITSAH_sp1 <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1_AHsp.csv"), sep = ",",  header = TRUE, na.strings = c("NR", "NA")) 
ITSAH_sp <- ITSAH_sp1[, -1]
richness <- rowSums(ITSAH_sp > 0) 
S<-specnumber(ITSAH_sp) 
H <- diversity(ITSAH_sp) #shannon diversity
simpson<-diversity(ITSAH_sp,index="simpson") #simpsons 1 - D
inversimp<-diversity(ITSAH_sp,index="invsimpson") #inverse simpsons 1/D
J<-H/log(S) #Pielou's evenness
evenness<-inversimp/S #simpson's evenness (1/D)/S
diversityITSAHsp <- data.frame(richness, S, H, simpson, inversimp, J, evenness)
diversityITSAHsp
ITSAH_sp_richness <- data.frame(ITSall_sampledata,diversityITSAHsp)
ITSAH_sp_richness
write.csv(ITSAH_sp_richness, file = "datafiles/ITSAH_sp_richness.csv")
ITSAHsp_richness <- group.CI(richness ~ fungicide_treatment + cycle, data=ITSAH_sp_richness, ci = 0.95)
ITSAHsp_richness

p4_AH_ITSspecies <- ggplot(ITSAHsp_richness, aes(x = fungicide_treatment)) + 
  geom_pointrange(aes(y = richness.mean, ymin = richness.lower, ymax = richness.upper, 
                      shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('Species richness')+ 
  xlab("")+ ylim(0,20)+ ggtitle("b) Aquatic Hyphomycetes (ITS species)")+ theme_classic()+ 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), 
        legend.position="none", strip.text.x = element_text(face = "bold" ), text = element_text(size=14))
p4_AH_ITSspecies

#1.5 - ASV richness 16S
ASV_16S1 <- read.table(file.path(prj, "datafiles/asv_table_DFG.16S.Fto_CT1.wo.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
ASV_16S <- ASV_16S1[, -1]
sampledata_16S <- read.table(file.path(prj, "datafiles/sample_data_DFG.16S.Fto_CT1.wo.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
richness <- rowSums(ASV_16S > 0) 
S<-specnumber(ASV_16S) 
H <- diversity(ASV_16S) #shannon diversity
simpson<-diversity(ASV_16S,index="simpson") #simpsons 1 - D
inversimp<-diversity(ASV_16S,index="invsimpson") #inverse simpsons 1/D
J<-H/log(S) #Pielou's evenness
evenness<-inversimp/S #simpson's evenness (1/D)/S
diversity <- data.frame(richness, S, H, simpson, inversimp, J, evenness)
diversity
ASV_16S_samp <- data.frame(sampledata_16S,diversity)
ASV_16S_samp
write.csv(ASV_16S_samp, file = "datafiles/16S_asv_richness.csv")
richness_16S <- group.CI(richness ~ fungicide_treatment + cycle, data=ASV_16S_samp, ci = 0.95)
richness_16S
richness_16S$fungicide_treatment <- gsub("TU-1", "fungicide", richness_16S$fungicide_treatment)
richness_16S

p5_richness16s <- ggplot(richness_16S, aes(x = fungicide_treatment)) + 
  geom_pointrange(aes(y = richness.mean, ymin = richness.lower, 
                      ymax = richness.upper, shape= fungicide_treatment)) + 
  facet_wrap(~cycle)+ ylab('ASV richness')+ xlab("")+ ylim(0,1600)+ ggtitle("f) Bacteria richness (16S)")+ 
  theme_classic()+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.y = element_text(colour="black", size = 12), 
                         legend.position = "none", 
                         strip.text.x = element_text(face = "bold" ), 
                         text = element_text(size=14), legend.title = element_blank())
p5_richness16s

##### 2. Statistical analysis of richness #####
#2.1 - analysis of Conidia species richness (morphology)
data_all <- read.table(file.path(prj, "datafiles/richness_conidia_DE.csv"), header = TRUE, sep = ",", na.strings = c("NR", "NA")) 
data_all
mod1_Spec <- glm(Richness ~ fungicide_treatment + cycle + fungicide_treatment:cycle, data = data_all, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod1_Spec)
drop1(mod1_Spec, test = 'Chisq') #No interaction found
anova(mod1_Spec)
summary(mod1_Spec)
mod2_Spec <- glm(Richness ~ fungicide_treatment + cycle, data = data_all, family = 'poisson') ##### Fit model without interaction
par(mfrow = c(2, 2))
plot(mod2_Spec)
drop1(mod2_Spec, test = 'Chisq') # not significant

#2.2 - analysis of ITS ASV richness - ALL
ITSall_ASV_richness$cycle <- factor(ITSall_ASV_richness$cycle)
mod1_ITS_all <- glm(richness ~ fungicide_treatment + cycle + fungicide_treatment:cycle, data = ITSall_ASV_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod1_ITS_all)
drop1(mod1_ITS_all, test = 'Chisq') #Significant interaction 
anova(mod1_ITS_all)
summary(mod1_ITS_all)
mod2_ITS_all <- glm(richness ~ fungicide_treatment + cycle, data = ITSall_ASV_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod2_ITS_all)
drop1(mod2_ITS_all, test = 'Chisq')

#2.3 - analysis of ITS ASV richness - Aquatic Hyphomycetes
ITSAH_ASV_richness$cycle <- factor(ITSAH_ASV_richness$cycle)
mod1_ITS_AH <- glm(richness ~ fungicide_treatment + cycle + fungicide_treatment:cycle, data = ITSAH_ASV_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod1_ITS_AH)
drop1(mod1_ITS_AH, test = 'Chisq') #No interaction
anova(mod1_ITS_AH)
summary(mod1_ITS_AH)
mod2_ITS_AH <- glm(richness ~ fungicide_treatment + cycle, data   = ITSAH_ASV_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod2_ITS_AH)
drop1(mod2_ITS_AH, test = 'Chisq')

#2.4 - analysis of ITS species richness - Aquatic Hyphomycetes
ITSAH_sp_richness$cycle <- factor(ITSAH_sp_richness$cycle)
mod1_ITS_AHsp <- glm(richness ~ fungicide_treatment + cycle + fungicide_treatment:cycle, data   = ITSAH_sp_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod1_ITS_AHsp)
drop1(mod1_ITS_AHsp, test = 'Chisq') #No interaction
anova(mod1_ITS_AHsp)
summary(mod1_ITS_AHsp)
mod2_ITS_AHsp <- glm(richness ~ fungicide_treatment + cycle,data   = ITSAH_sp_richness, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod2_ITS_AHsp)
drop1(mod2_ITS_AHsp, test = 'Chisq')

#2.5 - analysis of 16S ASV richness
ASV_16S_samp$cycle <- factor(ASV_16S_samp$cycle)
mod1_16S <- glm(richness ~ fungicide_treatment + cycle + fungicide_treatment:cycle, data = ASV_16S_samp, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod1_16S)
drop1(mod1_16S, test = 'Chisq') #sign interaction 
anova(mod1_16S)
summary(mod1_16S)
mod2_16S <- glm(richness ~ fungicide_treatment + cycle, data   = ASV_16S_samp, family = 'poisson')
par(mfrow = c(2, 2))
plot(mod2_16S)
drop1(mod2_16S, test = 'Chisq')


##### 3. RDA analyses #####
#3.1 - RDA for ITS all fungi (ASV level)
#read data
data_ITS <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1.ra_withsampledata.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_ITS$replicate_num <- factor(data_ITS$replicate_num)
data_ITS$cycle <- factor(data_ITS$cycle)
# recode fungicide order
data_ITS$fungicide_treatment <- factor(data_ITS$fungicide_treatment,  levels = c("control", "TU-1"))
dITS <- data_ITS[,c(13:1254)]
Fak_ITS <- data_ITS[,c(1:12)]
Fak_ITS <- droplevels(Fak_ITS)
#RDA
rda_ITS_1 <- rda(decostand(dITS, "hellinger") ~ cycle, data = Fak_ITS)
summary(rda_ITS_1, display = NULL)
anova(rda_ITS_1, by = "margin")
rda_ITS_2 <- rda(decostand(dITS, "hellinger") ~ fungicide_treatment, data = Fak_ITS)
summary(rda_ITS_2, display = NULL)
anova(rda_ITS_2, by = "margin")
rda_ITS_3 <- rda(decostand(dITS, "hellinger") ~ cycle + fungicide_treatment + fungicide_treatment:cycle, data = Fak_ITS)
summary(rda_ITS_3, display = NULL)
anova(rda_ITS_3, by = "margin")
rda_ITS_4 <- rda(decostand(dITS, "hellinger") ~ cycle + fungicide_treatment, data = Fak_ITS)
summary(rda_ITS_4, display = NULL)
anova(rda_ITS_4, by = "margin")
scl = 'symmetric'
si_sc_ITS <- scores(rda_ITS_4, choices = 1:2, scaling = scl, display = 'sites')
sp_sc_ITS <- scores(rda_ITS_4, choices = 1:2, scaling = scl, display = 'species')
va_sc_ITS <- scores(rda_ITS_4, choices = 1:2, scaling = scl, display = 'cn')

#3.2 - RDA - ITS Aquatic Hyphomycetes (ASV level) 
#read data
data_ITS_AH <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1_AH.ra_withsampledata.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_ITS_AH$replicate_num <- factor(data_ITS_AH$replicate_num)
data_ITS_AH$cycle <- factor(data_ITS_AH$cycle)
# recode fungicide order
data_ITS_AH$fungicide_treatment <- factor(data_ITS_AH$fungicide_treatment, levels = c("control", "TU-1"))
ITS_AH <- data_ITS_AH[,c(13:283)]
Fak_ITS_AH <- data_ITS_AH[,c(1:12)]
Fak_ITS_AH <- droplevels(Fak_ITS_AH)
#RDA
rda_ITSAH_1 <- rda(decostand(ITS_AH, "hellinger") ~ cycle, data = Fak_ITS_AH)
summary(rda_ITSAH_1, display = NULL)
anova(rda_ITSAH_1, by = "margin") #sign
rda_ITSAH_2 <- rda(decostand(ITS_AH, "hellinger") ~ fungicide_treatment, data = Fak_ITS_AH)
summary(rda_ITSAH_2, display = NULL)
anova(rda_ITSAH_2, by = "margin")#ns
rda_ITSAH_3 <- rda(decostand(ITS_AH, "hellinger") ~ cycle + fungicide_treatment + fungicide_treatment:cycle, data = Fak_ITS_AH)
summary(rda_ITSAH_3, display = NULL)
anova(rda_ITSAH_3, by = "margin")#sign
rda_ITSAH_4 <- rda(decostand(ITS_AH, "hellinger") ~ cycle + fungicide_treatment, data = Fak_ITS_AH)
summary(rda_ITSAH_4, display = NULL)
anova(rda_ITSAH_4, by = "margin")
scl = 'symmetric'
si_sc_AH <- scores(rda_ITSAH_4, choices = 1:2, scaling = scl, display = 'sites')
sp_sc_AH <- scores(rda_ITSAH_4, choices = 1:2, scaling = scl, display = 'species')
va_sc_AH <- scores(rda_ITSAH_4, choices = 1:2, scaling = scl, display = 'cn')

#3.3 - RDA - ITS Aquatic Hyphomycetes (species level) 
#read data
data_ITS_AHsp <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1_AHsp.ra_withsampledata_ID.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_ITS_AHsp$replicate_num <- factor(data_ITS_AHsp$replicate_num)
data_ITS_AHsp$cycle <- factor(data_ITS_AHsp$cycle)
# recode fungicide order
data_ITS_AHsp$fungicide_treatment <- factor(data_ITS_AHsp$fungicide_treatment, levels = c("control", "TU-1"))
ITS_AHsp <- data_ITS_AHsp[,c(13:30)]
Fak_ITS_AHsp <- data_ITS_AHsp[,c(1:12)]
Fak_ITS_AHsp <- droplevels(Fak_ITS_AHsp)
#RDA
rda_ITSAHsp_1 <- rda(decostand(ITS_AHsp, "hellinger") 
                     ~ cycle, data = Fak_ITS_AHsp)
summary(rda_ITSAHsp_1, display = NULL)
anova(rda_ITSAHsp_1, by = "margin") #sign
rda_ITSAHsp_2 <- rda(decostand(ITS_AHsp, "hellinger") ~ fungicide_treatment, data = Fak_ITS_AHsp)
summary(rda_ITSAHsp_2, display = NULL)
anova(rda_ITSAHsp_2, by = "margin")#ns
rda_ITSAHsp_3 <- rda(decostand(ITS_AHsp, "hellinger") ~ cycle + fungicide_treatment + fungicide_treatment:cycle, data = Fak_ITS_AHsp)
summary(rda_ITSAHsp_3, display = NULL)
anova(rda_ITSAHsp_3, by = "margin")#sign
rda_ITSAHsp_4 <- rda(decostand(ITS_AHsp, "hellinger") ~ cycle + fungicide_treatment, data = Fak_ITS_AHsp)
summary(rda_ITSAHsp_4, display = NULL)
anova(rda_ITSAHsp_4, by = "margin")
scl = 'symmetric'
si_sc_AHsp <- scores(rda_ITSAHsp_4, choices = 1:2, scaling = scl, display = 'sites')
sp_sc_AHsp <- scores(rda_ITSAHsp_4, choices = 1:2, scaling = scl, display = 'species')
va_sc_AHsp <- scores(rda_ITSAHsp_4, choices = 1:2, scaling = scl, display = 'cn')

#3.4 - RDA - Conidia Aquatic Hyphomycetes (morphology) cycles 2 and 3
#read data
data_conid_DE <- read.table(file.path(prj, "datafiles/Raw_Conidia_DE.csv"), header = TRUE, sep = ",", na.strings = c("NR", "NA")) 
rownames(data_conid_DE) <- NULL
data_conid_DE$cycle <- factor(data_conid_DE$cycle)
data_conid_DE$fungicide_treatment <- factor(data_conid_DE$fungicide_treatment, levels = c("control", "TU-1"))
conid_DE <- data_conid_DE[,c(5:25)]
Fak_conid_DE <- data_conid_DE[,c(1:4)]
Fak_conid_DE <- droplevels(Fak_conid_DE)
#RDA
rda_conid_DE_1 <- rda(decostand(conid_DE, "hellinger") ~ cycle, data = Fak_conid_DE)
summary(rda_conid_DE_1, display = NULL)
anova(rda_conid_DE_1, by = "margin")
rda_conid_DE_2 <- rda(decostand(conid_DE, "hellinger") ~ fungicide_treatment, data = Fak_conid_DE)
summary(rda_conid_DE_2, display = NULL)
anova(rda_conid_DE_2, by = "margin")
rda_conid_DE_3 <- rda(decostand(conid_DE, "hellinger") ~ cycle + fungicide_treatment + fungicide_treatment:cycle, data = Fak_conid_DE)
summary(rda_conid_DE_3, display = NULL)
anova(rda_conid_DE_3, by = "margin")
rda_conid_DE_4 <- rda(decostand(conid_DE, "hellinger") ~ cycle + fungicide_treatment, data = Fak_conid_DE)
summary(rda_conid_DE_4, display = NULL)
anova(rda_conid_DE_4, by = "margin")
scl = 'symmetric'
si_sc_conid <- scores(rda_conid_DE_4, choices = 1:2, scaling = scl, display = 'sites')
sp_sc_conid <- scores(rda_conid_DE_4, choices = 1:2, scaling = scl, display = 'species')
va_sc_conid <- scores(rda_conid_DE_4, choices = 1:2, scaling = scl, display = 'cn')

#3.5 - RDA - 16S Bacteria (ASV level) 
#read data
data_16S <- read.table(file.path(prj, "datafiles/asv_table_16S.Fto_CT1.ra_withsampledata.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_16S$replicate_num <- factor(data_16S$replicate_num)
data_16S$cycle <- factor(data_16S$cycle)
# recode fungicide order
data_16S$fungicide_treatment <- factor(data_16S$fungicide_treatment,  levels = c("control", "TU-1"))
d16S <- data_16S[,c(13:9171)]
Fak_16S <- data_16S[,c(1:12)]
Fak_16S <- droplevels(Fak_16S)
#RDA
rda_16S_1 <- rda(decostand(d16S, "hellinger") ~ cycle, data = Fak_16S)
summary(rda_16S_1, display = NULL)
anova(rda_16S_1, by = "margin")
rda_16S_2 <- rda(decostand(d16S, "hellinger") ~ fungicide_treatment, data = Fak_16S)
summary(rda_16S_2, display = NULL)
anova(rda_16S_2, by = "margin")
rda_16S_3 <- rda(decostand(d16S, "hellinger") ~ cycle + fungicide_treatment + fungicide_treatment:cycle, data = Fak_16S)
summary(rda_16S_3, display = NULL)
anova(rda_16S_3, by = "margin")
rda_16S_4 <- rda(decostand(d16S, "hellinger") ~ cycle + fungicide_treatment, data = Fak_16S)
summary(rda_16S_4, display = NULL)
anova(rda_16S_4, by = "margin")
scl = 'symmetric'
si_sc_16S <- scores(rda_16S_4, choices = 1:2, scaling = scl, display = 'sites')
sp_sc_16S <- scores(rda_16S_4, choices = 1:2, scaling = scl, display = 'species')
va_sc_16S <- scores(rda_16S_4, choices = 1:2, scaling = scl, display = 'cn')


##### 4. individual species  comparison #####
#Looking closer at following species:
#4.1 Alatospora acuminata
#4.2 Articulospora tetracladia
#4.3 Flagellospora curvula
#4.4 Lemonniera aquatica
#4.5 Tetrachaetum elegans
#4.6 Tetracladium marchalianum

data_AqH_conidia <- read.table(file.path(prj, "datafiles/Raw_Conidia_DE.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_AqH_conidia
data_AqH_conidia$Number.replicate <- factor(data_AqH_conidia$Number.replicate)
data_AqH_conidia$cycle <- factor(data_AqH_conidia$cycle)
# recode fungicide order
data_AqH_conidia$fungicide_treatment <- factor(data_AqH_conidia$fungicide_treatment,  levels = c("control", "TU-1"))

data_AqH_ITS <- read.table(file.path(prj, "datafiles/asv_table_ITS.Fto_CT1_AHsp.ra_withsampledata_IDs.csv"), sep = ",",  header = TRUE,  na.strings = c("NR", "NA")) 
data_AqH_ITS 
data_AqH_ITS$replicate_num <- factor(data_AqH_ITS$replicate_num)
data_AqH_ITS$cycle <- factor(data_AqH_ITS$cycle)
# recode fungicide order
data_AqH_ITS$fungicide_treatment <- factor(data_AqH_ITS$fungicide_treatment,  levels = c("control", "TU-1"))

#4.1 Alatospora acuminata
mod1_conidia_Aa <- lm(Alatospora.acuminata~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_Aa)
drop1(mod1_conidia_Aa, test = 'F') #not significant
mod2_conidia_Aa <- lm(Alatospora.acuminata ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_Aa)
drop1(mod2_conidia_Aa, test = 'F') #not significant
#plot
Alatospora.acuminata_data1 <-  group.CI(Alatospora.acuminata ~ fungicide_treatment + cycle, data=data_AqH_conidia, ci = 0.95)
p1_conidia_Aa <- ggplot(Alatospora.acuminata_data1, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Alatospora.acuminata.mean, 
  ymin = Alatospora.acuminata.lower, ymax = Alatospora.acuminata.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('Sporulation')+ 
  xlab("")+ ylim(-2000,6000)+ ggtitle("")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

mod1_ITS_Aa <- lm(Alatospora.acuminata~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_Aa)
drop1(mod1_ITS_Aa, test = 'F') #not significant
mod2_ITS_Aa <- lm(Alatospora.acuminata ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_Aa)
drop1(mod2_ITS_Aa, test = 'F')#significant
#plot
Alatospora.acuminata_data2 <-  group.CI(Alatospora.acuminata ~ fungicide_treatment + cycle, data=data_AqH_ITS, ci = 0.95)
p1_ITS_Aa <- ggplot(Alatospora.acuminata_data2, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Alatospora.acuminata.mean, 
  ymin = Alatospora.acuminata.lower, ymax = Alatospora.acuminata.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('Relative abundance')+ 
  xlab("")+ ylim(-0.2,1)+ ggtitle("Alatospora acuminata")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

#4.2 Articulospora tetracladia / Hymenoscyphus tetracladius
mod1_conidia_At <- lm(Articulospora.tetracladia~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_At)
drop1(mod1_conidia_At, test = 'F') #not significant
mod2_conidia_At <- lm(Articulospora.tetracladia ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_At)
drop1(mod2_conidia_At, test = 'F')#not significant
#plot
Articulospora.tetracladia_data1 <-  group.CI(Articulospora.tetracladia ~ fungicide_treatment + cycle, data=data_AqH_conidia, ci = 0.95)
p2_conidia_At <- ggplot(Articulospora.tetracladia_data1, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Articulospora.tetracladia.mean, 
  ymin = Articulospora.tetracladia.lower, ymax = Articulospora.tetracladia.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-40,120)+ ggtitle("")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

mod1_ITS_At <- lm(Hymenoscyphus.tetracladius~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_At)
drop1(mod1_ITS_At, test = 'F') #not significant
mod2_ITS_At <- lm(Hymenoscyphus.tetracladius ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_At)
drop1(mod2_ITS_At, test = 'F')#significant
#plot
Hymenoscyphus.tetracladius_data2 <-  group.CI(Hymenoscyphus.tetracladius ~ fungicide_treatment + cycle, data=data_AqH_ITS, ci = 0.95)
p2_ITS_At <- ggplot(Hymenoscyphus.tetracladius_data2, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Hymenoscyphus.tetracladius.mean, 
  ymin = Hymenoscyphus.tetracladius.lower, ymax = Hymenoscyphus.tetracladius.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-0.2,1)+ ggtitle("Articulospora tetracladia")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

#4.3 Flagellospora curvula
mod1_conidia_Fc <- lm(Flagellospora.curvula ~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_Fc)
drop1(mod1_conidia_Fc, test = 'F') #not significant
mod2_conidia_Fc <- lm(Flagellospora.curvula ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_Fc)
drop1(mod2_conidia_Fc, test = 'F') #not significant

mod1_ITS_Fc <- lm(Flagellospora.curvula~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_Fc)
drop1(mod1_ITS_Fc, test = 'F') #not significant
mod2_ITS_Fc <- lm(Flagellospora.curvula ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_Fc)
drop1(mod2_ITS_Fc, test = 'F') #not significant

#4.4 Lemonniera aquatica
mod1_conidia_La <- lm(Lemonniera.aquatica ~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_La)
drop1(mod1_conidia_La, test = 'F') #not significant
mod2_conidia_La <- lm(Lemonniera.aquatica ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_La)
drop1(mod2_conidia_La, test = 'F') #not significant

mod1_ITS_La <- lm(Lemonniera.aquatica ~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_La)
drop1(mod1_ITS_La, test = 'F') #not significant
mod2_ITS_La <- lm(Lemonniera.aquatica ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_La)
drop1(mod2_ITS_La, test = 'F') #not significant

#4.5 Tetrachaetum elegans
mod1_conidia_Te <- lm(Tetrachaetum.elegans~ fungicide_treatment + cycle  +  fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_Te)
drop1(mod1_conidia_Te, test = 'F') #not significant
mod2_conidia_Te <- lm(Tetrachaetum.elegans ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_Te)
drop1(mod2_conidia_Te, test = 'F')
#plot
Tetrachaetum.elegans_data1 <-  group.CI(Tetrachaetum.elegans ~ fungicide_treatment + cycle, data=data_AqH_conidia, ci = 0.95)
p5_conidia_Te <- ggplot(Tetrachaetum.elegans_data1, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Tetrachaetum.elegans.mean, 
  ymin = Tetrachaetum.elegans.lower, ymax = Tetrachaetum.elegans.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-200,1200)+ ggtitle("")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

mod1_ITS_Te <- lm(Tetrachaetum.elegans ~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_Te)
drop1(mod1_ITS_Te, test = 'F') #not significant
mod2_ITS_Te <- lm(Tetrachaetum.elegans ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_Te)
drop1(mod2_ITS_Te, test = 'F')
#plot
Tetrachaetum.elegans_data2 <-  group.CI(Tetrachaetum.elegans ~ fungicide_treatment + cycle, data=data_AqH_ITS, ci = 0.95)
p5_ITS_Te <- ggplot(Tetrachaetum.elegans_data2, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Tetrachaetum.elegans.mean, 
  ymin = Tetrachaetum.elegans.lower, ymax = Tetrachaetum.elegans.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-0.2,1)+ ggtitle("Tetrachaetum elegans")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

#4.6 Tetracladium marchalianum
mod1_conidia_Tm <- lm(Tetracladium.marchalianum~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod1_conidia_Tm)
drop1(mod1_conidia_Tm, test = 'F') #not significant
mod2_conidia_Tm <- lm(Tetracladium.marchalianum ~ fungicide_treatment + cycle, data = data_AqH_conidia)
par(mfrow = c(2, 2))
plot(mod2_conidia_Tm)
drop1(mod2_conidia_Tm, test = 'F')
#plot
Tetracladium.marchalianum_data1 <-  group.CI(Tetracladium.marchalianum ~ fungicide_treatment + cycle, data=data_AqH_conidia, ci = 0.95)
p6_conidia_Tm <- ggplot(Tetracladium.marchalianum_data1, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Tetracladium.marchalianum.mean, 
  ymin = Tetracladium.marchalianum.lower, ymax = Tetracladium.marchalianum.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-500,1000)+ ggtitle("")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

mod1_ITS_Tm <- lm(Tetracladium.marchalianum~ fungicide_treatment + cycle  + fungicide_treatment:cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod1_ITS_Tm)
drop1(mod1_ITS_Tm, test = 'F') #not significant
mod2_ITS_Tm <- lm(Tetracladium.marchalianum ~ fungicide_treatment + cycle, data = data_AqH_ITS)
par(mfrow = c(2, 2))
plot(mod2_ITS_Tm)
drop1(mod2_ITS_Tm, test = 'F') 
#plot
Tetracladium.marchalianum_data2 <-  group.CI(Tetracladium.marchalianum ~ fungicide_treatment + cycle, data=data_AqH_ITS, ci = 0.95)
p6_ITS_Tm <- ggplot(Tetracladium.marchalianum_data2, aes(x = fungicide_treatment)) + geom_pointrange(aes(y = Tetracladium.marchalianum.mean, 
  ymin = Tetracladium.marchalianum.lower, ymax = Tetracladium.marchalianum.upper, shape= fungicide_treatment))+ facet_wrap(~cycle)+ ylab('')+ 
  xlab("")+ ylim(-0.2,1)+ ggtitle("Tetracladium marchalianum")+ theme_classic()+ theme(axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", size = 12), legend.position="none", 
  strip.text.x = element_text(face = "bold" ), text = element_text(size=12))

##### 5. DESeq Analyses - ITS and 16S #####
#BiocManager::install("DESeq2")
library(DESeq2); packageVersion("DESeq2")
library(ggrepel)
#load phyloseq objects
ITS.Fto_CT1 <- readRDS(file.path(prj, "datafiles/ITS.Fto_CT1.rds"))
B16S.Fto_CT1 <- readRDS(file.path(prj, "datafiles/B16S.Fto_CT1.rds"))
#subset for DESeq pair-wise comparisons
cycle1.ITS  = subset_samples(ITS.Fto_CT1, cycle == "1")
cycle1.ITS
cycle2.ITS  = subset_samples(ITS.Fto_CT1, cycle == "2")
cycle2.ITS
cycle3.ITS  = subset_samples(ITS.Fto_CT1, cycle == "3")
cycle3.ITS
cycle1.16S  = subset_samples(B16S.Fto_CT1, cycle == "1")
cycle1.16S
cycle2.16S  = subset_samples(B16S.Fto_CT1, cycle == "2")
cycle2.16S
cycle3.16S  = subset_samples(B16S.Fto_CT1, cycle == "3")
cycle3.16S

#5.1 - ITS cycle 1 
diagdds = phyloseq_to_deseq2(cycle1.ITS, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(cycle1.ITS)[rownames(sigtab), ], "matrix"))
sigtab
write.csv(sigtab, file = "datafiles/ITS_cycle1_DESeq_sigtab.csv")
sigtab1a <- read.csv("datafiles/ITS_cycle1_DESeq_sigtab_ID.csv")
sigtab1 <- sigtab1a[, -1]
row.names(sigtab1)<- sigtab1a$ASV #sample data row names must align with dada2 rowname outputs
sigtab1 = as.data.frame(sigtab1)
head(sigtab1)
# ID order
x = tapply(sigtab1$log2FoldChange, sigtab1$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab1$ID = factor(as.character(sigtab1$ID), levels=names(x))
cycle1plot <- ggplot(sigtab1, aes(y=ID, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("ITS - cycle 1") #
### geom_label_repel
pITSc1 <- cycle1plot + 
  geom_label_repel(aes(label=rownames(sigtab1)),
                   box.padding   = 0.35, 
                   point.padding = 0.25,
                   segment.color = 'grey50') +
  theme_classic()
g1 <- ggplotGrob(pITSc1)
pITSc1

#5.2 - ITS cycle 2 
diagdds = phyloseq_to_deseq2(cycle2.ITS, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(cycle2.ITS)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, file = "datafiles/ITS_cycle2_DESeq_sigtab.csv")
sigtab2a <- read.csv("datafiles/ITS_cycle2_DESeq_sigtab_ID.csv")
sigtab2 <- sigtab2a[, -1]
row.names(sigtab2)<- sigtab2a$ASV #sample data row names must align with dada2 rowname outputs
sigtab2 = as.data.frame(sigtab2)
head(sigtab2)
# ID order
x = tapply(sigtab2$log2FoldChange, sigtab2$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab2$ID = factor(as.character(sigtab2$ID), levels=names(x))
cycle2plot <- ggplot(sigtab2, aes(y=ID, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("ITS - cycle 2") #
### geom_label_repel
pITSc2 <- cycle2plot + 
  geom_label_repel(aes(label=rownames(sigtab2)),
                   box.padding   = 0.35, 
                   point.padding = 0.25,
                   segment.color = 'grey50') +
  theme_classic()
g2 <- ggplotGrob(pITSc2)
pITSc2

#5.3 - ITS cycle 3 
diagdds = phyloseq_to_deseq2(cycle3.ITS, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(cycle3.ITS)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, file = "datafiles/ITS_cycle3_DESeq_sigtab.csv")
sigtab3a <- read.csv("datafiles/ITS_cycle3_DESeq_sigtab_ID.csv")
sigtab3 <- sigtab3a[, -1]
row.names(sigtab3)<- sigtab3a$ASV #sample data row names must align with dada2 rowname outputs
sigtab3 = as.data.frame(sigtab3)
head(sigtab3)
# ID order
x = tapply(sigtab3$log2FoldChange, sigtab3$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab3$ID = factor(as.character(sigtab3$ID), levels=names(x))
cycle3plot <- ggplot(sigtab3, aes(y=ID, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("ITS - cycle 3") 
### geom_label_repel
pITSc3 <- cycle3plot + 
  geom_label_repel(aes(label=rownames(sigtab3)),
                   box.padding   = 0.35, 
                   point.padding = 0.25,
                   segment.color = 'grey50') +
  theme_classic()
g3 <- ggplotGrob(pITSc3)
pITSc3

#5.4 - 16S cycle 1 
diagdds = phyloseq_to_deseq2(cycle1.16S, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab16Sc1 = res[(res$padj < alpha), ]
sigtab16Sc1 = cbind(as(sigtab16Sc1, "data.frame"), as(tax_table(cycle1.16S)[rownames(sigtab16Sc1), ], "matrix"))
head(sigtab16Sc1)
write.csv(sigtab16Sc1, file = "datafiles/16S_cycle1_DESeq_sigtab.csv")
sigtab16S1a <- read.csv("datafiles/16S_cycle1_DESeq_sigtab_ID.csv")
sigtab16S1 <- sigtab16S1a[, -1]
row.names(sigtab16S1)<- sigtab16S1a$ASV #sample data row names must align with dada2 rowname outputs
sigtab16S1 = as.data.frame(sigtab16S1)
head(sigtab16S1)
# Order order
x = tapply(sigtab16S1$log2FoldChange, sigtab16S1$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab16S1$Order = factor(as.character(sigtab16S1$Order), levels=names(x))
# ID order
x = tapply(sigtab16S1$log2FoldChange, sigtab16S1$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab16S1$ID = factor(as.character(sigtab16S1$ID), levels=names(x))

cycle1plot16s <- ggplot(sigtab16S1, aes(y=ID, x=log2FoldChange, color=Order)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("16S - cycle 1") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
### geom_label_repel
p16Sc1 <- cycle1plot16s + 
  geom_label_repel(aes(label=rownames(sigtab16S1)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()
p16Sc1

#5.5 - 16S cycle 2 
diagdds = phyloseq_to_deseq2(cycle2.16S, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab16S2 = res[(res$padj < alpha), ]
sigtab16S2 = cbind(as(sigtab16S2, "data.frame"), as(tax_table(cycle2.16S)[rownames(sigtab16S2), ], "matrix"))
head(sigtab16S2)
write.csv(sigtab16S2, file = "datafiles/16S_cycle2_DESeq_sigtab.csv")
sigtab16S2a <- read.csv("datafiles/16S_cycle2_DESeq_sigtab_ID.csv")
sigtab16S2 <- sigtab16S2a[, -1]
row.names(sigtab16S2)<- sigtab16S2a$ASV #sample data row names must align with dada2 rowname outputs
sigtab16S2 = as.data.frame(sigtab16S2)
head(sigtab16S2)
# Order order
x = tapply(sigtab16S2$log2FoldChange, sigtab16S2$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab16S2$Order = factor(as.character(sigtab16S2$Order), levels=names(x))
# ID order
x = tapply(sigtab16S2$log2FoldChange, sigtab16S2$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab16S2$ID = factor(as.character(sigtab16S2$ID), levels=names(x))

cycle2plot16s <- ggplot(sigtab16S2, aes(y=ID, x=log2FoldChange, color=Order)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("16S - cycle 2") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
### geom_label_repel
p16Sc2 <- cycle2plot16s + 
  geom_label_repel(aes(label=rownames(sigtab16S2)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()
p16Sc2

#5.6 - 16S cycle 3 
diagdds = phyloseq_to_deseq2(cycle3.16S, ~ fungicide_treatment)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
#results table
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab16S3 = res[(res$padj < alpha), ]
sigtab16S3 = cbind(as(sigtab16S3, "data.frame"), as(tax_table(cycle3.16S)[rownames(sigtab16S3), ], "matrix"))
head(sigtab16S3)
write.csv(sigtab16S3, file = "datafiles/16S_cycle3_DESeq_sigtab.csv")
sigtab16S3a <- read.csv("datafiles/16S_cycle3_DESeq_sigtab_ID.csv")
sigtab16S3 <- sigtab16S3a[, -1]
row.names(sigtab16S3)<- sigtab16S3a$ASV #sample data row names must align with dada2 rowname outputs
sigtab16S3 = as.data.frame(sigtab16S3)
head(sigtab16S3)
# Order order
x = tapply(sigtab16S3$log2FoldChange, sigtab16S3$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab16S3$Order = factor(as.character(sigtab16S3$Order), levels=names(x))
# ID order
x = tapply(sigtab16S3$log2FoldChange, sigtab16S3$ID, function(x) max(x))
x = sort(x, TRUE)
sigtab16S3$ID = factor(as.character(sigtab16S3$ID), levels=names(x))
cycle3plot16s <- ggplot(sigtab16S3, aes(y=ID, x=log2FoldChange, color=Order)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  ggtitle("16S - cycle 3") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
### geom_label_repel
p16Sc3 <- cycle3plot16s + 
  geom_label_repel(aes(label=rownames(sigtab16S3)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()
p16Sc3

#### 6 leaf decomposition + bact density - data VS ####
#### 6.1 - decomposed leaf mass ####
# read in data
data_all <- read.table(file.path(prj, "datafiles/rawdata_DE.csv"), 
                       header = TRUE, 
                       sep = ",", 
                       na.strings = c("NR", "NA")) 
# encode factors
data_all$Number.replicate <- factor(data_all$Number.replicate)
data_all$cycle <- factor(data_all$cycle)

# recode fungicide order
data_all$treatment <- factor(data_all$treatment, 
                             levels = c("control", "TU-1"))
# calculate decomposed leaf mass (DLM)
# The data from this endpint are % of lost leaf mass, divided with number degreedays 
# (to normalise for temperature differences)
data_all$DLM <- ifelse(data_all$cycle =="1", 
                       (((data_all$weight.start - data_all$weight.end)
                         / data_all$weight.start)
                        / ((1*data_all$Temp_field + 2*data_all$Temp_lab)/3)),
                       ((data_all$weight.start - data_all$weight.end) 
                        / data_all$weight.start) / (data_all$Temp_lab))
# fit model
mod1_DE_OMD <- lm(DLM ~ treatment + cycle + treatment:cycle,
                  data = data_all)
# diagnose model
par(mfrow = c(2, 2))
plot(mod1_DE_OMD)
# model fits
drop1(mod1_DE_OMD, test = 'F')
# significant interaction between fungicide and cycle.
mod2_DE_OMD <- lm(DLM ~ treatment + cycle ,
                  data = data_all)
drop1(mod2_DE_OMD, test = 'F')

# calculate CI of DLM
data_all$DLM <- data_all$DLM*1000
DLM_all <- group.CI(DLM ~ treatment + cycle,
                    data=data_all,
                    ci = 0.95)
# prepare figure
DLM <- ggplot(DLM_all, aes(x = treatment)) +
  geom_pointrange(aes(y = DLM.mean, ymin = DLM.lower, ymax = DLM.upper, shape= treatment))+
  facet_wrap(~cycle)+
  ylab('Decomposed leaf mass (%) per dday')+
  xlab("")+
  ylim(2,40)+
  ggtitle("e) Decomposed leaf mass")+ 
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", 
                                                                  size = 12),
        legend.position="none", strip.text.x = element_text(face = "bold" ), 
        text = element_text(size=14))

# svg(filename="DLM_absolute_values.svg", width=5, height=5, pointsize=12)
DLM
# dev.off()
# Colouring of the colonisation cycles (boxes on the top) afterwards in inkscape.

#### 6.2 - bacterial density ####
# fit model
mod1_DE_bac <- lm(Bact_density ~ treatment + cycle + treatment:cycle,
                  data = data_all)
# diagnose model
par(mfrow = c(2, 2))
plot(mod1_DE_bac)
# model fits
drop1(mod1_DE_bac, test = 'F')
# no significant interaction, simplify model
mod2_DE_bac <- lm(Bact_density ~ treatment + cycle ,
                  data = data_all)
drop1(mod2_DE_bac, test = 'F')
summary(mod1_DE_bac)


# calculate CI of bacterial density
data_all$log_bact <- log10(data_all$Bact_density)
Bact_dens <- group.CI(log_bact ~ treatment + cycle,
                      data=data_all,
                      ci = 0.95)
# prepare figure
Fig_bac <- ggplot(Bact_dens, aes(x = treatment)) +
  geom_pointrange(aes(y = log_bact.mean, ymin = log_bact.lower, ymax = log_bact.upper, shape= treatment))+
  facet_wrap(~cycle)+
  ylab('Density (log bacteria per cm²)')+
  xlab("")+
  ggtitle("g) Bacteria density")+ 
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_text(colour="black", 
                                                                  size = 12),
        legend.position="none", strip.text.x = element_text(face = "bold" ), 
        text = element_text(size=14))
# svg(filename="Bact_density.svg", width=5, height=5, pointsize=12)
Fig_bac
# dev.off()
# Colouring of the colonisation cycles (boxes on the top) afterwards in inkscape.

#7
# 7.1 - Figure 2: Species richness for fungi and bacteria, decomposed leaf mass and bacterial density across cycles and treatments.
ggdraw() +
  draw_plot(p1_AH_Morphspecies, x = 0, y = 0.66, width = 0.33, height = 0.33) +
  draw_plot(p4_AH_ITSspecies, x = 0.33, y = 0.66, width = 0.33, height = 0.33) +
  draw_plot(p3_AH_asv, x = 0.66, y = 0.66, width = 0.33, height = 0.33) +
  draw_plot(p2_all_asv, x = 0.17, y = 0.33, width = 0.33, height = 0.33) +
  draw_plot(DLM, x = 0.5, y = 0.33, width = 0.33, height = 0.33) +
  draw_plot(p5_richness16s, x = 0.17, y = 0, width = 0.33, height = 0.33)+
  draw_plot(Fig_bac, x = 0.5, y = 0, width = 0.33, height = 0.33)

# 7.2 - Figure 3: RDA plots
#7.2.1 - ITS all (ASV level) RDA plot
par(mfrow = c(3, 2))
p1 <- plot(rda_ITS_4, scaling = 'symmetric', type = 'n', xlab = paste(attributes(summary(rda_ITS_4)$cont$importance)$dimnames[[2]][1], " (", round(summary(rda_ITS_4)$cont$importance[2,1]*100), "% total variance)", sep=""), ylab = paste(attributes(summary(rda_ITS_4)$cont$importance)$dimnames[[2]][2], " (", round(summary(rda_ITS_4)$cont$importance[2,2]*100), "% total variance)", sep=""), main = "Entire ITS data set (ASVs)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
cols <- c("#19A1BF", "#F00303") 
pchs <- c(15, 16, 17) 
points(si_sc_ITS, col = cols[Fak_ITS$fungicide],  # colour by fungicide
       pch = pchs[Fak_ITS$cycle],      # shape by cycle
       cex = 1.5)
ordispider(p1, interaction(Fak_ITS$cycle, Fak_ITS$fungicide), col = rep(cols, each = 3), label = F) # connect singe points per fungicide treatment and cycle

#7.2.2 - ITS AH (ASV level) RDA plot
p2 <- plot(rda_ITSAH_4, scaling = 'symmetric', type = 'n', xlab = paste(attributes(summary(rda_ITSAH_4)$cont$importance)$dimnames[[2]][1], " (", round(summary(rda_ITSAH_4)$cont$importance[2,1]*100), "% total variance)", sep=""), ylab = paste(attributes(summary(rda_ITSAH_4)$cont$importance)$dimnames[[2]][2], " (", round(summary(rda_ITSAH_4)$cont$importance[2,2]*100), "% total variance)", sep=""), main = "Aquatic hyphomycetes (ITS ASVs)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
cols <- c("#19A1BF", "#F00303") 
pchs <- c(15, 16, 17)
points(si_sc_AH, col = cols[Fak_ITS_AH$fungicide_treatment],  # colour by fungicide
       pch = pchs[Fak_ITS_AH$cycle],      # shape by cycle
       cex = 1.5)
ordispider(p2, interaction(Fak_ITS_AH$cycle, Fak_ITS_AH$fungicide_treatment), col = rep(cols, each = 3), label = F) # connect singe points per fungicide treatment and cycle

#7.2.3 - ITS AH (species level) RDA plot
p3 <- plot(rda_ITSAHsp_4, scaling = 'symmetric', type = 'n', xlab = paste(attributes(summary(rda_ITSAHsp_4)$cont$importance)$dimnames[[2]][1], " (", round(summary(rda_ITSAHsp_4)$cont$importance[2,1]*100), "% total variance)", sep=""), ylab = paste(attributes(summary(rda_ITSAHsp_4)$cont$importance)$dimnames[[2]][2], " (", round(summary(rda_ITSAHsp_4)$cont$importance[2,2]*100), "% total variance)", sep=""), main = "Aquatic hyphomycetes (ITS species)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
cols <- c("#19A1BF", "#F00303") 
pchs <- c(15, 16, 17) 
points(si_sc_AHsp, col = cols[Fak_ITS_AHsp$fungicide_treatment],  # colour by fungicide
       pch = pchs[Fak_ITS_AHsp$cycle],      # shape by cycle
       cex = 1.5)
ordispider(p3, interaction(Fak_ITS_AHsp$cycle, Fak_ITS_AHsp$fungicide_treatment), col = rep(cols, each = 3), label = F) # connect singe points per fungicide treatment and cycle

#7.2.4 Conidia RDA plot
p4 <- plot(rda_conid_DE_4, scaling = 'symmetric', type = 'n', xlab = paste(attributes(summary(rda_conid_DE_4)$cont$importance)$dimnames[[2]][1], " (", round(summary(rda_conid_DE_4)$cont$importance[2,1]*100), "% total variance)", sep=""), ylab = paste(attributes(summary(rda_conid_DE_4)$cont$importance)$dimnames[[2]][2], " (", round(summary(rda_conid_DE_4)$cont$importance[2,2]*100), "% total variance)", sep=""), main = "Aquatic hyphomycetes (morph.)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
cols <- c("#19A1BF", "#F00303")
pchs <- c(16, 17)
points(si_sc_conid, col = cols[Fak_conid_DE$fungicide_treatment],  # colour by fungicide
       pch = pchs[Fak_conid_DE$cycle],      # shape by cycle
       cex = 1.5)
ordispider(p4, interaction(Fak_conid_DE$cycle, Fak_conid_DE$fungicide_treatment), col = rep(cols, each = 2), label = F) # connect singe points per fungicide treatment and cycle

#7.2.5 - bacteria 16S (ASV level) RDA plot
p5 <- plot(rda_16S_4, scaling = 'symmetric', type = 'n', xlab = paste(attributes(summary(rda_16S_4)$cont$importance)$dimnames[[2]][1], " (", round(summary(rda_16S_4)$cont$importance[2,1]*100), "% total variance)", sep=""), ylab = paste(attributes(summary(rda_16S_4)$cont$importance)$dimnames[[2]][2], " (", round(summary(rda_16S_4)$cont$importance[2,2]*100), "% total variance)", sep=""), main = "Bacteria (16S)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
cols <- c("#19A1BF", "#F00303")
pchs <- c(15, 16, 17)
points(si_sc_16S, col = cols[Fak_16S$fungicide],  # colour by fungicide
       pch = pchs[Fak_16S$cycle],      # shape by cycle
       cex = 1.5)
ordispider(p5, interaction(Fak_16S$cycle, Fak_16S$fungicide), col = rep(cols, each = 3), label = F) # connect singe points per fungicide treatment and cycle
legend('topleft', legend = c("control", "fungicide", paste0('cycle ', levels(Fak_16S$cycle))), col = c(cols, 'grey20', 'grey20', 'grey20'), pch = c(rep(16, 2), 15, 16, 17), bty = "n", cex=1.3)

# 7.3 - Figure S7: Aq. Hyphomycetes taxa plots
ggdraw() +
  draw_plot(p1_ITS_Aa, x = 0, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(p2_ITS_At, x = 0.25, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(p5_ITS_Te, x = 0.5, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(p6_ITS_Tm, x = 0.75, y = 0.5, width = 0.25, height = 0.5) +
  draw_plot(p1_conidia_Aa, x = 0, y = 0, width = 0.25, height = 0.5) +
  draw_plot(p2_conidia_At, x = 0.25, y = 0, width = 0.25, height = 0.5) +
  draw_plot(p5_conidia_Te, x = 0.5, y = 0, width = 0.25, height = 0.5) +
  draw_plot(p6_conidia_Tm, x = 0.75, y = 0, width = 0.25, height = 0.5)


