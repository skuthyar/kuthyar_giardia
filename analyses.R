# PERMANOVA
library(vegan)
library(data.table)
library(tidyverse)

# All sites, unweighted UniFrac
metadata = fread("allsites_metadata.txt", header = T)
unweighted_bothyears <- read.table("unweighted-both.txt", header = T)
metadata_both = inner_join(unweighted_bothyears, metadata, by = "sampleid")
unweighted_bothyears <- as.dist(unweighted_bothyears[,2:59])

adonis(formula = unweighted_bothyears ~ humaninteraction+group+year+giardiapresence, data = metadata_both, permutations = 5000)

# All sites, weighted Unifrac
weighted_bothyears <-read.table("weighted-bothyears-distance-matrix.tsv", header = T)
metadata_together = inner_join(weighted_bothyears, metadata, by="sampleid")
weighted_bothyears <- as.dist(weighted_bothyears[,2:59])

adonis(formula = weighted_bothyears~humaninteraction+group+year+giardiapresence, data=metadata_together, permutations = 5000)

# filtered unweighted Unifrac 
# Rural
giardiametadata = fread("giardia-metadata.txt", header = T)
unweighted_rural <- read.table("unweighted-rural-distance-matrix.tsv", header=T)
unweighted_rural_metadata = inner_join(unweighted_rural, giardiametadata, by = "sampleid")
unweighted_rural<-as.dist(unweighted_rural[,2:28])

adonis(formula = unweighted_rural ~ group+year+giardiapresence, data=unweighted_rural_metadata, permutations=5000)

# Village
unweighted_village <- read.table("unweighted-village-distance-matrix.tsv", header=T)
unweighted_village_metadata = inner_join(unweighted_village, giardiametadata, by = "sampleid")
unweighted_village<-as.dist(unweighted_village[,2:23])

adonis(formula = unweighted_village ~ group+year+giardiapresence, data=unweighted_village_metadata, permutations = 5000)

# Across habitats 
giardiametadata = fread("giardia-metadata
.txt", header = T)
unweighted_interaction <- read.table("unweighted-2017-distance-matrix.tsv", header = T)
giardia_metadata = inner_join(unweighted_interaction, giardiametadata, by = "sampleid")
unweighted_interaction <- as.dist(unweighted_interaction[,2:49])

adonis(formula = unweighted_interaction ~ humaninteraction+group+giardiapresence, data = giardia_metadata, permutations = 5000)

adonis(formula = unweighted_interaction ~ humaninteraction*giardiapresence+group, data = giardia_metadata, permutations = 5000)

# filtered weighted Unifrac
# Rural
weighted_rural <- read.table("weighted-rural-distance-matrix.tsv", header=T)
weighted_rural_metadata = inner_join(weighted_rural, giardiametadata, by = "sampleid")
weighted_rural<-as.dist(weighted_rural[,2:28])

adonis(formula = weighted_rural ~ group+year+giardiapresence, data = weighted_rural_metadata, permutations = 5000)

# Village
weighted_village <- read.table("weighted-village-distance-matrix.tsv", header=T)
weighted_village_metadata = inner_join(weighted_village, giardiametadata, by = "sampleid")
weighted_village<-as.dist(weighted_village[,2:23])

adonis(formula = weighted_village ~ group+year+giardiapresence, data = weighted_village_metadata, permutations = 5000)

# Across habitats
weighted_interaction <- read.table("weighted-2017-distance-matrix.tsv", header = T)
giardia_metadata = inner_join(weighted_interaction, giardiametadata, by = "sampleid")
weighted_interaction <- as.dist(weighted_interaction[,2:49])

adonis(formula = weighted_interaction~humaninteraction+group+giardiapresence, data = giardia_metadata, permutations = 5000)

adonis(formula = weighted_interaction ~ humaninteraction*giardiapresence+group, data = giardia_metadata, permutations = 5000)

# Pairwise comparisons
library(pairwiseAdonis)

# All sites, unweighted UniFrac
pairwise.adonis(unweighted_bothyears, factors=metadata_both$humaninteraction, perm = 5000, p.adjust.m = 'holm')

# All sites, weighted UniFrac
pairwise.adonis(weighted_bothyears, factors=metadata_together$humaninteraction, perm = 5000, p.adjust.m = 'holm')

# Linear mixed models
library(nlme)
library(multcomp)

# Alpha diversity
alpha<-read.csv("alpha-diversity-2017.csv")

faith_alpha=lme(fixed=faith_pd~giardiapresence, data=alpha, random=~1|sampleid)
summary(faith_alpha)
anova(faith_alpha)
summary(glht(faith_alpha,linfct=mcp(giardiapresence="Tukey")))

faith_group=lme(fixed=faith_pd~group, data=alpha, random=~1|sampleid)
summary(faith_group)
anova(faith_group)
summary(glht(faith_group,linfct=mcp(group="Tukey")))

faith_interaction_alpha=lme(fixed=faith_pd~humaninteraction, data=alpha, random=~1|sampleid)
summary(faith_interaction_alpha)
anova(faith_interaction_alpha)
summary(glht(faith_interaction_alpha,linfct=mcp(humaninteraction="Tukey")))

shannon_alpha=lme(fixed=shannon~giardiapresence, data=alpha, random=~1|sampleid)
summary(shannon_alpha)
anova(shannon_alpha)
summary(glht(shannon_alpha,linfct=mcp(giardiapresence="Tukey")))

shannon_group=lme(fixed=shannon~group, data = alpha, random = ~1|sampleid)
summary(shannon_group)
anova(shannon_group)
summary(glht(shannon_group,linfct=mcp(group="Tukey")))

shannon_interaction_alpha=lme(fixed=shannon~humaninteraction, data=alpha, random=~1|sampleid)
summary(shannon_interaction_alpha)
anova(shannon_interaction_alpha)
summary(glht(shannon_interaction_alpha,linfct=mcp(humaninteraction="Tukey")))

otu_alpha=lme(fixed=observed_otus~giardiapresence, data=alpha, random=~1|sampleid)
summary(otu_alpha)
anova(otu_alpha)
summary(glht(otu_alpha,linfct=mcp(giardiapresence="Tukey")))

otu_group=lme(fixed=observed_otus~group, data=alpha, random=~1|sampleid)
summary(otu_group)
anova(otu_group)
summary(glht(otu_group,linfct=mcp(group="Tukey")))

otu_interaction_alpha=lme(fixed=observed_otus~humaninteraction, data=alpha, random=~1|sampleid)
summary(otu_interaction_alpha)
anova(otu_interaction_alpha)
summary(glht(otu_interaction_alpha,linfct=mcp(humaninteraction="Tukey")))

# Taxa 
phyla=read.csv("giardia-2017-phyla.csv", header=TRUE)
family<-read.csv("giardia-2017-family.csv", header=T)
genus<-read.csv("giardia-2017-genus.csv", header=T)

acti=lme(fixed=Actinobacteria~giardiapresence, data=phyla, random=~1|sampleid)
summary(acti)
anova(acti)
summary(glht(acti,linfct=mcp(giardiapresence="Tukey")))

acti_humaninteraction=lme(fixed=Actinobacteria~humaninteraction, data=phyla, random=~1|sampleid)
summary(acti_humaninteraction)
anova(acti_humaninteraction)
summary(glht(acti_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

bact=lme(fixed=Bacteroidetes~giardiapresence, data=phyla, random=~1|sampleid)
summary(bact)
anova(bact)
summary(glht(bact,linfct=mcp(giardiapresence="Tukey")))

bact_humaninteraction=lme(fixed=Bacteroidetes~humaninteraction, data=phyla, random=~1|sampleid)
summary(bact_humaninteraction)
anova(bact_humaninteraction)
summary(glht(bact_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

prot=lme(fixed=Proteobacteria~giardiapresence, data=phyla, random=~1|sampleid)
summary(prot)
anova(prot)
summary(glht(prot,linfct=mcp(giardiapresence="Tukey")))

prot_humaninteraction=lme(fixed = Proteobacteria~humaninteraction, data=phyla, random=~1|sampleid)
summary(prot_humaninteraction)
anova(prot_humaninteraction)
summary(glht(prot_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

firm=lme(fixed=Firmicutes~giardiapresence, data=phyla, random=~1|sampleid)
summary(firm)
anova(firm)
summary(glht(firm,linfct=mcp(giardiapresence="Tukey")))

firm_humaninteraction=lme(fixed=Firmicutes~humaninteraction, data=phyla, random=~1|sampleid)
summary(firm_humaninteraction)
anova(firm_humaninteraction)
summary(glht(firm_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

prev_giardia=lme(fixed=Prevotellaceae~giardiapresence, data=family, random=~1|sampleid)
summary(prev_giardia)
anova(prev_giardia)
summary(glht(prev_giardia, linfct=mcp(giardiapresence="Tukey")))

prev_humaninteraction=lme(fixed=Prevotellaceae~humaninteraction, data=family, random=~1|sampleid)
summary(prev_humaninteraction)
anova(prev_humaninteraction)
summary(glht(prev_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

para_giardia=lme(fixed=Paraprevotellaceae~giardiapresence, data=family, random=~1|sampleid)
summary(para_giardia)
anova(para_giardia)
summary(glht(para_giardia, linfct=mcp(giardiapresence="Tukey")))

para_humaninteraction=lme(fixed=Paraprevotellaceae~humaninteraction, data=family, random=~1|sampleid)
summary(para_humaninteraction)
anova(para_humaninteraction)
summary(glht(para_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

lach_giardia=lme(fixed=Lachnospiraceae~giardiapresence, data=family, random=~1|sampleid)
summary(lach_giardia)
anova(lach_giardia)
summary(glht(lach_giardia, linfct=mcp(giardiapresence="Tukey")))

lach_humaninteraction=lme(fixed=Lachnospiraceae~humaninteraction, data=family, random=~1|sampleid)
summary(lach_humaninteraction)
anova(lach_humaninteraction)
summary(glht(lach_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

rumi_giardia=lme(fixed=Ruminococcaceae~giardiapresence, data=family, random=~1|sampleid)
summary(rumi_giardia)
anova(rumi_giardia)
summary(glht(rumi_giardia, linfct=mcp(giardiapresence="Tukey")))

rumi_humaninteraction=lme(fixed=Ruminococcaceae~humaninteraction, data=family, random=~1|sampleid)
summary(rumi_humaninteraction)
anova(rumi_humaninteraction)
summary(glht(rumi_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

prevo_giardia=lme(fixed=Prevotella~giardiapresence, data=genus, random=~1|sampleid)
summary(prevo_giardia)
anova(prevo_giardia)
summary(glht(prevo_giardia,linfct=mcp(giardiapresence="Tukey")))

prevo_humaninteraction=lme(fixed=Prevotella~humaninteraction, data=genus, random=~1|sampleid)
summary(prevo_humaninteraction)
anova(prevo_humaninteraction)
summary(glht(prevo_humaninteraction,linfct=mcp(humaninteraction="Tukey")))

