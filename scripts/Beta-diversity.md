# Beta diversity analysis

[TOC]

## Load libraries

```R
library(qiime2R)
library(lubridate)
library(openxlsx)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(metagenomeSeq)
library(vegan)
```


## Make a phyloseq object

```R
#Metadata with alpha metrics
metadata <- readRDS("metadata_sparrows.rds")

# Create phyloseq object
ps <- qza_to_phyloseq(
  features="filtered-table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "sparrows-meta.tsv")

#Extract taxonomy 
taxonomy <- as.data.frame(tax_table(ps))

#Edit taxonomy file (for some reason Kingdom name comes with "d_" before)
taxonomy$Kingdom <- gsub("d__","",as.character(taxonomy$Kingdom))
taxonomy <- as.matrix(taxonomy)

# Extract phylogeny file
tree <- phy_tree(ps)

```



## CSS transformation

```R
#Convert the phyloseq object to a metagenomeSeq object (MRexperiment)
meta.obj <- phyloseq_to_metagenomeSeq(ps)

#Normalise counts
meta.obj <- cumNorm(meta.obj, p = cumNormStatFast(meta.obj))

#Convert CSS data into data.frame-formatted OTU table (log transformed data)
asv_table_css <- MRcounts(meta.obj, norm = TRUE, log = TRUE)

## Make a new phyloseq object with with the new CSS transformed ASV table----
asv_table_css <- otu_table(asv_table_css, taxa_are_rows = TRUE)
taxonomy <- tax_table(taxonomy)
metadata <- sample_data(ps)
tree <- phy_tree(tree)

ps_css <- phyloseq(asv_table_css, taxonomy, metadata, tree)
otu <- as.data.frame(otu_table(ps_css))

#Save phyloseq object as rds 
saveRDS(ps_css, "phyloseq_css.rds")


# Explore ps object
summarize_phyloseq(ps_css)
sample_names(ps_css) # looks at the sample names on the phyloseq object
meta(ps_css) # retrieves the metadata file
sample_data(ps_css) # retrieves the metadata file
taxa(ps) # retrieves taxa name (ASV_1, ASV_2...etc)
abundances(ps_css) # retrieves ASV counts table
abundances(ps_css, "compositional") # computes relative abundaces
readcount(ps_css) # number of reads per sample
```



## Edit metadata file

```R
ps_css <- readRDS("phyloseq_css.rds")
metadata <- sample_data(ps_css)
metadata <- clean_names(metadata)
metadata$place <- gsub("IBISS, yard", "IBISS", metadata$place)
metadata$place <- gsub("Kalemegdan, ZOO", "ZOO", metadata$place)
metadata$identifier <- as.factor(metadata$identifier)
metadata$ring_number <- as.factor(metadata$ring_number)
metadata$place <- as.factor(metadata$place)
metadata$sex <- as.factor(metadata$sex)
metadata$m_infection <- as.factor(metadata$m_infection)
metadata$bmi <- as.numeric(metadata$bmi)
metadata$std_bmi <- as.numeric(metadata$std_bmi)
metadata$parasetimia <- as.numeric(metadata$parasetimia)

# re build the phyloseq object
asv_table <- otu_table(ps_css)
taxonomy <- tax_table(ps_css)
tree <- phy_tree(ps_css)
asv <- otu_table(ps_css)

ps_css <- phyloseq(asv, taxonomy, metadata, tree)
saveRDS(ps_css,"phyloseq_css.rds")

```



## Plot all variables with bray curtis

```R
#Calculate Bray-Curtis distance using the vegan package
dist_bc <- as.matrix(vegdist(otu_table(ps_css), method="bray"))

set.seed(1234)
data_pcoa_bray_css  <- ordinate(physeq = ps_css, method = "PCoA", distance = "bray") # compute pcoa

#bray vs place
bray_place <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "place",title = "PCOA")
bray_place + geom_point(size=4)+theme_classic()

### bray vs sex
bray_sex <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "sex",title = "PCOA")
bray_sex + geom_point(size=4)

### bray vs bmi
bray_bmi <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "std_bmi",title = "PCOA")
plot(bray_bmi)
bray_bmi + geom_point(size=4)

### bray vs m_infection
bray_infection <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "m_infection",title = "PCOA")
bray_infection + geom_point(size=4)

### bray vs parasetimia
bray_parasetimia <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "parasetimia",title = "PCOA")
bray_parasetimia + geom_point(size=4)

```


## PERMANOVA (bray curtis)
```R
asv_table <- as.data.frame(otu_table(ps_css))
asv_table <- t (asv_table)
metadata <- data.frame(sample_data(ps_css))

BC.dist=vegdist(asv_table, distance="bray")
perm1 <- how(nperm = 9999)
set.seed(1234)
```



### Model 1 place + sex

```R
#Run PERMANOVA
perm_place_sex <- adonis2(asv_table ~ place + sex  , data = metadata, permutations = perm1, method = "bray", by= "margin", diag=TRUE)
summary(perm_place_sex)
p1 <- bc_perm_place_sex

# Analysis of dispersion
betadisp_place <-betadisper(BC.dist, metadata$place) # analysis of dispersion (non significant)
anova(betadisp_place)
permutest(betadisp_place, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

betadisp_sex <-betadisper(BC.dist, metadata$sex) # analysis of dispersion (non significant)
anova(betadisp_sex)
permutest(betadisp_sex, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

plot(betadisp_place)

#Plot results 
my_cols <- c("#4CAF50", "#D32F2F")
plot(betadisp_place, col = my_cols)

par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
plot(betadisp_place, ellipse = FALSE, hull = FALSE, col = my_cols,  pch = c(16,17), cex = 2.2,seg.lty = "dashed",  main = "")# 1 sd data ellipse
```



### Model 2 place + bmi
```R
#Run PERMANOVA
bc_perm_place_bmi<- adonis2(asv_table ~ place + std_bmi  , data = metadata, permutations = perm1, method = "bray", by= "margin")
summary(perm_place_bmi)
p2 <- bc_perm_place_bmi

# Analysis of dispersion
betadisp_bmi <-betadisper(BC.dist, metadata$std_bmi) # analysis of dispersion (non significant)
anova(betadisp_bmi)
permutest(betadisp_bmi, pairwise = FALSE, permutations = 9999) # same as betadisp but with permutations

plot(betadisp_bmi)
```



### Model 3 place + infection

```R
bc_perm_place_infection<- adonis2(asv_table ~ place + m_infection  , data = metadata, 
                   permutations = perm1, method = "bray", by= "margin")
p3 <- bc_perm_place_infection

betadisp_infection <-betadisper(BC.dist, metadata$m_infection) # analysis of dispersion (non significant)
anova(betadisp_infection)
permutest(betadisp_infection, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations
plot(betadisp_infection)
```



### Adjust p values BH correction

```R
p_adjust_bray <- c(p1$`Pr(>F)`[1], p2$`Pr(>F)`[1], p3$`Pr(>F)`[1])

p_adjust_bray <- p.adjust(p_adjust_bray, method = "BH")
p_adjust_bray <- format(round(p_adjust_bray, digits = 3), scientific = FALSE)
p_adjust_bray
```



## Plot all variables with Weighted UniFrac

```R
#Calculate Wheighted UniFrac distance using the vegan package
dist_wu <- distance(ps_css, method = "wunifrac", type = "samples")

set.seed(1234)
data_pcoa_wu_css  <- ordinate(physeq = ps_css, method = "PCoA", distance = "wunifrac") # compute pcoa

#wu vs place
wu_place <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "place",title = "PCOA")
wu_place + geom_point(size=4)+theme_classic()

#wu vs sex
wu_sex <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "sex",title = "PCOA")
wu_sex + geom_point(size=4)

#wu vs bmi
wu_bmi <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "std_bmi",title = "PCOA")
plot(wu_bmi)
wu_bmi + geom_point(size=4)

#wu vs m_infection
wu_infection <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "m_infection",title = "PCOA")
wu_infection + geom_point(size=4)

#wu vs parasetimia
wu_parasetimia <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "parasetimia",title = "PCOA")
wu_parasetimia + geom_point(size=4)
```



## PERMANOVA (wheighted unifrac)

```R
asv_table <- as.data.frame(otu_table(ps_css))
asv_table <- t (asv_table)
metadata <- data.frame(sample_data(ps_css))

perm1 <- how(nperm = 9999)
set.seed(1234) 
```



### Model 1 place + sex

```R
#Run PERMANOVA
wu_perm_place_sex <- adonis2(dist_wu ~ place + sex  , data = metadata, permutations = perm1, by= "margin")
p1 <- wu_perm_place_sex

# Analysis of dispersion
wu_betadisp_place <-betadisper(dist_wu, metadata$place) # analysis of dispersion (non significant)
anova(wu_betadisp_place)
permutest(wu_betadisp_place, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

wu_betadisp_sex <-betadisper(dist_wu, metadata$sex) # analysis of dispersion (non significant)
anova(wu_betadisp_sex)
permutest(wu_betadisp_sex, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations


#Plot the resutls
my_cols <- c("#4CAF50", "#D32F2F")
plot(wu_betadisp_place, col = my_cols)

par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
plot(wu_betadisp_place, ellipse = FALSE, hull = FALSE, col = my_cols,  pch = c(16,17), cex = 2.2,seg.lty = "dashed",  main = "")
```



### Model 2 place + bmi
```R
#Run PERMANOVA
wu_perm_place_bmi<- adonis2(dist_wu ~ place + std_bmi  , data = metadata, permutations = perm1, by= "margin")
p2 <- wu_perm_place_bmi

# Analysis of dispersion
wu_betadisp_bmi <-betadisper(dist_wu, metadata$std_bmi) # analysis of dispersion (non significant)
anova(wu_betadisp_bmi)
permutest(wu_betadisp_bmi, pairwise = FALSE, permutations = 9999) # same as betadisp but with permutations
plot(wu_betadisp_bmi)
```



### Model 3 place + infection

```R
#Run PERMANOVA
wu_perm_place_infection<- adonis2(dist_wu ~ place + m_infection  , data = metadata, permutations = perm1, by= "margin")
p3 <- wu_perm_place_infection

# Analysis of dispersion
wu_betadisp_infection <-betadisper(dist_wu, metadata$m_infection) # analysis of dispersion (non significant)
anova(wu_betadisp_infection)
permutest(wu_betadisp_infection, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations
plot(wu_betadisp_infection)
```



### Adjust p values BH correction

```R
p_adjust_wu <- c(p1$`Pr(>F)`[1], p2$`Pr(>F)`[1], p3$`Pr(>F)`[1])

p_adjust_wu <- p.adjust(p_adjust_wu, method = "BH")
p_adjust_wu <- format(round(p_adjust_wu, digits = 3), scientific = FALSE)
p_adjust_wu
```



## Interaction models

### Model  place*bmi
```R
inter_place_bmi<- adonis2( asv_table ~ place*std_bmi  , data = metadata, permutations = perm1, method = "bray", by= "margin")
inter_place_bmi

wu_inter_place_bmi<- adonis2(dist_wu ~ place*std_bmi  , data = metadata, permutations = perm1, by= "margin")
wu_inter_place_bmi
```



### Model  place*infection

```R
inter_place_infection<- adonis2( asv_table ~ place*m_infection  , data = metadata, permutations = perm1, method = "bray", by= "margin")
inter_place_infection

wu_inter_place_infection<- adonis2(WU.dist ~ place*m_infection  , data = metadata, permutations = perm1, by= "margin")
wu_inter_place_infection
```

