# Beta diversity analysis

- [Load libraries](#load-libraries)
- [Make a phyloseq object](#make-a-phyloseq-object)
- [CSS transformation](#css-transformation)
- [Edit metadata file](#edit-metadata-file)
- [Plot all variables with bray curtis](#plot-all-variables-with-bray-curtis)
- [PERMANOVA (bray curtis)](#permanova--bray-curtis-)
  * [Model 1 location + sex](#model-1-location---sex)
  * [Model 2 location + smi](#model-2-location---smi)
  * [Model 3 location + infection](#model-3-location---infection)
  * [Adjust p values BH correction](#adjust-p-values-bh-correction)
- [Plot all variables with Weighted UniFrac](#plot-all-variables-with-weighted-unifrac)
- [PERMANOVA (wheighted unifrac)](#permanova--wheighted-unifrac-)
  * [Model 1 location + sex](#model-1-location---sex-1)
  * [Model 2 location + smi](#model-2-location---smi-1)
  * [Model 3 location + infection](#model-3-location---infection-1)
  * [Adjust p values BH correction](#adjust-p-values-bh-correction-1)
- [Interaction models](#interaction-models)
  * [Model  location*smi](#model--location-smi)
  * [Model  location*infection](#model--location-infection)

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
metadata$location <- gsub("IBISS, yard", "IBISS", metadata$location)
metadata$location <- gsub("Kalemegdan, ZOO", "ZOO", metadata$location)
metadata$identifier <- as.factor(metadata$identifier)
metadata$ring_number <- as.factor(metadata$ring_number)
metadata$location <- as.factor(metadata$location)
metadata$sex <- as.factor(metadata$sex)
metadata$m_infection <- as.factor(metadata$m_infection)
metadata$smi <- as.numeric(metadata$smi)
metadata$std_smi <- as.numeric(metadata$std_smi)
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

#bray vs location
bray_location <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "location",title = "PCOA")
bray_location + geom_point(size=4)+theme_classic()

### bray vs sex
bray_sex <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "sex",title = "PCOA")
bray_sex + geom_point(size=4)

### bray vs smi
bray_smi <- plot_ordination(physeq = ps_css, ordination = data_pcoa_bray_css, color = "std_smi",title = "PCOA")
plot(bray_smi)
bray_smi + geom_point(size=4)

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



### Model 1 location + sex

```R
#Run PERMANOVA
perm_location_sex <- adonis2(asv_table ~ location + sex  , data = metadata, permutations = perm1, method = "bray", by= "margin", diag=TRUE)
summary(perm_location_sex)
p1 <- bc_perm_location_sex

# Analysis of dispersion
betadisp_location <-betadisper(BC.dist, metadata$location) # analysis of dispersion (non significant)
anova(betadisp_location)
permutest(betadisp_location, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

betadisp_sex <-betadisper(BC.dist, metadata$sex) # analysis of dispersion (non significant)
anova(betadisp_sex)
permutest(betadisp_sex, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

plot(betadisp_location)

#Plot results 
my_cols <- c("#4CAF50", "#D32F2F")
plot(betadisp_location, col = my_cols)

par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
plot(betadisp_location, ellipse = FALSE, hull = FALSE, col = my_cols,  pch = c(16,17), cex = 2.2,seg.lty = "dashed",  main = "")# 1 sd data ellipse
```



### Model 2 location + smi
```R
#Run PERMANOVA
bc_perm_location_smi<- adonis2(asv_table ~ location + std_smi  , data = metadata, permutations = perm1, method = "bray", by= "margin")
summary(perm_location_smi)
p2 <- bc_perm_location_smi

# Analysis of dispersion
betadisp_smi <-betadisper(BC.dist, metadata$std_smi) # analysis of dispersion (non significant)
anova(betadisp_smi)
permutest(betadisp_smi, pairwise = FALSE, permutations = 9999) # same as betadisp but with permutations

plot(betadisp_smi)
```



### Model 3 location + infection

```R
bc_perm_location_infection<- adonis2(asv_table ~ location + m_infection  , data = metadata, 
                   permutations = perm1, method = "bray", by= "margin")
p3 <- bc_perm_location_infection

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

#wu vs location
wu_location <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "location",title = "PCOA")
wu_location + geom_point(size=4)+theme_classic()

#wu vs sex
wu_sex <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "sex",title = "PCOA")
wu_sex + geom_point(size=4)

#wu vs smi
wu_smi <- plot_ordination(physeq = ps_css, ordination = data_pcoa_wu_css, color = "std_smi",title = "PCOA")
plot(wu_smi)
wu_smi + geom_point(size=4)

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



### Model 1 location + sex

```R
#Run PERMANOVA
wu_perm_location_sex <- adonis2(dist_wu ~ location + sex  , data = metadata, permutations = perm1, by= "margin")
p1 <- wu_perm_location_sex

# Analysis of dispersion
wu_betadisp_location <-betadisper(dist_wu, metadata$location) # analysis of dispersion (non significant)
anova(wu_betadisp_location)
permutest(wu_betadisp_location, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations

wu_betadisp_sex <-betadisper(dist_wu, metadata$sex) # analysis of dispersion (non significant)
anova(wu_betadisp_sex)
permutest(wu_betadisp_sex, pairwise = TRUE, permutations = 9999) # same as betadisp but with permutations


#Plot the resutls
my_cols <- c("#4CAF50", "#D32F2F")
plot(wu_betadisp_location, col = my_cols)

par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
plot(wu_betadisp_location, ellipse = FALSE, hull = FALSE, col = my_cols,  pch = c(16,17), cex = 2.2,seg.lty = "dashed",  main = "")
```



### Model 2 location + smi
```R
#Run PERMANOVA
wu_perm_location_smi<- adonis2(dist_wu ~ location + std_smi  , data = metadata, permutations = perm1, by= "margin")
p2 <- wu_perm_location_smi

# Analysis of dispersion
wu_betadisp_smi <-betadisper(dist_wu, metadata$std_smi) # analysis of dispersion (non significant)
anova(wu_betadisp_smi)
permutest(wu_betadisp_smi, pairwise = FALSE, permutations = 9999) # same as betadisp but with permutations
plot(wu_betadisp_smi)
```



### Model 3 location + infection

```R
#Run PERMANOVA
wu_perm_location_infection<- adonis2(dist_wu ~ location + m_infection  , data = metadata, permutations = perm1, by= "margin")
p3 <- wu_perm_location_infection

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

### Model  location*smi
```R
inter_location_smi<- adonis2( asv_table ~ location*std_smi  , data = metadata, permutations = perm1, method = "bray", by= "margin")
inter_location_smi

wu_inter_location_smi<- adonis2(dist_wu ~ location*std_smi  , data = metadata, permutations = perm1, by= "margin")
wu_inter_location_smi
```



### Model  location*infection

```R
inter_location_infection<- adonis2( asv_table ~ location*m_infection  , data = metadata, permutations = perm1, method = "bray", by= "margin")
inter_location_infection

wu_inter_location_infection<- adonis2(WU.dist ~ location*m_infection  , data = metadata, permutations = perm1, by= "margin")
wu_inter_location_infection
```

