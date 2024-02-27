# Sequence data processing





## Table of Contents

[TOC]

------



## 1.  Activate Qiime2

```bash
# Activate base env
. ~/.bashrc

# Activate qiime2
conda activate qiime2-2022.11
```



## 2. Import sequences

```python
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path seqs --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza 
```



## 3. Visualize quality plots

```python
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-quality-plots.qzv 
```



## 4. Exit qiime

```python
conda deactivate
```



## 5. In R run cutadapt pipeline

```python
Rscript cutadapt.R
```



## 6.  Import trimmed sequences to qiime2

```bash
# Activate qiime2
conda activate qiime2-2022.11

# Import sequences
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path seqs/cutadapt --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path primer-trimmed-demux.qza

```



## 7. Generate quality plots 2nd time
```python
qiime demux summarize --i-data primer-trimmed-demux.qza --o-visualization primer-trimed-demux-quality-plots.qzv
```



## 8. Run dada2 

```python
qiime dada2 denoise-paired --i-demultiplexed-seqs primer-trimmed-demux.qza --p-trunc-len-f 253 --p-trunc-len-r 185  --p-trunc-q 2 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza
```



### 8.1 DADA2 results

```python
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv | 

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv | 

qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file buzzard_meta.tsv
```



## 9. Taxonomy assignment

```python
qiime feature-classifier classify-sklearn --i-classifier silva-138.1-SSU-nr99-515F-806R-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza 
```



### 10. Taxonomy visualisation

```python
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
```



## 11. Exit qiime2

```python
conda deactivate 
```



## 12. In R run decontam

### 12.1. Read in the data

```R
# Load Libraries
library(decontam)
library(qiime2R)
library(phyloseq)
library(biomformat)

# Make a phyloseq object
ps <- qza_to_phyloseq(features = "table.qza", taxonomy = "taxonomy.qza", metadata = "buzzard_meta.tsv")

# Choose which samples are the negative controls 
sample_data(ps)$is.neg <- sample_data(ps)$type == "negative"
```



### 12.2. Run decontam 

```R
# Identify contaminants based on prevalence method (treshold 0.1 is the standard)
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Remove contaminants from the phyloseq object
ps.nocontam <- prune_taxa(!contamdf.prev$contaminant,ps)
```



### 12.3. Export feature table as biom file

```R
# Extract asv table from the phyloseq object
table_nocontam <- as(otu_table(ps.nocontam),"matrix",) 

#'t' to transform if taxa_are_rows=FALSE table_nocontam<- t(as(otu_table(ps.nocontam),"matrix",))#if taxa_are_rows=TRUE

# Make a biom table
table_nocontam_biom <- make_biom(data=table_nocontam)
write_biom(table_nocontam_biom,"table-nocontam.biom")
```



## 13. Import biom table from R to qiime2

```python
conda activate qiime2-2022.11

qiime tools import --input-path table-nocontam.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path table-nocontam.qza
```



## 14. Remove control samples from dataset
```python
qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file sparrows_meta.tsv--p-where "type='positive'" --p-exclude-ids --o-filtered-table table-nocontam.qza

qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file sparrows_meta.tsv--p-where "type='negative'" --p-exclude-ids --o-filtered-table table-nocontam.qza
```



## 15. Taxonomy based filtering

```python
qiime taxa filter-table --i-table table-nocontam.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast,Unassigned,Vertebrata,Eukaryota --p-include p_ --o-filtered-table table-taxa-filter.qza

qiime feature-table summarize --i-table table-taxa-filter.qza --o-visualization table-taxa-filter.qzv --m-sample-metadata-file sparrows_meta.tsv
```



## 16. Filter unique features

```python
qiime feature-table filter-features --i-table table-taxa-filter.qza --p-min-samples 2 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file sparrows_meta.tsv
```



## 17. Filter samples with less than 500 reads

```python
qiime feature-table filter-samples --i-table table-taxa-filter-no_singles.qza --p-min-frequency 500 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file sparrows_meta.tsv
```



## 18. Filter representative sequences

```python
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table filtered-table.qza --o-filtered-data filter-seqs.qza
qiime feature-table tabulate-seqs --i-data filter-seqs.qza --o-visualization filter-seqs.qzv 
```



## 19. Building a phylogenetic tree

```python
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences filter-seqs.qza  --o-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-seqs.qza --o-tree unrooted-tree.qza  --o-rooted-tree rooted-tree.qza 
```



## 20. Rarefaction curves

```python
qiime diversity alpha-rarefaction --i-table filtered-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 14856 --m-metadata-file sparrows_meta.tsv --o-visualization alpha-rarefaction.qzv
```



## 21. Calculate alpha diversity metrics and rarefy the data-set

```python
# Calculates observed features, shannon diversity, Faith PD and retrieves rarefied table
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table filtered-table.qza --p-sampling-depth 14856 --m-metadata-file sparrows_meta.tsv--output-dir alpha-metrics-results

#Add alpha diversity metrics to the metadata
qiime metadata tabulate --m-input-file buzzard_metadata.tsv --m-input-file shannon_vector.qza --m-input-file observed_features.qza --m-input-file faith_pd_vector.qza --o-visualization sparrows_meta_alpha.qzv
```
