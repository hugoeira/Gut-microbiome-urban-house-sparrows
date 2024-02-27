# CUTADAPT PIPELINE



- [1. Load libraries](#1-load-libraries)
- [2. Define path to fastq files](#2-define-path-to-fastq-files)
- [3. Generate matched lists of the forward and reverse read files](#3-generate-matched-lists-of-the-forward-and-reverse-read-files)
- [4. Identify Primers](#4-identify-primers)
- [5. Verify the presence and orientation of the primers in the data](#5-verify-the-presence-and-orientation-of-the-primers-in-the-data)
  * [5.1 Create all orientantions of the primers](#51-create-all-orientantions-of-the-primers)
- [5.3. Count the number of times the primers appear in the forward and reverse reads](#53-count-the-number-of-times-the-primers-appear-in-the-forward-and-reverse-reads)
- [5.4. Remove primers](#54-remove-primers)
  * [5.4.1 Call cutadapt in R](#541-call-cutadapt-in-r)
  * [5.4.2 Create output filenames for the cutadapt-ed files. Define the parameters of cutadapt](#542-create-output-filenames-for-the-cutadapt-ed-files-define-the-parameters-of-cutadapt)



------



## 1. Load libraries

```R
library(dada2)
library(ShortRead)
library(Biostrings)
```



## 2. Define path to fastq files

```R
path_seqs <- "/path-to-seqs"  ## CHANGE to the directory containing the fastq files.
list.files(path_seqs) 
```



## 3. Generate matched lists of the forward and reverse read files

```R
fnFs <- sort(list.files(path_seqs, pattern = "_R1_001.fastq.gz", full.names = TRUE)) # Just select forward read files
fnRs <- sort(list.files(path_seqs, pattern = "_R2_001.fastq.gz", full.names = TRUE)) # Just select reverse read files 
```



## 4. Identify Primers

```R
FWD <- "GTGYCAGCMGCCGCGGTAA" ## INPUT forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ## INPUT reverse primer sequence
```



## 5. Verify the presence and orientation of the primers in the data



### 5.1 Create all orientantions of the primers

```R
allorients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allorients(FWD)
REV.orients <- allorients(REV)
FWD.orients
REV.orients
```




## 5.3. Count the number of times the primers appear in the forward and reverse reads
Considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient,assuming all the files were created using the same library preparation, so we’ll just process the first sample.

```R
primerHits <- function(primer, fn) {

  # Counts number of reads in which the primer is found

  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[4]]))
      
      
#Note: Orientation mixups are a common trip-up. If, for example, the REV primer is matching the Reverse reads 
#in its RevComp orientation, then replace REV with its reverse-complement orientation 
#(REV <- REV.orient[["RevComp"]]) before proceeding.
```



## 5.4. Remove primers

### 5.4.1 Call cutadapt in R
```R
cutadapt <- "/vol/cluster-data/hugoeira/miniconda3/bin/cutadapt" # CHANGE to cutadapt path 
system2(cutadapt, args = "--version") # Run shell commands from R
```



### 5.4.2 Create output filenames for the cutadapt-ed files. Define the parameters of cutadapt

```R
# Create output files
path.cut <- file.path(path_seqs, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Cutadapt parameters
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-m 50", "--match-read-wildcards","--trim-n", "--discard-untrimmed", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# Sanity check. Look for primers in cutadapt-ed files (should be 0)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[4]]))
```

