# pondeR

## Description

Nonsense-Mediated Decay (NMD) is an RNA surveillance mechanism which eliminates 
abberant messenger RNAs (mRNAs) harbouring a premature stop codon (PTC). pondeR
predicts the sensitivity of protein-coding mRNAs to NMD by scanning for NMD-inducing
features. pondeR can analyze mRNAs from any transcript annotation database or 
from high-throughput RNA-seq experiments

## Features
* Scans for NMD-inducing features on coding RNAs from various sources:
  * Most transcript annotation databases
  * Custom transcript annotation (GTF) generated from RNA-seq experiemnts
* Predict coding region (CDS) on custom-generated GTF annotation using reference CDS as guide
* Match chromosome levels and gene_id levels between query and reference GRanges object
* Resize start and end of GRanges object containing exons of transcripts
* *In development*
  * Analyze alternative segments between query and reference transcripts
  * Analyze uORFs and uATGs on a list of coding mRNAs



## Installation
```r
# install.packages("devtools")
devtools::install_github("fursham-h/ponder")
```

## Usage
pondeR is pre-loaded with GenomicRangesList objects containing sample
exon and cds ranges from 13 transcripts from a gene.
```r
library(pondeR)

names(query_exons)
#[1] "transcript1"  "transcript10" "transcript11" "transcript12" "transcript13" 
#[6] "transcript2"  "transcript3" "transcript4"  "transcript5"  "transcript6"  
#[11] "transcript7"  "transcript8"  "transcript9" 

names(query_cds)
#[1] "transcript1"  "transcript11" "transcript12" "transcript13" "transcript2"  
#[6] "transcript3"  "transcript4" "transcript9" 

query_exons[1]
# GRangesList object of length 1:
# $transcript1 
# GRanges object with 14 ranges and 0 metadata columns:
#        seqnames              ranges strand
#           <Rle>           <IRanges>  <Rle>
#    [1]     chr3 119783202-119783388      -
#    [2]     chr3 119781504-119781534      -
#    [3]     chr3 119761703-119761778      -
#    [4]     chr3 119752944-119753116      -
#    [5]     chr3 119751864-119752007      -
#    ...      ...                 ...    ...
#   [10]     chr3 119724612-119724645      -
#   [11]     chr3 119724094-119724186      -
#   [12]     chr3 119720789-119721005      -
#   [13]     chr3 119720607-119720684      -
#   [14]     chr3 119718742-119720460      -
# 
# -------
# seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

Transcripts without cds information can be part of a GRangesList but will not
be analysed for NMD-inducing features. Transcript architecture of these 
transcripts can be visualized using _wiggleplotr_ package
```r
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("wiggleplotr")

library(wiggleplotr)
plotTranscripts(query_exons, query_cds)
```
<img src="wiggleplot_query.png" width="800">

Run predictNMD with exon and cds information as input
```r
predictNMD(query_exons, query_cds)
## A tibble: 2 x 6
#  tx          is_NMD dist_to_lastEJ num_of_down_EJs dist_to_downEJs threeUTRlength
#  <chr>       <lgl>           <int>           <dbl> <chr>                    <dbl>
#1 transcript3 TRUE              361               3 66,283,361                 502
#2 transcript4 TRUE              361               3 66,283,361                 474
```

