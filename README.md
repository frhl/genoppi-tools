# pRoteomics

## Description
Downstream analysis of proteomics data for iTRAQ (labelled data). Modular workflow, i.e. preprocess and analysis. Some of the functions will attempt to guess column nanmes if nothing have been specified. These tools have only been tested on data generated from the Whitehed Institution, Cambridge MA.

## workflow
The following steps can be accomplished in pipeline:

### Preprocess
using prepare.R:
- log2 transformation
- median normalization
- filter out proteins from unwanted species
- filter out protein which fragments have only been observed once.
- calculate log fold change

### Analysis
using genoppi.R
- create replicate plots
- create volcano plots
- look for protein-protein interactions
- user hypergeometric (fisher exact test) to test for enrichment

## note
Tools have not yet been fully tested.



## example on how to setup personal pipeline
library(rProteomics)

infile = 'data/raw/frederik/frederik_BCAS3_proteins.csv'
prelim <- prepare(c("EC", "BCAS3"), infile = infile, impute = list(shift = -0.8, stdwidth = 0.3))
known.interactors <- interactors("BCAS3") # get known interactors

data <- prelim %>% mttest()
data %>% designate(FDR < 0.1) %>% plotScatter(bait, paste(bait, '[FDR < 0.1]'))
data %>% designate(FDR < 0.1) %>% plotVolcano(bait)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors)
