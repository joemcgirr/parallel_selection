######################
# parallel selection #
######################
# Use this script to identify and visualize genes shwoing
# signs of a hard selective sweep at user defined thresholds


# 1. set working directory containing:
#
#     SweeD_Report.SPECIES_NAME_grid1000_sweeps.txt
#          

required_packages <- c("BiocManager","tidyverse","plyranges","seqnir","UpSetR")

for (required_package in required_packages){

if (!requireNamespace(required_package, quietly = TRUE))
  install.packages(required_package)
}

library(tidyverse)
library(plyranges)
library(seqinr)
library(UpSetR)



for (species_name in species_names){
  sweed <- read.table(paste("SweeD_Report.",species_name,"_grid1000_sweeps.txt",sep = ""), header = FALSE, stringsAsFactors = FALSE)
  names(sweed) <- c("seqnames","start","CLR","alpha")
  sweed$start <- round(sweed$start)
  # set CLR_window_size to remove small scaffolds.
  # chroms split into 1000 equal chunks, so CLR_window_size > 100 filters chroms smaller than 100kb
  sweed <- sweed %>% 
    group_by(seqnames) %>%
    mutate(CLR_window_size = c(0,diff(start))) %>% filter(CLR_window_size>100) %>%
    as.data.frame()
  unique(paste(sweed$seqnames, sweed$CLR_window_size, sep = ":"))
  sweed$end <- sweed$start + sweed$CLR_window_size
  sweed <- sweed[c("seqnames","start","end","CLR","alpha","CLR_window_size")]
  
  
  # get best candidate for each chromosome
  #candidate_outliers <- sweed %>% 
  #  group_by(seqnames) %>%
  #  mutate(highest_CLR = max(CLR)) %>%
  #  group_by(seqnames) %>%
  #  filter(CLR == highest_CLR) %>%
  #  as.data.frame()
  
  # get top 90 percentile for each chromosome
  candidate_outliers <- sweed %>% 
    group_by(seqnames) %>%
    mutate(highest_CLR = quantile(CLR,.90)) %>%
    group_by(seqnames) %>%
    filter(CLR >= highest_CLR) %>%
    as.data.frame()
  
  head(sweed)
  head(candidate_outliers)
  
  prot <- read.delim(paste("ProteinTable_",species_name,".csv",sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = ",")
  # only ovelap with unique proteins (only appear in protein table once)
  prot <-prot[!duplicated(prot$Protein.name),]
  head(prot)
  prot$seqnames <- prot$SeqName2a
  prot$start <- prot$Start
  prot$end <- prot$Stop
  prot$strand <- prot$Strand
  prot <- prot[c("seqnames", "start","end","strand","Protein.product","Protein.name","GenBank.Accn","RefSeq.Accn")]
  
  prot <- prot %>% as_granges()
  candidate_outliers <- candidate_outliers %>% as_granges()
  
  candidate_genes <- join_overlap_intersect(prot,candidate_outliers) %>% as.data.frame()
  candidate_genes$species <- species_name
  head(candidate_genes) 
  
  write.table(candidate_genes, paste("C:/Users/jmcgirr/Documents/matute/candidates/",species_name,"_candidate_genes_top_90_CLR_for_scaffold_unique_proteins.txt",sep = ""), row.names = FALSE, quote= FALSE,sep="\t")
}

