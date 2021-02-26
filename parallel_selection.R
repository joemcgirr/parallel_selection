######################
# parallel selection -------------------------------------------------------------------------------
######################
#
# Use this script to identify and visualize genes shwoing
# signs of a hard selective sweep at user defined
# significance thresholds
#
####################################################
# 1. set working directory containing:
#
#     SweeD_Report.<SPECIES_NAME>_grid1000_sweeps.txt
#     ProteinTable_<SPECIES_NAME>.csv"
#
####################################################

setwd("/path/to/wrkdir")
#setwd("C:/Users/jmcgirr/Desktop/git_test/")

####################################################
# 2. Set SPECIES_NAMES (prefixes above)
####################################################

species_names <- c("species1","species2")
#species_names <- c("histoplasmaG186ar","Neoformans")

####################################################
# 3. Set CLR percentile threshold
####################################################

clr_thresh <- 0.9


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
  #  candidate_outliers <- sweed %>% 
  #  group_by(seqnames) %>%
  #  mutate(highest_CLR = max(CLR)) %>%
  #  group_by(seqnames) %>%
  #  filter(CLR == highest_CLR) %>%
  #  as.data.frame()
  
  # get top "clr_thresh" percentile for each chromosome
    candidate_outliers <- sweed %>% 
    group_by(seqnames) %>%
    mutate(highest_CLR = quantile(CLR,clr_thresh)) %>%
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
  
  write.table(candidate_genes, paste(species_name,"_candidate_genes.txt",sep = ""), row.names = FALSE, quote= FALSE,sep="\t")
}


#####
#######################################################
# 3. UpsetR plot showing shared genes under selection #
#######################################################

mybiglist <- list()
for (species_name in species_names){

  cans <- read.delim(paste(species_name,"_candidate_genes.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  cans <- cans$Protein.name
  cans <- str_remove(cans, "conserved hypothetical protein")
  cans <- str_remove(cans, "predicted protein")
  cans2 <- unique(cans)
  mybiglist[[species_name]] <- cans
}


#png("upsetr.png", width = 12, height = 7, units = 'in', res = 600)

upset(fromList(mybiglist), sets = species_names,order.by = "freq", point.size = 3.4, line.size = 1.2, 
      mainbar.y.label = "genes intersecting", sets.x.label = "genes", nsets = 6,
      text.scale = c(2, 1.6, 1.5, 1.4, 2,1.6))
#dev.off()

#####
##################################################
# 4. plot CLR for scaffolds with candidate genes #
##################################################

plot_candidate <- function(focal_species_name,protein_product){
  
  sweed <- read.table(paste("SweeD_Report.",focal_species_name,"_grid1000_sweeps.txt",sep = ""), header = FALSE, stringsAsFactors = FALSE)
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
  
  cans <- read.delim(paste(focal_species_name,"_candidate_genes.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  cans <- cans %>% filter(Protein.name == protein_product) 
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  sw <- filter(sweed, seqnames == cans$seqnames[1])
  smooth_sw <-   smooth.spline(sw$start, range01(sw$CLR), spar = .2)
  
  #png(paste(species_name,"_",protein_product,".png",sep = ""), width = 7, height = 4, units = 'in', res = 300)
  plot(sw$start, range01(sw$CLR), xlab = sw$seqnames[1], main = paste(paste(focal_species_name,cans$Protein.product[1],sep = " "),protein_product, sep = "\n"),ylab = "CLR")#,col = "white")
  lines(smooth_sw, lwd = 2, col = "#E8000B")
  sw <- sw[order(sw$CLR, decreasing = TRUE),]
  abline(v=sw$start[1] , col="#1f77b4", lwd=1, lty=2)
  #dev.off()
  
}

plot_candidate("Neoformans","streptomycin biosynthesis protein StrI")
plot_candidate("histoplasmaG186ar","streptomycin biosynthesis protein StrI")

