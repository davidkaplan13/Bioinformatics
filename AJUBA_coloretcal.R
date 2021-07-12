# Correlating AJUBA gene to all genes in colorectal cancer

setwd('/Users/davidkaplan/Desktop/Ewing/Ahood')

#install.packages(c('tibble','tidyr','dplyr'))

#Load packages
library(dplyr)
library(tidyr)
library(tibble)


sample_info <- read.csv('sample_info.csv')
CRISPR_gene_effect <- read.csv('/Users/davidkaplan/Desktop/Ewing/Ahood/CRISPR_gene_effect.csv')

#Filter by disease name
target_disease <- filter(
  .data=sample_info,
  primary_disease == 'Colon/Colorectal Cancer'
)

DepMapID <- CRISPR_gene_effect %>%
  filter(DepMap_ID %in% target_disease$DepMap_ID)


AJUBA <- read.csv('Ajuba_to_all_cell_lineages.csv')
Ajuba_colorectal <- as.data.frame(cor(DepMapID[,2:ncol(DepMapID)],DepMapID[['AJUBA..84962.']],use='pairwise.complete.obs'))
# Rename the column and tranpose
Ajuba_colorectal <- rownames_to_column(Ajuba_colorectal, var='gene')

# Order the correlation values in descending order
Ajuba_colorectal <- Ajuba_colorectal %>% arrange(-abs(V1))
# Remove the first value (Ajuba to Ajuba correlation)
Ajuba_colorectal <- Ajuba_colorectal[-1,]




write.csv(Ajuba_colorectal, file='Ajuba_to_colorectal.csv', row.names = FALSE)


