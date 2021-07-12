# Comparing AJUBA to all cell-lineages 
# and writing to a CSV file.

# Setting Work directory
setwd('/Users/davidkaplan/Desktop/Ewing/Ahood')

# Install packages
install.packages(c('tibble','tidyr'))

#Load packages
library(dplyr)
library(tidyr)
library(tibble)

# Load CSV file into dataframe
CRISPR_gene_effect <- read.csv('/Users/davidkaplan/Desktop/Ewing/Ahood/CRISPR_gene_effect.csv')


# Correlate AJUBA to all other genes in data-set
Cell_cor <- as.data.frame(cor(CRISPR_gene_effect[,2:ncol(CRISPR_gene_effect)],CRISPR_gene_effect[['AJUBA..84962.']],use='pairwise.complete.obs'))
# Rename the column and tranpose
Cell_cor <- rownames_to_column(Cell_cor, var = "gene")

# Order the correlation values in descending order
Cell_cor <- Cell_cor %>% arrange(-abs(V1))
# Remove the first value (Ajuba to Ajuba correlation)
Cell_cor <- Cell_cor[-1,]

# write CSV
write.csv(Cell_cor, file='Ajuba_to_all_cell_lineages_orderr.csv', row.names = FALSE)