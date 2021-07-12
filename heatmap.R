setwd('/Users/davidkaplan/Desktop/Ewing/Ahood')

install.packages(c('tibble','tidyr','dplyr'))
install.packages('gsubfn')

#Load packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(plotly)
library(cluster)
library(NbClust)
library(factoextra)
library(dendextend)
library(gsubfn)

CRISPR_gene_effect <- read.csv('CRISPR_gene_effect.csv')
AJUBA <- read.csv('Ajuba_to_colorectal.csv')

#Filter top 50 gene names (Gets rid of values)
AJUBA_top <- AJUBA[1:50,-2]
AJUBA_top <- as.data.frame(AJUBA_top)

#Filter only for colorectal cancer cell lines
sample_info <- read.csv('sample_info.csv')
target_disease <- filter(
  .data=sample_info,
  primary_disease == 'Colon/Colorectal Cancer'
)

DepMapID <- CRISPR_gene_effect %>%
  filter(DepMap_ID %in% target_disease$DepMap_ID)

#Filter for only the top 50 genes associated with AJUBA
Genes <- DepMapID[,which((names(DepMapID) %in% AJUBA_top$AJUBA_top)==TRUE)]

#Add Cell lines as index
rownames(Genes) <- DepMapID$DepMap_ID

#Transpose data
gene_t <- t(Genes)
gene_t[,0]
data_melt <- melt(gene_t)

# Plot heatmap
my_matrix <- as.matrix(gene_t)
# coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
Heatmap(my_matrix)

# Gene clustering

# Finding optimal clustering number
fviz_nbclust(gene_t,kmeans,method='wss')
fviz_nbclust(gene_t,kmeans,method='silhouette')

# Distance calculation
# Maximum - Euclidean.
# ward.D
distance <- dist(gene_t,method='euclidean')
clusterGenes <- hclust(distance, method='ward.D2')
plot(clusterGenes)

#Extract top cell lines. 


avg_dend_obj <- as.dendrogram(clusterGenes)
avg_col_dend <- color_branches(avg_dend_obj,h=4)
plot(avg_col_dend)
#
clustergroup = cutree(clusterGenes,k=4)
clustergroup <- as.data.frame(clustergroup)
cluster1 <- filter(
  .data=clustergroup,
  clustergroup == 1
)


cluster2 <- filter(
  .data=clustergroup,
  clustergroup == 2
)

cluster3 <- filter(
  .data=clustergroup,
  clustergroup == 3
)

cluster4 <- filter(
  .data=clustergroup,
  clustergroup == 4
)
gene_t <- as.data.frame(gene_t)
gene_t_cl <- mutate(gene_t,cluster=clustergroup)
count(gene_t_cl,cluster)


ggp <- ggplot(data_melt, aes(X2,X1))+
  geom_tile(aes(fill=value))+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=40,size=6),
        axis.text.y = element_text(size = 8))+
  

ggp







