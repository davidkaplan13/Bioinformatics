setwd('/Users/davidkaplan/Desktop/Ewing/Ahood')

library(dplyr)
library(tidyr)
library(tibble)
library(cluster)
library(NbClust)
library(factoextra)
library(dendextend)

mRNA <- read.csv('CCLE_expression.csv')
sample_info <- read.csv('sample_info.csv')

target_diseases <- filter(
  .data=sample_info,
  primary_disease == 'Colon/Colorectal Cancer'
)

target_diseases_select <- select(
  .data=target_diseases,
  DepMap_ID,
  stripped_cell_line_name
)

transformed_mRNA <- mRNA_colon[,names(mRNA_colon)]
rownames(transformed_mRNA) <- transformed_mRNA$X
transformed_mRNA <- transformed_mRNA[,-1]

transformed_mRNA <- t(transformed_mRNA)
distancemRNA <- dist(transformed_mRNA,method='maximum')
clustermRNA <- hclust(distancemRNA, method='ward.D')

clustermRNAs = cutree(clustermRNA,k=70)
clustermRNAs <- as.data.frame(clustermRNAs)

clustermRNAs <- cbind(gene = rownames(clustermRNAs), clustermRNAs)
clustermRNAs <- clustermRNAs[,-1]
clustermRNAs$gene <- sub("\\..*","",clustermRNAs$gene)

AJUBA_cluster_mRNA <- filter(
  .data=clustermRNAs,
  clustermRNAs == 1
)

USP7_cluster_mRNA <- filter(
  .data=clustermRNAs,
  clustermRNAs == 17
)

write.table(AJUBA_cluster_mRNA$gene,"ajuba_cluster_mRNA_70s_max_ward.txt",sep="\t",row.names = FALSE,quote=FALSE)
write.table(USP7_cluster_mRNA$gene,"USP7_cluster_mRNA_70s_max_ward.txt",sep='\t',row.names = FALSE,quote=FALSE)

USP7_mRNA <- as.data.frame(cor(mRNA_colon[,2:ncol(mRNA_colon)],mRNA_colon[['USP7..7874.']],use='pairwise.complete.obs'))
USP7_mRNA <- rownames_to_column(USP7_mRNA, var = "gene")
USP7_mRNA$gene <- sub("\\..*","",USP7_mRNA$gene)
USP7_mRNA <- USP7_mRNA %>% arrange(-abs(V1))

AJUBA_mRNA <- as.data.frame(cor(mRNA_colon[,2:ncol(mRNA_colon)],mRNA_colon[['AJUBA..84962.']],use='pairwise.complete.obs'))
AJUBA_mRNA <- rownames_to_column(AJUBA_mRNA, var = "gene")
AJUBA_mRNA$gene <- sub("\\..*","",AJUBA_mRNA$gene)
AJUBA_mRNA <- AJUBA_mRNA %>% arrange(-abs(V1))

both500_mRNA <- as.data.frame(intersect(AJUBA_mRNA$gene[1:500],USP7_mRNA$gene[1:500]),var='gene')
write.table(both500_mRNA,"both500_mRNA.txt",sep='\t',row.names=FALSE,quote=FALSE)

