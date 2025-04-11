---
title: Add full name and description for target genes
date: 2025-01-10 09:17:36
categories: Bioinformatic
statistics: true
sticky: 9994
---





### 	Add full name and description for target genes 

.

Some genes' functions may not be easily inferred from their abbreviated names, therefore, to be better understand the genes, we need to include their full names and other forms of ID numbers for the target genes.

```R
# ------
# function and detail description of target genes by gtf

# anno.
library(org.Dm.eg.db)
eg2Symbol=toTable(org.Dm.egSYMBOL)
eg2name=toTable(org.Dm.egGENENAME)
eg2ENSEMBL=toTable(org.Dm.egENSEMBL)
anno=merge(eg2ENSEMBL,merge(eg2Symbol,eg2name,by='gene_id'),by='gene_id')  
colnames(anno) <- c("ENTREID","GeneID","Symbol","Gene_name")
anno <- anno[!duplicated(anno$GeneID),]

# target genes
library(readxl)
setwd("E://220625_PC/R workplace/220320_SXL/240115_tranlation/RiboSeq/featureCount_CDS_multi/")
target_genes <- read_xlsx("riborex_translation.xlsx",sheet = "tranlation_S2",skip = 0)


# add detail description to target genes
new_target_genes <- merge(target_genes,anno, by="GeneID", all.x=TRUE)
length(intersect(anno$GeneID,target_genes$GeneID))

# save
setwd("E://220625_PC/R workplace/220320_SXL/240222_question.daily/17.detail.description_diff_translation.genes/")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "new_target_genes")
writeData(wb, "new_target_genes", new_target_genes,startCol = 1, startRow = 1, rowNames = T)
saveWorkbook(wb, file='target_genes_with_detial_description.xlsx')

```





output file :

the full names are indicated by the red box.

![target genes with their full names already added](./Add-Full-Name-and-Description-for-Target-Gene/target%20genes%20with%20their%20full%20names%20already%20added.png)

