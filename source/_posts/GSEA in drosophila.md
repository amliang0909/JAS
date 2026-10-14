---
title: GSEA in drosophila
date: 2024-12-31 14:56:23
categories: Bioinformatic
statistics: true
sticky: 9997
---



### GSEA in drosophila

.



This is a work log using R script to perform Gene Set Enrichment Analysis by differentialy expression genes in drosophila.

and this script also applies to other species, such as human and mouse.

human : org.Hs.eg.db

mouse : org.Mm.eg.db

```R

library(org.Dm.eg.db)
library(clusterProfiler)
library(DOSE)
library(topGO)
library(pathview)
library(KEGG.db)
library(enrichplot)


# RNA-Seq data
setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/240506")
library(readxl)
expr_gene <- read_xlsx("DEGs_240506.xlsx",sheet = "new.diff_gene_all_S2",skip = 0)


# transfer ENSEMBL into ENTREZID
geneList<-expr_gene[,c("GeneID","log2FoldChange")]
eg <- bitr(geneList$GeneID, 
           fromType="ENSEMBL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Dm.eg.db")
head(eg)


# sorted by log2FC
input_data <- merge(eg,geneList,by.x="ENSEMBL",by.y="GeneID")
input_data_sort <- input_data[order(input_data$log2FoldChange, decreasing = T),]
gene_fc = input_data_sort$log2FoldChange
names(gene_fc) <- input_data_sort$ENTREZID
head(gene_fc)


# gseGO
GO <- gseGO(
  gene_fc, #gene_fc
  ont = "ALL", # "BP","MF","CC" or "ALL"
  OrgDb = org.Dm.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 1,
  pAdjustMethod = "BH")
head(GO[,1:4])
go_terms <- as.data.frame(GO)


# diaw plot
if (!file.exists("E:/220625_PC/R workplace/220320_SXL/202404_Fig/241231")){
  dir.create("E:/220625_PC/R workplace/220320_SXL/202404_Fig/241231")
  setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/241231")
} else {
  setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/241231")
}

tiff("sex_differentiation.tiff",units = "in",width = 10,height = 7.2,res = 600)
gseaplot2(GO, "GO:0007548", color = "firebrick", 
          title = "sex differentiation",
          rel_heights=c(1, .2, .6))

dev.off()
```





output files :

I chose to show a term about sex determination pathway:

![Terms](./GSEA%20in%20drosophila/Terms.png)





