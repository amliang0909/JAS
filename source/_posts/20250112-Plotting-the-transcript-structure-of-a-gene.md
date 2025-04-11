---
title: Plotting the transcript structure of a gene
date: 2025-01-12 10:26:09
categories: R
statistics: true
sticky: 9992
---





### 	 Plotting the transcript structure of a gene

.

Sometimes, we need to display the different transcript structures of a gene.  And this R script is designed for quickly accomplishing this task.

First, we need an annotation file that contains information about the gene transcripts.

```R
library(rtracklayer)
library(transPlotR)

setwd("E:/220625_PC/R workplace/220910_annotation/Drosaphila/")
gene.anno <- import("Drosophila_melanogaster.BDGP6.32.105.gtf" ,format = "gtf") %>% data.frame()
colnames(gene.anno)


trancriptVis(gtfFile = gene.anno,
             gene = c('Rala',"ssx"), # Target gene to plot
             facetByGene = T,
             addNormalArrow = F,
             newStyleArrow = T,
             base_size =15, # Theme basesize, default(14).
             textLabel="transcript_name", # The text label aesthetic mappings, default('transcript_id')
             textLabelSize=5, # The text label size
             textLabelColor = "black", # The text label color
             intronSize=2, # Intron line size
             exonWidth=0.2,     # exon width 
             relTextDist=0.2, # marked_name relative to exon
             reverse.y=TRUE,
             xAxis.info=TRUE) # Whether retain X axis ticks and text



```

output：

![all transcripts for genes](./20250112-Plotting-the-transcript-structure-of-a-gene/all%20transcripts%20for%20genes.png)

.

.

If we only want to display the specific transcripts of a gene that we are interested in, we can use this script as follows:

```R
# select and plot only the target transcript by myTranscript ("transcript"_id)
trancriptVis(gtfFile = gene.anno,
             gene = c('ssx'), 
             myTranscript = c('FBtr0070215','FBtr0339773'),
             addNormalArrow = F,
             newStyleArrow = T,
             base_size =15, 
             textLabel="transcript_name", 
             textLabelSize=5, 
             textLabelColor = "black", 
             intronSize=2, 
             exonWidth=0.2,     
             relTextDist=0.2, 
             reverse.y=TRUE,
             xAxis.info=TRUE) 
```

output：

![only target transcripts](./20250112-Plotting-the-transcript-structure-of-a-gene/only%20target%20transcripts.png)

-------

There are better ways to solve this problem.

 link:

https://www.jingege.wang/2022/06/04/ggtranscript/

