---
title: Proportional-area Venn diagram
date: 2025-02-19 21:05:38
categories: R
statistics: true
sticky: 9983
---

### drawing a proportional area Venn diagram in R

.

Sometimes, when creating a Venn diagram, we need the size of each area to be proportional to the values it represents. In this case, we require a proportional-area Venn diagram.

.

the input data need to be a list

```R
setwd("E:/220625_PC/R workplace/250219_veen/")
library(readxl)
peak_3UTR <- read_xlsx("input_data.xlsx",sheet = "3UTR",skip = 0)
peak_5UTR <- read_xlsx("input_data.xlsx",sheet = "5UTR",skip = 0)
peak_Intron <- read_xlsx("input_data.xlsx",sheet = "Intron",skip = 0)


data_list <- list("3UTR" = peak_3UTR$geneID,
                  "5UTR" = peak_5UTR$geneID,
                  "Intron" = peak_Intron$geneID)
data_list <- lapply(data_list, function(x) x[x != ""]) 
# Remove empty strings ("") from each vector in the list
data_list  <- lapply(data_list , unique)  
# Remove duplicate elements from each vector 
```



.

plot

```R
library(eulerr)
plot(euler(
  data_list,
  shape = "circle"), 
  # Shape of the pattern: ellipse or circle
  quantities = list(type = c("percent","counts"),cex=1),
  # Display type: percentage and numbers, with the ability to adjust the size of the numbers.
  labels=list(cex=1),                   
  # Size of the group labels
  edges = list(col = "black", lex = 2), 
  # Edge color and size of the shapes
  fills = list(fill = c("#f18c8d","#8ec7ff","#bfff7f"),alpha=0.7), 
  # Fill color and transparency (opacity)
  legend = list(side = "right")       
  # Position of the legend
)

```

.



![Veen plot](./20250219-proportional-area-Venn-diagram/Veen%20plot.png)





