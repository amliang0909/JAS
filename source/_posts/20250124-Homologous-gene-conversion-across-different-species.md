---
title: Homologous gene conversion across different species
date: 2025-01-24 00:47:46
categories: Bioinformatic
statistics: true
sticky: 9990
---

### Homologous gene conversion across different species

.

 Here, we use three mouse genes as an example. The method is the same for more genes.

```R

library(homologene)
genelist<-c("Acadm","Eno2","Acadvl")
homologene(genelist, inTax = 10090, outTax = 9606)
# Use homologene function for gene conversion
# genelist: the gene list to be converted
# inTax: the taxonomic ID of the species for the input gene list, 10090 for mouse
# outTax: the taxonomic ID of the species to convert to, 9606 for human
```

.

result:

![convert mouse to human](./20250124-Homologous-gene-conversion-across-different-species/convert%20mouse%20to%20human.png)

.

The species IDs supported by Homologene 

```R
homologene::taxData
```

.

For commonly used organisms like mice and humans, there are even dedicated functions available

```R
mouse2human(genelist)

# and human to mouse
human2mouse(c("ACADM","ENO2","ACADVL"))
human2mouse(c("H1.2","TP53"))

# and dme to homo
genelist<-c("Sxl","msl-2")
homologene::taxData
homologene(genelist, inTax = 7227, outTax = 9606)
```



________

ref : https://cloud.tencent.com/developer/article/2115727

