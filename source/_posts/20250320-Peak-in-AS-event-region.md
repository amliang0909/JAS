---
title: Peak in AS-event region
date: 2025-03-20 16:18:10
categories: Bioinformatic
statistics: true
sticky: 9980
---

### Peak in AS-event region

.



**Algorithm:**

To quantify the number of splicing events **directly** influenced by RNA-binding proteins,

we consider the peak which locating in following region (yellow) as a regulator for alternative splicing.

![region of splicing event with peak](./20250320-Peak-in-AS-event-region/region%20of%20splicing%20event%20with%20peak.png)

.

**Scripts:**

Peak in AS-event-region 

the region of AS-event shows in "Figure peak_in_AS_events 2025-03-18.cvx"

```shell
# genome annotation for geneID transform
setwd("E:/220625_PC/R workplace/20200110_exp_count_WM/")
library(rtracklayer)
genome.anno = import("dmel-all-r6.44.gtf") %>% as.data.frame() 
genome.anno=genome.anno[,c("gene_id","gene_symbol","seqnames","start","end","width","type")] 
head(genome.anno)
genome.anno <- genome.anno[genome.anno$type=="gene",]
colnames(genome.anno) <- c("GeneID","gene_symbol","seqnames","start","end","width","type")

```

.

In light of the output file format from rMATS, the following diagram is instrumental in facilitating the extraction of each individual region within splicing events.



![postion of AS-events by rMATS](./20250320-Peak-in-AS-event-region/postion%20of%20AS-events%20by%20rMATS.png)



.

Extract the region of spliced sites of splcing events.

```shell
# ------ splice site
# SE
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/")
se_sxl_novelSS_S2 <- read.table("SE.MATS.JC.txt",header = T)
se_sxl_novelSS_S2 <- merge(se_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
se_sxl_novelSS_S2$chr <- gsub("chr", "", se_sxl_novelSS_S2$chr)

S2_se_upstreamEE <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE")]
S2_se_exonStart_0base <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","exonStart_0base")]
S2_se_exonEnd <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","exonEnd")]
S2_se_downstreamES <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
write.table(S2_se_upstreamEE,file = "S2_se_upstreamEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_se_exonStart_0base,file = "S2_se_exonStart_0base.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_se_exonEnd,file = "S2_se_exonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_se_downstreamES,file = "S2_se_downstreamES.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# RI
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/new/without_novel/")
ri_sxl_novelSS_S2 <- read.table("RI.MATS.JC.txt",header = T)
ri_sxl_novelSS_S2 <- merge(ri_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
ri_sxl_novelSS_S2$chr <- gsub("chr", "", ri_sxl_novelSS_S2$chr)

S2_ri_upstreamEE <- ri_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE")]
S2_ri_downstreamES <- ri_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
write.table(S2_ri_upstreamEE,file = "S2_ri_upstreamEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_ri_downstreamES,file = "S2_ri_downstreamES.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# A5SS
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/") 
a5ss_sxl_novelSS_S2 <- read.table("A5SS.MATS.JC.txt",header = T)
a5ss_sxl_novelSS_S2 <- merge(a5ss_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
a5ss_sxl_novelSS_S2$chr <- gsub("chr", "", a5ss_sxl_novelSS_S2$chr)
a5ss_sxl_novelSS_S2_p <- a5ss_sxl_novelSS_S2[a5ss_sxl_novelSS_S2$strand == "+",]
a5ss_sxl_novelSS_S2_m <- a5ss_sxl_novelSS_S2[a5ss_sxl_novelSS_S2$strand == "-",]

S2_a5ss_p_shortEE <- a5ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","shortEE")]
S2_a5ss_p_longExonEnd <- a5ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","longExonEnd")]
S2_a5ss_p_flankingES <- a5ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","flankingES")]

S2_a5ss_m_flankingEE <- a5ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","flankingEE")]
S2_a5ss_m_longExonStart_0base <- a5ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","longExonStart_0base")]
S2_a5ss_m_shortES <- a5ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","shortES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
write.table(S2_a5ss_p_shortEE,file = "S2_a5ss_p_shortEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_p_longExonEnd,file = "S2_a5ss_p_longExonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_p_flankingES,file = "S2_a5ss_p_flankingES.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_m_flankingEE,file = "S2_a5ss_m_flankingEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_m_longExonStart_0base,file = "S2_a5ss_m_longExonStart_0base.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_m_shortES,file = "S2_a5ss_m_shortES.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# A3SS
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/") 
a3ss_sxl_novelSS_S2 <- read.table("A3SS.MATS.JC.txt",header = T)
a3ss_sxl_novelSS_S2 <- merge(a3ss_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
a3ss_sxl_novelSS_S2$chr <- gsub("chr", "", a3ss_sxl_novelSS_S2$chr)
a3ss_sxl_novelSS_S2_p <- a3ss_sxl_novelSS_S2[a3ss_sxl_novelSS_S2$strand == "+",]
a3ss_sxl_novelSS_S2_m <- a3ss_sxl_novelSS_S2[a3ss_sxl_novelSS_S2$strand == "-",]

S2_a3ss_p_flankingEE <- a3ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","flankingEE")]
S2_a3ss_p_longExonStart_0base <- a3ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","longExonStart_0base")]
S2_a3ss_p_shortES <- a3ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","shortES")]

S2_a3ss_m_shortEE <- a3ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","shortEE")]
S2_a3ss_m_longExonEnd <- a3ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","longExonEnd")]
S2_a3ss_m_flankingES <- a3ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","flankingES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
write.table(S2_a3ss_p_flankingEE,file = "S2_a3ss_p_flankingEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_p_longExonStart_0base,file = "S2_a3ss_p_longExonStart_0base.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_p_shortES,file = "S2_a3ss_p_shortES.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_m_shortEE,file = "S2_a3ss_m_shortEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_m_longExonEnd,file = "S2_a3ss_m_longExonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_m_flankingES,file = "S2_a3ss_m_flankingES.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# MXE
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/")
mxe_sxl_novelSS_S2 <- read.table("MXE.MATS.JC.txt",header = T)
mxe_sxl_novelSS_S2 <- merge(mxe_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
mxe_sxl_novelSS_S2$chr <- gsub("chr", "", mxe_sxl_novelSS_S2$chr)

S2_mxe_upstreamEE <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE")]
S2_mxe_X1stExonStart_0base <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X1stExonStart_0base")]
S2_mxe_X1stExonEnd <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X1stExonEnd")]
S2_mxe_X2ndExonStart_0base <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X2ndExonStart_0base")]
S2_mxe_X2ndExonEnd <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X2ndExonEnd")]
S2_mxe_downstreamES <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
write.table(S2_mxe_upstreamEE,file = "S2_mxe_upstreamEE.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_X1stExonStart_0base,file = "S2_mxe_X1stExonStart_0base.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_X1stExonEnd,file = "S2_mxe_X1stExonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_X2ndExonStart_0base,file = "S2_mxe_X2ndExonStart_0base.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_X2ndExonEnd,file = "S2_mxe_X2ndExonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_X2ndExonEnd,file = "S2_mxe_X2ndExonEnd.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_downstreamES,file = "S2_mxe_downstreamES.txt",quote = F,row.names = F,col.names = F,sep = '\t')

```

​	.

Extract the region of intron related to splcing events.

```shell
# ------ splicing_related_intron
# SE
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/")
se_sxl_novelSS_S2 <- read.table("SE.MATS.JC.txt",header = T)
se_sxl_novelSS_S2 <- merge(se_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
se_sxl_novelSS_S2$chr <- gsub("chr", "", se_sxl_novelSS_S2$chr)

S2_se_intron1 <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE","exonStart_0base")]
S2_se_intron2 <- se_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","exonEnd","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
write.table(S2_se_intron1,file = "S2_se_intron1.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_se_intron2,file = "S2_se_intron2.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# RI
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/new/without_novel/")
ri_sxl_novelSS_S2 <- read.table("RI.MATS.JC.txt",header = T)
ri_sxl_novelSS_S2 <- merge(ri_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
ri_sxl_novelSS_S2$chr <- gsub("chr", "", ri_sxl_novelSS_S2$chr)

S2_ri_intron <- ri_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
write.table(S2_ri_intron,file = "S2_ri_intron.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# A5SS
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/") 
a5ss_sxl_novelSS_S2 <- read.table("A5SS.MATS.JC.txt",header = T)
a5ss_sxl_novelSS_S2 <- merge(a5ss_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
a5ss_sxl_novelSS_S2$chr <- gsub("chr", "", a5ss_sxl_novelSS_S2$chr)
a5ss_sxl_novelSS_S2_p <- a5ss_sxl_novelSS_S2[a5ss_sxl_novelSS_S2$strand == "+",]
a5ss_sxl_novelSS_S2_m <- a5ss_sxl_novelSS_S2[a5ss_sxl_novelSS_S2$strand == "-",]

S2_a5ss_p_intron <- a5ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","longExonEnd","flankingES")]
S2_a5ss_m_intron <- a5ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","flankingEE","longExonStart_0base")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
write.table(S2_a5ss_p_intron,file = "S2_a5ss_p_intron.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a5ss_m_intron,file = "S2_a5ss_m_intron.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# A3SS
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/") 
a3ss_sxl_novelSS_S2 <- read.table("A3SS.MATS.JC.txt",header = T)
a3ss_sxl_novelSS_S2 <- merge(a3ss_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
a3ss_sxl_novelSS_S2$chr <- gsub("chr", "", a3ss_sxl_novelSS_S2$chr)
a3ss_sxl_novelSS_S2_p <- a3ss_sxl_novelSS_S2[a3ss_sxl_novelSS_S2$strand == "+",]
a3ss_sxl_novelSS_S2_m <- a3ss_sxl_novelSS_S2[a3ss_sxl_novelSS_S2$strand == "-",]

S2_a3ss_p_intron <- a3ss_sxl_novelSS_S2_p[,c("GeneID","ID","gene_symbol","chr","strand","flankingEE","longExonStart_0base")]
S2_a3ss_m_intorn <- a3ss_sxl_novelSS_S2_m[,c("GeneID","ID","gene_symbol","chr","strand","longExonEnd","flankingES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
write.table(S2_a3ss_p_intron,file = "S2_a3ss_p_intron.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_a3ss_m_intorn,file = "S2_a3ss_m_intorn.txt",quote = F,row.names = F,col.names = F,sep = '\t')

# MXE
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/")
mxe_sxl_novelSS_S2 <- read.table("MXE.MATS.JC.txt",header = T)
mxe_sxl_novelSS_S2 <- merge(mxe_sxl_novelSS_S2,genome.anno, by="GeneID", all.x=TRUE)
mxe_sxl_novelSS_S2$chr <- gsub("chr", "", mxe_sxl_novelSS_S2$chr)

S2_mxe_intron1 <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","upstreamEE","X1stExonStart_0base")]
S2_mxe_intron2 <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X1stExonEnd","X2ndExonStart_0base")]
S2_mxe_intron3 <- mxe_sxl_novelSS_S2[,c("GeneID","ID","gene_symbol","chr","strand","X2ndExonEnd","downstreamES")]

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
write.table(S2_mxe_intron1,file = "S2_mxe_intron1.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_intron2,file = "S2_mxe_intron2.txt",quote = F,row.names = F,col.names = F,sep = '\t')
write.table(S2_mxe_intron3,file = "S2_mxe_intron3.txt",quote = F,row.names = F,col.names = F,sep = '\t')
```

​	.

Get the peak which located on the identified region by shell script

```shell
# ------ splice site(-50bp & +50bp) 

# to bed
ls *.txt |while read id; do(less $id | awk '{print $4"\t"$6-50"\t"$6+50"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'  > $(basename ${id} ".txt").bed );done

# get intersected peak
ls *bed | while read id; do(bedtools intersect -a ../../../F.peak.bed -b $id -wa -wb |sort -k4 > peak_in_$(basename ${id} ".bed"));done


# ------ intron
# Convert all irregular delimiters in the file to tab ("\t") separation 
ls *.txt |while read id; do(tr -d '\r' < $id | sed 's/ \+/\t/g'  > $(basename ${id} ".txt").TXT );done

# to bed
ls *.TXT |while read id; do(less $id | awk '{print $4"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$5}'  > $(basename ${id} ".TXT").bed );done

# get intersected peak
ls *bed | while read id; do(bedtools intersect -a ../../F.peak.bed -b $id -wa -wb |sort -k4 > peak_in_$(basename ${id} ".bed"));done
```

.

Peak in spliced-related region

```R
#------------------------------------ with peak 

# the shell script & data in : 
# /data1/amliang/projects/sxl_DC/220906_peak_in_AS.events/20250318_peak_in_S2

# ------ peak in splice site
setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/splice.site/")
peak_in_S2_se_upstreamEE <- read.table("peak_in_S2_se_upstreamEE",header = F)
peak_in_S2_se_exonStart_0base <- read.table("peak_in_S2_se_exonStart_0base",header = F)
peak_in_S2_se_exonEnd <- read.table("peak_in_S2_se_exonEnd",header = F)
peak_in_S2_se_downstreamES <- read.table("peak_in_S2_se_downstreamES",header = F)

peak_in_S2_ri_upstreamEE <- read.table("peak_in_S2_ri_upstreamEE",header = F)
peak_in_S2_ri_downstreamES <- read.table("peak_in_S2_ri_downstreamES",header = F)

peak_in_S2_a5ss_p_shortEE <- read.table("peak_in_S2_a5ss_p_shortEE",header = F)
peak_in_S2_a5ss_p_longExonEnd <- read.table("peak_in_S2_a5ss_p_longExonEnd",header = F)
peak_in_S2_a5ss_p_flankingES <- read.table("peak_in_S2_a5ss_p_flankingES",header = F)
peak_in_S2_a5ss_m_flankingEE <- read.table("peak_in_S2_a5ss_m_flankingEE",header = F)
peak_in_S2_a5ss_m_longExonStart_0base <- read.table("peak_in_S2_a5ss_m_longExonStart_0base",header = F)
peak_in_S2_a5ss_m_shortES <- read.table("peak_in_S2_a5ss_m_shortES",header = F)

peak_in_S2_a3ss_p_flankingEE <- read.table("peak_in_S2_a3ss_p_flankingEE",header = F)
peak_in_S2_a3ss_p_longExonStart_0base <- read.table("peak_in_S2_a3ss_p_longExonStart_0base",header = F)
peak_in_S2_a3ss_p_shortES <- read.table("peak_in_S2_a3ss_p_shortES",header = F)
peak_in_S2_a3ss_m_shortEE <- read.table("peak_in_S2_a3ss_m_shortEE",header = F)
peak_in_S2_a3ss_m_longExonEnd <- read.table("peak_in_S2_a3ss_m_longExonEnd",header = F)
peak_in_S2_a3ss_m_flankingES <- read.table("peak_in_S2_a3ss_m_flankingES",header = F)

peak_in_S2_mxe_upstreamEE<- read.table("peak_in_S2_mxe_upstreamEE",header = F)
peak_in_S2_mxe_X1stExonStart_0base <- read.table("peak_in_S2_mxe_X1stExonStart_0base",header = F)
peak_in_S2_mxe_X1stExonEnd<- read.table("peak_in_S2_mxe_X1stExonEnd",header = F)
peak_in_S2_mxe_X2ndExonStart_0base<- read.table("peak_in_S2_mxe_X2ndExonStart_0base",header = F)
peak_in_S2_mxe_X2ndExonEnd<- read.table("peak_in_S2_mxe_X2ndExonEnd",header = F)
peak_in_S2_mxe_downstreamES<- read.table("peak_in_S2_mxe_downstreamES",header = F)


# 
peak_in_S2_se_SS <- rbind(peak_in_S2_se_upstreamEE,
                       peak_in_S2_se_exonStart_0base,
                       peak_in_S2_se_exonEnd,
                       peak_in_S2_se_downstreamES)

peak_in_S2_ri_SS <- rbind(peak_in_S2_ri_upstreamEE,
                       peak_in_S2_ri_downstreamES)

peak_in_S2_a5ss_SS <- rbind(peak_in_S2_a5ss_p_shortEE,
                       peak_in_S2_a5ss_p_longExonEnd,
                       peak_in_S2_a5ss_p_flankingES,
                       peak_in_S2_a5ss_m_flankingEE,
                       peak_in_S2_a5ss_m_longExonStart_0base,
                       peak_in_S2_a5ss_m_shortES)

peak_in_S2_a3ss_SS <- rbind(peak_in_S2_a3ss_p_flankingEE,
                            peak_in_S2_a3ss_p_longExonStart_0base,
                            peak_in_S2_a3ss_p_shortES,
                            peak_in_S2_a3ss_m_shortEE,
                            peak_in_S2_a3ss_m_longExonEnd,
                            peak_in_S2_a3ss_m_flankingES)

peak_in_S2_mxe_SS <- rbind(peak_in_S2_mxe_upstreamEE,
                            peak_in_S2_mxe_X1stExonStart_0base,
                            peak_in_S2_mxe_X1stExonEnd,
                            peak_in_S2_mxe_X2ndExonStart_0base,
                            peak_in_S2_mxe_X2ndExonEnd,
                            peak_in_S2_mxe_downstreamES)

col_name_ss <- c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                     "chr","target_start","target_end","GeneID","ID","gene_symbol","chr","strand","pos")

colnames(peak_in_S2_se_SS) <- col_name_ss 
colnames(peak_in_S2_ri_SS) <- col_name_ss 
colnames(peak_in_S2_a5ss_SS) <- col_name_ss 
colnames(peak_in_S2_a3ss_SS) <- col_name_ss 
colnames(peak_in_S2_mxe_SS) <- col_name_ss 


# ------ peak in splice site
setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/250318/intron/")
peak_in_S2_se_intron1 <- read.table("peak_in_S2_se_intron1",header = F)
peak_in_S2_se_intron2 <- read.table("peak_in_S2_se_intron2",header = F)

peak_in_S2_ri_intron <- read.table("peak_in_S2_ri_intron",header = F)

peak_in_S2_a5ss_p_intron <- read.table("peak_in_S2_a5ss_p_intron",header = F)
peak_in_S2_a5ss_m_intron <- read.table("peak_in_S2_a5ss_m_intron",header = F)

peak_in_S2_a3ss_p_intron <- read.table("peak_in_S2_a3ss_p_intron",header = F)
peak_in_S2_a3ss_m_intron <- read.table("peak_in_S2_a3ss_m_intorn",header = F)

peak_in_S2_mxe_intron1 <- read.table("peak_in_S2_mxe_intron1",header = F)
peak_in_S2_mxe_intron2 <- read.table("peak_in_S2_mxe_intron2",header = F)
peak_in_S2_mxe_intron3 <- read.table("peak_in_S2_mxe_intron3",header = F)

peak_in_S2_se_intron <- rbind(peak_in_S2_se_intron1,
                              peak_in_S2_se_intron2)
peak_in_S2_a5ss_intron <- rbind(peak_in_S2_a5ss_p_intron,
                                peak_in_S2_a5ss_m_intron)
peak_in_S2_a3ss_intron <- rbind(peak_in_S2_a3ss_p_intron,
                                peak_in_S2_a3ss_m_intron)
peak_in_S2_mxe_intron <- rbind(peak_in_S2_mxe_intron1,
                               peak_in_S2_mxe_intron2,
                               peak_in_S2_mxe_intron3)

col_name_intron <- c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                 "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")

colnames(peak_in_S2_se_intron) <- col_name_intron 
colnames(peak_in_S2_ri_intron) <- col_name_intron 
colnames(peak_in_S2_a5ss_intron) <- col_name_intron 
colnames(peak_in_S2_a3ss_intron) <- col_name_intron 
colnames(peak_in_S2_mxe_intron) <- col_name_intron 


# ------ merge peak in splice site and intron
peak_in_S2_se_SS <- peak_in_S2_se_SS[,c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                        "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")]
peak_in_S2_ri_SS <- peak_in_S2_ri_SS[,c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                        "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")]
peak_in_S2_a5ss_SS <- peak_in_S2_a5ss_SS[,c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                        "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")]
peak_in_S2_a3ss_SS <- peak_in_S2_a3ss_SS[,c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                        "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")]
peak_in_S2_mxe_SS <- peak_in_S2_mxe_SS[,c("peak_chr","peak_start","peak_end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                        "chr","target_start","target_end","GeneID","ID","gene_symbol","strand")]

peak_SE_S2 <- rbind(peak_in_S2_se_SS,peak_in_S2_se_intron)
peak_RI_S2 <- rbind(peak_in_S2_ri_SS,peak_in_S2_ri_intron)
peak_A5SS_S2 <- rbind(peak_in_S2_a5ss_SS,peak_in_S2_a5ss_intron)
peak_A3SS_S2 <- rbind(peak_in_S2_a3ss_SS,peak_in_S2_a3ss_intron)
peak_MXE_S2 <- rbind(peak_in_S2_mxe_SS,peak_in_S2_mxe_intron)

# filter
peak_SE_S2  <- peak_SE_S2[peak_SE_S2$peak.strand == peak_SE_S2$strand,]
peak_RI_S2  <- peak_RI_S2[peak_RI_S2$peak.strand == peak_RI_S2$strand,]
peak_A5SS_S2  <- peak_A5SS_S2[peak_A5SS_S2$peak.strand == peak_A5SS_S2$strand,]
peak_A3SS_S2  <- peak_A3SS_S2[peak_A3SS_S2$peak.strand == peak_A3SS_S2$strand,]
peak_MXE_S2  <- peak_MXE_S2[peak_MXE_S2$peak.strand == peak_MXE_S2$strand,]

```

.

Significant alternative splicing events

```R
#------ AS-events

setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/raw_gtf/novel/") 
S2_se_novelSS <- read.table("SE.MATS.JC.txt",header = T)
S2_a3ss_novelSS <- read.table("A3SS.MATS.JC.txt",header = T)
S2_a5ss_novelSS <- read.table("A5SS.MATS.JC.txt",header = T)
S2_mxe_novelSS <- read.table("MXE.MATS.JC.txt",header = T)
setwd("E:/220625_PC/R workplace/220320_SXL/220421_SXL.S2/AS/230806_AS/new/without_novel/") 
S2_ri_novelSS <- read.table("RI.MATS.JC.txt",header = T)

library(stringr)
S2_se_novelSS$SJC_WT_1 <- str_split(S2_se_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,1]
S2_se_novelSS$SJC_WT_2 <- str_split(S2_se_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,2]
S2_se_novelSS$SJC_WT_3 <- str_split(S2_se_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,3]
S2_se_novelSS$SJC_Sxl_1 <- str_split(S2_se_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,1]
S2_se_novelSS$SJC_Sxl_2 <- str_split(S2_se_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,2]
S2_se_novelSS$SJC_Sxl_3 <- str_split(S2_se_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,3]

S2_ri_novelSS$SJC_WT_1 <- str_split(S2_ri_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,1]
S2_ri_novelSS$SJC_WT_2 <- str_split(S2_ri_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,2]
S2_ri_novelSS$SJC_WT_3 <- str_split(S2_ri_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,3]
S2_ri_novelSS$SJC_Sxl_1 <- str_split(S2_ri_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,1]
S2_ri_novelSS$SJC_Sxl_2 <- str_split(S2_ri_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,2]
S2_ri_novelSS$SJC_Sxl_3 <- str_split(S2_ri_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,3]

S2_a5ss_novelSS$SJC_WT_1 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,1]
S2_a5ss_novelSS$SJC_WT_2 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,2]
S2_a5ss_novelSS$SJC_WT_3 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,3]
S2_a5ss_novelSS$SJC_Sxl_1 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,1]
S2_a5ss_novelSS$SJC_Sxl_2 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,2]
S2_a5ss_novelSS$SJC_Sxl_3 <- str_split(S2_a5ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,3]

S2_a3ss_novelSS$SJC_WT_1 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,1]
S2_a3ss_novelSS$SJC_WT_2 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,2]
S2_a3ss_novelSS$SJC_WT_3 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,3]
S2_a3ss_novelSS$SJC_Sxl_1 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,1]
S2_a3ss_novelSS$SJC_Sxl_2 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,2]
S2_a3ss_novelSS$SJC_Sxl_3 <- str_split(S2_a3ss_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,3]

S2_mxe_novelSS$SJC_WT_1 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,1]
S2_mxe_novelSS$SJC_WT_2 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,2]
S2_mxe_novelSS$SJC_WT_3 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_1,"\\,",simplify = T)[,3]
S2_mxe_novelSS$SJC_Sxl_1 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,1]
S2_mxe_novelSS$SJC_Sxl_2 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,2]
S2_mxe_novelSS$SJC_Sxl_3 <- str_split(S2_mxe_novelSS$SJC_SAMPLE_2,"\\,",simplify = T)[,3]


# gene ID transform
S2_se_novelSS <- merge(S2_se_novelSS,genome.anno, by="GeneID", all.x=TRUE)
S2_ri_novelSS <- merge(S2_ri_novelSS,genome.anno, by="GeneID", all.x=TRUE)
S2_a5ss_novelSS <- merge(S2_a5ss_novelSS,genome.anno, by="GeneID", all.x=TRUE)
S2_a3ss_novelSS <- merge(S2_a3ss_novelSS,genome.anno, by="GeneID", all.x=TRUE)
S2_mxe_novelSS <- merge(S2_mxe_novelSS, genome.anno, by="GeneID", all.x=TRUE)

S2_diff_se_up <- S2_se_novelSS[S2_se_novelSS$FDR < 0.05 & S2_se_novelSS$IncLevelDifference > 0.05,]
S2_diff_se_down <- S2_se_novelSS[S2_se_novelSS$FDR < 0.05 & S2_se_novelSS$IncLevelDifference < -0.05,]
S2_diff_ri_up <- S2_ri_novelSS[S2_ri_novelSS$FDR < 0.05 & S2_ri_novelSS$IncLevelDifference > 0.05,]
S2_diff_ri_down <- S2_ri_novelSS[S2_ri_novelSS$FDR < 0.05 & S2_ri_novelSS$IncLevelDifference < -0.05,]
S2_diff_a3ss_up <- S2_a3ss_novelSS[S2_a3ss_novelSS$FDR < 0.05 & S2_a3ss_novelSS$IncLevelDifference > 0.05,]
S2_diff_a3ss_down <- S2_a3ss_novelSS[S2_a3ss_novelSS$FDR < 0.05 & S2_a3ss_novelSS$IncLevelDifference < -0.05,]
S2_diff_a5ss_up <- S2_a5ss_novelSS[S2_a5ss_novelSS$FDR < 0.05 & S2_a5ss_novelSS$IncLevelDifference > 0.05,]
S2_diff_a5ss_down <- S2_a5ss_novelSS[S2_a5ss_novelSS$FDR < 0.05 & S2_a5ss_novelSS$IncLevelDifference < -0.05,]
S2_diff_mxe_up <- S2_mxe_novelSS[S2_mxe_novelSS$FDR < 0.05 & S2_mxe_novelSS$IncLevelDifference > 0.05,]
S2_diff_mxe_down <- S2_mxe_novelSS[S2_mxe_novelSS$FDR < 0.05 & S2_mxe_novelSS$IncLevelDifference < -0.05,]

S2_diff_se_up_filter <- S2_diff_se_up[apply(data.frame(as.numeric(S2_diff_se_up$SJC_WT_1),
                                                       as.numeric(S2_diff_se_up$SJC_WT_2)),
                                            1,  mean) > 5 |
                                        apply(data.frame(as.numeric(S2_diff_se_up$SJC_Sxl_1),
                                                         as.numeric(S2_diff_se_up$SJC_Sxl_2)),
                                              1, mean) > 5,]
S2_diff_se_down_filter <- S2_diff_se_down[apply(data.frame(as.numeric(S2_diff_se_down$SJC_WT_1),
                                                           as.numeric(S2_diff_se_down$SJC_WT_2)),
                                                1, mean) > 5 |
                                            apply(data.frame(as.numeric(S2_diff_se_down$SJC_Sxl_1),
                                                             as.numeric(S2_diff_se_down$SJC_Sxl_2)),
                                                  1, mean) > 5,]

S2_diff_ri_up_filter <- S2_diff_ri_up[apply(data.frame(as.numeric(S2_diff_ri_up$SJC_WT_1),
                                                       as.numeric(S2_diff_ri_up$SJC_WT_2)),
                                            1,  mean) > 5 |
                                        apply(data.frame(as.numeric(S2_diff_ri_up$SJC_Sxl_1),
                                                         as.numeric(S2_diff_ri_up$SJC_Sxl_2)),
                                              1, mean) > 5,]
S2_diff_ri_down_filter <- S2_diff_ri_down[apply(data.frame(as.numeric(S2_diff_ri_down$SJC_WT_1),
                                                           as.numeric(S2_diff_ri_down$SJC_WT_2)),
                                                1, mean) > 5 |
                                            apply(data.frame(as.numeric(S2_diff_ri_down$SJC_Sxl_1),
                                                             as.numeric(S2_diff_ri_down$SJC_Sxl_2)),
                                                  1, mean) > 5,]

S2_diff_a3ss_up_filter <- S2_diff_a3ss_up[apply(data.frame(as.numeric(S2_diff_a3ss_up$SJC_WT_1),
                                                           as.numeric(S2_diff_a3ss_up$SJC_WT_2)),
                                                1,  mean) > 5 |
                                            apply(data.frame(as.numeric(S2_diff_a3ss_up$SJC_Sxl_1),
                                                             as.numeric(S2_diff_a3ss_up$SJC_Sxl_2)),
                                                  1, mean) > 5,]
S2_diff_a3ss_down_filter <- S2_diff_a3ss_down[apply(data.frame(as.numeric(S2_diff_a3ss_down$SJC_WT_1),
                                                               as.numeric(S2_diff_a3ss_down$SJC_WT_2)),
                                                    1, mean) > 5 |
                                                apply(data.frame(as.numeric(S2_diff_a3ss_down$SJC_Sxl_1),
                                                                 as.numeric(S2_diff_a3ss_down$SJC_Sxl_2)),
                                                      1, mean) > 5,]

S2_diff_a5ss_up_filter <- S2_diff_a5ss_up[apply(data.frame(as.numeric(S2_diff_a5ss_up$SJC_WT_1),
                                                           as.numeric(S2_diff_a5ss_up$SJC_WT_2)),
                                                1,  mean) > 5 |
                                            apply(data.frame(as.numeric(S2_diff_a5ss_up$SJC_Sxl_1),
                                                             as.numeric(S2_diff_a5ss_up$SJC_Sxl_2)),
                                                  1, mean) > 5,]
S2_diff_a5ss_down_filter <- S2_diff_a5ss_down[apply(data.frame(as.numeric(S2_diff_a5ss_down$SJC_WT_1),
                                                               as.numeric(S2_diff_a5ss_down$SJC_WT_2)),
                                                    1, mean) > 5 |
                                                apply(data.frame(as.numeric(S2_diff_a5ss_down$SJC_Sxl_1),
                                                                 as.numeric(S2_diff_a5ss_down$SJC_Sxl_2)),
                                                      1, mean) > 5,]

S2_diff_mxe_up_filter <- S2_diff_mxe_up[apply(data.frame(as.numeric(S2_diff_mxe_up$SJC_WT_1),
                                                         as.numeric(S2_diff_mxe_up$SJC_WT_2)),
                                              1,  mean) > 5 |
                                          apply(data.frame(as.numeric(S2_diff_mxe_up$SJC_Sxl_1),
                                                           as.numeric(S2_diff_mxe_up$SJC_Sxl_2)),
                                                1, mean) > 5,]
S2_diff_mxe_down_filter <- S2_diff_mxe_down[apply(data.frame(as.numeric(S2_diff_mxe_down$SJC_WT_1),
                                                             as.numeric(S2_diff_mxe_down$SJC_WT_2)),
                                                  1, mean) > 5 |
                                              apply(data.frame(as.numeric(S2_diff_mxe_down$SJC_Sxl_1),
                                                               as.numeric(S2_diff_mxe_down$SJC_Sxl_2)),
                                                    1, mean) > 5,]
```

.

Peak in AS-events-region

```R
# ------ peak in AS-events-region
S2_diff_se_up_peak <- S2_diff_se_up_filter[S2_diff_se_up_filter$ID %in% peak_SE_S2$ID,]
S2_diff_se_down_peak <- S2_diff_se_down_filter[S2_diff_se_down_filter$ID %in% peak_SE_S2$ID,]

S2_diff_ri_up_peak <- S2_diff_ri_up_filter[S2_diff_ri_up_filter$ID %in% peak_RI_S2$ID,]
S2_diff_ri_down_peak <- S2_diff_ri_down_filter[S2_diff_ri_down_filter$ID %in% peak_RI_S2$ID,]

S2_diff_a5ss_up_peak <- S2_diff_a5ss_up_filter[S2_diff_a5ss_up_filter$ID %in% peak_A5SS_S2$ID,]
S2_diff_a5ss_down_peak <- S2_diff_a5ss_down_filter[S2_diff_a5ss_down_filter$ID %in% peak_A5SS_S2$ID,]

S2_diff_a3ss_up_peak <- S2_diff_a3ss_up_filter[S2_diff_a3ss_up_filter$ID %in% peak_A3SS_S2$ID,]
S2_diff_a3ss_down_peak <- S2_diff_a3ss_down_filter[S2_diff_a3ss_down_filter$ID %in% peak_A3SS_S2$ID,]

S2_diff_mxe_up_peak <- S2_diff_mxe_up_filter[S2_diff_mxe_up_filter$ID %in% peak_MXE_S2$ID,]
S2_diff_mxe_down_peak <- S2_diff_mxe_down_filter[S2_diff_mxe_down_filter$ID %in% peak_MXE_S2$ID,]

```

.

Peak in gene with AS-events

```R
# ------ peak in AS-events-genes

# peak
setwd("E:/220625_PC/R workplace/220320_SXL/Fan.data/sxl_peak_bind_DEGs/20220328/")
peak_intersect_feature <- read.table("Sxl-bind-feature.txt",header = F)
head(peak_intersect_feature)
colnames(peak_intersect_feature) <- c("peak.chr","peak.strat","peak.end","peak.xuhao","log10Pv","peak.strand","siteScore","V8",
                                      "chr","start","end","V12","ID_feature","strand")
peak_feature <- peak_intersect_feature
library(stringr)
peak_feature$GeneID <- str_split(peak_feature$ID_feature,"\\_",simplify = T)[,1]
peak_feature$transID <- str_split(peak_feature$ID_feature,"\\_",simplify = T)[,2]
peak_feature$feature <- str_split(peak_feature$ID_feature,"\\_",simplify = T)[,3]

peak_feature <- peak_feature[peak_feature$strand == peak_feature$peak.strand,]

peak_feature$pos <- paste0(peak_feature$peak.chr,":",peak_feature$peak.strat)
peak_feature$pos <- paste0(peak_feature$pos,"-",peak_feature$peak.end)
peak_gene <- peak_feature[!duplicated(peak_feature$GeneID),]


# 
S2_diff_se_up_peak.gene <- S2_diff_se_up_filter[S2_diff_se_up_filter$GeneID %in% peak_gene$GeneID,]
S2_diff_se_down_peak.gene <- S2_diff_se_down_filter[S2_diff_se_down_filter$GeneID %in% peak_gene$GeneID,]

S2_diff_ri_up_peak.gene <- S2_diff_ri_up_filter[S2_diff_ri_up_filter$GeneID %in% peak_gene$GeneID,]
S2_diff_ri_down_peak.gene <- S2_diff_ri_down_filter[S2_diff_ri_down_filter$GeneID %in% peak_gene$GeneID,]

S2_diff_a3ss_up_peak.gene <- S2_diff_a3ss_up_filter[S2_diff_a3ss_up_filter$GeneID %in% peak_gene$GeneID,]
S2_diff_a3ss_down_peak.gene <- S2_diff_a3ss_down_filter[S2_diff_a3ss_down_filter$GeneID %in% peak_gene$GeneID,]

S2_diff_a5ss_up_peak.gene <- S2_diff_a5ss_up_filter[S2_diff_a5ss_up_filter$GeneID %in% peak_gene$GeneID,]
S2_diff_a5ss_down_peak.gene <- S2_diff_a5ss_down_filter[S2_diff_a5ss_down_filter$GeneID %in% peak_gene$GeneID,]

S2_diff_mxe_up_peak.gene <- S2_diff_mxe_up_filter[S2_diff_mxe_up_filter$GeneID %in% peak_gene$GeneID,]
S2_diff_mxe_down_peak.gene <- S2_diff_mxe_down_filter[S2_diff_mxe_down_filter$GeneID %in% peak_gene$GeneID,]
```

.

Count the number of splicing events associated with peak binding.

```R
# ------ STATICS

# get all var in environment
peak_in_ss_names <- ls(pattern = "^S2_diff_.*_peak$")  # Get all variables with named "S2_diff_xx_xx_peak"
peak_in_gene_names <- ls(pattern = "^S2_diff_.*_peak\\.gene$") #  "S2_diff_xx_xx_peak.gene"
diff_AS_names <- ls(pattern = "^S2_diff_.*_filter$")  # "d.mut_f.wt_XX_XX"

# sort and match names
peak_in_ss_names <- sort(peak_in_ss_names)
peak_in_gene_names <- sort(peak_in_gene_names)
diff_AS_names <- sort(diff_AS_names)

# statics the nrow
row_counts_df <- data.frame(
  Dataset = peak_var_names, 
  peak_in_ss_Counts = sapply(mget(peak_in_ss_names, ifnotfound = list(data.frame())), nrow),
  peak_in_gene_Counts = sapply(mget(peak_in_gene_names, ifnotfound = list(data.frame())), nrow),
  diff_AS_Counts = sapply(mget(diff_AS_names, ifnotfound = list(data.frame())), nrow)
)
```

.
