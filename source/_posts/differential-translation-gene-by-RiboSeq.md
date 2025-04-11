---
title: Differential translation gene by RiboSeq
date: 2025-01-02 21:38:15
categories: Bioinformatic
statistics: true
sticky: 9996
---



### Ribo-Seq analysis for mRNA translation

.



we use Ribo-Seq in combination with RNA-Seq data to investigate the sig. differential translation genes in transcriptome by condition.



#### 1.  Data prepare 

###### 	1.1  verify the integrity of raw sequence data

```shell
ls *fastq.gz | xargs  md5sum > md5.lam.txt
```



###### 	1.2  quality control

```shell
mkdir qc_raw && cd qc_raw
source activate rnaseq
ls ../*.gz | xargs fastqc -o ./
conda deactivate

```

##### 	

###### 	1.3  adapter trimming

```shell
mkdir clean.SE && cd clean.SE
source activate rnaseq
dir='/data1/amliang/projects/sxl_DC/2308_RiboSeq/231221_Ribo_S2/rawdata/clean.SE'
trim_galore -q 25 --phred33 -a AAAAAAAAAA  --clip_R1 4  --length 20  -e 0.1 --stringency 3  -o $dir   ../*.R1.fastq.gz 

# && and QC 
mkdir qc_clean.SE && cd qc_clean.SE
ls ../clean.SE/*.gz | xargs fastqc -o ./
conda deactivate

```



#### 2. remove rRNA

###### 2.1 download rRNA fasta

```shell
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz


# rRNA index building
less GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz | grep "^>" | grep "gbkey=rRNA" | awk '{print $1}'|sed 's/>//g' > id.rRNA.list
source activate py2
seqkit grep -f id.rRNA.list  GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz  > rRNA.fa

mkdir rRNA && cd rRNA 
bowtie2-build ../rRNA.fa  rRNA

# tRNA index building
less GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz | grep "^>" | grep "gbkey=tRNA" | awk '{print $1}'|sed 's/>//g' > id.tRNA.list
source activate py2
seqkit grep -f id.tRNA.list  GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz  > tRNA.fa

mkdir tRNA && cd tRNA 
bowtie2-build ../tRNA.fa  rRNA

```



###### 2.2 remove rRNA and tRNA

```shell
# remove rRNA by mapping reads to rRNA
# save reads which unmapped as non_rRNA reads
mkdir rm.rRNA && cd rm.rRNA
ln -s ../clean.SE/*gz ./

source activate rnaseq
ls ./*.gz | while read id; do (echo $id && bowtie2  -p 18  -x /data1/amliang/reference/rRNA/Drosophila.M/rRNA/rRNA --un-gz $(basename ${id} ".fq.gz").rm.rRNA.fq.gz  -U $id  -S $(basename ${id} ".fq.gz").rRNA.sam);done


# remove rRNA by mapping reads to rRNA
# save reads which unmapped as non_rRNA reads
mkdir rm.tRNA && cd rm.tRNA
ln -s ../rm.tRNA/*.rm.rRNA.fq.gz ./

ls ./*.rm.rRNA.fq.gz | while read id; do (echo $id && bowtie2  -p 18  -x /data1/amliang/reference/rRNA/Drosophila.M/tRNA/tRNA --un-gz $(basename ${id} ".rm.rRNA.fq.gz").rm.rRNA.tRNA.fq.gz  -U $id  -S $(basename ${id} ".rm.rRNA.fq.gz").rRNA.tRNA.sam);done
conda deactivate

```



#### 3. mapping to genome & read counting

###### 3.1 mapping to genome

```shell
# mapping
ln -s ../rm.tRNA/*rm.rRNA.tRNA.fq.gz
Ref_index=/data1/amliang/reference/index/star/flybase_m6/
Read_dir=./
Run_log=runlog.txt
GTF=/data1/amliang/annotation/fly/flyase/gtf/dmel-all-r6.31.gtf
#input/output file 
for id in {P1_L2_Q0023W0073,P2_L3_Q0013W0073,P3_L3_Q0028W0073,S4_L3_Q0016W0073,S5_L3_Q0014W0073,S6_L3_Q0021W0073};
do
Read=${Read_dir}/${id}.R1.rm.rRNA.tRNA.fq.gz
BAM=${id}.bam

#Cycle Running
STAR \
--runThreadN 10 \
--outFilterType BySJout \
--outFilterMismatchNmax 2 \
--outFilterMultimapNmax 1 \
--genomeDir $Ref_index \
--readFilesIn $Read \
--readFilesCommand zcat \
--outFileNamePrefix ${id}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All \
--outSAMattrRGline ID:1 LB:ribo_seq PL:ILLUMINA SM:${id} \
--outBAMcompression 6 \
--outReadsUnmapped Fastx
done
conda deactivate

# sort & index
source activate chipseq
ls *.out.bam|while read id ;do (samtools sort -O bam  -T sorted -o $(basename ${id} ".out.bam").sort.bam   ${id});done

ls *.sort.bam | while read id ; do (samtools index $id); done
conda deactivate
```



###### 3.2 CDS reads count

```shell
source activate rnaseq
featureCounts -T 16  -t CDS -g gene_id -a /data1/amliang/annotation/fly/Drosophila_melanogaster.BDGP6.22.42.gtf  -o  Ribo.CDS_count.a.txt  ./*.bam
conda deactivate

```



#### 4. corresponding RNA-Seq data

###### 4.1 mapping to genome

```shell
source activate rnaseq
Ref_index=/data1/amliang/reference/index/bowtie2/flybase_m6_bowtie2/
Read_dir=/data1/amliang/projects/sxl_DC/2308_RiboSeq/240108_RNASeq_S2/data
Run_log=runlog.txt

#input/output file 
for id in {P1_L2_UDI296,P2_L2_UDI297,P3_L2_UDI298,S4_L2_UDI299,S5_L2_UDI300,S6_L2_UDI301};
do
R1=${Read_dir}/${id}.R1.fastq.gz
R2=${Read_dir}/${id}.R2.fastq.gz

#Cycle Running
bowtie2  -p 18 -N 0 -x $Ref_index/fly -1 $R1 -2 $R2 -S ${id}.sam 2>>$Run_log 
 
done
conda deactivate
```



###### 4.1 mapping to genome

```shell
# Note: reads_count for CDS not exon

source activate rnaseq
featureCounts -T 16 -p -M -t CDS -g gene_id  -a /data1/amliang/annotation/fly/Drosophila_melanogaster.BDGP6.22.42.gtf  -o  RNASeq.CDS_count.txt  ../*.bam
conda deactivate
```





#### 5. differential translation

###### 5.1 expr. count

```R
#Ribo
setwd("E://220625_PC/R workplace/220320_SXL/240115_tranlation/RiboSeq/featureCount_CDS_multi/")
Ribo_count <- read.table("Ribo.CDS.multi_count.a.txt",header = T)
rownames(Ribo_count) <- Ribo_count$Geneid
colnames(Ribo_count)
colnames(Ribo_count) <- c( "Geneid","Chr","Start","End","Strand","Length",                          "WT_1","WT_2","WT_3","SXL_1","SXL_2","SXL_3")
Ribo_exp <- Ribo_count[,c(7:12)]
Ribo_exp <- Ribo_exp[rowSums(Ribo_exp)>0,]

# RNASeq
setwd("E://220625_PC/R workplace/220320_SXL/240115_tranlation/transla_by_Cyt.RNAseq/")
RNA_count <- read.table("RNASeq.CDS_count.txt",header = T)
rownames(RNA_count) <- RNA_count$Geneid
colnames(RNA_count)
colnames(RNA_count) <-c( "Geneid","Chr","Start","End","Strand","Length","wt_1","wt_2","wt_3","Sxl_1","Sxl_2","Sxl_3")
RNA_exp <- RNA_count[,c("C.wt_1","C.wt_2","C.wt_3","C.Sxl_1","C.Sxl_2","C.Sxl_3")]
colnames(RNA_exp) <- c("wt_1","wt_2","wt_3","Sxl_1","Sxl_2","Sxl_3")
RNA_exp <- RNA_exp[rownames(RNA_exp) %in% rownames(Ribo_exp),]

```



###### 5.2 PCA

```R
#normalized - Ribo
library(DESeq2)
Ribo_exp_read.count <- Ribo_count[,7:12]
Ribo_exp_read.count <- Ribo_exp_read.count[rowSums(Ribo_exp_read.count)>0,]
condition <-factor(c(rep("WT",3), rep("SXL",3)),
                   levels = c("WT","SXL"))
colData <- data.frame(row.names=colnames(Ribo_exp_read.count),condition)
dds <- DESeqDataSetFromMatrix(Ribo_exp_read.count, colData, design= ~ condition)
dds <- DESeq(dds)
Ribo_normalized_counts <- counts(dds, normalized=TRUE)
Ribo_normalized_counts <- as.data.frame(Ribo_normalized_counts)
Ribo_normalized_counts$GeneID <- rownames(Ribo_normalized_counts)
Ribo_rpkm<- Ribo_normalized_counts[,-7]

library(FactoMineR)
gene <- t(Ribo_rpkm)
gene.pca <- PCA(gene, ncp = 4, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)
# chose: WT_1,WT_3,SXL_1,SXL2


#normalized - RNAS
library(DESeq2)
RNA_exp_read.count <- RNA_count[,c(10:12,7:9)]
RNA_exp_read.count <- RNA_exp_read.count[rowSums(RNA_exp_read.count)>0,]
condition <-factor(c(rep("WT",3), rep("SXL",3)),
                   levels = c("WT","SXL"))
colData <- data.frame(row.names=colnames(RNA_exp_read.count),condition)
dds <- DESeqDataSetFromMatrix(RNA_exp_read.count, colData, design= ~ condition)
dds <- DESeq(dds)
RNA_normalized_counts <- counts(dds, normalized=TRUE)
RNA_normalized_counts <- as.data.frame(RNA_normalized_counts)
RNA_normalized_counts$GeneID <- rownames(RNA_normalized_counts)
RNA_rpkm<- RNA_normalized_counts[,-7]

library(FactoMineR)
gene <- t(RNA_rpkm)
gene.pca <- PCA(gene, ncp = 4, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)
# chose: WT_1,WT_2,SXL_2,SXL3
```



###### 5.3 calculating differential translation genes by riborex

```R
Ribo_exp <- na.omit(Ribo_exp)
RNA_exp <- na.omit(RNA_exp)

# riborex
library(fdrtool)
library(riborex)
library(DESeq2)

RNACntTable <- RNA_exp[,c(1,2,5,6)]
RiboCntTable <- Ribo_exp[,c(1,3,4,5)]
head(RNACntTable,3)
head(RiboCntTable,3)

rnaCond <- c("control", "control",  "treated", "treated")
riboCond <- c("control", "control",  "treated", "treated")
res.deseq2 <- riborex(RNACntTable, RiboCntTable, rnaCond, riboCond, engine = "DESeq2")
head(res.deseq2,3)
hist(res.deseq2$pvalue, main = 'DESeq2 unadjusted p-values', xlab='Unadjusted p-values',
     col = '#F4C7AB')

result <- as.data.frame(res.deseq2)
tranla_diff <- result[abs(result$log2FoldChange) >1 & result$padj< 0.05,]
tranla_diff <- tranla_diff[tranla_diff$baseMean > 10,]
tranla_diff <- na.omit(tranla_diff)
tranla_up <- tranla_diff[tranla_diff$log2FoldChange >0,]
tranla_down <- tranla_diff[tranla_diff$log2FoldChange < 0,]

setwd("E:/220625_PC/R workplace/220910_annotation/Drosaphila/")
library(rtracklayer)
library(dplyr)
genome.anno <- import("dmel-all-r6.44.gtf") %>% as.data.frame()
genome.anno=genome.anno[,c("gene_id","gene_symbol","seqnames","start","end","width",
                           "strand","type")] 
head(genome.anno)
genome.anno <- genome.anno[genome.anno$type=="gene",]
colnames(genome.anno) <- c("GeneID","gene_symbol","seqnames","start","end","length","strand","type")

result$GeneID <- rownames(result)
T_S2 <- merge(result,genome.anno, by="GeneID", all.x=TRUE)
Ribo_exp$GeneID <-  rownames(Ribo_exp)
T_S2 <- merge(T_S2,Ribo_exp,by="GeneID", all.x=TRUE)
RNA_exp$GeneID <- rownames(RNA_exp)
T_S2 <- merge(T_S2,RNA_exp,by="GeneID", all.x=TRUE)

T_diff <- T_S2[abs(T_S2$log2FoldChange) > 1 & T_S2$pvalue < 0.05,]
T_diff <- dplyr::filter(T_diff, !is.na(GeneID))
T_test <- T_diff[T_diff$baseMean <= 10,]
T_diff <- T_diff[T_diff$baseMean > 10,]
T_up <- T_diff[T_diff$log2FoldChange > 0,]
T_down <- T_diff[T_diff$log2FoldChange < 0,]
```

