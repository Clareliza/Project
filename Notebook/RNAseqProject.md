RNA seq project
#Vamos a trabajar con R para realizar los análisis estadísticos:
#1.-Instalamos los paquetes de bioconductor:

source("http://bioconductor.org/biocLite.R")
#Instalar DESeq2
biocLite("DESeq2")
#Una vez instalado, se instala library:
library(DESeq2)
#2.- Ahora vamos a bajar los datos 
data <- read.table("/Users/alumnomatlab/Desktop/Project/data/Cgigas-HS-count.txt", header = T, sep = "\t")
rownames(data) <- data$Feature
data <- data[,-1]
# Vamos a construir la estructura en columnas 
#Build Objects
# Specify which columns are in which groups
```
deseq2.colData <- data.frame(condition=factor(c(rep("Treated", 3), rep("Control", 3))), 
							 type=factor(rep("single-read", 6)))
rownames(deseq2.colData) <- colnames(data)
deseq2.dds <- DESeqDataSetFromMatrix(countData = data,
									 colData = deseq2.colData, 
									 design = ~ condition)
 ```
#Para correr el análisis:
```deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ] 
 ```
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#PARA VER LA INFORMACION:
 ```
head(deseq2.res)  
```				 
# Count number of hits with adjusted p-value less then 0.05
 ```
dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])
tmp <- deseq2.res
```
# The main plot
 ```
plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
	 main="DEG Crassostrea gigas Heat Shock  (pval <= 0.05)",
	 xlab="mean of normalized counts",
	 ylab="Log2 Fold Change")
```
# Getting the significant points and plotting them again so they're a different color
 ```
 tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="red")
```
# 2 FC lines
 ```
 abline(h=c(-1,1), col="blue")
 ```
#Guardar la informacion como tabla.
 
  ```
  write.table(tmp.sig, "/Users/alumnomatlab/Desktop/Project/data/Cgigas-HS-count.tab", row.names = T)
!plot(plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
	 main="DEG Crassostrea gigas Heat Shock  (pval <= 0.05)",
	 xlab="mean of normalized counts",
	 ylab="Log2 Fold Change")
 ```
 
# Getting the significant points and plotting them again so they're a different color
 
 ```
  tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="red")
```
# 2 FC lines
```
abline(h=c(-1,1), col="blue")
```
#Guardar la informacion como tabla.

```
write.table(tmp.sig, "/Users/alumnomatlab/Desktop/Project/data/Cgigas-HS-count.tab", row.names = T)
```
#Grafica
```
!plot
```