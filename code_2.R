
#1. Preparación del entorno
setwd("~/Desktop/Master UOC/Semestre 3/Asignaturas/157-Analisis de datos omicos/PEC2/ADO_PEC2")
workingDir <- getwd()
dir.create("data")
dir.create("results")
dataDir <- file.path(workingDir, "data/")
resultsDir <- file.path(workingDir, "results/")

#2. Lectura de datos

counts <- read.csv("./counts.csv", sep = ";", header = T)
targets <- read.csv("./targets.csv", header = T)

head(counts[,1:3],3)
head(targets, 3)

#3. Modificación de los objetos `counts` y `targets`

head(names(counts))
names(counts)[1] <- "Ensembl"
head(names(counts))
head(counts$Ensembl)
counts$Ensembl <- gsub("\\..*", "", counts$Ensembl, fixed = FALSE)
head(counts$Ensembl)

targets$Sample_Name <- gsub("-", ".", x = as.character(targets$Sample_Name))
head(targets, 3)




#4. Selección aleatoria de muestras

targets_NIT <- targets[targets$Group=="NIT",]
targets_ELI <- targets[targets$Group=="ELI",]
targets_SFI <- targets[targets$Group=="SFI",]

set.seed(2020)

targets_ELI_rndm <- targets_ELI[sample(nrow(targets_ELI), 10),]
targets_SFI_rndm <- targets_SFI[sample(nrow(targets_SFI), 10),]
targets_NIT_rndm <- targets_NIT[sample(nrow(targets_NIT), 10),]

targets_ELI_rndm
targets_SFI_rndm
targets_NIT_rndm

targets_ELI_vs_NIT <- do.call(what = rbind,
                              args =list(targets_ELI_rndm, targets_NIT_rndm))

targets_ELI_vs_SFI <- do.call(what = rbind,
                              args =list(targets_ELI_rndm, targets_SFI_rndm))

targets_NIT_vs_SFI <- do.call(what = rbind,
                              args =list(targets_NIT_rndm, targets_SFI_rndm))


#5. Modificación del archivo `counts` según el muestreo


counts_red_ELI_vs_NIT <- cbind(counts$Ensembl,
                        counts[colnames(counts) %in% targets_ELI_vs_NIT$Sample_Name])
counts_ELI_vs_NIT <- counts_red_ELI_vs_NIT[,-1]
rownames(counts_ELI_vs_NIT) <- counts_red_ELI_vs_NIT[,1]
all(rownames(targets_ELI_vs_NIT$Sample_Name) == colnames(counts_ELI_vs_NIT))

counts_red_ELI_vs_SFI <- cbind(counts$Ensembl,
                               counts[colnames(counts) %in% targets_ELI_vs_SFI$Sample_Name])
counts_ELI_vs_SFI <- counts_red_ELI_vs_SFI[,-1]
rownames(counts_ELI_vs_SFI) <- counts_red_ELI_vs_SFI[,1]
all(rownames(targets_ELI_vs_SFI$Sample_Name) == colnames(counts_ELI_vs_SFI))

counts_red_NIT_vs_SFI <- cbind(counts$Ensembl,
                               counts[colnames(counts) %in% targets_NIT_vs_SFI$Sample_Name])
counts_NIT_vs_SFI <- counts_red_NIT_vs_SFI[,-1]
rownames(counts_NIT_vs_SFI) <- counts_red_NIT_vs_SFI[,1]
all(rownames(targets_NIT_vs_SFI$Sample_Name) == colnames(counts_NIT_vs_SFI))


#6. Construcción de los objetos `DESeqDataSet` usando de datos sin normalizar


library("DESeq2")

dds_ELI_vs_NIT <- DESeqDataSetFromMatrix(countData = counts_ELI_vs_NIT,
                              colData = targets_ELI_vs_NIT,
                              design = ~ Group)
dds_ELI_vs_SFI <- DESeqDataSetFromMatrix(countData = counts_ELI_vs_SFI,
                              colData = targets_ELI_vs_SFI,
                              design = ~ Group)
dds_NIT_vs_SFI <- DESeqDataSetFromMatrix(countData = counts_NIT_vs_SFI,
                              colData = targets_NIT_vs_SFI,
                              design = ~ Group)


dds_ELI_vs_NIT
dds_ELI_vs_SFI
dds_NIT_vs_SFI

nrow(dds_ELI_vs_NIT)
keep_ELI_vs_NIT <- rowSums(counts(dds_ELI_vs_NIT)) > 1
dds_ELI_vs_NIT <- dds_ELI_vs_NIT[keep_ELI_vs_NIT,]
nrow(dds_ELI_vs_NIT)

nrow(dds_ELI_vs_SFI)
keep_ELI_vs_SFI <- rowSums(counts(dds_ELI_vs_SFI)) > 1
dds_ELI_vs_SFI <- dds_ELI_vs_SFI[keep_ELI_vs_SFI,]
nrow(dds_ELI_vs_SFI)

nrow(dds_NIT_vs_SFI)
keep_NIT_vs_SFI <- rowSums(counts(dds_NIT_vs_SFI)) > 1
dds_NIT_vs_SFI <- dds_NIT_vs_SFI[keep_NIT_vs_SFI,]
nrow(dds_NIT_vs_SFI)




dds_ELI_vs_NIT <- DESeq(dds_ELI_vs_NIT)
dds_ELI_vs_NIT
res_ELI_vs_NIT <- results(dds_ELI_vs_NIT)
res_ELI_vs_NIT


dds_ELI_vs_SFI <- DESeq(dds_ELI_vs_SFI)
dds_ELI_vs_SFI
res_ELI_vs_SFI <- results(dds_ELI_vs_SFI)
res_ELI_vs_SFI


dds_NIT_vs_SFI <- DESeq(dds_NIT_vs_SFI)
dds_NIT_vs_SFI
res_NIT_vs_SFI <- results(dds_NIT_vs_SFI)
res_NIT_vs_SFI


head(res_ELI_vs_NIT, 5)


resLFC_ELI_vs_NIT <- lfcShrink(dds_ELI_vs_NIT,
                               coef="Group_NIT_vs_ELI", type="apeglm")
resLFC_ELI_vs_NIT


resLFC_ELI_vs_SFI <- lfcShrink(dds_ELI_vs_SFI,
                               coef="Group_SFI_vs_ELI", type="apeglm")
resLFC_ELI_vs_SFI


resLFC_NIT_vs_SFI <- lfcShrink(dds_NIT_vs_SFI,
                               coef="Group_SFI_vs_NIT", type="apeglm")
resLFC_NIT_vs_SFI

#añadir despues de filtraje

vsd_ELI_vs_NIT <- vst(dds_ELI_vs_NIT, blind = F)
head(assay(vsd_ELI_vs_NIT), 3)

vsd_ELI_vs_SFI <- vst(dds_ELI_vs_SFI, blind = F)
head(assay(vsd_ELI_vs_SFI), 3)

vsd_NIT_vs_SFI <- vst(dds_NIT_vs_SFI, blind = F)
head(assay(vsd_NIT_vs_SFI), 3)




rld_ELI_vs_NIT <- rlog(dds_ELI_vs_NIT, blind = F)
head(assay(rld_ELI_vs_NIT), 3)

rld_ELI_vs_SFI <- rlog(dds_ELI_vs_SFI, blind = F)
head(assay(rld_ELI_vs_SFI), 3)

rld_NIT_vs_SFI <- rlog(dds_NIT_vs_SFI, blind = F)
head(assay(rld_NIT_vs_SFI), 3)




sampleDists_ELI_vs_NIT <- dist(t(assay(vsd_ELI_vs_NIT)))
sampleDists_ELI_vs_SFI <- dist(t(assay(vsd_ELI_vs_SFI)))
sampleDists_NIT_vs_SFI <- dist(t(assay(vsd_NIT_vs_SFI)))


library("pheatmap")
library("RColorBrewer")

sampleDistMatrix_ELI_vs_NIT <- as.matrix(sampleDists_ELI_vs_NIT)
rownames(sampleDistMatrix_ELI_vs_NIT) <- vsd_ELI_vs_NIT$Group
colnames(sampleDistMatrix_ELI_vs_NIT) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_ELI_vs_NIT,
         clustering_distance_rows = sampleDists_ELI_vs_NIT,
         clustering_distance_cols = sampleDists_ELI_vs_NIT,
         col = colors)

sampleDistMatrix_ELI_vs_SFI <- as.matrix(sampleDists_ELI_vs_SFI)
rownames(sampleDistMatrix_ELI_vs_SFI) <- vsd_ELI_vs_SFI$Group
colnames(sampleDistMatrix_ELI_vs_SFI) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_ELI_vs_SFI,
         clustering_distance_rows = sampleDists_ELI_vs_SFI,
         clustering_distance_cols = sampleDists_ELI_vs_SFI,
         col = colors)

sampleDistMatrix_NIT_vs_SFI <- as.matrix(sampleDists_NIT_vs_SFI)
rownames(sampleDistMatrix_NIT_vs_SFI) <- vsd_NIT_vs_SFI$Group
colnames(sampleDistMatrix_NIT_vs_SFI) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_NIT_vs_SFI,
         clustering_distance_rows = sampleDists_NIT_vs_SFI,
         clustering_distance_cols = sampleDists_NIT_vs_SFI,
         col = colors)


#plotting despues de DEA

par(mfrow=c(2,1))
plotMA(res_ELI_vs_NIT, ylim=c(-2.5,2.5), main="ELI_vs_NIT")
plotMA(resLFC_ELI_vs_NIT, ylim=c(-2.5,2.5), main="ELI_vs_NIT shrinkage")

plotMA(res_ELI_vs_SFI, ylim=c(-2.5,2.5), main="ELI_vs_SFI")
plotMA(resLFC_ELI_vs_SFI, ylim=c(-2.5,2.5), main="ELI_vs_SFI shrinkage")

plotMA(res_NIT_vs_SFI, ylim=c(-2.5,2.5), main="NIT_vs_SFI")
plotMA(resLFC_NIT_vs_SFI, ylim=c(-2.5,2.5), main="NIT_vs_SFI shrinkage")

dev.off()


plotCounts(dds_ELI_vs_NIT,
           gene=which.min(res_ELI_vs_NIT$padj), intgroup="Group")
plotCounts(dds_ELI_vs_SFI,
           gene=which.min(res_ELI_vs_SFI$padj), intgroup="Group")
plotCounts(dds_NIT_vs_SFI,
           gene=which.min(res_NIT_vs_SFI$padj), intgroup="Group")

plotPCA(vsd_ELI_vs_NIT, intgroup="Group")
plotPCA(vsd_ELI_vs_SFI, intgroup="Group")
plotPCA(vsd_NIT_vs_SFI, intgroup="Group")


library(vsn)

ntd_ELI_vs_NIT <- normTransform(dds_ELI_vs_NIT)
ntd_ELI_vs_SFI <- normTransform(dds_ELI_vs_SFI)
ntd_NIT_vs_SFI <- normTransform(dds_NIT_vs_SFI)


meanSdPlot(assay(ntd_ELI_vs_NIT))
meanSdPlot(assay(vsd_ELI_vs_NIT))
meanSdPlot(assay(rld_ELI_vs_NIT))


meanSdPlot(assay(ntd_ELI_vs_SFI))
meanSdPlot(assay(vsd_ELI_vs_SFI))
meanSdPlot(assay(rld_ELI_vs_SFI))


meanSdPlot(assay(ntd_NIT_vs_SFI))
meanSdPlot(assay(vsd_NIT_vs_SFI))
meanSdPlot(assay(rld_NIT_vs_SFI))


#QC

library("pheatmap")
library("RColorBrewer")
library("genefilter")

topVarGenes_ELI_vs_NIT <- head(order(rowVars(assay(vsd_ELI_vs_NIT)),
                                     decreasing = TRUE), 20)
mat_ELI_vs_NIT  <- assay(vsd_ELI_vs_NIT)[topVarGenes_ELI_vs_NIT, ]
mat_ELI_vs_NIT  <- mat_ELI_vs_NIT - rowMeans(mat_ELI_vs_NIT)
anno_ELI_vs_NIT <- as.data.frame(colData(vsd_ELI_vs_NIT)[,"Group"])
rownames(anno_ELI_vs_NIT) <- colnames(mat_ELI_vs_NIT)
pheatmap(mat_ELI_vs_NIT, annotation_col = anno_ELI_vs_NIT)


topVarGenes_ELI_vs_SFI <- head(order(rowVars(assay(vsd_ELI_vs_SFI)),
                                     decreasing = TRUE), 20)
mat_ELI_vs_SFI  <- assay(vsd_ELI_vs_SFI)[topVarGenes_ELI_vs_SFI, ]
mat_ELI_vs_SFI  <- mat_ELI_vs_SFI - rowMeans(mat_ELI_vs_SFI)
anno_ELI_vs_SFI <- as.data.frame(colData(vsd_ELI_vs_SFI)[,"Group"])
rownames(anno_ELI_vs_SFI) <- colnames(mat_ELI_vs_SFI)
pheatmap(mat_ELI_vs_SFI, annotation_col = anno_ELI_vs_SFI)

topVarGenes_NIT_vs_SFI <- head(order(rowVars(assay(vsd_NIT_vs_SFI)),
                                     decreasing = TRUE), 20)
mat_NIT_vs_SFI  <- assay(vsd_NIT_vs_SFI)[topVarGenes_NIT_vs_SFI, ]
mat_NIT_vs_SFI  <- mat_NIT_vs_SFI - rowMeans(mat_NIT_vs_SFI)
anno_NIT_vs_SFI <- as.data.frame(colData(vsd_NIT_vs_SFI)[,"Group"])
rownames(anno_NIT_vs_SFI) <- colnames(mat_NIT_vs_SFI)
pheatmap(mat_NIT_vs_SFI, annotation_col = anno_NIT_vs_SFI)





#heatmap matriz contajes
select_ELI_vs_NIT <- order(rowMeans(counts(dds_ELI_vs_NIT,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df_ELI_vs_NIT <- as.data.frame(colData(dds_ELI_vs_NIT)[,"Group"])
rownames(df_ELI_vs_NIT) <- colnames(ntd_ELI_vs_NIT)

pheatmap(assay(ntd_ELI_vs_NIT)[select_ELI_vs_NIT,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_NIT)

pheatmap(assay(vsd_ELI_vs_NIT)[select_ELI_vs_NIT,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_NIT)

pheatmap(assay(rld_ELI_vs_NIT)[select_ELI_vs_NIT,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_NIT)


select_ELI_vs_SFI <- order(rowMeans(counts(dds_ELI_vs_SFI,normalized=TRUE)),
                           decreasing=TRUE)[1:20]
df_ELI_vs_SFI <- as.data.frame(colData(dds_ELI_vs_SFI)[,"Group"])
rownames(df_ELI_vs_SFI) <- colnames(ntd_ELI_vs_SFI)

pheatmap(assay(ntd_ELI_vs_SFI)[select_ELI_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_SFI)

pheatmap(assay(vsd_ELI_vs_SFI)[select_ELI_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_SFI)

pheatmap(assay(rld_ELI_vs_SFI)[select_ELI_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_ELI_vs_SFI)


select_NIT_vs_SFI <- order(rowMeans(counts(dds_NIT_vs_SFI,normalized=TRUE)),
                           decreasing=TRUE)[1:20]
df_NIT_vs_SFI <- as.data.frame(colData(dds_NIT_vs_SFI)[,"Group"])
rownames(df_NIT_vs_SFI) <- colnames(ntd_NIT_vs_SFI)

pheatmap(assay(ntd_NIT_vs_SFI)[select_NIT_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_NIT_vs_SFI)

pheatmap(assay(vsd_NIT_vs_SFI)[select_NIT_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_NIT_vs_SFI)

pheatmap(assay(rld_NIT_vs_SFI)[select_NIT_vs_SFI,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_NIT_vs_SFI)


#Annotation

library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

res_ELI_vs_NIT$symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res_ELI_vs_NIT),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res_ELI_vs_NIT$entrez <- mapIds(org.Hs.eg.db,
                     keys = row.names(res_ELI_vs_NIT),
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
resOrdered_ELI_vs_NIT <- res_ELI_vs_NIT[order(res_ELI_vs_NIT$pvalue),]
head(resOrdered_ELI_vs_NIT)


res_ELI_vs_SFI$symbol <- mapIds(org.Hs.eg.db,
                                keys = row.names(res_ELI_vs_SFI),
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")
res_ELI_vs_SFI$entrez <- mapIds(org.Hs.eg.db,
                                keys = row.names(res_ELI_vs_SFI),
                                column = "ENTREZID",
                                keytype = "ENSEMBL",
                                multiVals = "first")
resOrdered_ELI_vs_SFI <- res_ELI_vs_SFI[order(res_ELI_vs_SFI$pvalue),]
head(resOrdered_ELI_vs_SFI)

res_NIT_vs_SFI$symbol <- mapIds(org.Hs.eg.db,
                                keys = row.names(res_NIT_vs_SFI),
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")
res_NIT_vs_SFI$entrez <- mapIds(org.Hs.eg.db,
                                keys = row.names(res_NIT_vs_SFI),
                                column = "ENTREZID",
                                keytype = "ENSEMBL",
                                multiVals = "first")
resOrdered_NIT_vs_SFI <- res_NIT_vs_SFI[order(res_NIT_vs_SFI$pvalue),]
head(resOrdered_NIT_vs_SFI)



                     
                     
library("clusterProfiler")

kk <- enrichKEGG(gene         = colnames(dds_ELI_vs_NIT),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

colnames(dds_ELI_vs_NIT)




#preprocesado de datos

counts <- read.csv("./counts.csv", sep = ";", header = T)
targets <- read.csv("./targets.csv", header = T)

names(targets)


targets_NIT <- targets[targets$Group=="NIT",]
targets_ELI <- targets[targets$Group=="ELI",]
targets_SFI <- targets[targets$Group=="SFI",]

nrow(targets_NIT)
nrow(targets_ELI)
nrow(targets_SFI)

set.seed(2020)

targets_ELI_rndm <- targets_ELI[sample(nrow(targets_ELI), 10),]
targets_SFI_rndm <- targets_SFI[sample(nrow(targets_SFI), 10),]
targets_NIT_rndm <- targets_NIT[sample(nrow(targets_NIT), 10),]

targets_ELI_rndm$Sample_Name <- gsub("-", ".", x = as.character(targets_ELI_rndm$Sample_Name))
counts_ELI <- counts[colnames(counts) %in% targets_ELI_rndm$Sample_Name]
targets_ELI_rndm$Sample_Name %in% colnames(counts_ELI)
targets_ELI_rndm$Sample_Name == colnames(counts_ELI)

counts_ELI <- counts_ELI[, rownames(targets_ELI_rndm$Sample_Name)]
all(rownames(targets_ELI_rndm$Sample_Name) == colnames(counts_ELI))
targets_ELI_rndm$Sample_Name %in% colnames(counts_ELI)
targets_ELI_rndm$Sample_Name == colnames(counts_ELI)

targets_SFI_rndm$Sample_Name <- gsub("-", ".", x = as.character(targets_SFI_rndm$Sample_Name))
counts_SFI <- counts[colnames(counts) %in% targets_SFI_rndm$Sample_Name]
targets_SFI_rndm$Sample_Name %in% colnames(counts_SFI)
all(rownames(targets_SFI_rndm$Sample_Name) == colnames(counts_SFI))

targets_NIT_rndm$Sample_Name <- gsub("-", ".", x = as.character(targets_NIT_rndm$Sample_Name))
counts_NIT <- counts[colnames(counts) %in% targets_NIT_rndm$Sample_Name]
targets_NIT_rndm$Sample_Name %in% colnames(counts_NIT)
all(rownames(targets_NIT_rndm$Sample_Name) == colnames(counts_NIT))


write.csv(targets_rndm, file = "targets_rndm.csv")
write.csv(counts_red, file = "counts_red.csv")

targets_rndm <- do.call(what = rbind,
                        args =list(targets_ELI_rndm, targets_NIT_rndm, targets_SFI_rndm))

all(rownames(targets_rndm$Sample_Name) %in% colnames(counts))
all(rownames(targets_rndm$Sample_Name) == colnames(counts))



targets_rndm$Sample_Name <- gsub("-", ".", x = as.character(targets_rndm$Sample_Name))
counts_red <- counts[colnames(counts) %in% targets_rndm$Sample_Name]
counts_red <- cbind(counts$X, counts_red)
dim(counts_red)
colnames(counts_red)
head(rownames(counts_red))

names(counts_red)[1] <- "X"

targets_rndm$Sample_Name %in% colnames(counts_red)

#filtrado y normalizado

targets_rndm[order(targets_rndm$Group),]

BiocManager::install("DESeq2")
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts_red[,2:31],
                              colData = targets_rndm,
                              design = ~ Group)

nrow(dds)
head(rownames(dds))
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

vsd <- vst(dds, blind = F)
head(assay(vsd), 3)

vsd <- vst(dds, blind = T)
head(assay(vsd), 3)
colData(vsd)

rld <- rlog(dds, blind = F)
head(assay(rld), 3)

rld <- rlog(dds, blind = T)
head(assay(rld), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDists

BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
library("pheatmap")
library("RColorBrewer")

vsd$Group

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Group
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

BiocManager::install("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds$Group
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = "Group")

pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

BiocManager::install("glmpca")
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Group <- dds$Group
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Group)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


BiocManager::install("magrittr")
library("magrittr")

mds <- as.data.frame(colData(vsd)) %>% cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Group)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = Group)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")


#DEA

dds <- DESeq(dds)
res_ELI_SFI <- results(dds, contrast = c("Group", "ELI", "SFI"))
res_ELI_NIT <- results(dds, contrast = c("Group", "ELI", "NIT"))
res_SFI_NIT <- results(dds, contrast = c("Group", "SFI", "NIT"))

head(rownames(dds))
head(rownames(res_ELI_SFI))

res_SFI_NIT
res_ELI_NIT
res_ELI_SFI

mcols(res_SFI_NIT, use.names = TRUE)
mcols(res_ELI_NIT, use.names = TRUE)
mcols(res_ELI_SFI, use.names = TRUE)

lapply(c(res_ELI_SFI, res_ELI_NIT,res_SFI_NIT), FUN = summary)
summary(res_ELI_SFI)
summary(res_ELI_NIT)
summary(res_SFI_NIT)

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

sum(res_ELI_NIT$padj < 0.1, na.rm=TRUE)
sum(res_ELI_SFI$padj < 0.1, na.rm=TRUE)
sum(res_SFI_NIT$padj < 0.1, na.rm=TRUE)


resSig_ELI_NIT <- subset(res_ELI_NIT, padj < 0.1)
head(resSig_ELI_NIT[ order(resSig_ELI_NIT$log2FoldChange, decreasing = F), ])

resSig_ELI_SFI <- subset(res_ELI_SFI, padj < 0.1)
head(resSig_ELI_SFI[ order(resSig_ELI_SFI$log2FoldChange, decreasing = F), ])

resSig_SFI_NIT <- subset(res_SFI_NIT, padj < 0.1)
head(resSig_SFI_NIT[ order(resSig_SFI_NIT$log2FoldChange, decreasing = F), ])

#Plotting


topGene_ELI_NIT <- rownames(res_ELI_NIT)[which.min(res_ELI_NIT$padj)]
plotCounts(dds, gene = topGene_ELI_NIT, intgroup=c("Group"))

topGene_ELI_SFI <- rownames(res_ELI_SFI)[which.min(res_ELI_SFI$padj)]
plotCounts(dds, gene = topGene_ELI_SFI, intgroup=c("Group"))

topGene_SFI_NIT <- rownames(res_SFI_NIT)[which.min(res_SFI_NIT$padj)]
plotCounts(dds, gene = topGene_SFI_NIT, intgroup=c("Group"))


BiocManager::install("ggbeeswarm")
library("ggbeeswarm")




geneCounts_ELI_NIT <- plotCounts(dds, gene = topGene_ELI_NIT, intgroup = "Group",
                                 returnData = TRUE)

ggplot(geneCounts_ELI_NIT, aes(x = Group, y = count, color = Group)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts_ELI_NIT, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


#si

BiocManager::install("BiocManager")

BiocManager::install("apeglm")
BiocManager::install("mvtnorm")
library("mvtnorm")
library("apeglm")
resultsNames(dds)

res <- lfcShrink(dds, coef="Group_NIT_vs_ELI", type="apeglm")
plotMA(res, ylim = c(-5, 5))
res.noshr <- results(dds, name="Group_NIT_vs_ELI")
plotMA(res.noshr, ylim = c(-5, 5))

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#gene clustering

BiocManager::install("genefilter")
library("genefilter")

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno)



qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

#Annotation

BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results2.csv")


#no
resGR <- lfcShrink(ddsMat, coef="Group_NIT_vs_ELI", type="apeglm", format="GRanges")
resGR
#si

BiocManager::install("edgeR")
library(edgeR)


BiocManager::install(version = "3.11")



BiocManager::install("DelayedArray")
library("DelayedArray")

BiocManager::install("clusterProfiler")
library("clusterProfiler")
BiocManager::install("magrittr")
library("magrittr")
BiocManager::install("DOSE")
library("DOSE")
BiocManager::install("enrichplot")
library("enrichplot")














