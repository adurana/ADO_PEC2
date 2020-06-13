
#1. Preparación del entorno

setwd("~/Desktop/Master UOC/Semestre 3/Asignaturas/
      157-Analisis de datos omicos/PEC2/ADO_PEC2")
workingDir <- getwd()
dir.create("data")
dir.create("results")
dataDir <- file.path(workingDir, "data/")
resultsDir <- file.path(workingDir, "results/")

#2. Lectura de datos

counts <- read.csv("./data/counts.csv", sep = ";", header = T)
targets <- read.csv("./data/targets.csv", header = T)

head(counts[,1:3],3)
head(targets, 3)

#3. Modificación de los objetos `counts` y `targets`

names(counts)[1] <- "Ensembl"
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
                               counts[colnames(counts) %in%
                                        targets_ELI_vs_NIT$Sample_Name])
counts_ELI_vs_NIT <- counts_red_ELI_vs_NIT[,-1]
rownames(counts_ELI_vs_NIT) <- counts_red_ELI_vs_NIT[,1]
all(rownames(targets_ELI_vs_NIT$Sample_Name) == colnames(counts_ELI_vs_NIT))

counts_red_ELI_vs_SFI <- cbind(counts$Ensembl,
                               counts[colnames(counts) %in%
                                        targets_ELI_vs_SFI$Sample_Name])
counts_ELI_vs_SFI <- counts_red_ELI_vs_SFI[,-1]
rownames(counts_ELI_vs_SFI) <- counts_red_ELI_vs_SFI[,1]
all(rownames(targets_ELI_vs_SFI$Sample_Name) == colnames(counts_ELI_vs_SFI))

counts_red_NIT_vs_SFI <- cbind(counts$Ensembl,
                               counts[colnames(counts) %in%
                                        targets_NIT_vs_SFI$Sample_Name])
counts_NIT_vs_SFI <- counts_red_NIT_vs_SFI[,-1]
rownames(counts_NIT_vs_SFI) <- counts_red_NIT_vs_SFI[,1]
all(rownames(targets_NIT_vs_SFI$Sample_Name) == colnames(counts_NIT_vs_SFI))

head(counts_ELI_vs_NIT[,1:3], 3)
head(counts_ELI_vs_SFI[,1:3], 3)
head(counts_NIT_vs_SFI[,1:3], 3)


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

#7. Filtraje de los datos

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

#8. Estabilización de la varianza

vsd_ELI_vs_NIT <- vst(dds_ELI_vs_NIT, blind = F)
vsd_ELI_vs_SFI <- vst(dds_ELI_vs_SFI, blind = F)
vsd_NIT_vs_SFI <- vst(dds_NIT_vs_SFI, blind = F)

rld_ELI_vs_NIT <- rlog(dds_ELI_vs_NIT, blind = F)
rld_ELI_vs_SFI <- rlog(dds_ELI_vs_SFI, blind = F)
rld_NIT_vs_SFI <- rlog(dds_NIT_vs_SFI, blind = F)

#9. Distancia entre las muestras

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

#10. Análisis de componentes principales (PCA)

plotPCA(vsd_ELI_vs_NIT, intgroup="Group")

plotPCA(vsd_ELI_vs_SFI, intgroup="Group")

plotPCA(vsd_NIT_vs_SFI, intgroup="Group")

#11. Análisis de expresión diferencial

dds_ELI_vs_NIT <- DESeq(dds_ELI_vs_NIT)
res_ELI_vs_NIT <- results(dds_ELI_vs_NIT)

dds_ELI_vs_SFI <- DESeq(dds_ELI_vs_SFI)
res_ELI_vs_SFI <- results(dds_ELI_vs_SFI)

dds_NIT_vs_SFI <- DESeq(dds_NIT_vs_SFI)
res_NIT_vs_SFI <- results(dds_NIT_vs_SFI)

dds_ELI_vs_NIT
head(res_ELI_vs_NIT, 5)

dds_ELI_vs_SFI
head(res_ELI_vs_SFI,5)

dds_NIT_vs_SFI
head(res_NIT_vs_SFI,5)

resultsNames(dds_ELI_vs_NIT)
resultsNames(dds_ELI_vs_SFI)
resultsNames(dds_NIT_vs_SFI)

#12. Contracción del cambio de log fold

library("apeglm")

resLFC_ELI_vs_NIT <- lfcShrink(dds_ELI_vs_NIT,
                               coef="Group_NIT_vs_ELI", type = "apeglm")
resLFC_ELI_vs_SFI <- lfcShrink(dds_ELI_vs_SFI,
                               coef="Group_SFI_vs_ELI", type = "apeglm")
resLFC_NIT_vs_SFI <- lfcShrink(dds_NIT_vs_SFI,
                               coef="Group_SFI_vs_NIT", type = "apeglm")

head(resLFC_ELI_vs_NIT, 5)
head(resLFC_ELI_vs_SFI, 5)
head(resLFC_NIT_vs_SFI, 5)

#13. Counts plot

plotCounts(dds_ELI_vs_NIT,
           gene=which.min(res_ELI_vs_NIT$padj), intgroup="Group")

plotCounts(dds_ELI_vs_SFI,
           gene=which.min(res_ELI_vs_SFI$padj), intgroup="Group")

plotCounts(dds_NIT_vs_SFI,
           gene=which.min(res_NIT_vs_SFI$padj), intgroup="Group")

#14. MA plot

par(mfrow=c(2,1))

plotMA(res_ELI_vs_NIT, ylim=c(-2.5,2.5), main="ELI_vs_NIT")
plotMA(resLFC_ELI_vs_NIT, ylim=c(-2.5,2.5), main="ELI_vs_NIT shrinkage")

par(mfrow=c(2,1))

plotMA(res_ELI_vs_SFI, ylim=c(-2.5,2.5), main="ELI_vs_SFI")
plotMA(resLFC_ELI_vs_SFI, ylim=c(-2.5,2.5), main="ELI_vs_SFI shrinkage")

par(mfrow=c(2,1))

plotMA(res_NIT_vs_SFI, ylim=c(-2.5,2.5), main="NIT_vs_SFI")
plotMA(resLFC_NIT_vs_SFI, ylim=c(-2.5,2.5), main="NIT_vs_SFI shrinkage")

#15. Gene clustering

library("genefilter")

topVarGenes_ELI_vs_NIT <- head(order(rowVars(assay(vsd_ELI_vs_NIT)),
                                     decreasing = TRUE), 20)
mat_ELI_vs_NIT  <- assay(vsd_ELI_vs_NIT)[topVarGenes_ELI_vs_NIT, ]
mat_ELI_vs_NIT  <- mat_ELI_vs_NIT - rowMeans(mat_ELI_vs_NIT)
Group <- colData(vsd_ELI_vs_NIT)[,"Group"]
anno_ELI_vs_NIT <- as.data.frame(Group)
rownames(anno_ELI_vs_NIT) <- colnames(mat_ELI_vs_NIT)
pheatmap(mat_ELI_vs_NIT, annotation_col = anno_ELI_vs_NIT,
         main = "Agrupamiento de genes para ELI_vs_NIT")

topVarGenes_ELI_vs_SFI <- head(order(rowVars(assay(vsd_ELI_vs_SFI)),
                                     decreasing = TRUE), 20)
mat_ELI_vs_SFI  <- assay(vsd_ELI_vs_SFI)[topVarGenes_ELI_vs_SFI, ]
mat_ELI_vs_SFI  <- mat_ELI_vs_SFI - rowMeans(mat_ELI_vs_SFI)
Group <- colData(vsd_ELI_vs_SFI)[,"Group"]
anno_ELI_vs_SFI <- as.data.frame(Group)
rownames(anno_ELI_vs_SFI) <- colnames(mat_ELI_vs_SFI)
pheatmap(mat_ELI_vs_SFI, annotation_col = anno_ELI_vs_SFI,
         main = "Agrupamiento de genes para ELI_vs_SFI")

topVarGenes_NIT_vs_SFI <- head(order(rowVars(assay(vsd_NIT_vs_SFI)),
                                     decreasing = TRUE), 20)
mat_NIT_vs_SFI  <- assay(vsd_NIT_vs_SFI)[topVarGenes_NIT_vs_SFI, ]
mat_NIT_vs_SFI  <- mat_NIT_vs_SFI - rowMeans(mat_NIT_vs_SFI)
Group <- colData(vsd_NIT_vs_SFI)[,"Group"]
anno_NIT_vs_SFI <- as.data.frame(Group)
rownames(anno_NIT_vs_SFI) <- colnames(mat_NIT_vs_SFI)
pheatmap(mat_NIT_vs_SFI, annotation_col = anno_NIT_vs_SFI,
         main = "Agrupamiento de genes para NIT_vs_SFI")

#16. Anotación

library("AnnotationDbi")
library("org.Hs.eg.db")

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
res_ord_ELI_vs_NIT <- res_ELI_vs_NIT[order(res_ELI_vs_NIT$pvalue),]
head(res_ord_ELI_vs_NIT)

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
res_ord_ELI_vs_SFI <- res_ELI_vs_SFI[order(res_ELI_vs_SFI$pvalue),]
head(res_ord_ELI_vs_SFI)

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
res_ord_NIT_vs_SFI <- res_NIT_vs_SFI[order(res_NIT_vs_SFI$pvalue),]
head(res_ord_NIT_vs_SFI)

#17. Resumen de resultados

write.csv(counts_ELI_vs_NIT, file = "./results/counts_ELI_vs_NIT.csv")
write.csv(counts_ELI_vs_SFI, file = "./results/counts_ELI_vs_SFI.csv")
write.csv(counts_NIT_vs_SFI, file = "./results/counts_NIT_vs_SFI.csv")

write.csv(targets_ELI_vs_NIT, file = "./results/targets_ELI_vs_NIT.csv")
write.csv(targets_ELI_vs_SFI, file = "./results/targets_ELI_vs_SFI.csv")
write.csv(targets_NIT_vs_SFI, file = "./results/targets_NIT_vs_SFI.csv")

write.csv(as.data.frame(res_ELI_vs_NIT), file = "./results/res_ELI_vs_NIT.csv")
write.csv(as.data.frame(res_ELI_vs_SFI), file = "./results/res_ELI_vs_SFI.csv")
write.csv(as.data.frame(res_NIT_vs_SFI), file = "./results/res_NIT_vs_SFI.csv")

write.csv(as.data.frame(res_ord_ELI_vs_NIT), file = "./results/res_ord_ELI_vs_NIT.csv")
write.csv(as.data.frame(res_ord_ELI_vs_SFI), file = "./results/res_ord_ELI_vs_SFI.csv")
write.csv(as.data.frame(res_ord_NIT_vs_SFI), file = "./results/res_ord_NIT_vs_SFI.csv")

listOfFiles <- dir("./results/") 
knitr::kable(
  listOfFiles, booktabs = TRUE,
  col.names="ListadoFicheros"
)



