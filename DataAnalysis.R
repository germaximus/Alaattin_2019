library(magrittr)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(rstudioapi)
library(dplyr)
library(stringr)
library(pheatmap)
library(data.table)
setwd(dirname(getActiveDocumentContext()$path))

# SK-MEL30 cell line - male
# HT-29 colon cancer cell line - female
# BJ fibroblast male - male

#------------------------------------------- Import data ----------------------------------------------------
SampleTable <- read.table(file = "./GeneCounts/SampleTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

GeneCounts <- lapply(SampleTable$file, function(x) {
              data <- read.table(file = paste0("./GeneCounts/", x), skip=4, row.names = 1, stringsAsFactors = FALSE) %>% select(1)
}) %>% as.data.frame() %>% setNames(SampleTable$file)

#------------------------------------------ Matrix plot -----------------------------------------------------
png('MatrixPlot.png', width = 1200, height = 1200, units = 'px', res = 300)
pairs(GeneCounts, log = "yx", pch = '.', gap = 0.25)
dev.off()

#-------------------------------------- PCA plot, all cell lines --------------------------------------------
PlotData <- GeneCounts[rowMeans(GeneCounts) >= 10, ] %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = paste0(SampleTable$CellLine, '-', SampleTable$TreatmentTime),
                   name  = SampleTable$file
)

png(file = "PCAplot.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
   geom_point(size = 6) + 
   xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
   ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
   theme_bw() +
   theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
   theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
   theme(axis.title = element_text(size = rel(1.5)))+
   theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()

#--------------------------- PCA plot, HT-29 colon cancer cell line separately -------------------------------
PlotData <- GeneCounts %>% select(grep('HT', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = factor(paste0(SampleTable$CellLine[10:18], '-', SampleTable$TreatmentTime[10:18]), 
                                        levels = unique(paste0(SampleTable$CellLine[10:18], '-', SampleTable$TreatmentTime[1:9]))),
                         name  = SampleTable$file[10:18]
)

png(file = "PCAplot-HT29.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()


#-------------------------------- PCA plot, BJ fibroblast cell line separately --------------------------------
PlotData <- GeneCounts %>% select(grep('BJ', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = factor(paste0(SampleTable$CellLine[1:9], '-', SampleTable$TreatmentTime[1:9]), 
                                        levels = unique(paste0(SampleTable$CellLine[1:9], '-', SampleTable$TreatmentTime[1:9]))),
                         name  = SampleTable$file[1:9]
)

png(file = "PCAplot-BJ.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()

#----------------------------------- Correlation Heatmap, all cell lines ---------------------------------------
PlotData  <- GeneCounts[rowMeans(GeneCounts) >= 10, ] %>% as.matrix() %>% rlog()
colors    <- colorRampPalette(brewer.pal(9, "GnBu"))(100) %>% rev()
distance  <- t(PlotData) %>% dist() %>% as.matrix()

png("CorrelatioHeatmap.png", width = 2400, height = 2400, res = 300, unit = 'px')
pheatmap(distance, color = colors, border_color = 'white')
dev.off()

#------------------------- Correlation Heatmap, HT-29 colon cancer cell line separately ------------------------
PlotData <- GeneCounts %>% select(grep('HT', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()
colors    <- colorRampPalette(brewer.pal(9, "GnBu"))(100) %>% rev()
distance  <- t(PlotData) %>% dist() %>% as.matrix()

png("CorrelatioHeatmap-HT29.png", width = 2400, height = 2400, res = 300, unit = 'px')
pheatmap(distance, color = colors, border_color = 'white')
dev.off()

#--------------------------- Correlation Heatmap, BJ fibroblast cell line separately --------------------------
PlotData <- GeneCounts %>% select(grep('BJ', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()
colors    <- colorRampPalette(brewer.pal(9, "GnBu"))(100) %>% rev()
distance  <- t(PlotData) %>% dist() %>% as.matrix()

png("CorrelatioHeatmap-BJ.png", width = 2400, height = 2400, res = 300, unit = 'px')
pheatmap(distance, color = colors, border_color = 'white')
dev.off()



#--------------------------- Heatmap genes vs samples, BJ fibroblast cell line separately ---------------------
PlotData <- GeneCounts %>% select(grep('BJ', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()

png("GenesHeatmap-BJ.png", width = 2400, height = 2400, res = 300, unit = 'px')
colors    <- colorRampPalette(c("blue","black","yellow"))(256)
annotation_col <- data.frame(TreatmentTime = c(rep('untreated', 3), rep('30 min', 3), rep('120 min', 3)))
row.names(annotation_col) <- colnames(PlotData)

pheatmap(PlotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = F,
         annotation_col = annotation_col, annotation_names_col = FALSE,
         color = colors,  border_color = NA, show_colnames = FALSE, show_rownames = FALSE)
dev.off()


#--------------------------- Heatmap genes vs samples, HT fibroblast cell line separately ---------------------
PlotData <- GeneCounts %>% select(grep('HT', colnames(.))) %>% filter(rowMeans(.) >= 10) %>% as.matrix() %>% rlog()

png("GenesHeatmap-HT.png", width = 2400, height = 2400, res = 300, unit = 'px')
colors    <- colorRampPalette(c("blue","black","yellow"))(256)
annotation_col <- data.frame(TreatmentTime = c(rep('untreated', 3), rep('30 min', 3), rep('120 min', 3)))
row.names(annotation_col) <- colnames(PlotData)

pheatmap(PlotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = F,
         annotation_col = annotation_col, annotation_names_col = FALSE,
         color = colors,  border_color = NA, show_colnames = FALSE, show_rownames = FALSE)
dev.off()



###############################################################################################################################
###################################   Differential Gene Expression Analysis   #################################################
###############################################################################################################################

# Create a relational table between annotation gene_id and NCBI gene name
annotation <- fread(file="GRCh38.p12.Refseq.coding.gff", 
                    skip = 8, 
                    stringsAsFactors = F, 
                    header = F, fill = T, 
                    na.strings = c("", "NA"), 
                    sep="\t") %>% 
              na.omit() %>%
              filter(V3 == 'gene')

local_gene_id   <- str_match(annotation$V9, 'ID=gene([^;]+)')
Entrez_gene_id   <- str_match(annotation$V9, 'Dbxref=GeneID:([^,;]+)')
geneName <- str_match(annotation$V9, 'Name=([^;]+)')
RelationTable <- data.frame(local_gene_id = local_gene_id[,2], Entrez_gene_id = Entrez_gene_id[,2], GeneName = geneName[,2])
write.table(RelationTable, './GeneCounts/GeneNames.txt', sep = '\t', row.names = F) #save this file for future use


# Import relational table between NCBI gene IDs and NCBI gene names
RelationTable <- read.table('./GeneCounts/GeneNames.txt', sep = '\t', header = T)


#----------------------------------  HT-29 colon cancer cell line  ------------------------------------
colData <- data.frame(sample = colnames(GeneCounts)[10:18],
                      TreatmentTime = factor(c(rep("0min", 3), rep("30min",3), rep("120min",3)), levels = c("0min","30min","120min")),
                      stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(countData = GeneCounts[, 10:18], colData = colData, design =~ TreatmentTime)
dds <- dds[rowMeans(counts(dds)) >= 10, ] 
dds <- DESeq(dds)    
result <- results(dds, contrast = c("TreatmentTime", "30min",  "0min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
summary(result)
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result, file="HT-29_30min_vs_zero.csv")

result <- results(dds, contrast = c("TreatmentTime", "120min",  "0min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
summary(result)
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result,file="HT-29_120min_vs_zero.csv")

result <- results(dds, contrast = c("TreatmentTime", "120min",  "30min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
summary(result)
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result,file="HT-29_120min_vs_30min.csv")

#----------------------------------  BJ fibroblast cell line  ------------------------------------
colData <- data.frame(sample = colnames(GeneCounts)[1:9],
                      TreatmentTime = factor(c(rep("0min", 3), rep("30min",3), rep("120min",3)), levels = c("0min","30min","120min")),
                      stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(countData = GeneCounts[, 1:9], colData = colData, design =~ TreatmentTime)
dds <- dds[rowMeans(counts(dds)) >= 10, ] 
dds <- DESeq(dds)    
result <- results(dds, contrast = c("TreatmentTime", "30min",  "0min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
summary(result)
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result,file="BJ_30min_vs_zero.csv")

result <- results(dds, contrast = c("TreatmentTime", "120min",  "0min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result,file="BJ_120min_vs_zero.csv")

result <- results(dds, contrast = c("TreatmentTime", "120min",  "30min"), cooksCutoff = FALSE, independentFiltering = FALSE)
result <- result[order(result$padj),]
summary(result)
result <- as.data.frame(result)
names <- row.names(result)
ids <- RelationTable[match(substring(row.names(result), 5), RelationTable$local_gene_id), 2:3]
result <- cbind(ids, result) %>% set_rownames(names)
write.csv(result,file="BJ_120min_vs_30min.csv")






















