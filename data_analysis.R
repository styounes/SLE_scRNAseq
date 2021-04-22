library(tidyr)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(AnnotationHub)
library(scran)
library(edgeR)
library(PCAtools)
library(pheatmap)
library(viridis)
library(scater)
library(BiocParallel)
library(doParallel)
library(SingleR)
library(monocle3)
library(AUCell)
library(scDblFinder)
library(topGO)
library(EnhancedVolcano)

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

# During our actual workflow, it makes little sense to re-run the analysis every time we boot RStudio
# Consequently, this block of code loads our saved data analysis when needed and is commented out when running the analysis de novo.
# load("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/base.analysis.RData")

sce <- read10xCounts("filtered_feature_bc_matrix")

# Sample and Feature Annotation -------------------------------------------
# This block of code appropriately annotates each sample and adds the genome annotation from AnnotationHub.
# In addition, mitochondrial genes are labeled for subsequent quality control.

working_strain <- character()
batch <- character()
mouse <- character()

for(strain in sce$Barcode){
  if(endsWith(strain, "1") | endsWith(strain, "2") | endsWith(strain, "3")) {
    working_strain <- append(working_strain, "NZB")
  } else if(endsWith(strain, "4") | endsWith(strain, "5") | endsWith(strain, "6")) {
    working_strain <- append(working_strain, "NZW")
  } else if(endsWith(strain, "7") | endsWith(strain, "8") | endsWith(strain, "9")) {
    working_strain <- append(working_strain, "NZBWF1")
  }
  
  if(endsWith(strain, "1") | endsWith(strain, "4") | endsWith(strain, "7")) {
    batch <- append(batch, "1")
  } else if(endsWith(strain, "2") | endsWith(strain, "5") | endsWith(strain, "8")) {
    batch <- append(batch, "2")
  } else if(endsWith(strain, "3") | endsWith(strain, "6") | endsWith(strain, "9")) {
    batch <- append(batch, "3")
  }
  
  if(endsWith(strain, "1")) {
    mouse <- append(mouse, "NZB1")
  } else if(endsWith(strain, "2")){
    mouse <- append(mouse, "NZB2")
  } else if(endsWith(strain, "3")){
    mouse <- append(mouse, "NZB3")
  } else if(endsWith(strain, "4")){
    mouse <- append(mouse, "NZW1")
  } else if(endsWith(strain, "5")){
    mouse <- append(mouse, "NZW2")
  } else if(endsWith(strain, "6")){
    mouse <- append(mouse, "NZW3")
  } else if(endsWith(strain, "7")){
    mouse <- append(mouse, "SLE1")
  } else if(endsWith(strain, "8")){
    mouse <- append(mouse, "SLE2")
  } else if(endsWith(strain, "9")){
    mouse <- append(mouse, "SLE3")
  }
  
}

sce$Strain <- working_strain
sce$Batch <- batch
sce$Mouse <- mouse
sce$Sample <- NULL
rm(working_strain)
rm(batch)
rm(mouse)
rm(strain)

mouse_annotation <- AnnotationHub()[["AH63799"]]
mouse_annotation <- subset(mouse_annotation, mouse_annotation$type == "gene")
mouse_annotation <- subset(mouse_annotation, mouse_annotation$gene_id %in% rownames(sce))
names(mouse_annotation) <- mouse_annotation$gene_name
SummarizedExperiment::rowRanges(sce) <- mouse_annotation[,c(5, 7, 9)]
rm(mouse_annotation)

location <- rowRanges(sce)
names(location) <- location$gene_id
location <- as.data.frame(location)
is.mito <- location$seqnames == "MT"

# Quality Control ---------------------------------------------------------

qc <- perCellQCMetrics(sce, subsets = list(Mito = is.mito), BPPARAM = MulticoreParam())
colData(sce) <- cbind(colData(sce), qc)

qc.lib <- qc$sum < 1000
qc.nexprs <- qc$detected < 500
qc.mito <- qc$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.mito

DataFrame(LibSize = sum(qc.lib), NExprs = sum(qc.nexprs), MitoProp = sum(qc.mito), Total = sum(discard))
sce$discard <- discard

# These plots allow us to visualize the data based on each metric and evaluate visually for batch effects
gridExtra::grid.arrange(
  plotColData(sce, x = "Mouse", y = "sum", colour_by = "discard", other_fields = "Strain") +
    scale_y_log10() + ggtitle("Total Count"),
  plotColData(sce, x = "Mouse", y = "detected", colour_by = "discard", other_fields = "Strain") +
    scale_y_log10() + ggtitle("Detected Features"),
  plotColData(sce, x = "Mouse", y = "subsets_Mito_percent", colour_by = "discard", other_fields = "Strain") 
    + ggtitle("Mitochondrial Percent"),
  ncol = 1
)

gridExtra::grid.arrange(
  plotColData(sce, x = "Batch", y = "sum", colour_by = "discard", other_fields = "Strain") +
    scale_y_log10() + ggtitle("Total Count"),
  plotColData(sce, x = "Batch", y = "detected", colour_by = "discard", other_fields = "Strain") +
    scale_y_log10() + ggtitle("Detected Features"),
  plotColData(sce, x = "Batch", y = "subsets_Mito_percent", colour_by = "discard", other_fields = "Strain") 
  + ggtitle("Mitochondrial Percent"),
  ncol = 1
)

plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "discard")
plot(qc$sum, qc$subsets_MT_percent, log="x", xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(is.mito, "thresholds")["higher"], col="red")

# This section of code allows us to evaluate whether we are inadvertently eliminating a cell type with our quality control
# It does so by evaluating whether the average 
lost <- calculateAverage(counts(sce)[,!sce$discard])
kept <- calculateAverage(counts(sce)[,sce$discard])
logged <- edgeR::cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)
plot(abundance, logFC, xlab = "Average Count", ylab = "Log-FC (lost/kept)", pch = 16)
points(abundance[is.mito], logFC[is.mito], col = "blue", pch = 16)

# Our system has limited memory; a little bit of memory management goes a long way
rm(lost)
rm(kept)
rm(abundance)
rm(logged)
rm(logFC)
rm(qc)
rm(location)
rm(is.mito)
rm(qc.lib)
rm(qc.mito)
rm(qc.nexprs)
rm(discard)

sce <- sce[,!sce$discard]

# Top Features and Dimensionality Reduction -------------------------------

# Selecting top differentially expressed features
set.seed(42)
clust.sce <- scran::quickCluster(sce, BPPARAM = MulticoreParam(), block.BPPARAM = MulticoreParam())
table(clust.sce)

set.seed(42)
sce <- scran::computeSumFactors(sce, cluster = clust.sce, min.mean = 0.1, BPPARAM = MulticoreParam())
sce <- logNormCounts(sce)

set.seed(42)
dec.log.sce <- modelGeneVar(sce, BPPARAM = MulticoreParam())
fit.log.sce <- metadata(dec.log.sce)
plot(fit.log.sce$mean, fit.log.sce$var, xlab = "Mean of Log-Expression", 
     ylab = "Variance of Log-Expression", title("Log-Count Feature Selection"))
curve(fit.log.sce$trend(x), col = "blue", add = TRUE, lwd = 2)
dec.log.sce <- dec.log.sce[order(dec.log.sce$bio, decreasing = TRUE),]
head(dec.log.sce)
rm(fit.log.sce)

hvg <- getTopHVGs(dec.log.sce, prop = 0.20)
rm(dec.log.sce)

set.seed(42)
sce <- runPCA(sce, subset_row = hvg)

set.seed(42)
sce <- runUMAP(sce, subset_row = hvg)

set.seed(42)
sce <- runTSNE(sce, subset_row = hvg, BPPARAM = MulticoreParam())

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab = "PC", ylab = "Variance Explained (%)")
rm(percent.var)

# Finally, we remove putative doublets (i.e. two cells which were in one droplet)
set.seed(42)
sce <- scDblFinder(sce, samples = "Batch", nfeatures = 750, propRandom = 1, verbose = TRUE, BPPARAM = MulticoreParam())
sce[,sce$scDblFinder.class == "doublet"]
sce <- sce[,sce$scDblFinder.class == "singlet"]

# Cell Type Annotation ----------------------------------------------------
colnames(sce) <- paste("V", as.character(1:ncol(counts(sce))), sep = "")

# Cell Marker Reference
cm.lthpsc <- c("Esam", "Cish", "Procr", "F11r", "Scarf1", "Mmrn1", "Tgm2", "Tinagl1", "Tie1", "Trpc6")
cm.hpsc <- c("Kit", "Cd34", "Ly6a")
cm.gmp <- c("Gfi1", "Elane", "Mpo")
cm.clp <- c("Igll1", "Rag1", "Rag2")
cm.mgk <- c("Pf4", "Mef2c", "Itga2b", "Pbx1")
cm.monocyte <- c("Csf1r", "F13a1", "Ly86")
cm.macrophage <- c("Adgre1", "Itgam", "Mrc1", "Arg1")
cm.neutrophil <- c("Ltf", "S100a8", "S100a9", "Lcn2", "Cebpe", "Ly6g")
cm.dendritic <- c("Itgax", "Siglech")
cm.eosinophil <- c("Prg2", "Il5ra")
cm.basophil <- c("Prss34", "Mcpt8")
cm.mastcell <- c("Cpa3", "Fcer1a")
cm.erythroid <- c("Car1", "Car2", "Hba-a1")
cm.tcell <- c("Cd2", "Cd5", "Cd3e")
cm.nkcell <- c("Klrb1c", "Gzma", "Ccl5")
cm.bcell <- c("Vpreb1", "Vpreb3", "Ebf1", "Igll1", "Pax5", "Cd79a", "Cd19", "Ms4a1")
cm.plasma <- c("Jchain")
cm.endothelial <- c("Cdh5", "Emcn", "Fabp4", "Flk1", "Ednrb")

cm.list <- list('Long-term Hematopoietic Stem Cell' = cm.lthpsc,
                'Hematopoietic Stem Cell' = cm.hpsc,
                'Granulocyte-Monocyte Precursor' = cm.gmp,
                'Common Lymphoid Progenitor' = cm.clp,
                'Megakaryocyte' = cm.mgk,
                'Monocyte' = cm.monocyte,
                'Macrophage' = cm.macrophage,
                'Neutrophil' = cm.neutrophil,
                'Dendritic Cell' = cm.dendritic,
                'Eosinophil' = cm.eosinophil,
                'Basophil' = cm.basophil,
                'Mast Cell' = cm.mastcell,
                'Erythroid Cell' = cm.erythroid,
                'NK Cell' = cm.nkcell,
                'T-Cell' = cm.tcell,
                'B-Cell' = cm.bcell,
                'Plasma Cell' = cm.plasma)

cell.rankings <- AUCell_buildRankings(counts(sce), plotStats = TRUE)
cell.AUC <- AUCell_calcAUC(cm.list, cell.rankings, verbose = TRUE, aucMaxRank = 750)
cell.assignments <- AUCell_exploreThresholds(cell.AUC, plotHist = TRUE, assignCells = TRUE)

cell.assignments$`B-Cell`$aucThr$selected <- 0.20
cell.assignments$`B-Cell`$assignment <- names(which(getAUC(cell.AUC)['B-Cell',] > 0.20))

cell.assignments$`Erythroid Cell`$aucThr$selected <- 0.40
cell.assignments$`Erythroid Cell`$assignment <- names(which(getAUC(cell.AUC)['Erythroid Cell',] > 0.40))

cell.assignments$`Eosinophil`$aucThr$selected <- 0.20
cell.assignments$`Eosinophil`$assignment <- names(which(getAUC(cell.AUC)['Eosinophil',] > 0.20))

cell.assignments$`Neutrophil`$aucThr$selected <- 0.52
cell.assignments$`Neutrophil`$assignment <- names(which(getAUC(cell.AUC)['Neutrophil',] > 0.52))

cell.assignments$`Monocyte`$aucThr$selected <- 0.10
cell.assignments$`Monocyte`$assignment <- names(which(getAUC(cell.AUC)['Monocyte',] > 0.10))

cell.assignments$`Macrophage`$aucThr$selected <- 0.05
cell.assignments$`Macrophage`$assignment <- names(which(getAUC(cell.AUC)['Macrophage',] > 0.05))

assignments <- matrix(nrow = length(colnames(sce)), ncol = 1) 
rownames(assignments) <- colnames(sce)
for(i in names(cell.assignments)){
  cells <- cell.assignments[[i]]$assignment
  assignments[cells,] <- i
}
assignments <- replace_na(assignments, "Unassigned")
sce$cm.type <- as.character(assignments)
plotReducedDim(sce, "UMAP", colour_by = "cm.type")

rm(cell.assignments)
rm(cell.rankings)
rm(cells)
rm(assignments)
rm(cm.list)

# Differential Analysis ---------------------------------------------------

# Differential Expression
sce.aggregate <- aggregateAcrossCells(sce, ids = colData(sce)[,c("cm.type", "Mouse")])
sce.aggregate <- sce.aggregate[,sce.aggregate$ncells >= 10]
de.results <- pseudoBulkDGE(sce.aggregate, label = sce.aggregate$cm.type, design = ~0 + Strain,
                            condition = sce.aggregate$Strain, contrast = c(-0.5, 1, -0.5))
is.de <- decideTestsPerLabel(de.results, threshold = 0.05)
summarizeTestsPerLabel(is.de)

up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing = TRUE), 10)

down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing = TRUE), 10)

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/nzbwf1_de/")
for(cell.type in names(de.results)){
    x <- as.data.frame(de.results[[cell.type]])
    x <- x[order(x$FDR),]
    x <- subset(x, x$FDR < 0.05)
    x <- cbind(gene = rownames(x), x)
    write.table(x, file = paste(as.character(cell.type), " - NZBWF1 Differential Gene Expression.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  }
setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/nzbwf1_de/")
for(cell.type in names(de.results)){
  x <- as.data.frame(de.results[[cell.type]])
  graph <- EnhancedVolcano(x, lab = rownames(x), x = 'logFC', y = 'PValue', FCcutoff = 0.5, pCutoff = 0.0001,
                  xlab = "Log Fold Change", title = cell.type, subtitle = "NZBWF1 vs. NZB/NZW")
  tiff(filename = paste(cell.type, " - Volcano Plot.tiff"), width = 850, height = 1000)
  print(graph)
  dev.off()
  ggplot(x) + geom_point(aes(x = logFC, y= -log10(PValue), color = threshold), size = 3) +
    ggtitle(cell.type) + xlab("Log Fold Change") + ylab("-log10 PValue") +
    scale_color_manual(values = c("#A9A9A9", "#009051")) +
    geom_vline(xintercept = 0.5, col = "black", linetype = "dotted", size = 1) +
    geom_vline(xintercept = -0.5, col = "black", linetype = "dotted", size = 1) +
    geom_hline(yintercept = 4, col = "black", linetype = "dotted", size = 1) +
    theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5), axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 15), legend.position = "none") +
    geom_label_repel(aes(x = logFC, y = -log10(PValue), label = threshold))
}
setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

go.list <- list()
kegg.list <- list()
for(cell.type in names(de.results)){
  genelist <- as.data.frame(de.results[[cell.type]])
  genelist <- cbind(genelist, gene_id = rowData(sce)$gene_id)
  significant <- subset(genelist, genelist$FDR < 0.05)
  go.list[[cell.type]] <- goana(significant$gene_id, species = "Mm")
  go.list[[cell.type]] <- go.list[[cell.type]][order(go.list[[cell.type]]$P.DE),]
  kegg.list[[cell.type]] <- kegga(significant$gene_id, species.KEGG = "mmu")
  kegg.list[[cell.type]] <- kegg.list[[cell.type]][order(kegg.list[[cell.type]]$P.DE),]
}

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/gene_ontology/")
for(cell.type in names(go.list)){
  x <- go.list[[cell.type]]
  print(range(x$P.DE))
  write.table(x, file = paste(as.character(cell.type), " - NZBWF1 Differential Gene Ontology.tsv"), quote = FALSE, sep = "\t")
  x <- kegg.list[[cell.type]]
  print(range(x$P.DE))
  write.table(x, file = paste(as.character(cell.type), " - NZBWF1 Differential KEGG.tsv"), quote = FALSE, sep = "\t")
}
setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

# Differential Abundance
abundances <- table(sce$cm.type, sce$Mouse) 
abundances <- unclass(abundances)  
extra.info <- colData(sce)[match(colnames(abundances), sce$Mouse), 2:4]
y.ab <- DGEList(abundances, samples = extra.info)
keep <- filterByExpr(y.ab, group = y.ab$samples$Strain)
y.ab <- y.ab[keep,]
summary(keep)

da.design <- model.matrix(~0 + factor(Strain), y.ab$samples)
y.ab <- estimateDisp(y.ab, da.design, trend = "none")
fit.ab <- glmQLFit(y.ab, da.design, robust = TRUE, abundance.trend = FALSE)
res <- glmQLFTest(fit.ab, coef = ncol(da.design), contrast = c(-0.5, 1, -0.5))
summary(decideTests(res))
topTags(res)

res.table <- as.data.frame(topTags(res))
threshold <- res.table$FDR < 0.05
ggplot(res.table) + geom_point(aes(x = logFC, y = -FDR, color = threshold))

rm(keep)
rm(da.design)
rm(y.ab)
rm(fit.ab)
rm(abundances)
rm(extra.info)

# Cell Cycling and Trajectory Analysis ------------------------------------

#Cell Cycling
set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names = rowData(sce)$gene_id, BPPARAM = MulticoreParam(workers = 10))
colData(sce)$phase <- assignments$phases

rm(mm.pairs)
rm(assignments)

#Trajectory Analysis
expression_matrix <- counts(sce)
colnames(expression_matrix) <- 1:ncol(counts(sce))
cell_metadata <- colData(sce)
rownames(cell_metadata) <- 1:nrow(colData(sce))
gene_metadata <- rowData(sce)
colnames(gene_metadata) <- c("gene_id", "gene_short_name", "gene_biotype")

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
rm(expression_matrix)
rm(cell_metadata)
rm(gene_metadata)
reducedDim(cds, "PCA") <- reducedDim(sce, "PCA")

set.seed(42)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", cores = 22, verbose = FALSE)
cds <- cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
cds <- learn_graph(cds, verbose = FALSE)

cds_3d <- reduce_dimension(cds, max_components = 3, cores = 22, verbose = FALSE, reduction_method = "UMAP")
cds_3d <- reduce_dimension(cds_3d, max_components = 3, cores = 22, verbose = FALSE, reduction_method = "tSNE")
cds_3d <- cluster_cells(cds_3d, reduction_method = "UMAP", verbose = FALSE)
cds_3d <- learn_graph(cds_3d, verbose = FALSE)

# B Cell Subclustering and Analysis ----------------------------------------------------

sce.bcell <- sce[,sce$cm.type == "B-Cell"]

set.seed(42)
sce.bcell <- runPCA(sce.bcell)

set.seed(42)
sce.bcell <- runUMAP(sce.bcell)

set.seed(42)
sce.bcell <- runTSNE(sce.bcell, BPPARAM = MulticoreParam())

pcs <- reducedDim(sce.bcell, "PCA")
choices <- getClusteredPCs(pcs, BPPARAM = MulticoreParam())
metadata(choices)$chosen
plot(choices$n.pcs, choices$n.clusters, xlab = "Number of PCs", ylab = "Number of clusters")
abline(a = 1, b = 1, col = "red")
abline(v = metadata(choices)$chosen, col = "blue", lty = 2)

reducedDim(sce.bcell, "PCA_Reduced") <- reducedDim(sce.bcell, "PCA")[,1:metadata(choices)$chosen]

g <- buildSNNGraph(sce.bcell, use.dimred = "PCA_Reduced")
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
colLabels(sce.bcell) <- factor(clust)
plotReducedDim(sce.bcell, "TSNE", colour_by = "label")

ratio <- clusterModularity(g, clust, as.ratio=TRUE)
dim(ratio)

pheatmap::pheatmap(log2(ratio + 1), cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100))

cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio + 1), mode = "upper", weighted = TRUE, diag = FALSE)

set.seed(42)
plot(cluster.gr, edge.width = igraph::E(cluster.gr)$weight * 5, layout = igraph::layout_with_lgl)

ClusterFUN <- function(x) {
  graph <- buildSNNGraph(x, use.dimred = "PCA_Reduced")
  igraph::cluster_walktrap(graph)$membership
}

originals <- ClusterFUN(sce.bcell)

set.seed(42)
coassign <- bluster::bootstrapStability(sce.bcell, FUN = ClusterFUN, clusters = originals)
pheatmap(coassign, cluster_row = FALSE, cluster_col = FALSE, color = rev(viridis::plasma(100)))
markers.wmw <- findMarkers(sce.bcell, test = "wilcox", direction = "up", lfc = 1, pval.type = "some")

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/bcell_subclusters/")
for(i in 1:length(markers.wmw)){
  x <- markers.wmw[[i]]
  x <- subset(x, x$p.value < 0.05)
  names <- rownames(x)
  rownames(x) <- NULL
  x <- cbind(names, x)
  write.table(x, file = paste("Cluster ", as.character(i), ": Wilcoxon-Mann-Whitney Determined Markers, pval.type = some.tsv"),
              quote = FALSE, row.names = FALSE, sep = "\t")
}
setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

abundances <- table(sce.bcell$label, sce.bcell$Mouse) 
abundances <- unclass(abundances)  
extra.info <- colData(sce.bcell)[match(colnames(abundances), sce.bcell$Mouse), 2:4]
y.ab <- DGEList(abundances, samples = extra.info)
keep <- filterByExpr(y.ab, group = y.ab$samples$Strain)
y.ab <- y.ab[keep,]
summary(keep)

da.design <- model.matrix(~0 + factor(Strain), y.ab$samples)
y.ab <- estimateDisp(y.ab, da.design, trend = "none")
fit.ab <- glmQLFit(y.ab, da.design, robust = TRUE, abundance.trend = FALSE)
res <- glmQLFTest(fit.ab, coef = ncol(da.design), contrast = c(-0.5, 1, -0.5))
summary(decideTests(res))
topTags(res)

rm(keep)
rm(da.design)
rm(y.ab)
rm(fit.ab)
rm(abundances)
rm(extra.info)

sce.bcell.aggregate <- aggregateAcrossCells(sce.bcell, ids = colData(sce.bcell)[,c("label", "Mouse")])
de.bcell <- pseudoBulkDGE(sce.bcell.aggregate, label = sce.bcell.aggregate$label, design = ~0 + Strain, condition = sce.bcell.aggregate$Strain, contrast = c(-0.5, 1, -0.5))
is.de.bcell <- decideTestsPerLabel(de.bcell, threshold = 0.05)
summarizeTestsPerLabel(is.de.bcell)
up.de.bcell <- is.de.bcell > 0 & !is.na(is.de.bcell)
head(sort(rowMeans(up.de.bcell), decreasing = TRUE), 10)
down.de.bcell <- is.de.bcell < 0 & !is.na(is.de.bcell)
head(sort(rowMeans(down.de.bcell), decreasing = TRUE), 10)

setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/analysis/bcell_de/")
for(cell.type in names(de.bcell)){
  x <- as.data.frame(de.results[[cell.type]])
  x <- x[order(x$FDR),]
  x <- subset(x, x$FDR < 0.05)
  write.table(x, file = paste(as.character(cell.type), " - NZBWF1 Differential Gene Expression.tsv"), quote = FALSE, sep = "\t")
}

for(cell.type in names(de.bcell)){
  x <- as.data.frame(de.bcell[[cell.type]])
  graph <- EnhancedVolcano(x, lab = rownames(x), x = 'logFC', y = 'PValue', FCcutoff = 0.5, pCutoff = 0.0001,
                           xlab = "Log Fold Change", title = paste('B-Cell Subcluster #', cell.type), subtitle = "NZBWF1 vs. NZB/NZW")
  tiff(filename = paste(cell.type, " - Volcano Plot.tiff"), width = 850, height = 1000)
  print(graph)
  dev.off()
}
setwd("~/MEGA/Bioinformatics/SLE_scRNAseq/data")

expression_matrix <- counts(sce.bcell)
colnames(expression_matrix) <- 1:ncol(counts(sce.bcell))
cell_metadata <- colData(sce.bcell)
rownames(cell_metadata) <- 1:nrow(colData(sce.bcell))
gene_metadata <- rowData(sce.bcell)
colnames(gene_metadata) <- c("gene_id", "gene_short_name", "gene_biotype")

cds.bcell <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
rm(expression_matrix)
rm(cell_metadata)
rm(gene_metadata)
reducedDim(cds.bcell, "PCA") <- reducedDim(sce.bcell, "PCA")

set.seed(42)
cds.bcell <- reduce_dimension(cds.bcell, reduction_method = "UMAP", preprocess_method = "PCA", cores = 22, verbose = FALSE)
