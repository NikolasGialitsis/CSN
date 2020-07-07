# buettner sc rna seq differential analysis

library(scRNAseq)
sce=BuettnerESCData()
#mouse embryonic stem cells

# Quality control.
library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection.
#BiocManager::install("scran")
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
library(tables)
sce$label <- factor(igraph::cluster_louvain(g)$membership)

# Visualization.
plotUMAP(sce, colour_by="label")



### Differential expression analysis using DEsingle
#BiocManager::install("DEsingle")
library(DEsingle)
library(SingleCellExperiment)
#gene expression matrix
gem=counts(sce)
subset2= gem[sample(nrow(gem), 2000), sample(ncol(gem), 100)]
subset= gem[sample(nrow(gem), 1000), sample(ncol(gem), 100)]
group <- factor(c(rep(1,60), rep(2,40)))
results_subset2 <- DEsingle(counts = subset2, group = group)
results <- DEsingle(counts = subset, group = group)

# Dividing the DE genes into 3 categories at threshold of FDR < 0.05
results.classified <- DEtype(results =results, threshold = 0.05)
# Extract DE genes at threshold of FDR < 0.05
results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]
# Extract three types of DE genes separately
results.DEs <- results.sig[results.sig$Type == "DEs", ]
results.DEa <- results.sig[results.sig$Type == "DEa", ]
results.DEg <- results.sig[results.sig$Type == "DEg", ]


###DE using limma
library(edgeR)
library(limma)
library(Glimma)#
library(gplots)#
library(org.Mm.eg.db)#
library(RColorBrewer)


myCPM <- cpm(subset2)
thresh <- myCPM > 0.7
keep <- rowSums(thresh) >= 70
counts.keep <- subset2[keep,]
rownames(counts.keep) <- substr(rownames(counts.keep),start=14,stop=18)
colnames(counts.keep) <- substr(colnames(counts.keep),start=1,stop=9)

names=rownames(counts.keep)
write.csv(names, "DEgenes.txt")
heatmap(counts.keep)
y <- DGEList(counts.keep)
y$samples$lib.size
barplot(y$samples$lib.size,names=colnames(y),las=2)
logcounts <- cpm(y,log=TRUE)
heatmap(logcounts)

