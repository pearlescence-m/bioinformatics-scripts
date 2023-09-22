# RNA-seq Analysis

> Exercises and solutions for Chapter 8 of <https://compgenomr.github.io/book/>

Here we will use a subset of the RNA-seq count table from a colorectal cancer study. We have filtered the original count table for only protein-coding genes (to improve the speed of calculation) and also selected only five metastasized colorectal cancer samples along with five normal colon samples. There is an additional column width that contains the length of the corresponding gene in the unit of base pairs. The length of the genes are important to compute RPKM and TPM values. The original count tables can be found from the [recount2 database](https://jhubiostatistics.shinyapps.io/recount/) using the SRA project code SRP029880, and the experimental setup along with other accessory information can be found from the NCBI Trace archive using the SRA project code [SRP029880](https://trace.ncbi.nlm.nih.gov/Traces/index.html?study=SRP029880).

### 1. Exploring the count tables

Here, import an example count table and do some exploration of the expression data.

``` {.r echo="FALSE," eval="FALSE"}
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")
```

1.  Normalize the counts using the TPM approach.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))

geneLengths <- as.vector(subset(counts, select = c(width)))

#find gene length normalized values 
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))

#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
colSums(tpm)
```

2.  Plot a heatmap of the top 500 most variable genes. Compare with the heatmap obtained using the 100 most variable genes.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#compute the variance of each gene across samples
V <- apply(tpm, 1, var)

#and select the top genes 
selectedGenes500 <- names(V[order(V, decreasing = T)][1:500])
selectedGenes100 <- names(V[order(V, decreasing = T)][1:100])

library(pheatmap)
pheatmap(tpm[selectedGenes500,], scale = 'row', show_rownames = FALSE)
pheatmap(tpm[selectedGenes100,], scale = 'row', show_rownames = FALSE)
```

3.  Re-do the heatmaps setting the `scale` argument to `none`, and `column`. Compare the results with `scale = 'row'`.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
pheatmap(tpm[selectedGenes500,], scale = 'none', show_rownames = FALSE)
pheatmap(tpm[selectedGenes500,], scale = 'column', show_rownames = FALSE)
pheatmap(tpm[selectedGenes100,], scale = 'column', show_rownames = FALSE)
pheatmap(tpm[selectedGenes100,], scale = 'none', show_rownames = FALSE)
```

4.  Draw a correlation plot for the samples depicting the sample differences as 'ellipses', drawing only the upper end of the matrix, and order samples by hierarchical clustering results based on `average` linkage clustering method.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(stats)
correlationMatrix <- cor(tpm)

library(corrplot)
corrplot(correlationMatrix, method='ellipse', type='upper', order = 'hclust', 
         hclust.method = c("average")) 
```

5.  How else could the count matrix be subsetted to obtain quick and accurate clusters? Try selecting the top 100 genes that have the highest total expression in all samples and re-draw the cluster heatmaps and PCA plots.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#total expression in top 100 genes
TE <- rowSums(tpm)
top_100_genes <- names(TE[order(TE, decreasing = TRUE)][1:100])

#re-draw heatmaps
colData <- read.table(coldata_file, header = T, sep = '\t', stringsAsFactors = TRUE)

pheatmap(tpm[top_100_genes,], scale = 'row', show_rownames = FALSE, annotation_col = colData)

#PCA plots
library(stats)
library(ggplot2)
library(ggfortify)

#transpose the matrix 
M <- t(tpm[top_100_genes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

autoplot(pcaResults, data = colData, colour = 'group')
summary(pcaResults)
```

6.  Add an additional column to the annotation data.frame object to annotate the samples and use the updated annotation data.frame to plot the heatmaps. (Hint: Assign different batch values to CASE and CTRL samples). Make a PCA plot and color samples by the added variable (e.g. batch).

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#extract batches
colData$batch <- as.numeric(sub("\\D+", "", rownames(colData)))

M <- t(tpm[top_100_genes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

autoplot(pcaResults, data = colData, colour = 'batch')
summary(pcaResults)
```

7.  Try making the heatmaps using all the genes in the count table, rather than sub-selecting.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}

pheatmap(tpm, scale = 'row', show_rownames = FALSE, annotation_col = colData)
```

8.  Use the [`Rtsne` package](https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf) to draw a t-SNE plot of the expression values. Color the points by sample group. Compare the results with the PCA plots.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(Rtsne)

tpm_unique <- unique(tpm) # Remove duplicates
tpm_matrix <- as.matrix(tpm_unique)
X <- normalize_input(tpm_matrix)
colMeans(X)
range(X)

tpm_df <- as.data.frame(tpm)
tpm_df$geneName <- rownames(tpm_df)

# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(tpm_matrix,pca=FALSE,theta=0.0) # Run TSNE

plot(tsne_out$Y,col=tpm_df$geneName, asp=1)
```

### 2. Differential expression analysis

Firstly, carry out a differential expression analysis starting from raw counts. Use the following datasets:

``` {.r echo="FALSE,eval=FALSE"}
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                         package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                          package = "compGenomRData")
```

-   Import the read counts and colData tables.
-   Set up a DESeqDataSet object.
-   Filter out genes with low counts.
-   Run DESeq2 contrasting the `CASE` sample with `CONTROL` samples.

``` {.r echo="FALSE,eval=FALSE"}
#remove the 'width' column
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
countData <- as.matrix(subset(counts, select = c(-width)))
#define the experimental setup 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
#define the design formula
designFormula <- "~ group"

library(DESeq2)
library(stats)
#create a DESeq dataset object from the count matrix and the colData 
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))
#print dds object to see the contents
print(dds)
```

Now, you are ready to do the following exercises:

1.  Make a volcano plot using the differential expression analysis results. (Hint: x-axis denotes the log2FoldChange and the y-axis represents the -log10(pvalue)).

**solution:**

``` {.r echo="FALSE,eval=FALSE"}

#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

DEresults$log10_pvalue <- -log10(DEresults$pvalue)

# Create a volcano plot
library(ggplot2)

# Customize the plot
ggplot(as.data.frame(DEresults), aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 1 & log10_pvalue > 1.3, "red", "black")), alpha = 0.6, size = 2) +
  scale_color_identity() +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none") 
```

2.  Use DESeq2::plotDispEsts to make a dispersion plot and find out the meaning of this plot. (Hint: Type ?DESeq2::plotDispEsts)

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
plotDispEsts(dds)
```

3.  Explore `lfcThreshold` argument of the `DESeq2::results` function. What is its default value? What does it mean to change the default value to, for instance, `1`?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
# Extract results with a log2 fold change threshold of 1
result_1 <- results(dds, lfcThreshold = 1)

# Extract results with the default log2 fold change threshold (0)
result_default <- results(dds)
```

In the example above, result_1 will contain only genes with an absolute log2 fold change of 1 or higher, while result_default will include all genes regardless of their log2 fold change.

Changing the lfcThreshold allows to focus on genes that exhibit larger changes in expression, which can be useful when you are primarily interested in highly differentially expressed genes and want to filter out genes with smaller changes that may not be biologically significant for analysis.

4.  What is independent filtering? What happens if we don't use it? Google `independent filtering statquest` and watch the online video about independent filtering.

Independent filtering is a statistical technique used to improve the accuracy of results and reduce the number of false positives while still maintaining reasonable statistical power. Here's a simplified overview of how independent filtering works:

```         
Calculate p-values for each gene based on statistical tests (e.g., negative binomial tests for RNA-seq data).

Rank genes by p-value.

Consider additional information, such as the fold change and gene expression levels. Genes with very low expression levels or very small fold changes may be less likely to be biologically meaningful.

Apply a more relaxed p-value threshold (higher than the initial threshold) to genes that meet certain criteria (e.g., larger fold change or higher expression levels).

Retain the genes that pass the adjusted threshold as differentially expressed.
```

5.  Re-do the differential expression analysis using the `edgeR` package. Find out how much DESeq2 and edgeR agree on the list of differentially expressed genes.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(edgeR)

# Create a DGEList object
dge <- DGEList(counts = countData, group = colData$group)

# Perform TMM normalization
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~group, data = colData)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit a generalized linear model
fit <- glmQLFit(dge, design)

# Perform likelihood ratio tests for differential expression
results_edgeR <- glmQLFTest(fit, coef = 2)

# Extract the differentially expressed genes
DE_genes_edgeR <- rownames(topTags(results_edgeR, n = Inf))
DE_genes_DESeq2 <- rownames(res)

# Find the overlapping genes between DESeq2 and edgeR results
common_genes <- intersect(DE_genes_DESeq2, DE_genes_edgeR)

# Calculate the level of agreement
agreement <- length(common_genes) / length(unique(c(DE_genes_DESeq2, DE_genes_edgeR)))

# Print the number of common genes and agreement percentage
cat("Number of common differentially expressed genes:", length(common_genes), "\n")
# Number of common differentially expressed genes: 19719 
cat("Percentage agreement between DESeq2 and edgeR:", agreement * 100, "%\n")
# Percentage agreement between DESeq2 and edgeR: 100 %
 
```

6.  Use the `compcodeR` package to run the differential expression analysis using at least three different tools and compare and contrast the results following the `compcodeR` vignette.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(compcodeR)
colnames(colData)<-c("condition", "group")
comp_data <- compData(count.matrix=countData, sample.annotations=colData, 
  info.parameters = list(dataset = "mydata", uID = "SRP029880"))

tmpdir <- normalizePath(tempdir(), winslash = "/")
mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, samples.per.cond = 5,   n.diffexp = 100, output.file = file.path(tmpdir, "mydata.rds"))

runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.exact",
Rmdfunction = "edgeR.exact.createRmd",
output.directory = tmpdir, norm.method = "TMM",
trend.method = "movingave", disp.type = "tagwise")

runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2",
Rmdfunction = "DESeq2.createRmd",
output.directory = tmpdir, fit.type = "parametric",
test = "Wald", beta.prior = TRUE,
independent.filtering = TRUE, cooks.cutoff = TRUE,
impute.outliers = TRUE)

runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma",
Rmdfunction = "voom.limma.createRmd", output.directory = tmpdir,
norm.method = "TMM")

file.table <- data.frame(input.files = file.path(tmpdir, c("mydata_voom.limma.rds", "mydata_edgeR.exact.rds", "mydata_DESeq2.rds")), stringsAsFactors = FALSE)

parameters <- list(incl.nbr.samples = 5, incl.replicates = 1, incl.dataset = "mydata",
incl.de.methods = NULL,
fdr.threshold = 0.05, tpr.threshold = 0.05, typeI.threshold = 0.05,
ma.threshold = 0.05, fdc.maxvar = 1500, overlap.threshold = 0.05,
fracsign.threshold = 0.05, mcc.threshold = 0.05,
nbrtpfp.threshold = 0.05,
comparisons = c("auc", "fdr", "tpr", "ma", "correlation"))

runComparison(file.table = file.table, parameters = parameters, output.directory = tmpdir)
 
```

### 3. Functional enrichment analysis

1.  Re-run gProfileR, this time using pathway annotations such as KEGG, REACTOME, and protein complex databases such as CORUM, in addition to the GO terms. Sort the resulting tables by columns `precision` and/or `recall`. How do the top GO terms change when sorted for `precision`, `recall`, or `p.value`?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(DESeq2)
library(gProfileR)
library(knitr)
library(stats)
#create a DESeq dataset object from the count matrix and the colData 
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))
#print dds object to see the contents
print(dds)

dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
dds <- DESeq(dds)

# Extract differential expression results
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))

# Remove genes with NA values
DE <- DEresults[!is.na(DEresults$padj),]

# Select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]

# Select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

# Get the list of genes of interest
genesOfInterest <- rownames(DE)

# Calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest, 
                       organism = 'hsapiens', 
                       src_filter = c('GO', 'KEGG', 'REACTOME', 'CORUM'), 
                       hier_filtering = 'moderate')

# Order the results by precision, recall, and p.value
goResults_by_precision <- goResults[order(goResults$precision),]
goResults_by_recall <- goResults[order(goResults$recall),]
goResults_by_pvalue <- goResults[order(goResults$p.value),]
 
```

2.  Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage` function. How do the results differ?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#gprofiler doesn't work - coming soon
```

3.  Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#gprofiler doesn't work - coming soon
 
```

4.  Do a gene-set enrichment analysis using gene sets from top 10 GO terms.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```

5.  What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```

6.  Use the topGO package (<https://bioconductor.org/packages/release/bioc/html/topGO.html>) to re-do the GO term analysis. Compare and contrast the results with what has been obtained using the `gProfileR` package. Which tool is faster, `gProfileR` or topGO? Why?

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```

7.  Given a gene set annotated for human, how can it be utilized to work on *C. elegans* data? (Hint: See `biomaRt::getLDS`).

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```

8.  Import curated pathway gene sets with Entrez identifiers from the [MSIGDB database](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) and re-do the GSEA for all curated gene sets.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```

### 4. Removing unwanted variation from the expression data

For the exercises below, use the datasets at:

```         
counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                           package = 'compGenomRData')
```

1.  Run RUVSeq using multiple values of `k` from 1 to 10 and compare and contrast the PCA plots obtained from the normalized counts of each RUVSeq run.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(EDASeq)
library(RUVSeq)
library(ggplot2)

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, 
                      sep = '\t', stringsAsFactors = TRUE)
# simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')
                              
# remove 'width' column from counts
countData <- as.matrix(subset(counts, select = c(-width)))

# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)
                           
#source for house-keeping genes collection:
#https://m.tau.ac.il/~elieis/HKG/HK_genes.txt
HK_genes <- read.table(file = system.file("extdata/rna-seq/HK_genes.txt", 
                                          package = 'compGenomRData'), 
                       header = FALSE)
# let's take an intersection of the house-keeping genes with the genes available
# in the count table
house_keeping_genes <- intersect(rownames(set), HK_genes$V1)

# now, we use these genes as the empirical set of genes as input to RUVg.
# we try different values of k and see how the PCA plots look 
 
par(mfrow = c(3, 3))
for(k in 1:10) {
  set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = k)
  plotPCA(set_g, col=as.numeric(colData$group), cex = 0.9, adj = 0.5, 
          main = paste0('with RUVg, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1), )
}
 
```

2.  Re-run RUVSeq using the `RUVr()` function. Compare PCA plots from `RUVs`, `RUVg` and `RUVr` using the same `k` values and find out which one performs the best.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
library(edgeR)

#RUVs

differences <- makeGroups(colData$group)
par(mfrow = c(3, 3))
for(k in 1:10) {
  set_s <- RUVs(set, unique(rownames(set)), k=k, differences)
  plotPCA(set_s, col=as.numeric(colData$group), 
          cex = 0.9, adj = 0.5, 
        main = paste0('with RUVs, k = ',k), 
        ylim = c(-1, 1), xlim = c(-0.6, 0.6))
}

#RUVr

# Residuals from negative binomial GLM regression of UQ-normalized
# counts on covariates of interest, with edgeR
# Create a DGEList object
dge <- DGEList(counts = countData, group = colData$group)

# Perform TMM normalization
dge <- calcNormFactors(dge, method="upperquartile")

# Design matrix
design <- model.matrix(~group, data = colData)

y <- estimateGLMCommonDisp(dge, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# RUVr normalization (after UQ)
seqUQ <- betweenLaneNormalization(set, which="upper")
controls <- rownames(set)

par(mfrow = c(3, 3))
for(k in 1:10) {
  set_r <- RUVr(seqUQ, controls, k=k, res)
  plotPCA(set_r, col=as.numeric(colData$group), cex = 0.9, adj = 0.5, 
          main = paste0('with RUVr, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1), )
}
```

3.  Do the necessary diagnostic plots using the differential expression results from the EHF count table.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
# RLE plots

par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'without RUVg')
plotRLE(set_g, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'with RUVg')
        
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'without RUVs')
plotRLE(set_s, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'with RUVs')
        
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'without RUVr')
plotRLE(set_r, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'with RUVr')
        
# RLE plots

par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVg', 
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
plotPCA(set_g, col=as.numeric(colData$group), adj = 0.5, 
        main = 'with RUVg',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
        
par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVs', 
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
plotPCA(set_s, col=as.numeric(colData$group), adj = 0.5, 
        main = 'with RUVs',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
        
par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVr', 
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
plotPCA(set_r, col=as.numeric(colData$group), adj = 0.5, 
        main = 'with RUVr',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
 
```

4.  Use the `sva` package to discover sources of unwanted variation and re-do the differential expression analysis using variables from the output of `sva` and compare the results with `DESeq2` results using `RUVSeq` corrected normalization counts.

**solution:**

``` {.r echo="FALSE,eval=FALSE"}
#coming soon
 
```
