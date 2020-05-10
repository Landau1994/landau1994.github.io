---
title: LearnSeurat_PBMC3k
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - R
  - scRNA-seq
description: LearnSeurat_1
date: 2020-05-10 19:55:56
keywords:
photos:
---


本系列假定读者对于单细胞测序的数据分析和Seurat的官方教程有所了解。

本篇研究最基础的PBMC3k。其实这里只有2700个外周血的细胞。注意到，由于取样是外周血，没有干细胞的存在，所以可以认为样品处于稳态。这个教程就是讲稳态下的单细胞测序分析是如何进行的。Seurat的官方教程的缺点之一就是没有涉及动态过程的单细胞分析如何进行。

如无特殊说明，本系列的代码均可以在自己的笔记本电脑上运行；

### 1. 构建Seurat object

使用作者已经构建好的数据进行构建。关于Seurat更详细的文档可见[satijalab的wiki](https://github.com/satijalab/seurat/wiki)

``` r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(igraph)
library(tidyverse)
library(patchwork)
### AvailableData() check avaliable data: we choose cbmc
### InstallData('pbmc3k')
library(pbmc3k.SeuratData)

### how this dataset generate?
# ## Not run: 
# if (requireNamespace(Seurat, quietly = TRUE)) {
#   url <- 'http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
#   curl::curl_download(url = url, destfile = basename(path = url))
#   untar(tarfile = basename(path = url))
#   pbmc.data <- Seurat::Read10X(data.dir = 'filtered_gene_bc_matrices/hg19/')
#   pbmc3k <- Seurat::CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)
#   # Annotations come from Seurat's PBMC3k Guided Clustering Tutorial
#   # https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
#   annotations <- readRDS(file = system.file('extdata/annotations/annotations.Rds', package = 'pbmc3k.SeuratData'))
#   pbmc3k <- Seurat::AddMetaData(object = pbmc3k, metadata = annotations)
#   # Clean up downloaded files
#   file.remove(basename(path = url))
#   unlink(x = 'filtered_gene_bc_matrices/', recursive = TRUE)
# }
# 
# ## End(Not run)

### how Create Seurat object work?
### by run `Seurat::CreateSeuratObject` you can get the source function
# ## Not run:
# Seurat::CreateSeuratObject
# function (counts, project = "SeuratProject", assay = "RNA", 
#     min.cells = 0, min.features = 0, names.field = 1, names.delim = "_", 
#     meta.data = NULL) 
# {
#     if (!is.null(x = meta.data)) {
#         if (is.null(x = rownames(x = meta.data))) {
#             stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
#         }
#         if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = counts)))) {
#             warning("Some cells in meta.data not present in provided counts matrix.")
#             meta.data <- meta.data[intersect(x = rownames(x = meta.data), 
#                 y = colnames(x = counts)), ]
#         }
#         if (is.data.frame(x = meta.data)) {
#             new.meta.data <- data.frame(row.names = colnames(x = counts))
#             for (ii in 1:ncol(x = meta.data)) {
#                 new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, 
#                   ii, drop = FALSE]
#             }
#             meta.data <- new.meta.data
#         }
#     }
#     assay.data <- CreateAssayObject(counts = counts, min.cells = min.cells, 
#         min.features = min.features)
#     Key(object = assay.data) <- paste0(tolower(x = assay), "_")
#     assay.list <- list(assay.data)
#     names(x = assay.list) <- assay
#     init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
#     idents <- factor(x = unlist(x = lapply(X = colnames(x = assay.data), 
#         FUN = ExtractField, field = names.field, delim = names.delim)))
#     if (any(is.na(x = idents))) {
#         warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name")
#     }
#     ident.levels <- length(x = unique(x = idents))
#     if (ident.levels > 100 || ident.levels == 0 || ident.levels == 
#         length(x = idents)) {
#         idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
#     }
#     names(x = idents) <- colnames(x = assay.data)
#     object <- new(Class = "Seurat", assays = assay.list, 
#         meta.data = init.meta.data, active.assay = assay, active.ident = idents, 
#         project.name = project, version = packageVersion(pkg = "Seurat"))
#     object[["orig.ident"]] <- idents
#     n.calc <- CalcN(object = assay.data)
#     if (!is.null(x = n.calc)) {
#         names(x = n.calc) <- paste(names(x = n.calc), assay, 
#             sep = "_")
#         object[[names(x = n.calc)]] <- n.calc
#     }
#     if (!is.null(x = meta.data)) {
#         object <- AddMetaData(object = object, metadata = meta.data)
#     }
#     return(object)
# }
# 
# Seurat:: CreateAssayObject
# function (counts, data, min.cells = 0, min.features = 0) 
# {
#     if (missing(x = counts) && missing(x = data)) {
#         stop("Must provide either 'counts' or 'data'")
#     }
#     else if (!missing(x = counts) && !missing(x = data)) {
#         stop("Either 'counts' or 'data' must be missing; both cannot be provided")
#     }
#     else if (!missing(x = counts)) {
#         if (anyDuplicated(rownames(x = counts))) {
#             warning("Non-unique features (rownames) present in the input matrix, making unique", 
#                 call. = FALSE, immediate. = TRUE)
#             rownames(x = counts) <- make.unique(names = rownames(x = counts))
#         }
#         if (anyDuplicated(colnames(x = counts))) {
#             warning("Non-unique cell names (colnames) present in the input matrix, making unique", 
#                 call. = FALSE, immediate. = TRUE)
#             colnames(x = counts) <- make.unique(names = colnames(x = counts))
#         }
#         if (is.null(x = colnames(x = counts))) {
#             stop("No cell names (colnames) names present in the input matrix")
#         }
#         if (any(rownames(x = counts) == "")) {
#             stop("Feature names of counts matrix cannot be empty", 
#                 call. = FALSE)
#         }
#         if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
#             stop("No feature names (rownames) names present in the input matrix")
#         }
#         if (!inherits(x = counts, what = "dgCMatrix")) {
#             counts <- as(object = as.matrix(x = counts), Class = "dgCMatrix")
#         }
#         if (min.features > 0) {
#             nfeatures <- Matrix::colSums(x = counts > 0)
#             counts <- counts[, which(x = nfeatures >= min.features)]
#         }
#         if (min.cells > 0) {
#             num.cells <- Matrix::rowSums(x = counts > 0)
#             counts <- counts[which(x = num.cells >= min.cells), 
#                 ]
#         }
#         data <- counts
#     }
#     else if (!missing(x = data)) {
#         if (anyDuplicated(rownames(x = data))) {
#             warning("Non-unique features (rownames) present in the input matrix, making unique", 
#                 call. = FALSE, immediate. = TRUE)
#             rownames(x = data) <- make.unique(names = rownames(x = data))
#         }
#         if (anyDuplicated(colnames(x = data))) {
#             warning("Non-unique cell names (colnames) present in the input matrix, making unique", 
#                 call. = FALSE, immediate. = TRUE)
#             colnames(x = data) <- make.unique(names = colnames(x = data))
#         }
#         if (is.null(x = colnames(x = data))) {
#             stop("No cell names (colnames) names present in the input matrix")
#         }
#         if (any(rownames(x = data) == "")) {
#             stop("Feature names of data matrix cannot be empty", 
#                 call. = FALSE)
#         }
#         if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
#             stop("No feature names (rownames) names present in the input matrix")
#         }
#         if (min.cells != 0 | min.features != 0) {
#             warning("No filtering performed if passing to data rather than counts", 
#                 call. = FALSE, immediate. = TRUE)
#         }
#         counts <- new(Class = "matrix")
#     }
#     if (!is.vector(x = rownames(x = counts))) {
#         rownames(x = counts) <- as.vector(x = rownames(x = counts))
#     }
#     if (!is.vector(x = colnames(x = counts))) {
#         colnames(x = counts) <- as.vector(x = colnames(x = counts))
#     }
#     if (!is.vector(x = rownames(x = data))) {
#         rownames(x = data) <- as.vector(x = rownames(x = data))
#     }
#     if (!is.vector(x = colnames(x = data))) {
#         colnames(x = data) <- as.vector(x = colnames(x = data))
#     }
#     if (any(grepl(pattern = "_", x = rownames(x = counts))) || 
#         any(grepl(pattern = "_", x = rownames(x = data)))) {
#         warning("Feature names cannot have underscores ('_'), replacing with dashes ('-')", 
#             call. = FALSE, immediate. = TRUE)
#         rownames(x = counts) <- gsub(pattern = "_", replacement = "-", 
#             x = rownames(x = counts))
#         rownames(x = data) <- gsub(pattern = "_", replacement = "-", 
#             x = rownames(x = data))
#     }
#     if (any(grepl(pattern = "|", x = rownames(x = counts), 
#         fixed = TRUE)) || any(grepl(pattern = "|", x = rownames(x = data), 
#         fixed = TRUE))) {
#         warning("Feature names cannot have pipe characters ('|'), replacing with dashes ('-')", 
#             call. = FALSE, immediate. = TRUE)
#         rownames(x = counts) <- gsub(pattern = "|", replacement = "-", 
#             x = rownames(x = counts), fixed = TRUE)
#         rownames(x = data) <- gsub(pattern = "|", replacement = "-", 
#             x = rownames(x = data), fixed = TRUE)
#     }
#     init.meta.features <- data.frame(row.names = rownames(x = data))
#     assay <- new(Class = "Assay", counts = counts, data = data, 
#         scale.data = new(Class = "matrix"), meta.features = init.meta.features)
#     return(assay)
# }

###update object to avoid warning.
data("pbmc3k")
pbmc <- UpdateSeuratObject(pbmc3k)
rm(pbmc3k)
pbmc
```

    ## An object of class Seurat 
    ## 13714 features across 2700 samples within 1 assay 
    ## Active assay: RNA (13714 features)

### 2. 基本预处理

作者在原教程说： 

> The steps below encompass the standard pre-processing
workflow for scRNA-seq data in Seurat. These represent the selection and
filtration of cells based on QC metrics, data normalization and scaling,
and the detection of highly variable features.

#### 2.1 细胞质控

三种基本的QC metrics

> 1.  The number of unique genes detected in each cell.
>       - Low-quality cells or empty droplets will often have very few
>         genes
>       - Cell doublets or multiplets may exhibit an aberrantly high
>         gene count
> 2.  Similarly, the total number of molecules detected within a cell
>     (correlates strongly with unique genes)
> 3.  The percentage of reads that map to the mitochondrial genome
>       - Low-quality / dying cells often exhibit extensive
>         mitochondrial contamination
>       - We calculate mitochondrial QC metrics with the
>         `PercentageFeatureSet` function, which calculates the
>         percentage of counts originating from a set of features
>       - We use the set of all genes starting with `MT-` as a set of
>         mitochondrial genes

注意，人的线粒体基因是“MT-”开头，而小鼠的线粒体基因是“mt-”开头

``` r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

### how does PercentageFeatureSet work
# PercentageFeatureSet
# function (object, pattern = NULL, features = NULL, col.name = NULL, 
#     assay = NULL) 
# {
#     assay <- assay %||% DefaultAssay(object = object)
#     if (!is.null(x = features) && !is.null(x = pattern)) {
#         warning("Both pattern and features provided. Pattern is being ignored.")
#     }
#     features <- features %||% grep(pattern = pattern, x = rownames(x = object[[assay]]), 
#         value = TRUE)
#     percent.featureset <- colSums(x = GetAssayData(object = object, 
#         assay = assay, slot = "counts")[features, , drop = FALSE])/object[[paste0("nCount_", 
#         assay)]] * 100
#     if (!is.null(x = col.name)) {
#         object <- AddMetaData(object = object, metadata = percent.featureset, 
#             col.name = col.name)
#         return(object)
#     }
#     return(percent.featureset)
# }

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
```

    ##                orig.ident nCount_RNA nFeature_RNA seurat_annotations percent.mt
    ## AAACATACAACCAC     pbmc3k       2419          779       Memory CD4 T  3.0177759
    ## AAACATTGAGCTAC     pbmc3k       4903         1352                  B  3.7935958
    ## AAACATTGATCAGC     pbmc3k       3147         1129       Memory CD4 T  0.8897363
    ## AAACCGTGCTTCCG     pbmc3k       2639          960         CD14+ Mono  1.7430845
    ## AAACCGTGTATGCG     pbmc3k        980          521                 NK  1.2244898

由于我们用的是作者给了metadata的数据，里面已经出现了细胞类型的注释，见`seurat_annotation`这一项；

QC metric的可视化：

``` r
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

``` r
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

最终选择的质控标准为：

>   - We filter cells that have unique feature counts over 2,500 or less
>     than 200
>   - We filter cells that have \>5% mitochondrial counts

``` r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

#### 2.2 标准化

> After removing unwanted cells from the dataset, the next step is to
> normalize the data. By default, we employ a global-scaling
> normalization method “LogNormalize” that normalizes the feature
> expression measurements for each cell by the total expression,
> multiplies this by a scale factor (10,000 by default), and
> log-transforms the result. Normalized values are stored in
> pbmc\[\[“RNA”\]\]@data.

``` r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

#### 2.3 特征选择

哪些基因能反应不同细胞之间的异质性？是那些表达差异大的基因；

注意`FindVariableFeatures`是S3 generic，泛型函数。
如何查看一个泛型函数的源代码呢，我们先用`methods`函数匹配该范型函数的名字：

``` r
methods(FindVariableFeatures)
```

    ## [1] FindVariableFeatures.Assay*   FindVariableFeatures.default*
    ## [3] FindVariableFeatures.Seurat* 
    ## see '?methods' for accessing help and source code

星号表明我们不能直接通过运行函数名字来查看其源代码，但是我们可以通过运行 **getAnywhere**函数来获取这个函数，

``` r
getAnywhere(FindVariableFeatures.Seurat)
```

    ## A single object matching 'FindVariableFeatures.Seurat' was found
    ## It was found in the following places
    ##   registered S3 method for FindVariableFeatures from namespace Seurat
    ##   namespace:Seurat
    ## with value
    ## 
    ## function (object, assay = NULL, selection.method = "vst", loess.span = 0.3, 
    ##     clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, 
    ##     num.bin = 20, binning.method = "equal_width", nfeatures = 2000, 
    ##     mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), verbose = TRUE, 
    ##     ...) 
    ## {
    ##     assay <- assay %||% DefaultAssay(object = object)
    ##     assay.data <- GetAssay(object = object, assay = assay)
    ##     assay.data <- FindVariableFeatures(object = assay.data, selection.method = selection.method, 
    ##         loess.span = loess.span, clip.max = clip.max, mean.function = mean.function, 
    ##         dispersion.function = dispersion.function, num.bin = num.bin, 
    ##         binning.method = binning.method, nfeatures = nfeatures, 
    ##         mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff, 
    ##         verbose = verbose, ...)
    ##     object[[assay]] <- assay.data
    ##     object <- LogSeuratCommand(object = object)
    ##     return(object)
    ## }
    ## <bytecode: 0x0000000022fa2410>
    ## <environment: namespace:Seurat>

我们可以发现，默认的**FindVariableFeatures.Seurat**method调用了**FindVariableFeatures.Assay**：

``` r
getAnywhere(FindVariableFeatures.Assay)
```

    ## A single object matching 'FindVariableFeatures.Assay' was found
    ## It was found in the following places
    ##   registered S3 method for FindVariableFeatures from namespace Seurat
    ##   namespace:Seurat
    ## with value
    ## 
    ## function (object, selection.method = "vst", loess.span = 0.3, 
    ##     clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, 
    ##     num.bin = 20, binning.method = "equal_width", nfeatures = 2000, 
    ##     mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), verbose = TRUE, 
    ##     ...) 
    ## {
    ##     if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 
    ##         2) {
    ##         stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
    ##     }
    ##     if (selection.method == "vst") {
    ##         data <- GetAssayData(object = object, slot = "counts")
    ##         if (IsMatrixEmpty(x = data)) {
    ##             warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
    ##             data <- GetAssayData(object = object, slot = "data")
    ##         }
    ##     }
    ##     else {
    ##         data <- GetAssayData(object = object, slot = "data")
    ##     }
    ##     hvf.info <- FindVariableFeatures(object = data, selection.method = selection.method, 
    ##         loess.span = loess.span, clip.max = clip.max, mean.function = mean.function, 
    ##         dispersion.function = dispersion.function, num.bin = num.bin, 
    ##         binning.method = binning.method, verbose = verbose, ...)
    ##     object[[names(x = hvf.info)]] <- hvf.info
    ##     hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 
    ##         0), ]
    ##     if (selection.method == "vst") {
    ##         hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, 
    ##             decreasing = TRUE), , drop = FALSE]
    ##     }
    ##     else {
    ##         hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), 
    ##             , drop = FALSE]
    ##     }
    ##     selection.method <- switch(EXPR = selection.method, mvp = "mean.var.plot", 
    ##         disp = "dispersion", selection.method)
    ##     top.features <- switch(EXPR = selection.method, mean.var.plot = {
    ##         means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 
    ##             1] < mean.cutoff[2])
    ##         dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & 
    ##             (hvf.info[, 3] < dispersion.cutoff[2])
    ##         rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    ##     }, dispersion = head(x = rownames(x = hvf.info), n = nfeatures), 
    ##         vst = head(x = rownames(x = hvf.info), n = nfeatures), 
    ##         stop("Unkown selection method: ", selection.method))
    ##     VariableFeatures(object = object) <- top.features
    ##     vf.name <- ifelse(test = selection.method == "vst", yes = "vst", 
    ##         no = "mvp")
    ##     vf.name <- paste0(vf.name, ".variable")
    ##     object[[vf.name]] <- rownames(x = object[[]]) %in% top.features
    ##     return(object)
    ## }
    ## <bytecode: 0x000000002dc5d050>
    ## <environment: namespace:Seurat>

千层饼的最后一层；

``` r
getAnywhere(FindVariableFeatures.default)
```

    ## A single object matching 'FindVariableFeatures.default' was found
    ## It was found in the following places
    ##   registered S3 method for FindVariableFeatures from namespace Seurat
    ##   namespace:Seurat
    ## with value
    ## 
    ## function (object, selection.method = "vst", loess.span = 0.3, 
    ##     clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, 
    ##     num.bin = 20, binning.method = "equal_width", verbose = TRUE, 
    ##     ...) 
    ## {
    ##     CheckDots(...)
    ##     if (!inherits(x = object, "Matrix")) {
    ##         object <- as(object = as.matrix(x = object), Class = "Matrix")
    ##     }
    ##     if (!inherits(x = object, what = "dgCMatrix")) {
    ##         object <- as(object = object, Class = "dgCMatrix")
    ##     }
    ##     if (selection.method == "vst") {
    ##         if (clip.max == "auto") {
    ##             clip.max <- sqrt(x = ncol(x = object))
    ##         }
    ##         hvf.info <- data.frame(mean = rowMeans(x = object))
    ##         hvf.info$variance <- SparseRowVar2(mat = object, mu = hvf.info$mean, 
    ##             display_progress = verbose)
    ##         hvf.info$variance.expected <- 0
    ##         hvf.info$variance.standardized <- 0
    ##         not.const <- hvf.info$variance > 0
    ##         fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
    ##             data = hvf.info[not.const, ], span = loess.span)
    ##         hvf.info$variance.expected[not.const] <- 10^fit$fitted
    ##         hvf.info$variance.standardized <- SparseRowVarStd(mat = object, 
    ##             mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
    ##             vmax = clip.max, display_progress = verbose)
    ##         colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
    ##     }
    ##     else {
    ##         if (!inherits(x = mean.function, what = "function")) {
    ##             stop("'mean.function' must be a function")
    ##         }
    ##         if (!inherits(x = dispersion.function, what = "function")) {
    ##             stop("'dispersion.function' must be a function")
    ##         }
    ##         feature.mean <- mean.function(object, verbose)
    ##         feature.dispersion <- dispersion.function(object, verbose)
    ##         names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
    ##         feature.dispersion[is.na(x = feature.dispersion)] <- 0
    ##         feature.mean[is.na(x = feature.mean)] <- 0
    ##         data.x.breaks <- switch(EXPR = binning.method, equal_width = num.bin, 
    ##             equal_frequency = c(-1, quantile(x = feature.mean[feature.mean > 
    ##                 0], probs = seq.int(from = 0, to = 1, length.out = num.bin))), 
    ##             stop("Unknown binning method: ", binning.method))
    ##         data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    ##         names(x = data.x.bin) <- names(x = feature.mean)
    ##         mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, 
    ##             FUN = mean)
    ##         sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, 
    ##             FUN = sd)
    ##         feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)])/sd.y[as.numeric(x = data.x.bin)]
    ##         names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    ##         hvf.info <- data.frame(feature.mean, feature.dispersion, 
    ##             feature.dispersion.scaled)
    ##         rownames(x = hvf.info) <- rownames(x = object)
    ##         colnames(x = hvf.info) <- paste0("mvp.", c("mean", "dispersion", 
    ##             "dispersion.scaled"))
    ##     }
    ##     return(hvf.info)
    ## }
    ## <bytecode: 0x0000000023f13fd0>
    ## <environment: namespace:Seurat>

忽略这些技术细节，进行特征选择；

``` r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

#### 2.4 Scaling the data

这步的目的是为了后续的PCA：

> Next, we apply a linear transformation (‘scaling’) that is a standard
> pre-processing step prior to dimensional reduction techniques like
> PCA. The **ScaleData** function: + Shifts the expression of each gene,
> so that the mean expression across cells is 0 + Scales the expression
> of each gene, so that the variance across cells is 1 + This step gives
> equal weight in downstream analyses, so that highly-expressed genes do
> not dominate + The results of this are stored in
> **pbmc\[\[“RNA”\]\]@scale.data**

回归掉percent.mt对于PCA的影响。这步是一步限速步骤；

``` r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes,vars.to.regress = "percent.mt")
```

有一个问题后面的marker基因一定是HVG吗？

#### 2.5 线性降维(PCA)

> Next we perform PCA on the scaled data. By default, only the
> previously determined variable features are used as input, but can be
> defined using features argument if you wish to choose a different
> subset.

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  CST3, TYROBP, LST1, AIF1, FTL 
    ## Negative:  MALAT1, LTB, IL32, IL7R, CD2 
    ## PC_ 2 
    ## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
    ## Negative:  NKG7, PRF1, CST7, GZMA, GZMB 
    ## PC_ 3 
    ## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPA1 
    ## Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
    ## PC_ 4 
    ## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
    ## Negative:  VIM, IL7R, S100A6, S100A8, IL32 
    ## PC_ 5 
    ## Positive:  GZMB, FGFBP2, S100A8, NKG7, GNLY 
    ## Negative:  LTB, IL7R, CKB, MS4A7, RP11-290F20.3

``` r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

``` r
DimPlot(pbmc, reduction = "pca")
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

> In particular **DimHeatmap** allows for easy exploration of the
> primary sources of heterogeneity in a dataset, and can be useful when
> trying to decide which PCs to include for further downstream analyses.
> Both cells and features are ordered according to their PCA scores.
> Setting cells to a number plots the ‘extreme’ cells on both ends of
> the spectrum, which dramatically speeds plotting for large datasets.
> Though clearly a supervised analysis, we find this to be a valuable
> tool for exploring correlated feature sets.

``` r
###     Plot an equal number of genes with both + and - scores.
mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE,fast = F)+scale_fill_gradientn(colors  = mypal)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

``` r
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

#### 2.6 Determine the ‘dimensionality’ of the dataset

> To overcome the extensive technical noise in any single feature for
> scRNA-seq data, Seurat clusters cells based on their PCA scores, with
> each PC essentially representing a ‘metafeature’ that combines
> information across a correlated feature set. The top principal
> components therefore represent a robust compression of the dataset.
> However, how many componenets should we choose to include? 10? 20?
> 100?

两种统计方法，`JackStraw`和`ElbowPlot`，前者比较耗时，不再展示了，用后者

``` r
ElbowPlot(pbmc)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

作者给出了更进一步的解释

> Identifying the true dimensionality of a dataset – can be
> challenging/uncertain for the user. We therefore suggest these three
> approaches to consider. The first is more supervised, exploring PCs to
> determine relevant sources of heterogeneity, and could be used in
> conjunction with GSEA for example. The second implements a statistical
> test based on a random null model, but is time-consuming for large
> datasets, and may not return a clear PC cutoff. The third is a
> heuristic that is commonly used, and can be calculated instantly. In
> this example, all three approaches yielded similar results, but we
> might have been justified in choosing anything between PC 7-12 as a
> cutoff.

> We chose 10 here, but encourage users to consider the following:

>   - Dendritic cell and NK aficionados may recognize that genes
>     strongly associated with PCs 12 and 13 define rare immune subsets
>     (i.e. MZB1 is a marker for plasmacytoid DCs). However, these
>     groups are so rare, they are difficult to distinguish from
>     background noise for a dataset of this size without prior
>     knowledge.
>   - We encourage users to repeat downstream analyses with a different
>     number of PCs (10, 15, or even 50\!). As you will observe, the
>     results often do not differ dramatically.
>   - We advise users to err on the higher side when choosing this
>     parameter. For example, performing downstream analyses with only 5
>     PCs does signifcanltly and adversely affect results.

### 3. 后续分析

#### 3.1 聚类

**FindNeighbors**构建构建SNN-graph, 而**FindClusters**用来实现Louvain
algorithm，进行图聚类；

``` r
methods(FindNeighbors)
```

    ## [1] FindNeighbors.Assay*   FindNeighbors.default* FindNeighbors.dist*   
    ## [4] FindNeighbors.Seurat* 
    ## see '?methods' for accessing help and source code

``` r
getAnywhere(FindNeighbors.Seurat)
```

    ## A single object matching 'FindNeighbors.Seurat' was found
    ## It was found in the following places
    ##   registered S3 method for FindNeighbors from namespace Seurat
    ##   namespace:Seurat
    ## with value
    ## 
    ## function (object, reduction = "pca", dims = 1:10, assay = NULL, 
    ##     features = NULL, k.param = 20, compute.SNN = TRUE, prune.SNN = 1/15, 
    ##     nn.method = "rann", annoy.metric = "euclidean", nn.eps = 0, 
    ##     verbose = TRUE, force.recalc = FALSE, do.plot = FALSE, graph.name = NULL, 
    ##     ...) 
    ## {
    ##     CheckDots(...)
    ##     if (!is.null(x = dims)) {
    ##         assay <- DefaultAssay(object = object[[reduction]])
    ##         data.use <- Embeddings(object = object[[reduction]])
    ##         if (max(dims) > ncol(x = data.use)) {
    ##             stop("More dimensions specified in dims than have been computed")
    ##         }
    ##         data.use <- data.use[, dims]
    ##         neighbor.graphs <- FindNeighbors(object = data.use, k.param = k.param, 
    ##             compute.SNN = compute.SNN, prune.SNN = prune.SNN, 
    ##             nn.method = nn.method, annoy.metric = annoy.metric, 
    ##             nn.eps = nn.eps, verbose = verbose, force.recalc = force.recalc, 
    ##             ...)
    ##     }
    ##     else {
    ##         assay <- assay %||% DefaultAssay(object = object)
    ##         data.use <- GetAssay(object = object, assay = assay)
    ##         neighbor.graphs <- FindNeighbors(object = data.use, features = features, 
    ##             k.param = k.param, compute.SNN = compute.SNN, prune.SNN = prune.SNN, 
    ##             nn.method = nn.method, annoy.metric = annoy.metric, 
    ##             nn.eps = nn.eps, verbose = verbose, force.recalc = force.recalc, 
    ##             ...)
    ##     }
    ##     graph.name <- graph.name %||% paste0(assay, "_", names(x = neighbor.graphs))
    ##     for (ii in 1:length(x = graph.name)) {
    ##         DefaultAssay(object = neighbor.graphs[[ii]]) <- assay
    ##         object[[graph.name[[ii]]]] <- neighbor.graphs[[ii]]
    ##     }
    ##     if (do.plot) {
    ##         if (!"tsne" %in% names(x = object@reductions)) {
    ##             warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    ##         }
    ##         else {
    ##             if (nrow(x = Embeddings(object = object[["tsne"]])) != 
    ##                 ncol(x = object)) {
    ##                 warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
    ##             }
    ##             else {
    ##                 net <- graph.adjacency(adjmatrix = as.matrix(x = neighbor.graphs[[2]]), 
    ##                   mode = "undirected", weighted = TRUE, diag = FALSE)
    ##                 plot.igraph(x = net, layout = as.matrix(x = Embeddings(object = object[["tsne"]])), 
    ##                   edge.width = E(graph = net)$weight, vertex.label = NA, 
    ##                   vertex.size = 0)
    ##             }
    ##         }
    ##     }
    ##     object <- LogSeuratCommand(object = object)
    ##     return(object)
    ## }
    ## <bytecode: 0x000000001e49f8e0>
    ## <environment: namespace:Seurat>

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2638
    ## Number of edges: 95930
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8737
    ## Number of communities: 9
    ## Elapsed time: 0 seconds

#### 3.2 Run UMAP/tsne

run tsne

``` r
pbmc <- RunTSNE(pbmc,dims = 1:10)
DimPlot(pbmc,label = T, reduction = "tsne")
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

draw snn graph on tsne-embeding

``` r
test <- pbmc[["RNA_snn"]]


net <- graph.adjacency(adjmatrix = as.matrix(x = test), 
                  mode = "undirected", weighted = TRUE, 
                  diag = FALSE)
plot.igraph(x = net, 
            layout = as.matrix(x = Embeddings(object = pbmc[["tsne"]])),
            edge.width = E(graph = net)$weight, vertex.label = NA, 
                  vertex.size = 0)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

run umap

``` r
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc,umap.method = "umap-learn", dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc,label = T, reduction = "umap")
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

``` r
test <- pbmc[["RNA_snn"]]

net <- graph.adjacency(adjmatrix = as.matrix(x = test), 
                  mode = "undirected", weighted = TRUE, 
                  diag = FALSE)
plot.igraph(x = net, 
            layout = as.matrix(x = Embeddings(object = pbmc[["umap"]])),
            edge.width = E(graph = net)$weight, vertex.label = NA, 
                  vertex.size = 0)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

#### 3.3 Finding differentially expressed features (cluster biomarkers)

之前分群结果做差异表达；

> Seurat can help you find markers that define clusters via differential
> expression. By default, it identifes positive and negative markers of
> a single cluster (specified in **ident.1**), compared to all other
> cells. **FindAllMarkers** automates this process for all clusters, but
> you can also test groups of clusters vs. each other, or against all
> cells.

> The **min.pct** argument requires a feature to be detected at a
> minimum percentage in either of the two groups of cells, and the
> thresh.test argument requires a feature to be differentially expressed
> (on average) by some amount between the two groups. You can set both
> of these to 0, but with a dramatic increase in time - since this will
> test a large number of features that are unlikely to be highly
> discriminatory. As another option to speed up these computations,
> **max.cells.per.ident** can be set. This will downsample each identity
> class to have no more cells than whatever this is set to. While there
> is generally going to be a loss in power, the speed increases can be
> significiant and the most highly differentially expressed features
> will likely still rise to the top.

``` r
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
```

    ## # A tibble: 18 x 7
    ## # Groups:   cluster [9]
    ##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene    
    ##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
    ##  1 1.88e-117     0.748 0.913 0.588 2.57e-113 0       LDHB    
    ##  2 5.01e- 85     0.931 0.437 0.108 6.88e- 81 0       CCR7    
    ##  3 0.            3.86  0.996 0.215 0.        1       S100A9  
    ##  4 0.            3.80  0.975 0.121 0.        1       S100A8  
    ##  5 2.61e- 81     0.886 0.981 0.65  3.58e- 77 2       LTB     
    ##  6 1.22e- 59     0.886 0.669 0.25  1.68e- 55 2       CD2     
    ##  7 0.            2.99  0.939 0.042 0.        3       CD79A   
    ##  8 1.06e-269     2.49  0.623 0.022 1.45e-265 3       TCL1A   
    ##  9 5.98e-221     2.23  0.987 0.226 8.20e-217 4       CCL5    
    ## 10 1.42e-173     2.08  0.572 0.051 1.94e-169 4       GZMK    
    ## 11 3.51e-184     2.30  0.975 0.134 4.82e-180 5       FCGR3A  
    ## 12 2.03e-125     2.14  1     0.315 2.78e-121 5       LST1    
    ## 13 3.17e-267     3.35  0.961 0.068 4.35e-263 6       GZMB    
    ## 14 1.04e-189     3.66  0.961 0.132 1.43e-185 6       GNLY    
    ## 15 1.48e-220     2.68  0.812 0.011 2.03e-216 7       FCER1A  
    ## 16 1.67e- 21     1.99  1     0.513 2.28e- 17 7       HLA-DPB1
    ## 17 7.73e-200     5.02  1     0.01  1.06e-195 8       PF4     
    ## 18 3.68e-110     5.94  1     0.024 5.05e-106 8       PPBP

可视化：

``` r
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

``` r
# you can plot raw counts as well
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), slot = "counts", log = TRUE)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

使用**FeatureScatter**获得和流式图一样的效果；

``` r
FeatureScatter(object = pbmc,
               feature1 = "MS4A1",
               feature2 = "CD79A")+
  ggtitle(label = NULL)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

用**FeaturePlot**在Embeding上展示表达量；

``` r
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
    "CD8A"))
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

气泡图**DotPlot**

``` r
DotPlot(object = pbmc,
        features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",  "CD8A"))+
  coord_flip()
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-31-1.png" style="display: block; margin: auto;" />

**RidgePlot**

``` r
RidgePlot(object = pbmc,
          features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",  "CD8A"))
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

热图**DoHeatmap**

``` r
top10 <- pbmc.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(pbmc, features = top10$gene) + 
  scale_fill_gradientn(colors  = mypal)
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

#### 3.4 Assigning cell type identity to clusters

``` r
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono", "Memory CD4 T",  "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

<img src="/figure/posts/LearnSeurat_PBMC3k_files/figure-gfm/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />

### 4. Seurat object 详解

这一部分来自wiki

#### 4.1 The Seurat object

一个Seurat对象有如下的`slots`:

| Slot           | Function                                                                        |
| -------------- | ------------------------------------------------------------------------------- |
| `assays`       | A list of assays within this object                                             |
| `meta.data`    | Cell-level meta data                                                            |
| `active.assay` | Name of active, or default, assay                                               |
| `active.ident` | Identity classes for the current object                                         |
| `graphs`       | A list of nearest neighbor graphs                                               |
| `reductions`   | A list of DimReduc objects                                                      |
| `project.name` | User-defined project name (optional)                                            |
| `tools`        | Empty list. Tool developers can store any internal data from their methods here |
| `misc`         | Empty slot. User can store additional information here                          |
| `version`      | Seurat version used when creating the object                                    |

这个对象把单细胞数据的所有的基本信息都包含进去了，可以用基本的一些函数去获取这些信息。例如，我们想要知道这个数据对应多少细胞，多少基因，可以用`dim`;`ncol`;`nrow`;细胞或者feature的名字，可以用`rownames`;`colnames`;
我们也可以通过`names`知道里面存储的如原始表达矩阵，或者降维后对象的名字。

``` r
names(x = pbmc)
```

    ## [1] "RNA"     "RNA_nn"  "RNA_snn" "pca"     "tsne"    "umap"

``` r
rna <- pbmc[['RNA']]
```

对于Seurat对象，有一系列的函数可以对其进行操作。这些函数可以称为其所属的**methods**。多说一句，Seurat采取的是S3对象的面向对象的数据结构。

可以使用如下命令访问与Seurat对象相关的操作。

``` r
utils::methods(class = 'Seurat')
```

    ##  [1] $                       $<-                     [                      
    ##  [4] [[                      [[<-                    AddMetaData            
    ##  [7] as.CellDataSet          as.loom                 as.SingleCellExperiment
    ## [10] Command                 DefaultAssay            DefaultAssay<-         
    ## [13] dim                     dimnames                droplevels             
    ## [16] Embeddings              FindClusters            FindMarkers            
    ## [19] FindNeighbors           FindVariableFeatures    GetAssay               
    ## [22] GetAssayData            HVFInfo                 Idents                 
    ## [25] Idents<-                Key                     levels                 
    ## [28] levels<-                Loadings                merge                  
    ## [31] Misc                    Misc<-                  names                  
    ## [34] NormalizeData           OldWhichCells           Project                
    ## [37] Project<-               RenameCells             RenameIdents           
    ## [40] ReorderIdent            RunALRA                 RunCCA                 
    ## [43] RunICA                  RunLSI                  RunPCA                 
    ## [46] RunTSNE                 RunUMAP                 ScaleData              
    ## [49] ScoreJackStraw          SetAssayData            SetIdent               
    ## [52] show                    StashIdent              Stdev                  
    ## [55] subset                  SubsetData              Tool                   
    ## [58] Tool<-                  VariableFeatures        VariableFeatures<-     
    ## [61] WhichCells             
    ## see '?methods' for accessing help and source code

#### 4.2 Assay

> The `Assay` class stores single cell data.

> For typical scRNA-seq experiments, a Seurat object will have a single
> Assay (“RNA”). This assay will also store multiple ‘transformations’
> of the data, including raw counts (@counts slot), normalized data
> (@data slot), and scaled data for dimensional reduction (@scale.data
> slot).

> For more complex experiments, an object could contain multiple assays.
> These could include multi-modal data types (CITE-seq antibody-derived
> tags, ADTs), or imputed/batch-corrected measurements. Each of those
> assays has the option to store the same data transformations as well.

一个**Assay** 所含有的Slots

| Slot            | Function                                                                     |
| --------------- | ---------------------------------------------------------------------------- |
| `counts`        | Stores unnormalized data such as raw counts or TPMs                          |
| `data`          | Normalized data matrix                                                       |
| `scale.data`    | Scaled data matrix                                                           |
| `key`           | A character string to facilitate looking up features from a specific `Assay` |
| `var.features`  | A vector of features identified as variable                                  |
| `meta.features` | Feature-level meta data                                                      |

**Assay**对象也可以使用以下方法

Summary information about `Assay` objects can be had quickly and easily
using standard R functions. Object shape/dimensions can be found using
the `dim`, `ncol`, and `nrow` functions; cell and feature names can be
found using the `colnames` and `rownames` functions, respectively, or
the `dimnames` function.

更多的方法见

``` r
utils::methods(class = 'Assay')
```

    ##  [1] [                    [[                   [[<-                
    ##  [4] AddMetaData          DefaultAssay         DefaultAssay<-      
    ##  [7] dim                  dimnames             FindNeighbors       
    ## [10] FindVariableFeatures GetAssayData         HVFInfo             
    ## [13] Key                  Key<-                merge               
    ## [16] Misc                 Misc<-               NormalizeData       
    ## [19] OldWhichCells        RenameCells          RunICA              
    ## [22] RunLSI               RunPCA               ScaleData           
    ## [25] SetAssayData         show                 subset              
    ## [28] SubsetData           VariableFeatures     VariableFeatures<-  
    ## [31] WhichCells          
    ## see '?methods' for accessing help and source code

Data Access

``` r
# GetAssayData allows pulling from a specific slot rather than just data
GetAssayData(object = rna, slot = 'scale.data')[1:3, 1:3]
```

    ##               AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC
    ## AL627309.1       -0.06433822    -0.06968772    -0.04966479
    ## AP006222.2       -0.02663018    -0.02065038    -0.04303249
    ## RP11-206L10.2    -0.03015459    -0.02024084    -0.05734758

``` r
head(x = HVFInfo(object = rna,selection.method = "vst"))
```

    ##                      mean    variance variance.standardized
    ## AL627309.1    0.003411676 0.003401325             0.9330441
    ## AP006222.2    0.001137225 0.001136363             0.9924937
    ## RP11-206L10.2 0.001895375 0.001892500             0.9627290
    ## RP11-206L10.9 0.001137225 0.001136363             0.9924937
    ## LINC00115     0.006823351 0.006779363             0.9062135
    ## NOC2L         0.107278241 0.159514698             0.7849309

The key

``` r
# Key both accesses and sets the key slot for an Assay object
> Key(object = rna)
"rna_"
> Key(object = rna) <- 'myRNA_'
> Key(object = rna)
"myRNA_"
# Pull a feature from the RNA assay on the Seurat level
> head(x = FetchData(object = pbmc, vars.fetch = 'rna_MS4A1'))
               rna_MS4A1
AAACATACAACCAC  0.000000
AAACATTGAGCTAC  2.583047
AAACATTGATCAGC  0.000000
AAACCGTGCTTCCG  0.000000
AAACCGTGTATGCG  0.000000
AAACGCACTGGTAC  0.000000
```

The `DimReduc` object represents a dimensional reduction taken upon the
Seurat object.

#### 4.3 The `DimReduc` object

The `DimReduc` object represents a dimensional reduction taken upon the
Seurat object.

| Slot                         | Function                                                                        |
| ---------------------------- | ------------------------------------------------------------------------------- |
| `cell.embeddings`            | A matrix with cell embeddings                                                   |
| `feature.loadings`           | A matrix with feature loadings                                                  |
| `feature.loadings.projected` | A matrix with projected feature loadings                                        |
| `assay.used`                 | Assay used to calculate this dimensional reduction                              |
| `stdev`                      | Standard deviation for the dimensional reduction                                |
| `key`                        | A character string to facilitate looking up features from a specific `DimReduc` |
| `jackstraw`                  | Results from the `JackStraw` function                                           |
| `misc`                       | …                                                                               |

和之前的很类似

``` r
pca <- pbmc[["pca"]]
# The following examples use the PCA dimensional reduction from the PBMC 3k dataset
> pca
A dimensional reduction object with key PC
 Number of dimensions: 20
 Projected dimensional reduction calculated: FALSE
 Jackstraw run: FALSE
# nrow and ncol provide the number of features and cells, respectively
# dim provides both nrow and ncol at the same time
> dim(x = pca)
[1] 1838 2638
# length provides the number of dimensions calculated
> length(x = pca)
[1] 20
# In addtion to rownames and colnames, one can use dimnames
# which provides a two-length list with both rownames and colnames
> head(x = rownames(x = rna))
[1] "TNFRSF4"  "CPSF3L"   "ATAD3C"   "C1orf86"  "RER1"     "TNFRSF25"
> head(x = colnames(x = rna))
[1] "AAACATACAACCAC" "AAACATTGAGCTAC" "AAACATTGATCAGC" "AAACCGTGCTTCCG"
[5] "AAACCGTGTATGCG" "AAACGCACTGGTAC"
```

Access data

``` r
# The key can be used to pull cell embeddings for specific dimensions from the Seurat level
> Key(object = pca)
"PC"
> head(x = FetchData(object = pbmc, vars.fetch = 'PC1'))
                      PC1
AAACATACAACCAC   5.569384
AAACATTGAGCTAC   7.216456
AAACATTGATCAGC   2.706629
AAACCGTGCTTCCG -10.134042
AAACCGTGTATGCG  -1.099311
AAACGCACTGGTAC   1.455335
# DefaultAssay gets the name of the Assay object used to calculate the DimReduc
> DefaultAssay(object = pca)
[1] "RNA"
# Stdev gets the vector of standard deviations for each dimension embedded.
Stdev(object = pca)
 [1] 5.666584 4.326466 3.952192 3.638124 2.191529 1.996551 1.877891 1.798251
 [9] 1.766873 1.753684 1.731568 1.720525 1.718079 1.715879 1.707009 1.702660
[17] 1.697318 1.692549 1.686149 1.683967
```

在其上可以执行的**method**有

``` r
utils::methods(class = "DimReduc")
```

    ##  [1] [              [[             Cells          DefaultAssay   DefaultAssay<-
    ##  [6] dim            dimnames       Embeddings     IsGlobal       JS            
    ## [11] JS<-           Key            Key<-          length         Loadings      
    ## [16] Loadings<-     names          print          RenameCells    RunTSNE       
    ## [21] ScoreJackStraw show           Stdev          subset        
    ## see '?methods' for accessing help and source code

#### 4.4 R面向对象编程的更多细节；

关于面向对象，以及S3对象的教程，更多可见：

1.  [R深入|面向对象——泛型函数](https://zhuanlan.zhihu.com/p/31160374)
2.  [OO field guide](http://adv-r.had.co.nz/OO-essentials.html)
3.  [R语言面向对象编程](https://dataxujing.github.io/R_oop/)
