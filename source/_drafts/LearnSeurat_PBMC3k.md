本系列假定读者对于单细胞测序的数据分析和Seurat的官方教程有所了解。

本篇研究最基础的PBMC3k。其实这里只有2700个外周血的细胞。注意到，由于取样是外周血，没有干细胞的存在，所以可以认为样品处于稳态。这个教程就是讲稳态下的单细胞测序分析是如何进行的。Seurat的官方教程的缺点之一就是没有涉及动态过程的单细胞分析如何进行。

### 1\. 导入数据

使用作者已经构建好的数据进行分析。

``` r
library(Seurat)
library(SeuratData)
library(ggplot2)
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

data("pbmc3k")
pbmc3k <- UpdateSeuratObject(pbmc3k)
pbmc3k
```

    ## An object of class Seurat 
    ## 13714 features across 2700 samples within 1 assay 
    ## Active assay: RNA (13714 features)
