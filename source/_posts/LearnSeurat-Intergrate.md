---
title: LearnSeurat_Intergrate
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
keywords:
description: Intergration in terminnal
date: 2020-05-27 11:44:09
photos:
---


### 0. 说明

更为详细的使用，可以参考Seurat官方教程，或者[Seurat进行单细胞RNA-seq数据整合](https://blog.mrdoge.cool/index.php/2020/02/28/seurat_integration/)。

### 1. 快速进行

实际分析中，如果细胞数目很多的话，进行整合分析的时候就会很慢，一种方式是将任务提交通过命令行提交到服务器或者集群上运行，并输出。

传递命令行的参数的脚本可以这么写：

```R
    #### reference base Integration
    .libPaths("/home/user/lib/R/library")
    library(Seurat)
    
    args = commandArgs(trailingOnly=TRUE)
    # test if there is at least one argument: if not, return an error
    if (length(args)==0) {
      stop("At least one argument must be supplied (input file).n", call.=FALSE)
    } 
    
    ###args[1] seu.list.rds
    ###args[2] seu.intergrated.rds
    ###args[3] UMAPplot.pdf
    
    ## program...
    seu.list <- readRDS(file = args[1])
    
    for (i in names(seu.list)) {
      seu.list[[i]] <- SCTransform(seu.list[[i]], verbose = FALSE)
    }
    seu.list.features <- SelectIntegrationFeatures(object.list = seu.list, 
                                                   nfeatures = 3000)
    seu.list <- PrepSCTIntegration(object.list = seu.list, 
                                   anchor.features = seu.list.features)
    
    reference_dataset <- 1
    names(seu.list)[1]
    
    seu.list.anchors <- FindIntegrationAnchors(object.list = seu.list, 
                                               normalization.method = "SCT", 
                                               anchor.features = seu.list.features, 
                                               reference = reference_dataset)
    
    seu.list.integrated <- IntegrateData(anchorset = seu.list.anchors, 
                                         normalization.method = "SCT")
    
    seu.list.integrated <- RunPCA(object = seu.list.integrated, verbose = FALSE)
    seu.list.integrated <- RunUMAP(object = seu.list.integrated, umap.method = "umap-learn",
                                   dims = 1:30)
    saveRDS(seu.list.integrated,file = args[2])
    
    #### show Integrated result
    plots <- DimPlot(seu.list.integrated, 
                     group.by = c("orig.ident", "cell.type"))
    plots & theme(legend.position = "top") & 
      guides(color = guide_legend(nrow = 4, byrow = TRUE, 
                                  override.aes = list(size = 2.5)))
    ggsave(filename = args[3],width = 12,height = 6)
```

将其命名为`Integrated_ref_based.R`

然后在终端输入`nohup Rscript Integrated_ref_based.R seu.list.rds
seu.intergrated.rds UMAPplot.pdf > log.txt &`提交即可；

### 2. 评论

当然，这个脚本读者还可以根据自己需求加以完善。
