---
title: Combine pheatmap
author:
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - R
date:
keywords:
  description: show how to combine pheatmap and test Rmd
photos:
---


Talk is cheap, this is code:

``` r
library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(colormap)
items=names(colormaps)
plot_list=list()
for (a in items[1:8]){
  x= pheatmap(volcano,
              cluster_rows = F,
              cluster_cols = F,
              main = a,
              height = 3,
              width = 3,
              color = colormap_pal(colormap = colormaps[[a]])(100),silent = T)
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}

cowplot::plot_grid(plotlist = plot_list[1:8],ncol = 2,nrow = 4)
```

<img src="/figure/posts/combine_pheatmap_files/figure-markdown_github/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

test equation: $E=mc^2$, math formula need second editing?
