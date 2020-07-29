---
title: Learn-igraph-More about igraph
author: Landau1994
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - Graph
  - Network
description: 更多关于igraph的知识
date: 2020-07-29 23:23:30
keywords:
photos:
---


#### 0. 说明

我们接着讲更多关于对igraph对象的操作，参考[Statistical Network Analysis with
igraph](https://sites.fas.harvard.edu/~airoldi/pub/books/BookDraft-CsardiNepuszAiroldi2016.pdf)第一章。

#### 1. 创建igraph 对象

使用管道

``` r
library(igraph)
library(igraphdata)
library(magrittr)
library(tidyverse)
library(ggraph)
library(ggnetwork)


# Notable graphs
# make_graph can create some notable graphs. The name of the graph (case insensitive), a character scalar must be suppliced as the edges argument, and other arguments are ignored. (A warning is given is they are specified.)
# eg.
# Cubical
# The Platonic graph of the cube. A convex regular polyhedron with 8 vertices and 12 edges.
g <- make_graph("Cubical") %>%
  set_vertex_attr("name",value = LETTERS[1:4])

g  %>%
  add_layout_(with_fr()) %>%
  plot()
```

<img src="/figure/posts/Learnigraph-2-MoreOnIgraph_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

补充：更为实际的案例中，需要使用数据集来创建图。**igraph**作者提供了一些根据数据集创建好的igraph对象：

``` r
library(igraphdata)
### data:Loads specified data sets, or list the available data sets.
data(package="igraphdata")
# Data sets in package ‘igraphdata’:
# 
# Koenigsberg                 Bridges of Koenigsberg from Euler's times
# UKfaculty                   Friendship network of a UK university faculty
# USairports                  US airport network, 2010 December
# enron                       Enron Email Network
# foodwebs                    A collection of food webs
# immuno                      Immunoglobulin interaction network
# karate                      Zachary's karate club network
# kite                        Krackhardt's kite
# macaque                     Visuotactile brain areas and connections
# rfid                        Hospital encounter network data
# yeast                       Yeast protein interaction network
```

#### 2. 使用iraph对象查看边和点的信息

``` r
### already
data("macaque")
macaque
```

    ## IGRAPH f7130f3 DN-- 45 463 -- 
    ## + attr: Citation (g/c), Author (g/c), shape (v/c), name (v/c)
    ## + edges from f7130f3 (vertex names):
    ##  [1] V1 ->V2     V1 ->V3     V1 ->V3A    V1 ->V4     V1 ->V4t    V1 ->MT    
    ##  [7] V1 ->PO     V1 ->PIP    V2 ->V1     V2 ->V3     V2 ->V3A    V2 ->V4    
    ## [13] V2 ->V4t    V2 ->VOT    V2 ->VP     V2 ->MT     V2 ->MSTd/p V2 ->MSTl  
    ## [19] V2 ->PO     V2 ->PIP    V2 ->VIP    V2 ->FST    V2 ->FEF    V3 ->V1    
    ## [25] V3 ->V2     V3 ->V3A    V3 ->V4     V3 ->V4t    V3 ->MT     V3 ->MSTd/p
    ## [31] V3 ->PO     V3 ->LIP    V3 ->PIP    V3 ->VIP    V3 ->FST    V3 ->TF    
    ## [37] V3 ->FEF    V3A->V1     V3A->V2     V3A->V3     V3A->V4     V3A->VP    
    ## [43] V3A->MT     V3A->MSTd/p V3A->MSTl   V3A->PO     V3A->LIP    V3A->DP    
    ## + ... omitted several edges

原作者是这么解释的：

``` 

This is the standard way of showing (printing) an igraph graph object on
the screen. The top line of the output declares that the object is an igraph
graph, and also lists its most important properties. A four-character long
code is printed first:

‘D/U’ The first character is either ‘D’ or ‘U’ and encodes whether the graph
is directed or undireted.
‘N’ The second letter is ‘N’ for named graphs (see Section 1.2.5). A dash
here means that the graph is not named.
‘W’ The third letter is ‘W’ if the graph is weighted (in other words, if the
graph is a valued graph, Section 2.4). Unweighted graphs have a dash in
this position.
‘B’ Finally, the fourth is ‘B’ if the graph is bipartite (two-mode, Section ??).
For unipartite (one-mode) graphs a dash is printed here.

This notation might seem quite dense at first, but it is easy to get used to and
conveys much information in a small space. Then two numbers are printed,
these are the number of vertices and the number of edges in the graph, 45
and 463 in our case. At the end of the line the name of the graph is printed,
if there is any. The next line(s) list attributes, meta-data that belong to the
vertices, edges or the graph itself. Finally, the edges of the graph are listed.
Except for very small graphs, this list is truncated, so that it fits to the screen.
```

一些基本量的展示，之前讲过，此外，还有更多关于边的操作：

``` r
###|V|
gorder(macaque)
###[1] 45
###|E|
gsize(macaque)
###[1] 463

V(macaque)
# + 45/45 vertices, named, from f7130f3:
#  [1] V1     V2     V3     V3A    V4     V4t    VOT    VP     MT     MSTd/p MSTl  
# [12] PO     LIP    PIP    VIP    DP     7a     FST    PITd   PITv   CITd   CITv  
# [23] AITd   AITv   STPp   STPa   TF     TH     FEF    46     3a     3b     1     
# [34] 2      5      Ri     SII    7b     4      6      SMA    Ig     Id     35    
# [45] 36    
E(macaque)
# + 463/463 edges from f7130f3 (vertex names):
#  [1] V1 ->V2     V1 ->V3     V1 ->V3A    V1 ->V4     V1 ->V4t    V1 ->MT    
#  [7] V1 ->PO     V1 ->PIP    V2 ->V1     V2 ->V3     V2 ->V3A    V2 ->V4    
# [13] V2 ->V4t    V2 ->VOT    V2 ->VP     V2 ->MT     V2 ->MSTd/p V2 ->MSTl  
# [19] V2 ->PO     V2 ->PIP    V2 ->VIP    V2 ->FST    V2 ->FEF    V3 ->V1    
# [25] V3 ->V2     V3 ->V3A    V3 ->V4     V3 ->V4t    V3 ->MT     V3 ->MSTd/p
# [31] V3 ->PO     V3 ->LIP    V3 ->PIP    V3 ->VIP    V3 ->FST    V3 ->TF    
# [37] V3 ->FEF    V3A->V1     V3A->V2     V3A->V3     V3A->V4     V3A->VP    
# [43] V3A->MT     V3A->MSTd/p V3A->MSTl   V3A->PO     V3A->LIP    V3A->DP    
# [49] V3A->FST    V3A->FEF    V4 ->V1     V4 ->V2     V4 ->V3     V4 ->V3A   
# [55] V4 ->V4t    V4 ->VOT    V4 ->VP     V4 ->MT     V4 ->LIP    V4 ->PIP   
# + ... omitted several edges

macaque %>% ends("V1|V2")
# 
#      [,1] [,2]
# [1,] "V1" "V2"

macaque %>% tail_of("V1|V2")
# + 1/45 vertex, named, from f7130f3:
# [1] V1

macaque %>% head_of("V1|V2")
# + 1/45 vertex, named, from f7130f3:
# [1] V2

macaque %>% neighbors("V1",mode = "out")

# + 8/45 vertices, named, from f7130f3:
# [1] V2  V3  V3A V4  V4t MT  PO  PIP

macaque %>% neighbors("V1",mode = "in")

# + 8/45 vertices, named, from f7130f3:
# [1] V2  V3  V3A V4  V4t MT  PO  PIP

E(macaque)[.from("V1")]
```

#### 3\. 子图

创建子图

``` r
V(macaque)["V1","V2",.nei("V1"),.nei("V2")] %>%
  induced_subgraph(graph = macaque) %>%
  summary()
```

    ## IGRAPH cb88d15 DN-- 16 156 -- 
    ## + attr: Citation (g/c), Author (g/c), shape (v/c), name (v/c)

连通

``` r
is_connected(macaque,mode = "weak")
```

    ## [1] TRUE

``` r
is_connected(macaque,mode = "strong")
```

    ## [1] TRUE

边和点的筛选：

``` r
V(macaque)[1:4]
# + 4/45 vertices, named, from f7130f3:
# [1] V1  V2  V3  V3A
V(macaque)[c("V1","V2","V3","V3A")]
# + 4/45 vertices, named, from f7130f3:
# [1] V1  V2  V3  V3
```

建立边或者点的索引向量：

``` r
E(macaque)[1:10] %>% as_ids()
# [1] "V1|V2"  "V1|V3"  "V1|V3A" "V1|V4"  "V1|V4t" "V1|MT"  "V1|PO"  "V1|PIP"
#  [9] "V2|V1"  "V2|V3" 
V(macaque)[1:10] %>% as_ids()
 # [1] "V1"     "V2"     "V3"     "V3A"    "V4"     "V4t"    "VOT"    "VP"    
 # [9] "MT"     "MSTd/p"
```

类似于算数操作，关于点的操作汇总：

``` r
data("kite")
V(kite)
# + 10/10 vertices, named, from 6b7ddad:
#  [1] A B C D E F G H I J
V(kite)[1:3,7:10]
# + 7/10 vertices, named, from 6b7ddad:
# [1] A B C G H I J
V(kite)[degree(kite) < 2]
# + 1/10 vertex, named, from 6b7ddad:
# [1] J
V(kite)[.nei("D")]
# + 6/10 vertices, named, from 6b7ddad:
# [1] A B C E F G
V(kite)[.innei("D")]
# + 6/10 vertices, named, from 6b7ddad:
# [1] A B C E F G
V(kite)[.outnei("D")]
# + 6/10 vertices, named, from 6b7ddad:
# [1] A B C E F G
V(kite)[.inc("A|D")]
# + 2/10 vertices, named, from 6b7ddad:
# [1] A D
c(V(kite)["A"],V(kite)["D"])
# + 2/10 vertices, named, from 6b7ddad:
# [1] A D
rev(V(kite))
# + 10/10 vertices, named, from 6b7ddad:
#  [1] J I H G F E D C B A
unique(V(kite)["A","A","C","C"])
# + 2/10 vertices, named, from 6b7ddad:
# [1] A C

### Set operation

union(V(kite)[1:5],v(kite)[6:10])
# + 2/10 vertices, named, from 6b7ddad:
# [1] A C
intersection(V(kite)[1:7],V(kite)[5:10])
# + 3/10 vertices, named, from 6b7ddad:
# [1] E F G
difference(V(kite),V(kite)[1:5])
# + 5/10 vertices, named, from 6b7ddad:
# [1] F G H I J


E(kite)
# + 18/18 edges from 6b7ddad (vertex names):
#  [1] A--B A--C A--D A--F B--D B--E B--G C--D C--F D--E D--F D--G E--G F--G F--H G--H
# [17] H--I I--J

E(kite,path = c("A","D","C"))
# + 2/18 edges from 6b7ddad (vertex names):
# [1] A--D C--D

E(kite)[ V(kite)[1:2] %--%  V(kite)[3:4] ]
# + 3/18 edges from 6b7ddad (vertex names):
# [1] A--C A--D B--D

E(kite)[1:3,7:10]
# + 7/18 edges from 6b7ddad (vertex names):
# [1] A--B A--C A--D B--G C--D C--F D--E

E(kite)[seq_len(gsize(kite))[seq_len(gsize(kite)) %%2 == 0]]
# + 9/18 edges from 6b7ddad (vertex names):
# [1] A--C A--F B--E C--D D--E D--G F--G G--H I--J

E(kite)[seq_len(gsize(kite)) %%2 == 0]
# + 9/18 edges from 6b7ddad (vertex names):
# [1] A--C A--F B--E C--D D--E D--G F--G G--H I--J

E(kite)[seq_len(gsize(kite)) %%2]

# + 9/18 edges from 6b7ddad (vertex names):
# [1] A--B A--B A--B A--B A--B A--B A--B A--B A--B

E(kite)[.inc("D")]
# + 6/18 edges from 6b7ddad (vertex names):
# [1] A--D B--D C--D D--E D--F D--G

E(macaque)[.from("V1")]
# + 8/463 edges from f7130f3 (vertex names):
# [1] V1->V2  V1->V3  V1->V3A V1->V4  V1->V4t V1->MT  V1->PO  V1->PIP
E(macaque)[.to("V1")]
# + 8/463 edges from f7130f3 (vertex names):
# [1] V2 ->V1 V3 ->V1 V3A->V1 V4 ->V1 V4t->V1 MT ->V1 PO ->V1 PIP->V1

### The remains are same as Vertices operations
```
