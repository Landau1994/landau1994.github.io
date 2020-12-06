---
title: Learn-igraph-Basic
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
date: 2020-04-20 
keywords: 
description: Networkdata manipulation and Baisc visualization
photos:
---

*Learn-igraph*系列是对[Statistical Analysis of Network Data with
R](https://link.springer.com/book/10.1007/978-1-4939-0983-4)一书的学习笔记，介绍如何使用R进行网络数据分析，网络数据的处理主要是基于`igraph`包，可视化用的是`ggnet`

### 0. 基本概念

一些需要知道的基本概念；

  - Network;

  - Graph;

  - Order of a graph;

  - Size of a graph;

  - directed graph;

  - undirected graph;

  - subgraph;

### 1. 创建igraph class

#### 1.1 无向图

igraph包处理网络图的数据结构为igraph class, 最基础的创建方式如下：

``` r
library(igraph)
library(ggraph)
library(ggnetwork)
g <- graph.formula(1-2,1-3,2-3,2-4,3-5,4-5,4-6,4-7,5-6,6-7)
l <- layout.auto(g)
plot(g, layout=l, vertex.color="skyblue")
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

该网络的基本信息可以通过如下方式获得：

``` r
V(g)
###+ 7/7 vertices, named, from 27d8280:
###[1] 1 2 3 4 5 6 7
E(g)
###+ 10/10 edges from 27d8280 (vertex names):
###[1] 1--2 1--3 2--3 2--4 3--5 4--5 4--6 4--7 5--6 6--7
###str(g)
get.adjedgelist(g)

# $`1`
# + 2/10 edges from f3f6e64 (vertex names):
#   [1] 1--2 1--3
# 
# $`2`
# + 3/10 edges from f3f6e64 (vertex names):
#   [1] 1--2 2--3 2--4
# 
# $`3`
# + 3/10 edges from f3f6e64 (vertex names):
#   [1] 1--3 2--3 3--5
# 
# $`4`
# + 4/10 edges from f3f6e64 (vertex names):
#   [1] 2--4 4--5 4--6 4--7
# 
# $`5`
# + 3/10 edges from f3f6e64 (vertex names):
#   [1] 3--5 4--5 5--6
# 
# $`6`
# + 3/10 edges from f3f6e64 (vertex names):
#   [1] 4--6 5--6 6--7
# 
# $`7`
# + 2/10 edges from f3f6e64 (vertex names):
#   [1] 4--7 6--7

get.edgelist(g)
# [,1] [,2]
# [1,] "1"  "2" 
# [2,] "1"  "3" 
# [3,] "2"  "3" 
# [4,] "2"  "4" 
# [5,] "3"  "5" 
# [6,] "4"  "5" 
# [7,] "4"  "6" 
# [8,] "4"  "7" 
# [9,] "5"  "6" 
# [10,] "6"  "7" 
print(g, e=TRUE, v=TRUE)
# IGRAPH f673c51 UN-- 7 10 -- 
#   + attr: name (v/c)
# + edges from f673c51 (vertex names):
#   [1] 1--2 1--3 2--3 2--4 3--5 4--5 4--6 4--7 5--6 6--7
get.adjacency(g)
# 7 x 7 sparse Matrix of class "dgCMatrix"
# 1 2 3 4 5 6 7
# 1 . 1 1 . . . .
# 2 1 . 1 1 . . .
# 3 1 1 . . 1 . .
# 4 . 1 . . 1 1 1
# 5 . . 1 1 . 1 .
# 6 . . . 1 1 . 1
# 7 . . . 1 . 1 .
```

#### 1.2 有向图

同样的方法，也可以用来创建有向图；

``` r
dg <- graph.formula(1-+2,1-+3,2++3)
op <- par(mfrow=c(1,2))
plot(g, vertex.size=10,layout=l, vertex.color="skyblue")
plot(dg,vertex.size=10,vertex.color="skyblue")
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

``` r
par(op)
```

#### 1.3 从邻接矩阵导入图；

我们选择一个神奇的数据Arecibo\_message\[<https://en.wikipedia.org/wiki/Arecibo_message>\],
来说明,有时候,信息所对应的矩阵，可能就是一张图片，而不是一个图。

``` r
###python command comes from
###https://codegolf.stackexchange.com/questions/182924/output-the-arecibo-message
mat <- reticulate::py_eval("''.join(bin(i)[3:]for i in b'`UP@JB`IDQKJjjd`@@@@@L@@Ah@@CP@@J`@@_@@@@@LNLLP@FPtXpu}}}|@@@@`@@`@@@A@@A~@@~@@@CCCcDA@DMCGM____@@@@HF@H@L@@PX@_`pO`A`@HA@HHF@`LLB@FHX@@s@@Xa`CC@`HD@``L@b@XAD@PDDA@PD@C@F@X@ck@A@P@BCx@DKi[@gI\x7f\\NC\\@TGY@hOrAPXDFp@@@@@\\D@@zbjipAU@@B`@Gp@@\x7fx@G@\\@X@LAh@lFXCLHhJHQHdPBJH@DHP@H@`@Dh@OOix')[1:]")
mat <- as.integer(unlist(strsplit(mat,split = "")))
mat <- matrix(data = mat,nrow = 23,ncol = 73)

expand.matrix <- function(A){
  m <- nrow(A)
  n <- ncol(A)
  B <- matrix(0,nrow = m, ncol = m)
  C <- matrix(0,nrow = n, ncol = n)
  cbind(rbind(B,t(A)),rbind(A,C))
}
g1 <- graph_from_adjacency_matrix(expand.matrix(mat),mode = "undirected")
plot(g1,vertex.size=10,edge.width=2,layout=layout.circle,vertex.color="coral")
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

如果直接可视化这个图，我们什么也看不出来，然而，如果我们用将原数据视为栅格数据，那么，我们能看出这个数据的内涵是很丰富的

``` r
dat_long <- reshape2::melt(mat)
dat_long$value <- as.factor(dat_long$value)
colnames(dat_long) <- c("V1","V2","value")
### plot
gg <- ggplot(dat_long)+
  geom_tile(aes(V1,V2,fill=value), color="#7f7f7f")+
  scale_fill_manual(values=c("black", "white"))+
  coord_equal()+
  labs(x=NULL, y=NULL)+
  scale_x_continuous(breaks = 1:6)+
  scale_y_reverse(breaks=1:6)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
gg
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#### 1.4 从data.frame中创建图

需要两个输入，一个是边的信息，一个是节点的信息

``` r
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)


## The opposite operation
as_data_frame(g, what="vertices")
```

    ##                name age gender
    ## Alice         Alice  48      F
    ## Bob             Bob  33      M
    ## Cecil         Cecil  45      F
    ## David         David  34      M
    ## Esmeralda Esmeralda  21      F

``` r
as_data_frame(g, what="edges")
```

    ##        from    to same.dept friendship advice
    ## 1       Bob Alice     FALSE          4      4
    ## 2     Cecil   Bob     FALSE          5      5
    ## 3     Cecil Alice      TRUE          5      5
    ## 4     David Alice     FALSE          2      4
    ## 5     David   Bob     FALSE          1      2
    ## 6 Esmeralda Alice      TRUE          1      3

可视化，

``` r
plot(g,vertex.size=10,vertex.color="skyblue")
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

#### 1.5 用预定义的函数生成

`igraph`里有很多带make的函数，是可以生成图的

``` r
# ls.str and lsf.str return an object of class "ls_str", basically the character vector of matching names (functions only for lsf.str), similarly to ls, with a print() method that calls str() on each object.
###head(lsf.str("package:igraph"))
grep(pattern = "^make",x=ls("package:igraph"),value = T)
```

    ##  [1] "make_"                     "make_bipartite_graph"     
    ##  [3] "make_chordal_ring"         "make_clusters"            
    ##  [5] "make_de_bruijn_graph"      "make_directed_graph"      
    ##  [7] "make_ego_graph"            "make_empty_graph"         
    ##  [9] "make_full_bipartite_graph" "make_full_citation_graph" 
    ## [11] "make_full_graph"           "make_graph"               
    ## [13] "make_kautz_graph"          "make_lattice"             
    ## [15] "make_line_graph"           "make_ring"                
    ## [17] "make_star"                 "make_tree"                
    ## [19] "make_undirected_graph"

我们展示其中的一些图：

``` r
g1 <- make_tree(10, 2)
g2 <- make_bipartite_graph( rep(0:1,length=10), c(1:10))
g3 <- make_star(10, mode = "out")
g4 <- make_star(10, mode = "in")
op <- par(mfrow=c(2,2))
plot(g1,vertex.size=20,vertex.color="skyblue")
plot(g2,vertex.size=20,vertex.color="skyblue")
plot(g3,vertex.size=20,vertex.color="skyblue")
plot(g4,vertex.size=20,vertex.color="skyblue")
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
par(op)
```

### 2. 基本操作

诱导子图

``` r
g <- graph.formula(1-2,1-3,2-3,2-4,3-5,4-5,4-6,4-7,5-6,6-7)
h <- induced.subgraph(g,1:5)
print(h)
```

    ## IGRAPH d91ee38 UN-- 5 6 -- 
    ## + attr: name (v/c)
    ## + edges from d91ee38 (vertex names):
    ## [1] 1--2 1--3 2--3 2--4 3--5 4--5

Exclusion：

``` r
h <- g - vertices(c(6,7))
print(h)
```

    ## IGRAPH d923ec9 UN-- 5 6 -- 
    ## + attr: name (v/c)
    ## + edges from d923ec9 (vertex names):
    ## [1] 1--2 1--3 2--3 2--4 3--5 4--5

Inclusion:

``` r
h <- h + vertices(c(6,7))
g <- h + edges(c(4,6),c(4,7),c(5,6),c(6,7))
print(g)
```

    ## IGRAPH d928f5d UN-- 7 10 -- 
    ## + attr: name (v/c)
    ## + edges from d928f5d (vertex names):
    ##  [1] 1--2 1--3 2--3 2--4 3--5 4--5 4--6 4--7 5--6 6--7

union:

``` r
h1 <- h
h2 <- graph.formula(4-6,4-7,5-6,6-7)
g <- graph.union(h1,h2)
print(g)
```

    ## IGRAPH d92f82f UN-- 7 10 -- 
    ## + attr: name (v/c)
    ## + edges from d92f82f (vertex names):
    ##  [1] 6--7 5--6 4--7 4--6 4--5 3--5 2--4 2--3 1--3 1--2

### 3. 查看/添加/修改 属性

首先创建一个示例的图，

``` r
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
```

我们可以通过`$`运算符来查看，添加，修改属性

``` r
###check edge attribute
names(edge_attr(g))
###[1] "same.dept"  "friendship" "advice" 
###vertext
names(vertex_attr(g))
###[1] "name"   "age"    "gender"
###Vertex
# list.vertex.attributes(g)
# list.edge.attributes(g)
V(g)$name
###[1] "Alice"     "Bob"       "Cecil"     "David"     "Esmeralda"
edge_attr(g)$same.dept
###[1] FALSE FALSE  TRUE FALSE FALSE  TRUE
edge_attr(g)$friendship
###[1] 4 5 5 2 1 1
```

可视化如下：

``` r
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)

V(g)$gender <- plyr::revalue(x=V(g)$gender,
                            replace=c("F"="Female","M"="Male"))
V(g)$gender
```

    ## [1] "Female" "Male"   "Female" "Male"   "Female"

``` r
g$name <- "Toy Graph"
set.seed(42)
tmp.df <- layout.graphopt(g)
V(g)$color <- plyr::revalue(x=V(g)$gender,
                            replace=c("Female"="skyblue",
                                      "Male"="coral"))
plot(g,layout=tmp.df,vertex.size=20,
     vertex.color=V(g)$color,main="Toy Graph")
legend('right',legend=unique(V(g)$gender),pch=c(19,19),col = c("skyblue","coral"))
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

``` r
set.seed(42)
tmp.df <- layout.graphopt(g)
gg.net = ggnetwork(g,
                   arrow.gap = 0.05, 
                   layout = tmp.df)
ggplot(gg.net, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "black", 
               alpha = 0.5, curvature = 0,
               arrow = arrow(length = unit(6, "pt"), 
                             type = "closed")) +
    geom_nodes(aes(color = gender), size = 10) +
  geom_nodetext(aes(label = name))+
  scale_color_manual(values = c("skyblue","coral"))+
    ggtitle("Toy Graph")+
    theme_blank()
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

### 4. 更多关于图的概念和术语

#### 4.1 概念

下述概念不搬运书里的定义；忘记就查书。后面的章节会再用到这些概念，进行图的可视化与统计分析。

  - multi-graph

  - simple-graph:
    可以用`is.simple()`判定，可以用`simplify()`将`multi-graph`转换为`simple-graph`.

  - neighbors

  - degree: The degree of a vertex v defined as the number of edges
    incident on v;

  - in-degree

  - out-degree

  - walk

  - trails

  - circuit & cylce;

  - reachable

  - graph connected

  - component of a graph

  - strong connected

  - weak connected

  - distance/geodesic distance

  - diameter

#### 4.2 一些特殊的图

与第一节有重叠

  - complet graph

  - clique

  - regular graph

  - tree

  - forest

  - root

  - ancestor

  - descendant

  - parents, children

  - k-star

  - dirrected acyclic graph(DAG)

  - bipartite graph

<!-- end list -->

``` r
g.bip <- graph.formula(actor1:actor2:actor3,
                       movie1:movie2,
                       actor1:actor2 - movie1,
                       actor2:actor3 - movie2)

V(graph = g.bip)$type <- grepl(pattern = "^movie",V(graph = g.bip)$name)

V(g.bip)$category <- ifelse(V(graph = g.bip)$type,"Movie","Actor")
V(g.bip)$category
```

    ## [1] "Actor" "Actor" "Actor" "Movie" "Movie"

``` r
g <- g.bip
set.seed(42)
### using matrxi product to do layout rotate 3/2pi
tmp.df <- layout.bipartite(g) %*% matrix(data = c(0,-1,1,0),nrow = 2)

gg.net = ggnetwork(g,
                   arrow.gap = 0.05, 
                   layout = tmp.df)
head(gg.net)
```

    ##   x   y   name  type category      xend      yend
    ## 1 0 0.0 actor1 FALSE    Actor 0.9514929 0.2378732
    ## 2 0 0.5 actor2 FALSE    Actor 0.9514929 0.2621268
    ## 3 0 0.5 actor2 FALSE    Actor 0.9514929 0.7378732
    ## 4 0 1.0 actor3 FALSE    Actor 0.9514929 0.7621268
    ## 5 0 0.0 actor1 FALSE    Actor 0.0000000 0.0000000
    ## 6 0 0.5 actor2 FALSE    Actor 0.0000000 0.5000000

``` r
ggplot(gg.net, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "black", 
               alpha = 0.5, curvature = 0
               # ,arrow = arrow(length = unit(6, "pt"), 
               #               type = "closed")
               ) +
    geom_nodes(aes(color = category), size = 16) +
  geom_nodetext(aes(label = name))+
  scale_color_manual(values = c("skyblue","coral"))+
    ggtitle("bipartite graph example")+
    theme_blank()
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

`igraph`自带的例子：

``` r
# Random bipartite graph
inc <- matrix(sample(0:1, 50, replace = TRUE, prob=c(2,1)), 10, 5)
g <- graph_from_incidence_matrix(inc)
plot(g, layout = layout_as_bipartite,vertex.size=20,
     vertex.color=c("skyblue","coral")[V(g)$type+1])
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

### 附录：R配色

基本颜色：

``` r
#### code provided by
####http://bc.bojanorama.pl/2013/04/r-color-reference-sheet/
m <- matrix(1:660, 60, 11)
kol <- colors()[m]
#op <- par(mar=c(.1, .1, 2, .1))
image(1:11, 1:60, t(m), col=kol, axes=FALSE, ann=FALSE)
txtcol <- ifelse( apply(col2rgb(kol), 2, mean) < 70, "white", "black")
text( as.numeric(col(m)), as.numeric(row(m)), kol, cex=.8, col=txtcol)
mtext("grDevices::colors", 3, cex=2)
```

<img src="/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

调色版

``` r
RColorBrewer::display.brewer.all()
mtext("RColorBrewer", 3, cex=2)
```

![](/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

渐变色

``` r
library(RColorBrewer)
library(colorRamps)
library(viridis)
### manu
rdylbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
navy <- colorRampPalette(c("navy", "white", "firebrick3"))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
warmcold <- colorRampPalette(c(rev(cold(21)), warm(20)))


### add manu with package function

N <- 100 # ramp length
funnames <- rev(c("manu::rdylbu","manu::rdbu","manu::navy","manu::jet.colors","manu::warmcold",
              "viridis::viridis",
              "grDevices::rainbow", "grDevices::heat.colors",
              "grDevices::terrain.colors", "grDevices::topo.colors",
              "grDevices::cm.colors", 
              "colorRamps::blue2red",
              "colorRamps::blue2green", "colorRamps::green2red",
              "colorRamps::blue2yellow", "colorRamps::cyan2yellow",
              "colorRamps::magenta2green", "colorRamps::matlab.like",
              "colorRamps::matlab.like2", "colorRamps::primary.colors",
              "colorRamps::ygobb"))
spl <- strsplit(funnames, "::")
pkgs <- sapply(spl, "[", 1)
funs <- sapply(spl, "[", 2)
kolmat <- sapply(funs, do.call, list(N))
mat <- matrix( seq(1, length(kolmat)), nrow(kolmat), ncol(kolmat))


image(seq(1, nrow(mat)), seq(1, ncol(mat)), mat, col=kolmat,
      axes=FALSE, ann=FALSE)
text( nrow(mat)/2, seq(1, ncol(mat)), funnames)
mtext("Color Ramps function", 3, cex=2)
```

![](/figure/posts/Learnigraph-1-ManuNetworkDataBasic_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->
