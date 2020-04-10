---
title: pheatmap_advanced
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
  description: 热图的breaks, dendrogram
photo:
---



### Case 1

The datasets were provided by
[data-to-viz](https://www.data-to-viz.com/story/SevCatOneNumNestedOneObsPerGroup.html)

``` r
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(viridis)
library(kableExtra)
### dataset 1
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)
# show data
data %>% head(3) %>% select(1:3) %>% kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Africa

</th>

<th style="text-align:right;">

East.Asia

</th>

<th style="text-align:right;">

Europe

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Africa

</td>

<td style="text-align:right;">

3.142471

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

2.107883

</td>

</tr>

<tr>

<td style="text-align:left;">

East Asia

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

1.630997

</td>

<td style="text-align:right;">

0.601265

</td>

</tr>

<tr>

<td style="text-align:left;">

Europe

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

2.401476

</td>

</tr>

</tbody>

</table>

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

### case 2

the codes were adapted from
[slowkow](https://slowkow.com/notes/pheatmap-tutorial/) Sort dendrogram
is very important

``` r
set.seed(42)
random_string <- function(n) {
  substr(paste(sample(letters), collapse = ""), 1, n)
}

mat <- matrix(rgamma(1000, shape = 1) * 5, ncol = 50)

colnames(mat) <- paste(
  rep(1:3, each = ncol(mat) / 3),
  replicate(ncol(mat), random_string(5)),
  sep = ""
)
rownames(mat) <- replicate(nrow(mat), random_string(3))

mat %>% as.data.frame %>% head(3) %>% select(1:3) %>% kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

1jrqxa

</th>

<th style="text-align:right;">

1pskvw

</th>

<th style="text-align:right;">

1ojvwz

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

abv

</td>

<td style="text-align:right;">

9.6964789

</td>

<td style="text-align:right;">

9.172811

</td>

<td style="text-align:right;">

2.827695

</td>

</tr>

<tr>

<td style="text-align:left;">

nft

</td>

<td style="text-align:right;">

0.9020955

</td>

<td style="text-align:right;">

15.575853

</td>

<td style="text-align:right;">

4.328376

</td>

</tr>

<tr>

<td style="text-align:left;">

xha

</td>

<td style="text-align:right;">

2.6721643

</td>

<td style="text-align:right;">

3.127039

</td>

<td style="text-align:right;">

1.765077

</td>

</tr>

</tbody>

</table>

split data into 3 groups, and increase the values in group1

``` r
col_groups <- substr(colnames(mat), 1, 1)
mat[,col_groups == "1"] <- mat[,col_groups == "1"] * 5
```

making the heatmap

``` r
# install.packages("pheatmap", "RColorBrewer", "viridis")
library(pheatmap)
library(RColorBrewer)
library(viridis)

# Data frame with column annotations.
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(mat)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(col_groups)

pheatmap(
  mat               = mat,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Default Heatmap"
)
```

![](/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The default color breaks in pheatmap are uniformly distributed across
the range of the data.

We can see that values in group 1 are larger than values in groups 2 and
3. However, we can’t distinguish different values within groups 2 and 3.

``` r
## ----uniform-color-breaks------------------------------------------------

mat_breaks <- seq(min(mat), max(mat), length.out = 10)

dat <- data.frame(values = as.numeric(mat))

## ----uniform-color-breaks-detail, fig.height=2, echo=FALSE---------------
dat_colors <- data.frame(
  xmin = mat_breaks[1:(length(mat_breaks)-1)],
  xmax = mat_breaks[2:length(mat_breaks)],
  ymin = 0,
  ymax = max(density(mat, bw = "SJ")$y),
  fill = rev(inferno(length(mat_breaks) - 1)),
  stringsAsFactors = FALSE
)
ggplot() +
  geom_rect(
    data = dat_colors,
    mapping = aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
    )
  ) +
  geom_density(
    data = dat,
    mapping = aes(values),
    bw = "SJ", color = "cyan"
  ) +
  scale_fill_manual(values = dat_colors$fill) +
  cowplot::theme_cowplot()+
  theme(legend.position = "none") +
  labs(title = "Uniform breaks")
```

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

there are 6 data points greater than or equal to 100 are represented
with 4 different colors.

``` r
dat2 <- as.data.frame(table(cut(
  mat, mat_breaks
)))
dat2$fill <- inferno(nrow(dat2))
ggplot() +
  geom_bar(
    data = dat2,
    mapping = aes(x = Var1, weight = Freq, fill = Var1),
    color = "black", size = 0.1
  ) +
  coord_flip() +
  scale_fill_manual(values = dat2$fill) +
  cowplot::theme_cowplot()+
  theme(legend.position = "none") +
  labs(y = "data points", x = "breaks",
       title = "Number of data points per color")
```

![](/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

If we reposition the breaks at the quantiles of the data, then each
color will represent an equal proportion of the data:

``` r
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat, n = 11)
```

lets see

``` r
dat_colors <- data.frame(
  xmin = mat_breaks[1:(length(mat_breaks)-1)],
  xmax = mat_breaks[2:length(mat_breaks)],
  ymin = 0,
  ymax = max(density(mat, bw = "SJ")$y),
  fill = rev(inferno(length(mat_breaks) - 1)),
  stringsAsFactors = FALSE
)
ggplot() +
  geom_rect(
    data = dat_colors,
    mapping = aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
    )
  ) +
  geom_density(
    data = dat,
    mapping = aes(values),
    bw = "SJ", color = "cyan"
  ) +
  scale_fill_manual(values = dat_colors$fill) +
  theme(legend.position = "none") +
  labs(title = "Quantile breaks")
```

![](/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
dat2 <- as.data.frame(table(cut(
  mat, mat_breaks
)))
dat2$fill <- inferno(nrow(dat2))
ggplot() +
  geom_bar(
    data = dat2,
    mapping = aes(x = Var1, weight = Freq, fill = Var1),
    color = "black", size = 0.1
  ) +
  coord_flip() +
  scale_fill_manual(values = dat2$fill) +
  theme(legend.position = "none") +
  labs(y = "data points", x = "breaks",
       title = "Number of data points per color")
```

![](/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

When we use quantile breaks in the heatmap, we can clearly see that
group 1 values are much larger than values in groups 2 and 3, and we can
also distinguish different values within groups 2 and 3:

``` r
pheatmap(
  mat               = mat,
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Quantile Color Scale"
)
```

![](/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

We can also transform data

``` r
pheatmap(
  mat               = log10(mat),
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Log10 Transformed Values"
)
```

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

sort dendrograms

``` r
library(dendsort)

mat_cluster_cols <- hclust(dist(t(mat)))


sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
```

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

sort Dendrogram heatmap

``` r
mat_cluster_rows <- sort_hclust(hclust(dist(mat)))
pheatmap(
  mat               = mat,
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  cluster_cols      = mat_cluster_cols,
  cluster_rows      = mat_cluster_rows,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Sorted Dendrograms"
)
```

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

change colnames angle

``` r
pheatmap(
  mat               = mat,
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  cluster_cols      = mat_cluster_cols,
  cluster_rows      = mat_cluster_rows,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  angle_col = 90,
  fontsize_col  = 8,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 10,
  main              = "Sorted Dendrograms"
)
```

<img src="/figure/posts/pheatmap_advanced_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />
