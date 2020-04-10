---
title: Rmarkdown_test
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
date: 2020-04-09
keywords:
description: Add rmarkdown function test
photos:
output: 
  md_document:
    variant: gfm
---

This post test blogdown, reference this [repo](https://github.com/yihui/blogdown-hexo)
This post generate by *blogdown::new_post(title = "Rmarkdown_test",ext=".Rmd")*







# R Markdown

This is an R Markdown document. Please note this page was **not** rendered using the [**rmarkdown**]( http://rmarkdown.rstudio.com) package or [Pandoc](http://pandoc.org). The R Markdown document is compiled to Markdown through **knitr**, and the Markdown document is rendered to HTML through [Hexo's Markdown renderer](https://github.com/hexojs/hexo-renderer-marked).

You can embed an R code chunk like this:


```r
summary(cars)
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
fit <- lm(dist ~ speed, data = cars)
fit
## 
## Call:
## lm(formula = dist ~ speed, data = cars)
## 
## Coefficients:
## (Intercept)        speed  
##     -17.579        3.932
```

# Including Plots

You can also embed R plots:


```r
par(mar = c(0, 1, 0, 1))
pie(
  c(280, 60, 20),
  c('Sky', 'Sunny side of pyramid', 'Shady side of pyramid'),
  col = c('#0292D8', '#F7EA39', '#C4B632'),
  init.angle = -50, border = NA
)
```

![](figure/pie-1.png)
