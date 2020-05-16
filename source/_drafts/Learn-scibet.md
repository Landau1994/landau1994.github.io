---
title: Learn-SciBet
author: Landau1994
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - scRNA-seq
  - R
  - sc-seq
date: 2020-05-15 22:16:52
keywords: 细胞类型注释，迁移学习
description: 介绍SciBet原理以及使用方法
photos:
---

### 1. 背景介绍

细胞类型注释是scRNA-seq里非常基础的一步，常见的策略是将基于无监督分群结果的基因差异表达分析结果作为marker,结合先验知识，判断该分群为何种类型。

近年来，随着scRNA-seq的数据集的逐渐积累，也有很多研究组提出了一些基于已有分群结果进行有监督的细胞类型注释的方法，比如`scmap`，或者Seurat里的`TransferData`。此外，也有[评估相关方法的benchmark研究](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z)：

![](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-019-1795-z/MediaObjects/13059_2019_1795_Fig1_HTML.png)

[SciBet](http://scibet.cancer-pku.cn/)是张泽民老师组最近新发表的一个快捷（实测真的比Seurat的`TransferData`快，而且效果确实差不多），方便的进行有监督的细胞类型注释的方法。论文见[SciBet as a portable and fast single cell type identifier](https://www.nature.com/articles/s41467-020-15523-2)，相关报道见[Nature Communications | 张泽民课题组发表单细胞转录组数据快速注释新方法](https://mp.weixin.qq.com/s/5Nmvzyk3_t9-eOawHZ-uGQ)

### 2. 方法原理

### 3. 软件使用

#### 3.1 在线版

在线版使用请按照[官网](http://scibet.cancer-pku.cn/)的[Online classification](http://scibet.cancer-pku.cn/download_references.html)教程。

值得一提的是，作者提供了很从不同组织，不同实验条件的单细胞测序数据中训练好的Signature：

![](https://imgkr.cn-bj.ufileos.com/bf56b36c-df96-49bb-af81-2e9118331e96.png)

这些不同的signature可以供不同研究者结合自己的兴趣使用。

#### 3.2 本地版

##### 3.2.1 安装

```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("PaulingLiu/scibet")
```
如果出现错误，请看相关[issue](https://github.com/PaulingLiu/scibet/issues/1)

##### 3.2.2 作者提供的测试数据

按照官网的教程[E-test and SciBet](http://scibet.cancer-pku.cn/document.html)；下载所需的数据；

+ TEST data: 不同类型的T细胞表达谱(TPM)，作者所在的组之前发表的数据;
+ 


```R
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
```
