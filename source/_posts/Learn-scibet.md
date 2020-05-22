---
title: Learn-SciBet
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - genomics
tags:
  - note
  - scRNA-seq
keywords: 细胞类型注释，迁移学习
description: 介绍SciBet原理以及使用方法
date: 2020-05-15 22:16:52
photos:
---


### 1. 背景介绍

细胞类型注释是scRNA-seq里非常基础的一步，常见的策略是将基于无监督分群结果的基因差异表达分析结果作为marker,结合先验知识，判断该分群为何种类型。

近年来，随着scRNA-seq的数据集的逐渐积累，也有很多研究组提出了一些基于已有分群结果进行有监督的细胞类型注释的方法，比如`scmap`，或者Seurat里的`TransferData`。此外，也有[评估相关方法的benchmark研究](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z)：

![](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-019-1795-z/MediaObjects/13059_2019_1795_Fig1_HTML.png)

[SciBet](http://scibet.cancer-pku.cn/)是张泽民老师组最近新发表的一个快捷（实测真的比Seurat的`TransferData`快，而且效果确实差不多），方便的进行有监督的细胞类型注释的方法。论文见[SciBet as a portable and fast single cell type identifier](https://www.nature.com/articles/s41467-020-15523-2)，相关报道见[Nature Communications | 张泽民课题组发表单细胞转录组数据快速注释新方法](https://mp.weixin.qq.com/s/5Nmvzyk3_t9-eOawHZ-uGQ)

### 2. 方法原理简介

如果读者只对软件使用感兴趣的本节可以略过。

根据上面提到的那篇[报道](https://mp.weixin.qq.com/s/5Nmvzyk3_t9-eOawHZ-uGQ)，如果用一句话概括SciBet的原理，应该是这样的：

> 张泽民实验室的博士生李辰威、刘宝琳联合任仙文副研究员开发的SciBet则有效地解决了这一问题：他们从“同一类型的单细胞表达谱服从同一多项分布”这一基本假设出发，对训练集数据中不同细胞类型分别进行建模，进而通过极大似然估计来对测试集细胞进行有监督注释。

接下来我们根据作者论文的方法部分和补充材料，学习一下这个巧妙地思路。分为如下几个部分：

+ 预备知识
+ 单细胞表达谱的统计模型
+ 根据表达谱计算信息熵
+ E-test
+ 构建有监督的细胞类型分类模型

我们假定读者已经修过概率论等相关课程。

#### 2.1 预备知识


预备知识1：

我们需要向读者回顾关于负二项分布的知识，如果一个离散随机变量$Y$的pmf(probability mass function)为:

$$ P(Y=y) = {y+r-1 \choose y}p^r(1-p)^{y},y=0,1,\dots,r=1,2,\dots$$

那么我们称其服从负二项分布。记为$Y \sim NB(r,p)$。

预备知识2：

对于两个随机变量$X$和$Y$，在给定$X=x$下，Y的条件pmf为：

$$P(Y=y|X=x)=\frac{P(X=x,Y=y)}{P(X=x)}$$

预备知识3：

关于负二项分布，有如下结论：

> 对于随机变量$Y,\Lambda$，若$Y|\Lambda \sim Poisson(\Lambda), \Lambda \sim Gamma(\alpha,\beta)$, 则 $Y \sim NB(r=\alpha,p=\frac{\beta}{1+\beta})$
其中
> $$ Poisson(Y=y;\Lambda=\lambda) = \frac{e^{-\lambda}\lambda^y}{y!},y \ge 0,  y\in N $$
> $$ Gamma(\Lambda=\lambda;a,b)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}\lambda^{\alpha-1}e^{-\beta\lambda},\lambda>0,\alpha >0,\beta >0 $$
> $$ \int_{0}^{\infty}\lambda^{\alpha-1}e^{-\beta\lambda}d\lambda=\frac{\Gamma(\alpha)}{\beta^{\alpha}} $$

证明：

因为：

$$\begin{aligned}
  P(Y=y) 
  &=\int_{0}^{\infty}P(Y=y,\Lambda=\lambda)d\lambda \\
  &=\int_{0}^{\infty}P(Y=y|\Lambda=\lambda)P(\Lambda=\lambda)d\lambda \\
  &=\int_{0}^{\infty}\frac{e^{-\lambda}\lambda^y}{y!}\frac{\beta^{\alpha}}{\Gamma(\alpha)}\lambda^{\alpha-1}e^{-\beta\lambda}d\lambda \\
  &=\frac{\beta^{\alpha}}{y!\Gamma(\alpha)}\int_{0}^{\infty}\lambda^{(y+\alpha)-1}e^{-(1+\beta)\lambda}d\lambda \\
  &= \frac{\beta^{\alpha}}{y!\Gamma(\alpha)}\frac{\Gamma(y+\alpha)}{(1+\beta)^{y+\alpha}}\\
  &={y+\alpha-1 \choose y}(\frac{\beta}{1+\beta})^{\alpha}(\frac{1}{1+\beta})^{y}
\end{aligned}
$$
所以 $Y \sim NB(r=\alpha,p=\frac{\beta}{1+\beta})$

根据预备知识3，我们有：

预备知识4：

> 负二项分布是泊松分布和伽马分布的混合分布。

预备知识5(该结论为概率论和数理统计的常识)：

> 对于随机变量的pdf为$f_X(x)$，做变换$Y=g(X)$,其中$g$为单调函数，X和Y各自样本空间的$\mathcal{X}$和$\mathcal{Y}$分别满足：
> $$\mathcal{X}=\{x|f_X(x)>0\},\mathcal{Y}=\{y|y=g(x),x\in\mathcal{X} \}$$
> 若$f_X(x)$在$\mathcal{X}$中连续且逆变换$g^{-1}(x)$在$\mathcal{Y}$中连续可微，则新的随机变量$Y$的pdf为：
> $$ f_Y(y) =\begin{cases} f_X(g^{-1}(y))|\frac{d}{dy}g^{-1}(y)|,y\in\mathcal{Y} \\ 0, otherwise \end{cases} $$ 

该结论为Statistical Inference(Casella Berger著)中的定理2.1.5，证明见该书。

预备知识6（只有这个结论作者在补充材料里给了证明）：

> **The scaling property of the Gamma distribution**
> 
> 若 $X \sim Gamma(\alpha,\beta)$, 则$Y=kX \sim Gamma(\alpha,\frac{\beta}{k}),k>0$

证明: 由题设：

$$ f_X(x;\alpha,\beta)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},\alpha >0,\beta >0 , x>0$$

根据题设，由预备知识5,有：

$$ \begin{aligned}f_Y(y) &= f_X(g^{-1}(y))|\frac{d}{dy}g^{-1}(y)| \\ &= f_X(y/k)|\frac{d}{dy}(y/k)| \\ &=\frac{1}{\Gamma(\alpha)\beta^{\alpha}}(\frac{y}{k})^{\alpha-1}e^{-\frac{\beta}{k}y}\frac{1}{k} \\ &= \frac{1}{\Gamma(\alpha)(k\beta)^{\alpha}}x^{\alpha-1}e^{-\frac{\beta}{k}y} \end{aligned} $$

故$Y \sim  Gamma(\alpha,\frac{\beta}{k})$

预备知识6：

> 若独立同分布样本$X_1,\dots,X_n$服从$X_i \sim Poisson(\lambda)$, 且观测值分别为$x_1,\dots,x_n$，则参数$\lambda$的矩估计量和极大似然估计量相等，均为$\bar{x}=\frac{1}{n}\sum_{i=1}^nx_i$

证明很简单，在很多教材也能找到，故从略。

预备知识7：

> 若独立同分布样本$X_1,\dots,X_n$服从$X_i \sim Gamma(\alpha,\beta)$, 且观测值分别为$x_1,\dots,x_n$，则参数$\alpha,\beta$的极大似然估计满足：
> $$ \hat{\beta}=\frac{\hat{\alpha}}{\bar{x}},\bar{x}=\frac{1}{n}\sum_{i=1}^nx_i $$



证明：

由题设：

$$f_{X_i}(x_i;\alpha,\beta)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x_{i}^{\alpha-1}e^{-\beta x_i}$$

故取对数后极大似然函数为：

$$ \begin{aligned}l(\alpha,\beta) &=\ln L(\alpha,\beta) \\ &= \ln \prod_{i=1}^{n} \frac{\beta^{\alpha}}{\Gamma(\alpha)}x_{i}^{\alpha-1}e^{-\beta x_i} \\ &=n(\alpha\ln\beta-\ln\Gamma(\alpha))+(\alpha-1)\sum_{i=1}^{n}\ln x_i - \beta\sum_{i=1}^{n}x_i \end{aligned}$$

由 
$$ \frac{\partial}{\partial\beta}l(\hat{\alpha},\beta) = 0$$

解得 $\hat{\beta}=\frac{\hat{\alpha}}{\bar{x}}$

预备知识8：

> 若随机变量$X$的pdf为$f_X(x)$,样本空间为$\mathcal{X}$, 定义`information generating function`:
> $$ T(u)=\int_{\mathcal{X}}f_{X}^{u}(x)dx $$
> 则：
> $$H(X)=-T^\prime(1)$$

证明：(严格证明的话，得要考虑函数$f_X(x)$和$\mathcal{X}$本身的性质，我们这里假定它们可以满足积分号下求导，如果要严格的话，请参考测度论相关教材)

$$ \begin{aligned}\frac{d}{du}T(u &)=\frac{d}{du}\int_{\mathcal{X}}f_{X}^{u}(x)dx \\ &=\int_{\mathcal{X}}\frac{d}{du}e^{u\ln f_X(x) }dx \\ &=\int_{\mathcal{X}}f_{X}^{u}(x)\ln f_{X}(x)dx  \end{aligned}$$

而随机变量$X$的信息熵的定义为：

$$ H(X)=-\int_{\mathcal{X}} f_{X}(x)\ln f_{X}(x)dx$$

故$H(X)=-T^\prime(1)$

预备知识9：

> 若$X\sim Gamma(\alpha,\beta)$，则其信息熵$S$为
> 
> $S=\alpha-\ln\beta+\ln\Gamma(\alpha)+(1-\alpha)\psi(\alpha),\psi=\frac{\Gamma^{\prime}(a)}{\Gamma(a)}$

由预备知识8经过计算即可证明，从略。

预备知识10：

> [多项分布](https://zhuanlan.zhihu.com/p/52481385)



注：更多关于伽马分布的知识可见：[理解Gamma分布、Beta分布与Dirichlet分布](https://zhuanlan.zhihu.com/p/37976562)


#### 2.2 单细胞表达谱的统计模型

经过大量的统计分析和后续的实验验证（相关证据可参考这篇[文献](https://www.nature.com/articles/nmeth.2930)）有这样一个经验性的结论：

**观察到单细胞基因表达的count(比如UMI count)的分布可以用负二项分布很好的拟合,且相同细胞类型的单细胞表达谱服从同一个分布。**

结合2.1中的预备知识，我们可以将单细胞基因表达的count表示为泊松分布的伽马分布的混合分布。所以作者参考了[SAVER](https://www.nature.com/articles/s41592-018-0033-z)可以进行如下建模：

假设我们现在有$C$个细胞，$m$个基因的原始表达谱数据，里面的数值为reads count或者是UMI count。如果我们记观察得到细胞c的某个基因i的UMI count为$Z_{ic},i=1,\dots,m,c=1,\dots,C$，那么对于同一类型的细胞而言，有：

$$ Z_{ic} \sim Poisson(\lambda_is_c), \lambda_i \sim Gamma(\alpha_i,\beta_i) $$

其中$\lambda_i$表示基因$i$在细胞$c$中的真实表达量，$s_c=\sum_{c}Z_{ic}$表示这个细胞中的UMI总数，与测序深度有关。而$\alpha_i,\beta_i$则是两个参数表征某个细胞类型中的基因的真实表达分布的参数。

接下来的可以根据预备知识里的结论进行参数估计：

由预备知识6，我们很容易得到$\lambda_{ic}=\frac{Z_ic}{s_c}$这也为我们常用的进行normalized的策略提供了一种依据。而由预备知识7，有$\hat{\beta_i}=\frac{\alpha_i}{E(\lambda_i)},E(\lambda_i)=\frac{1}{C}\sum_{c}\lambda_{ic}$。


#### 2.3 根据表达数据计算信息熵

不是所有的基因都是对后续的统计学习有用的，需要进行特征选择，也就是说，挑选出那些能表示不同群细胞之间表达差异的基因。本文的新意是基于信息熵(也就是香农熵)的概念引入了新的进行特征选择的方法：E-test。在讲E-test之前，我们需要看看作者是如何实现从表达量中计算信息熵的。

根据2.2我们知道可以将观测得到单细胞表达gene count$Z_{ic}$表示成$Z_{ic}$与真实表达量$\lambda_{ic}$的混合分布，真正能反映不同细胞之间表达差异的是$\lambda_{ic}$的分布。所以接下来要计算$\lambda_{ic}$的分布的信息熵。

对于相同类型的细胞而言，根据2.1中的结论9，真实表达量$\lambda_i$的信息熵为：

$$ S_i = \alpha_i - \ln\beta_i + \ln\Gamma(\alpha_i)+(1-\alpha_i)\psi(\alpha_i) \tag{1}$$

用$\hat{\beta_i}=\frac{\alpha_i}{E(\lambda_i)}$代入（1）可得：

$$ S_i = \alpha_i - \ln \alpha_i + \ln E(\lambda_i) + \ln\Gamma(\alpha_i)+(1-\alpha_i)\psi(\alpha_i) \tag{2} $$

记：
$$ h_i(\alpha_i) = \alpha_i - \ln \alpha_i  + \ln\Gamma(\alpha_i)+(1-\alpha_i)\psi(\alpha_i) \tag{3} $$


并且记$C$个细胞的平均normalized表达量为$X_i$, 显然有$X_i=E(\lambda_i)=\frac{1}{C}\sum_{c}\lambda_{ic}$

根据上述记号，（1）最终化简为：

$$ S_i = \ln E(\lambda_i)+h_i = \ln X_i + h_i \tag{4} $$

接下来，我们考虑道不同的细胞类型。若细胞$c,c=1,\dots,C_j$属于细胞类型$j$，定义所有属于$j$的细胞的平均标准化的表达量为$X_{ij}=\frac{1}{C_j}\sum_{c\in j}\lambda_{ic}$，根据上面的结论，可得：

$$ S_{ij} = \ln X_{ij} + h_{ij} \tag{5}$$

其中

$$ h_{ij}=\alpha_{ij}-\ln\alpha_{ij}+\ln\Gamma(\alpha_{ij})+(1-\alpha_{ij})\psi(\alpha_{ij}) \tag{6} $$

$X_{ij}$是直接可以通过实验数据计算的，而$h_{ij}$则需要估计。作者假设$h_{ij}$是只是基因特异的参数，也就是$h_{ij}=h_i$。理由如下：

若$\lambda_{i,c \in j}$为随机变量$\lambda_{ij}\sim\Gamma(\alpha_{ij},\beta_{ij})$的观测值，如果我们记基因$i$从细胞类型$j$到细胞类型$j^\prime$的表达量的fold change为$F_{i,j\rightarrow j^{\prime}}$, 并且假设$\lambda_{ij^{\prime}}=F_{i,j\rightarrow j^{\prime}}\lambda_{ij}$，那么由2.1中的预备知识6，可知 $\lambda_{ij^{\prime}}\sim\Gamma(\alpha_{ij},\frac{\beta_{ij}}{F_{i,j\rightarrow j^{\prime}}})$。
故$\alpha_{ij}=\alpha_i$, 结合（6），$h_{ij}=h_i$

最终，我们可以得到

 $$ S_{ij}=\ln X_{ij}+h_{i} \tag{7} $$

#### 2.4 E-test

首先，零假设为所有不同类型细胞都是从同一个细胞类型（记为 group 0）中均匀随机采样，那么基因的平均表达量为$X_{i0}=\frac{1}{n}\sum_{j=1}^nX_{ij}= AM_i$, 则group 0 的信息熵可以计算为：

$$ S_{i0} = \ln X_{i0}+h_i  \tag{8}$$

接着计算基因$i$在所有细胞类型$j$中与group 0 的信息熵的差之和：

(8)-(7)并求和，得：

$$ \begin{aligned}
  \Delta S_i &= \sum_{j=1}^{n}(S_{i0}-S_{ij}) \\
  & = \sum_{j=1}^n(\ln X_{i0}-\ln X_{ij}) \\
  & = n(\ln X_{i0}-\frac{1}{n}\ln \prod_{j=n}^nX_{ij}) \\
  & = n\ln\frac{AM_i}{GM_i}
\end{aligned} $$

其中$GM_i=(\prod_{j=n}^nX_{ij})^{\frac{1}{n}}$

利用Jesen不等式可以证明$GM_i\le AM_i$，故$\Delta S \ge 0$

若要进行假设检验，还需要计算$\Delta S$的显著性，作者的策略是基于置换检验的：

若所有预先定义分群的细胞类型的细胞均来自同一个样本，则对任意的标签为细胞类型$j$的细胞的size-factor normalized的表达量$X_{ij}=\frac{1}{C_j}\sum_{c\in j}\lambda_{ic}$，根据中心极限定理，单细胞数目足够多的时候，$X_{ij} \sim N(\mu_i,\sigma_i)$，很容易得到参数的无偏估计$\hat{\mu_i}=X_{i0},\hat{\sigma}_i=\frac{1}{n-1}\sum_{j=1}^n(X_{ij}-\hat{\mu_i})^2$，所以置换被简化成了每次从分布$N(\hat{\mu_i},\hat{\sigma_i})$n个不同的随机数$X_{i}$。接下来就是这么生成一个$\Delta S$的分布，然后计算这个分布中大于从真实数据中测得的$\Delta S$比例，作为显著性。

默认的Feature为500个基因。

#### 2.5 构建有监督学习的模型

根据2.4，在训练集中，我们可以进行特征选择选出一些“informative gene”进行模型训练。

作者假设从相同的转录出的mRNA是不可区分的，而且每个mRNA的产生是相互独立的，记录对于细胞类型为$j$的细胞，基因$i$产生一个mRNA的概率为$p_{ij}$,若有$m$个informative gene 则对细胞类型$j$,我们得到了一个随机向量$\boldsymbol{p}_j= (p_{1j},\cdots,p_{mj} )$其中$\sum_{i}p_{ij}=1$且其服从多项式分布,则 对于属于细胞类型$j$的细胞$y$，其后验表达谱$y=(y_1,\dots,y_n)$可以计算为：

$$ P(y|j)=\frac{(\sum_{i}y_i)!}{\prod_{i}(y_i!)}\prod_i(p_{ij}^{y_i}) $$

其中概率$p_{ij}$可估计为$\hat{p_{ij}}=\frac{1+X_{ij}}{\sum_{i}(1+X_{ij})}$

同样的，在测试集中，未知细胞类型的细胞$y$属于细胞类型$j$的概率也为

$$ P(y|j)=\frac{(\sum_{i}y_i)!}{\prod_{i}(y_i!)}\prod_i(p_{ij}^{y_i}) $$

其中$p_{ij}$是训练集中学习到的参数。

如果，我们引入$q_i=\frac{y_i}{\sum_i y_i}$,由极大似然的原则可知，细胞$y$最有可能的细胞类型$\hat{j}$为

$$ \begin{aligned}\hat{j} &= \argmax_{j}(P(y|j)) \\ &= \argmax_{j}(\frac{(\sum_{i}y_i)!}{\prod_{i}(y_i!)}\prod_{i}(p_{ij}^{y_i})) \\
&=\argmax_{j}(\prod_{i}(p_{ij}^{y_i}) \\
&=\argmax_{j}(y_i\ln p_{ij}) \\
&=\argmax_{j}(\frac{y_i}{\sum_{i}y_i}\ln p_{ij}-q_i\ln q_i) \\
&=\argmax_{j}(p_i\ln p_{ij}-q_i\ln q_i) \\
&=\argmax_{j}D_{KL}(q||p_j) \end{aligned}$$

#### 2.6 方法总结

可以用原文献中的流图对SciBet进行总结：

![](https://imgkr.cn-bj.ufileos.com/6c7712b3-9058-4a8d-a4ed-f32414601e4b.png)



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

按照官网的教程[E-test and SciBet](http://scibet.cancer-pku.cn/document.html)；下载所需的数据；然后后按照其说明文档进行即可。

接下来，我们看如何用使用SciBet结合Seurat进行单细胞分析。

##### 3.2.3 SciBet结合Seurat

为了说明问题，我们选择Seurat自带的pbmcsca数据集, 该数据集已经提供了预先定义好的细胞标签，代码如下：


```r
library(Seurat)
library(pbmcsca.SeuratData)
library(ggplot2)
library(scibet)
library(tidyverse)
library(viridis)

####0.---define utilized function--------
####---Export expr data from 10x to tibble----
#' @param seuv3 a seuv3 object 
#' 
myGetExpr <- function(seuv3,...){
	expr <- GetAssayData(object = seuv3, slot = "data")
	expr <- as_tibble(t(as.matrix(expr)),rownames = NA)
	return(expr)
}

###1.load data and preprocessing
data("pbmcsca")
###---avoid warning-----
pbmcsca <- UpdateSeuratObject(pbmcsca)
###---show predifined cell type------
table(pbmcsca$CellType)
```

```
# B cell              CD14+ monocyte 
# 5020                        5550 
# CD16+ monocyte                 CD4+ T cell 
# 804                        7391 
# Cytotoxic T cell              Dendritic cell 
# 9071                         433 
# Megakaryocyte         Natural killer cell 
# 977                        1565 
# Plasmacytoid dendritic cell                  Unassigned 
# 164                          46 
```


```r
###---split data-----
pbmc.list <- SplitObject(pbmcsca, split.by = "Method")
###---normalize-------
for (i in names(pbmc.list)) {
	pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
}
names(pbmc.list)
```

```
[1] "Smart-seq2"          "CEL-Seq2"            "10x Chromium (v2) A"
[4] "10x Chromium (v2) B" "10x Chromium (v3)"   "Drop-seq"           
[7] "Seq-Well"            "inDrops"             "10x Chromium (v2)"  
```
测试 reference base模式的效果

```r
###----2. test scibet----
###----define reference and query-----
reference <- pbmc.list[[1]]
query <- pbmc.list[[2]]


###----test query reference mode----
reference.expr <- myGetExpr(reference)
query.expr <- myGetExpr(query)


reference.label <- as.character(reference$CellType)
test.label <- as.character(query$CellType)

reference.expr <- cbind(reference.expr,label=reference.label)

prd.label <- SciBet(train = reference.expr, test = query.expr)
Confusion_heatmap(test.label, prd.label)
ggsave(filename = "res/fig/learn_scibet_confusionheatmap_refmode.pdf",
	   width = 6,height = 6)
```

![](https://imgkr.cn-bj.ufileos.com/a81d21b1-edf6-46a9-b46b-26afc649fc11.png)

准确率为

```r
num1 <- length(test.label)
num2 <- tibble(
	ori = test.label,
	prd = prd.label
) %>%
	dplyr::filter(ori == prd) %>%
	nrow(.)

num2/num1
```

```
0.851711
```

测试用作者提供的训练好的模型

```r
###----test load_model mode-----

###using 30 major cell types signature----

model <- readr::read_csv(file = "http://scibet.cancer-pku.cn/major_human_cell_types.csv")
model <- pro.core(model)

prd <- LoadModel(model)
prd.label <- prd(query.expr)

Confusion_heatmap(test.label,prd.label)
ggsave(filename = "res/fig/learn_scibet_confusionheatmap_signaturemode.pdf",
	   width = 6,height = 6)

```

![](https://imgkr.cn-bj.ufileos.com/ff705079-0b71-4114-8e2e-c8508fc02adb.png)