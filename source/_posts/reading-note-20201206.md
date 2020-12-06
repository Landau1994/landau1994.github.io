---
title: reading_note_20201206
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
date: 2020-12-06 20:49:41
keywords:
description:
photos:
---
	
  题目：[3D] bioRxiv 2020 3D Genome Contributes to Protein-Protein Interactome
  - 生物问题：3D 基因组的与蛋白质互作是否存在关系？
	- 实验设计：
		○ PPI数据：收集不同数据库蛋白质互作数据，作为正样本；从不同亚细胞定位的非PPI互作蛋白中抽样，作为负样本；
		○ HiC数据：不同细胞系的HiC数据，采用相同方法重新处理，使得互作矩阵分辨率相同；
		○ 从3D基因组信息中，重建基因的空间信息，采用不同的机器学习方法预测PPI数据，
	- 分析思路：
		○ 分组比较PPI gene counterparts 组和 Non-PPI组的分布图，HiCHeatmap，Gene-gene pair projections of PPIs overlaid on Hi-C heatmaps.，发现了PPI gene counterparts 在空间更为邻近；（Figure2,3,4）；
		○ 分析至少一个蛋白有信号肽的PPI（SigPep PPI）与无信号肽PPI(Non-SigPep PPI)的关系，得到结论是：“This can be explained that for the interacting proteins that are brought together by signal peptides, their gene counterparts can be more freely located on the 3D genome, with larger spatial distances“
		○ 用3维基因组的信息在不同的模型中预测PPI，发现引入3维基因组的信息之后，”
		the prediction accuracy in terms of AUC can be significantly improved if 3D genome information is employed“（Table1)，由于是预印本，特征工程的那一部分作者并没有详细写；
	- 评论：
		○ 目前有很多研究组致力于解析3D基因组结构和疾病的关系，报道了一些3D基因组结构改变，影响转录，从而影响疾病的案例；本文的分析虽然简单，但是给出了3D基因组结构对于更为下游的蛋白质互作有影响的可能性。但是这种可能性，还需要更为solid的技术，数据和分析方法来证明；
		○ alphafold是基于序列信息来预测蛋白质结构，目前已经取得重大突破。目前也陆续有基于机器学习的方法，整合不同基因组学的数据进行功能基因组学研究的报道。可能在未来的研究中，结合多组学数据，像alphafold这样的AI框架，才能实现更为深刻的蛋白动态功能的预测；
