---
title: mmet 2020 Age-related loss of gene-to-gene transcriptional coordination among single-cells
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
date: 2020-12-06 11:07:29
keywords:
description:
photos:
---

  + 在HSC中没有观察到cell-cell variation随着衰老的改变；For example, despite observations of genetic and epigenteic damage in ageing haematopoietic stem cells(HSCs); cell to cell transcriptional variability is not observed.
	
  + 基因共表达方法的缺陷：
  
    + First, co-expression networks estimates the direct or indirect correlations between pairs of genes, while an individual gene may be controlled by multiple regulators.
  
    + Second, each co-expression measure is designed to capture a specific feature that is not necessarily optimal for depicting all types of gene-to-gene transcriptional interrelations ( PCC, linear relatiosnhips)
  
    + Third, large calculated coexpresion matrices contain a considerable amount of noise, which raises an additional difficulty in explporing their differentiation across cohorts
	
  + Gcl 本质上是对bcdcorr的bootstrap.
	
  + 一个重要的观察和建设：
	  As an illustration of the coordination measured by the GCL, consider the expression profiles of cells with N genes, that are represented as points in an N-dimensional space. If the gene expression levels are not independent, the set of points do not fill the N- dimensional space but are rather located near lower-dimensional manifold.
	
  + The GCL measures has two main advantages.
    + The dependency level is defined with respect to a general dependency form, not specific relations ( such as linear)
    + It can include high-order dependencies between multiple variables
