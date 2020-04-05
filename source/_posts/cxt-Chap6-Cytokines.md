---
title: cxt_Chap6_Cytokines
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
date: 2020-04-05 21:44:02
keywords:
description: Notes on cytokine
photos:
---

如无特别说明，引用部分出自《医学免疫学》（第七版，曹雪涛主编，人卫社出版）第六章，细胞因子。

### 1. 细胞因子的定义
> 细胞因子是由免疫细胞及组织细胞分泌的在细胞间发挥相互调控作用的的一类小分子可溶性蛋白质，通过结合响应受体调节细胞生长分化和效应，调控免疫应答，在一定条件下也参与炎症等多种疾病的发生。

>作用方式：自分泌方式，旁分泌方式，内分泌方式；

>功能特点：多效性，重叠性，协同性，拮抗性；

### 2. 细胞因子的种类
>根据结构和功能可以分为如下六大类：
- 白细胞介素(interleukin, IL)，IL1-IL38。
- 集落刺激因子(colony-stimulating factor, CSF)，是指能够刺激多能造血干细胞和不同分化阶段的造血祖细胞分化和增殖的细胞因子。主要包括粒细胞-巨噬细胞集落刺激因子(GM-CSF), 巨噬细胞集落刺激因子（M-CSF)，红细胞生成素(EPO)，干细胞因子(SCF)和血小板生成素(TPO)等。分别诱导造血干细胞或祖细胞分化增殖为相应的细胞。
- 干扰素（interferon, IFN), 因具有干扰病毒复制的功能而得名。IFN根据其结构特征及生物学活性可分为I型、II型和III型。I型IFN主要包括IFN-$\alpha$、IFN-$\beta$, 主要由病毒感染的细胞、pDC细胞等产生；II型IFN即IFN-$\gamma$，主要由活化T细胞和NK细胞产生。III型IFN包括IFN-$\lambda 1$(IL-29)，IFN-$\lambda 2$(IL-28A)和IFN-$\lambda 3$(IL-28B)，主要由DC细胞产生。IFN具有抗病毒、抗细胞增殖、抗肿瘤和免疫调节等作用。
-  肿瘤坏死因子(tumor necrosis factor, TNF)家族。肿瘤坏死因子因最初被发现其能造成肿瘤组织坏死而得名，包括TNF-$\alpha$和TNF-$\beta$，前者主要由活化的单核/巨噬细胞产生，后者主要由活化的T细胞产生，又称淋巴毒素(lymphotoxin, LT)。TNF家族目前已经发现TRAIL(TNF related apoptosis-inducing ligand)、FasL、CD40L等30余种细胞因子。TNF家族成员在调节免疫应答、杀伤靶细胞和诱导细胞凋亡等过程中发挥重要作用。
-  生长因子（growth factor，GF)泛指一类可促进相应细胞生长和分化的细胞因子。其种类较多，包括转化生长因子-$\beta$(transforming growth factor-$\beta$,TGF-$\beta$)、血管内皮细胞生长因子(VEGF)、表皮生长因子(EGF)、成纤维细胞生长因子(FGF)、神经生长因子(NGF)、血小板生长因子(PDGF)等。
-  趋化因子(chemokine) 是一类结构相似，分子量约8~12kD，具有趋化功能的细胞因子。几乎所有的趋化因子都含有由2对或一对保守的半胱氨酸残基(C)形成的分子内二硫化键。可以根据靠近氨基端的C的个数以及排列顺序将趋化因子分为四个家族：1）C亚家族：氨基酸端只有1个C，该分子内只有一个分子内二硫化键；2）CC亚家族：氨基端2个C相邻；3）CXC亚家族：氨基酸2个C被1个氨基酸残基隔开；4）CX3C亚家族：氨基端2个C被3个氨基酸残基隔开，羧基端跨细胞膜。
已经发现的趋化因子有，CXCL1~16，CCL1~28，XCL1~2，CX3CL1.

### 3. 细胞因子受体

> 细胞因子受体可以根据其结构特点被分为如下六个家族：
- I型细胞因子受体家族，也称为血细胞生辰素受体家族(hematopoietin receptor family)， 通过 JAK-STAT通路转导信号；
- II型细胞因子受体家族，也称为干扰素受体家族(interferon receptor family)，也是通过JAK-STAT通路转导信号；
- 肿瘤坏死因子受体家族(tumor necrosis factor family)，主要通过TRAF-NF-kB，TRAF-AP-1 通路转导信号；
- 免疫球蛋白超家族受体(Ig superfamily receptor, IgSFR)，会结合集落刺激因子；
- IL-17受体家族(IL-17 receptor family)，主要通过TRAF-NF-kB通路转导信号；
- 趋化因子受体家族(chemokine receptor family) ，属于GPCR中的一员。

### 4.remark

可否从基因表达调控出发，来描述细胞因子的功能？