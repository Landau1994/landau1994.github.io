---
title: AlphaGenome 简介
author: 研精极锐
avatar: /images/faceicon.png
authorLink: https://github.com/Landau1994
authorAbout: https://github.com/Landau1994
authorDesc: An apprentice in bioinformatics
mathjax: true
categories:
  - genomics
tags:
  - note
date: 2026-02-25 13:43:12
keywords: 
 - AlphageGenome
 - Unet
description:
photos:
---



# Nature 2026 Advancing regulatory variant effect prediction with AlphaGenome

## Quick Summary

> AlphaGenome 是由 Google DeepMind 开发并发表于 _Nature_ 的一种统一的深度学习模型，旨在解决基因组序列建模中“长距离上下文”与“高分辨率预测”之间的权衡难题。该模型接受 **1 Mb** 的 DNA 序列作为输入，能够以 **单碱基分辨率** 同时预测 5,930 个（人类）和 1,128 个（小鼠）功能基因组轨道，涵盖基因表达、剪接（位点、使用率及连接）、染色质可及性、组蛋白修饰、转录因子结合以及 3D 染色体接触图谱。在 26 项变异效应预测（VEP）基准测试中，AlphaGenome 在 **25 项** 上达到或超过了现有的最佳模型（SOTA），并展示了在解释罕见病和癌症（如 T-ALL 中的 _TAL1_ 增强子突变）非编码区变异机制方面的强大能力。

## Key Points

- **统一的全能架构**：在一个模型中同时实现了长序列建模（1 Mb 上下文）和单碱基输出精度，打破了以往模型（如 Enformer 牺牲分辨率，SpliceAI 牺牲长度）的局限。
    
- **SOTA 性能**：在 22/24 个基因组轨道预测任务和 25/26 个变异效应预测任务中表现优于现有最佳模型（包括 Borzoi, Orca, Pangolin 等）。
    
- **剪接预测的突破**：不仅预测剪接位点，还引入了直接预测 **剪接连接（Splice Junctions）** 的机制，能更准确地识别外显子跳跃等复杂剪接变异。
    
- **高效蒸馏**：通过“预训练+蒸馏”策略，将大型集成模型的知识压缩到单个学生模型中，在单个 GPU 上推理时间小于 1 秒，大幅降低了应用门槛。
    

## Methods

### Data

- **训练数据**：来自 ENCODE, GTEx, FANTOM5, 4D Nucleome 等项目的公开数据。
    
- **涵盖模态**：RNA-seq, CAGE, PRO-cap, DNase-seq, ATAC-seq, ChIP-seq (TF & Histone), Hi-C/Micro-C。
    
- **物种**：人类 (hg38) 和 小鼠 (mm10)。
    

### Model Architecture

- **主干**：U-Net 风格的编码器-解码器架构 。
    
    - **Encoder**：卷积层将 1 Mb 序列下采样。
        
    - **Transformer Tower**：处理长距离依赖关系，并生成用于预测接触图谱的 2D 嵌入。
        
    - **Decoder**：通过跳跃连接（Skip connections）将特征上采样回 1 bp 分辨率。
        
- **输出头**：针对不同模态在不同分辨率（1 bp, 128 bp, 2048 bp）输出预测结果。特别设计了能够捕捉供体-受体相互作用的剪接连接预测头 。
    

### Training Strategy

- **两阶段训练**：
    
    1. **预训练 (Pre-training)**：在 8 个 TPUv3 核心上进行序列并行训练，学习预测实验轨道数据 。
        
    2. **蒸馏 (Distillation)**：训练一个学生模型来模仿“全折叠（all-fold）”教师模型集成的预测。此阶段引入了更强的增强策略（如随机突变），显著提升了模型的鲁棒性和变异效应预测能力 。
        

## Results

|**Metric**|**Value (AlphaGenome)**|**Baseline (SOTA)**|**Description**|
|---|---|---|---|
|VEP Benchmarks|**25/26**|-|在 26 项变异效应预测基准中，25 项达到或超越 SOTA|
|Gene Expr LFC (Pearson r)|**+14.7%**|(vs Borzoi)|细胞类型特异性基因表达预测的相对提升|
|Contact Map (Pearson r)|**+6.3%**|(vs Orca)|3D 基因组接触图谱预测的相对提升|
|eQTL Sign Prediction (auROC)|**0.80**|0.75 (Borzoi)|预测 eQTL 效应方向的准确性|
|ClinVar (Deep Intronic)|**0.66**|0.64 (Pangolin)|深层内含子致病变异的分类精度 (auPRC)|

## Figures

|**Figure**|**Description**|
|---|---|
|Fig 1|模型架构、训练流程（预训练与蒸馏）概览，以及与 SOTA 模型在各项基准上的性能对比摘要 。|
|Fig 2|展示了模型在 LDLR 等基因座上的高精度轨道预测（包括剪接连接），以及各模态预测值与观测值的高相关性 。|
|Fig 3|剪接变异预测深入分析：展示了模型如何准确预测导致外显子跳跃的罕见变异，并在 ClinVar 基准测试中领先 。|
|Fig 4|基因表达变异预测（eQTL）：展示了在 GTEx eQTL 效应值和方向预测上的显著提升，以及在 GWAS 信号解读中的应用 。|
|Fig 5|染色质状态变异预测（caQTL/bQTL）：展示了对染色质可及性和转录因子结合变异的精准预测 。|
|Fig 6|跨模态案例分析：解析 T-ALL 中 TAL1 癌基因的非编码突变机制，展示模型如何通过多模态输出揭示致病机理 。|

## Critical Analysis

### Strengths

- **多模态整合**：成功将此前分裂的多个领域（如剪接、表达、3D结构）整合到一个统一框架中，且在各子任务上均不输于专用模型。
    
- **分辨率优势**：相比 Borzoi/Enformer，提供了单碱基分辨率的输出，这对于识别精细的剪接位点和 TF 结合位点至关重要。
    
- **机制可解释性**：不仅能预测“变异致病”，还能通过多模态输出（如“该变异破坏了 CTCF 结合并改变了 3D 结构从而影响表达”）提供机制解释 。
    

### Weaknesses

- **远端调控局限**：尽管输入长达 1 Mb，但在预测距离超过 100 kb 的远端调控元件影响时，性能仍有下降，仍是一个挑战 。
    
- **组织特异性**：虽然表现优异，但在跨细胞背景下精确复现某些组织特异性模式方面仍有提升空间 。
    
- **非编码基因覆盖**：目前的训练和评估主要集中在蛋白编码基因，对 microRNAs 等非编码基因的覆盖不足 。
    

### Questions

- 如何进一步扩展模型以处理更长范围（>1 Mb）的相互作用，以捕捉超长距离的增强子-启动子互作？
    
- 未来的版本是否会整合单细胞数据，以提高在复杂组织中细胞类型特异性的分辨率？
    

## Connections

### Related Papers

- **Borzoi (Linder et al., 2025)**: 上一代多模态 SOTA，AlphaGenome 的主要对比基线。
    
- **Enformer (Avsec et al., 2021)**: 本文第一作者之前的开创性工作，引入了 Transformer 处理长序列。
    
- **Orca (Zhou, 2022)**: 3D 基因组预测的专用模型，被 AlphaGenome 在接触图谱任务上超越。
    
- **Pangolin / SpliceAI**: 剪接预测领域的标杆，AlphaGenome 在剪接任务上的对比对象。
    

### Related Concepts

- **Sequence-to-Function**: 从序列直接预测功能的建模范式。
    
- **In Silico Mutagenesis (ISM)**: 计算机模拟诱变，用于解析模型预测背后的序列基序（Motif）。
    
- **Knowledge Distillation**: 知识蒸馏，本文用于提升模型推理效率和鲁棒性的关键技术。
    

### Potential Applications

- **全基因组关联分析 (GWAS) 解读**：为非编码区 GWAS 信号提供因果变异和分子机制的假设（如确定 eQTL 的方向）。
    
- **临床变异解读**：辅助诊断罕见病，特别是针对那些现有工具难以解释的深层内含子变异或非编码调控变异。
    
- **序列设计**：用于设计具有特定组织特异性的增强子或优化反义寡核苷酸（ASO）疗法 。
    

## Notes

- 这是 AlphaGenome 的正式发表版本（Nature 2026），与之前的 bioRxiv 版本相比，基准测试结果更加完善（如 VEP SOTA 数从 24/26 更新为 25/26）。
    
- 论文明确指出，该模型并未在“个人基因组预测（personal genome prediction）”任务上进行基准测试，这是该领域模型的一个已知弱点 。