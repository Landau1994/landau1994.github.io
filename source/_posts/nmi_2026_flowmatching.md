---
title: NMI 2026 Flow matching for generative modelling in bioinformatics and computational biology
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
date: 2026-04-26 14:25:34
keywords:
description:
photos:
---


# Flow matching for generative modelling in bioinformatics and computational biology

> [!abstract] TL;DR (Executive Summary)
> **One-sentence summary**: 本文综述了流匹配（Flow Matching, FM）作为一种新兴的生成式人工智能范式，如何通过高效学习高维生物数据分布间的映射，推动分子建模、单细胞/多细胞分析及虚拟细胞（AI-based virtual cell）的发展。

## 1. Context & Rationale
- **Research Question**: 如何在生物信息学和计算生物学中高效地学习从一个生物状态（如患病状态）到另一个状态（如健康状态）的非平凡映射，并解决传统生成模型（如扩散模型）采样速度慢和约束处理难的问题？
- **The Gap**: 手动推导生物映射需要极高的专业知识；现有的生成模型（如 GANs, VAEs, Diffusion）在处理生物分子的几何约束（如 SE(3) 不变性）、离散序列数据以及大规模采样效率方面存在局限。
- **Hypothesis**: 流匹配（FM）作为一种基于最优传输（Optimal Transport）理论的无模拟（Simulation-free）训练范式，可以提供比扩散模型更简洁、更高效且更具可解释性的生物数据建模手段。

## 2. Methodology (Bio + Comp)

### 💻 Dry Lab (In Silico)
- **Algorithm/Model**: 
    - **核心理论**: 基于连续正则化流（CNFs），通过回归向量场（Vector Field）来学习概率路径（Probability Paths）。
    - **变体**: 条件流匹配（CFM）、修正流（Rectified Flow）、离散流匹配（Discrete FM，利用连续时间马尔可夫链）。
    - **架构**: 结合了 Transformer, 图神经网络（GNNs）和 SE(3) 等变网络（SE(3)-equivariant networks）。
- **Datasets**: 
    - 分子层面：PDB, GEOM-Drugs, Swiss-Prot, SAbDab。
    - 细胞层面：HLCA (Human Lung Cell Atlas), MERFISH, PBMC。
    - 影像层面：EMDB, fastMRI, JUMP。
- **Code Availability**: 
    - [Awesome Generative Flows (GitHub)](https://github.com/amorehead/awesome-generative-flows)
    - [Meta AI Flow Matching Guide](https://ai.meta.com/research/publications/flow-matching-guide-and-code)
    - [torchcfm](https://github.com/atong01/conditional-flow-matching)

## 3. Key Results

### Main Findings
1. **分子建模（Molecular Modelling）**: FM 在蛋白质折叠（如 FoldFlow）、小分子构象生成（如 FlowMol）中表现卓越，采样速度比扩散模型快 100 倍，且能更好地处理手性及几何约束。
2. **细胞建模（Cellular Modelling）**: FM 能够模拟单细胞在药物扰动下的表型演变（如 CellFlow），通过学习连续向量场预测细胞状态的轨迹。
3. **技术优势**: 相比扩散模型，FM 需要更少的推理步数（Inference Steps），实现更简单，且支持非高斯先验（Non-Gaussian priors），极大地增强了生物学合理性。
4. **虚拟细胞愿景**: FM 被视为构建“人工智能虚拟细胞”的关键脚手架，能够跨越分子、结构和表型尺度进行端到端的模拟。

### Key Figures
| Figure | Description & Takeaway |
| :--- | :--- |
| **Fig 1** | **FM 及其应用概览**: 展示了从细胞群体 A 到 B 的速度场映射，以及在蛋白质和分子建模中的生成策略。 |
| **Fig 2** | **发展时间线**: 描绘了从 2015 年正则化流起步，到 2022 年 FM 模型提出，再到 2024-2025 年在生物信息领域爆发式增长的过程。 |
| **Fig 3** | **应用分类法**: 将 FM 应用分为分子建模（基础）、单/多细胞建模（中间）和虚拟细胞建模（最高层）三个层次。 |

## 4. Critical Analysis (For Grants/Projects)

> [!success] Innovation 
> 1. **Conceptual**: 将复杂的生物过程建模为高维概率分布间的最优传输路径，提供了统一的数学框架。
> 2. **Technical**: 实现了“无模拟”训练，显著降低了计算成本；通过引入离散流匹配解决了生物序列（DNA/RNA/蛋白质）的生成难题。

> [!failure] Limitations & Critiques
> - **模式崩溃（Mode Collapse）**: 在特定路径选择下（如修正流），FM 仍可能面临稳定性挑战。
> - **数据稀疏性**: 对于极高维且稀疏的生物数据集，FM 的收敛性及表达能力仍需进一步理论验证。
> - **端点不稳定性**: 向量场在时间端点附近可能存在梯度方差不稳定的问题。

> [!tip] Relevance to My Research
> :
> - *Method to adapt*: 借鉴 CellFlow 的向量场学习方法，用于单细胞蛋白质组学轨迹推断。
> - *Idea to steal*: 利用最优传输路径（Optimal-transport paths）来减少生物样本模拟过程中的采样步数，提升实时分析效率。

## 5. Manuscript Snippets
- "Flow matching is a powerful and principled, data-driven framework for efficiently learning a mapping between arbitrary pairs of high-dimensional data distributions..." (描述 FM 的核心优势)
- "The rapid emergence of FM as a unifying paradigm for generative modelling marks a pivotal moment in the computational life sciences." (强调领域地位)

## 6. References
- **Related Papers**:
    - [[Lipman et al., 2023]] (Flow matching for generative modeling)
    - [[Tong et al., 2024]] (Improving and generalizing flow-based generative models with minibatch optimal transport)
- **Cited By**:
    - (该综述发表较新，预计将成为 2026 年后生成式生物 AI 的重要引文)