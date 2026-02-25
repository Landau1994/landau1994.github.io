---
title: 我花了三年时间，整理了这份生物信息学资源清单（2026最新版）
author: 研精极锐
avatar: /images/faceicon.png
authorLink: https://github.com/Landau1994
authorAbout: https://github.com/Landau1994
authorDesc: An apprentice in bioinformatics
mathjax: true
date: 2026-02-25
categories:
  - note
tags:
  - 生物信息学
  - 资源导航
  - AI工具
  - 学习路径
comments: true
---

## 写在前面

作为一名生物信息学 PhD，这三年来最常被师弟师妹们问到的问题就是：

> "学长，有什么好的生信学习资源推荐吗？"
> "AlphaFold 怎么用啊？"
> "哪个大模型适合写代码？"
> "R 和 Python 用什么 IDE 比较好？"

每次都要翻聊天记录找链接，实在太麻烦了。于是，我决定把这些年收藏的优质资源**系统性地整理出来**，做成一个完整的资源导航页面。

今天就把这份清单分享给大家。**全文 7000+ 字，建议先收藏再看！**

<!-- more -->

## 为什么要做这个清单？

说实话，生物信息学的学习曲线真的很陡峭：

- 生物学背景的同学，**编程是拦路虎**
- 计算机背景的同学，**生物知识是天书**
- 两边都懂一点的，**不知道从哪里开始系统学习**

市面上的资源虽然很多，但**良莠不齐**。我见过太多人因为找不到合适的资源，在入门阶段就放弃了。

所以这个清单的目标很明确：

✅ **只收录真正有价值的资源**
✅ **覆盖从入门到进阶的完整路径**
✅ **紧跟最新技术发展**（尤其是 AI 工具）

---

## 📋 清单概览

这份清单包含 **13 个主题板块**，**190+ 个精选资源**，涵盖：

| 板块 | 资源数 | 适合人群 |
|------|--------|----------|
| 🔬 基本资源 | 15+ | 所有人 |
| 📚 生物信息学资源 | 20+ | 入门&进阶 |
| 🏛️ 研究机构与大型项目 | 18+ | 科研工作者 |
| 🛠️ 工具与教程 | 12+ | 实践者 |
| 🧬 AI 生物学工具 | 25+ | 前沿研究者 |
| 🤖 AI 大语言模型 | 40+ | 效率提升 |
| 📊 可视化工具 | 5+ | 数据分析 |
| 💻 数据分析工具 | 25+ | 日常开发 |
| 🎓 课程与视频 | 5+ | 系统学习 |

完整清单在这里 👉 [**资源导航页面**](/links/)

下面我挑几个**最值得关注的板块**详细说说。

---

## 🧬 重点板块一：AI 生物学工具

这是 2024-2026 年变化最大的领域！**不夸张地说，AI 正在重新定义生物学研究。**

### 1. 蛋白质结构预测 - 已经是基操了

还记得 2020 年 AlphaFold 2 横空出世的震撼吗？现在：

- **[AlphaFold 3](https://alphafoldserver.com/)** - 最新版本，精度更高，速度更快
  - 数据库已经有 **2 亿+** 蛋白质结构
  - 基本上你想要的蛋白，都能找到预测结构

- **[RoseTTAFold](https://robetta.bakerlab.org/)** - Baker Lab 出品，开源实现
  - 在某些蛋白质家族上比 AlphaFold 还准

- **[ESMFold](https://esmatlas.com/)** - Meta 出品，**最快**
  - 宏基因组蛋白结构图谱达到 **6 亿+**
  - 适合大规模筛选

**实用建议**：现在做蛋白相关研究，第一步就是去 AlphaFold Database 查结构。实验验证之前先预测一下，能省好多时间和经费！

### 2. 蛋白质设计 - 从预测到创造

更激进的是，AI 现在不只能预测结构，还能**从头设计蛋白质**：

- **[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)** - 给定结构，生成序列
- **[RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)** - 扩散模型设计蛋白
- **[Chroma](https://generatebiomedicines.com/chroma)** - Generate Biomedicines 的商业化产品

这个方向太火了，已经有好几个 Biotech 公司拿到巨额融资。**未来的药物设计，可能不需要从自然界筛选，而是直接 AI 生成。**

### 3. 分子对接 - 药物发现的加速器

- **[DiffDock](https://github.com/gcorso/DiffDock)** - 扩散模型做分子对接
  - 比传统方法快 100 倍
  - 精度还更高

**我的观察**：传统的虚拟筛选管线正在被 AI 快速替代。如果你做药物发现，不学这些工具就真的落伍了。

---

## 🤖 重点板块二：AI 大语言模型

大模型不只是聊天工具，**用对了能让科研效率翻倍！**

### 1. 通用大模型 - 日常必备

**国际版：**
- **ChatGPT (GPT-5.3)** - 最强综合能力，但需要魔法上网
- **Claude 4** - **长文本之王**（可以处理整本书）
  - 我现在写论文都先让 Claude 帮忙润色
- **Gemini** - Google 出品，多模态能力强

**国内版：**
- **Kimi** - 长文本处理很不错，免费！
- **DeepSeek** - 推理能力强，适合做数学题和代码
- **豆包** - 字节出品，速度快

**实用技巧**：
- 写代码用 **Claude** 或 **DeepSeek**
- 读文献用 **Kimi**（一次能塞好几篇论文）
- 画图解用 **GPT-5.3** 或 **Gemini**

### 2. 生物医学专用模型 - 专业问题找专家

- **[BioGPT](https://github.com/microsoft/BioGPT)** - 微软训练的生物医学文本生成模型
- **[Med-PaLM](https://sites.research.google/med-palm/)** - Google 的医疗专用模型
- **[GeneGPT](https://github.com/ncbi/GeneGPT)** - NCBI 出品，专门回答基因问题

**举个例子**：我最近在研究一个不熟悉的基因，直接问 GeneGPT，它给出的解释比我自己查文献快多了。

### 3. AI 辅助工具 - 效率加倍器

- **[Cursor](https://cursor.sh/)** - AI 代码编辑器
  - **强烈推荐！** 写代码效率至少提升 50%
  - 直接在编辑器里问 AI，不用来回切换

- **[GitHub Copilot](https://github.com/features/copilot)** - AI 编程助手
  - 代码补全非常智能
  - 学生可以免费用

- **[Consensus](https://consensus.app/)** - AI 文献搜索
  - 输入问题，直接给你文献总结
  - 比自己翻 PubMed 快太多

- **[Elicit](https://elicit.org/)** - AI 研究助手
  - 帮你做系统性文献综述
  - 自动提取论文关键信息

**真心话**：如果你还在纯手工写代码、读文献，真的该试试这些工具了。**工具用得好，下班下得早！**

---

## 💻 重点板块三：数据分析工具链

工欲善其事，必先利其器。这个板块我花了很多心思整理。

### 1. IDE 选择 - 找到最适合你的

**R 用户：**
- **[RStudio](https://posit.co/)** - 老牌IDE，功能完善
- **[Positron](https://github.com/posit-dev/positron)** - Posit 新产品
  - **同时支持 Python 和 R**
  - 界面现代化，推荐尝鲜

**Python 用户：**
- **[VS Code](https://code.visualstudio.com/)** - 万能编辑器
  - 轻量、插件丰富
  - 我现在 80% 的代码都在 VS Code 写

- **[PyCharm](https://www.jetbrains.com/pycharm/)** - 专业 Python IDE
  - 功能强大，但比较重

- **[Jupyter Lab](https://jupyter.org/)** - 交互式分析
  - 数据探索必备
  - 出图特别方便

### 2. 包管理工具 - 依赖管理不再头疼

**Python 环境：**
- **[uv](https://github.com/astral-sh/uv)** - **新星工具！**
  - 用 Rust 写的，比 pip 快 **10-100 倍**
  - 强烈推荐，真的快到飞起

- **[Poetry](https://python-poetry.org/)** - 现代化依赖管理
  - 比 pip + requirements.txt 好用太多

- **[Conda](https://docs.conda.io/)** - 跨语言包管理
  - 生信领域的标配

**R 环境：**
- **[rig](https://github.com/r-lib/rig)** - R 版本管理
  - 多版本 R 切换很方便

- **[renv](https://rstudio.github.io/renv/)** - 项目环境管理
  - 保证项目的可复现性

**避坑经验**：很多人遇到的 "我电脑上能跑，你电脑上不能跑" 问题，90% 是因为环境管理不当。**学会用这些工具，能省无数时间！**

---

## 📚 重点板块四：学习路径推荐

资源虽好，但**如何系统学习**才是关键。根据我的经验：

### 新手路径（0-6 个月）

1. **编程基础**
   - Python：[Rosalind](http://rosalind.info/) - 生信编程练习，边玩边学
   - R：直接上手 Bioconductor，边用边学

2. **生信基础**
   - [生信技能树](http://www.biotrainee.com/) - 中文教程丰富
   - [The Carpentries](https://carpentries.org/) - 系统化培训

3. **实战项目**
   - 跟着 [NGS Analysis Tutorial](https://github.com/griffithlab/rnaseq_tutorial) 做一遍
   - 去 [Kaggle](https://www.kaggle.com/) 找生信竞赛练手

### 进阶路径（6-18 个月）

1. **深入特定方向**
   - 单细胞：从 [Seurat](https://satijalab.org/seurat/) 入手
   - 蛋白质：学 AlphaFold 和分子对接

2. **关注前沿动态**
   - 定期看 [bioRxiv](https://www.biorxiv.org/) 预印本
   - 关注几个大牛的 Twitter/X

3. **参与开源项目**
   - 去 GitHub 找感兴趣的项目贡献代码
   - 建立自己的作品集

### 高级路径（18 个月+）

1. **发表研究**
   - 找到细分领域的 gap
   - 开发自己的工具/算法

2. **建立影响力**
   - 写技术博客分享经验
   - 在社区里活跃

**核心建议**：**不要囤积资源！** 选 2-3 个最适合你的，深入学习，比收藏 100 个资源更有用。

---

## 🏛️ 值得关注的研究机构

如果你想了解最前沿的研究，关注这些机构准没错：

### 国际顶尖机构

- **[Wellcome Sanger Institute](https://www.sanger.ac.uk/)** - 基因组学研究先驱
  - 人类基因组计划的主力军

- **[Broad Institute](https://www.broadinstitute.org/)** - 哈佛-MIT 联合研究所
  - 单细胞测序技术的引领者

- **[Arc Institute](https://arcinstitute.org/)** - 新兴研究机构
  - Patrick Hsu 创立，专注前沿生物技术

### 大型科研计划

- **[Human Cell Atlas](https://www.humancellatlas.org/)** - 人类细胞图谱计划
  - 想了解单细胞领域，这是必看的

- **[Chan Zuckerberg Initiative](https://chanzuckerberg.com/)** - CZI
  - [CELLxGENE](https://cellxgene.cziscience.com/) 单细胞数据门户超好用

**实用建议**：这些机构经常有公开的讲座和培训，可以关注他们的 YouTube 频道。

---

## 💡 我的使用心得

### 1. 建立自己的工具箱

不要什么都想学，根据你的研究方向，**精选 10-15 个常用工具深度掌握**。

我自己的常用工具箱：
- IDE: VS Code + Jupyter Lab
- 包管理: uv + Conda
- 代码助手: Cursor + GitHub Copilot
- 文献管理: Zotero
- 大模型: Claude + Kimi
- 绘图: ggplot2 + BioRender

### 2. 善用 AI 提效

**能让 AI 做的，就不要手工做：**
- 读文献 → 用 Consensus 或 Elicit
- 写代码 → 用 Cursor 或 Copilot
- 画图 → 用 BioRender 或 GPT-4
- 润色文章 → 用 Claude

但记住：**AI 是助手，不是替代品。** 批判性思维永远不能丢。

### 3. 持续学习的节奏

**不要试图一次性学完所有东西！** 建议：

- 每周固定 1-2 个小时学新工具
- 每月读 2-3 篇前沿论文
- 每季度学一个新技能

**关键是持续，而不是爆发。**

---

## 📍 常见问题

### Q1: 这么多资源，从哪里开始？

**A**: 根据你的目标选择：
- 想快速上手 → 先看 **工具与教程** 板块
- 想系统学习 → 先看 **课程与视频** 板块
- 想了解前沿 → 先看 **AI 生物学工具** 板块

### Q2: 免费资源够用吗？要不要付费？

**A**: **90% 的情况下免费资源就够了。** 付费工具我只推荐：
- GitHub Copilot（学生免费）
- ChatGPT Plus（如果你高频使用）

其他的等你真正有需求了再考虑。

### Q3: 英文不好怎么办？

**A**:
1. **优先用国内资源**：生信技能树、生信菜鸟团的中文教程很全面
2. **借助翻译工具**：DeepL 的翻译质量很高
3. **边用边学**：实践中积累专业词汇，比死记硬背快

### Q4: 资源会持续更新吗？

**A**: 会的！我会保持每月更新，添加新的优质资源。

完整清单在这里 👉 [**资源导航页面**](/links/)

---

## 写在最后

这份清单耗费了我大量心血，但如果能帮到一个人少走弯路，就值了。

**生物信息学的学习之路很长，但不孤单。** 我们都是在这条路上摸索前进的。

几句肺腑之言：

1. **不要被资源的数量吓到** - 选择适合你的，深入学习
2. **不要闭门造车** - 多参与社区讨论，分享你的经验
3. **不要停止学习** - 这个领域变化太快，保持好奇心
4. **不要忘记初心** - 我们做这些，是为了推动科学进步

---

## 🙏 致谢

这份清单的诞生离不开：
- 生信社区的前辈们无私分享
- 师兄师姐们的经验传授
- 开源社区的贡献者们

**如果你觉得有帮助，欢迎分享给更多人！**

有任何建议或补充，欢迎在评论区留言，或者直接访问 [资源导航页面](/links/) 查看完整内容。

---

**相关阅读：**
- [计算生物学成长路线图](/links/) （附路线图）

---

*本文约 7000 字，阅读时间约 15 分钟*
*最后更新：2026-02-25*
*转载请注明出处*
