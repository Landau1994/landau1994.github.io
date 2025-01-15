---
title: udl-chap3-4
author: 研精极锐
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: An apprentice in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - AI
  - note
date: 2025-01-15 19:29:22
keywords:
description:
photos:
---

# chapter 3: Shallow neural networks

如上一章末尾所言，第3章讨论了浅层神经网络，它比线性回归稍复杂，但能描述更广泛的输入/输出关系。这章主要以浅层神经网络为例，初步介绍神经网络的基本概念。浅层神经网络的架构，激活函数（主要介绍ReLU），以及universal approximation theorem

> 浅层神经网络有一个隐藏层。它们（i）计算输入的多个线性函数，（ii）将每个结果通过一个激活函数，然后（iii）取这些激活的线性组合以形成输出。浅层神经网络通过将输入空间划分为连续的分段线性区域来基于输入 \( x \) 进行预测 \( y \)。只要有足够的隐藏单元（神经元），浅层神经网络可以任意精确地近似任何连续函数。

> 第4章讨论了深度神经网络，它通过增加更多的隐藏层来扩展本章中的模型。第5-7章描述了如何训练这些模型。


# chapter 4: deep neural networks

这更进一步，初步介绍了深层神经网络。


$h_1 = a[\beta_0 + \Omega_1x]$
$h_2 = a[\beta_1 + \Omega_2h_1]$
$h_3 = a[\beta_2 + \Omega_2h_2]$
$\vdots$
$h_K = a[\beta_{K-1} + \Omega_{K}h_{K-1}]$
$y = \beta_K + \Omega_Kh_K$

```Mermaid
graph LR
    x1((x₁)) --> h11((h₁₁))
    x1 --> h12((h₁₂))
    x1 --> h13((h₁₃))
    x2((x₂)) --> h11
    x2 --> h12
    x2 --> h13
    h11 --> h21((h₂₁))
    h12 --> h21
    h13 --> h21
    h11 --> h22((h₂₂))
    h12 --> h22
    h13 --> h22
    h21 --> y((y))
    h22 --> y
    
    style x1 fill:#f9f,stroke:#333,stroke-width:2px
    style x2 fill:#f9f,stroke:#333,stroke-width:2px
    style h11 fill:#bbf,stroke:#333,stroke-width:2px
    style h12 fill:#bbf,stroke:#333,stroke-width:2px
    style h13 fill:#bbf,stroke:#333,stroke-width:2px
    style h21 fill:#fbb,stroke:#333,stroke-width:2px
    style h22 fill:#fbb,stroke:#333,stroke-width:2px
    style y fill:#bfb,stroke:#333,stroke-width:2px
```

最重要的是末尾的：

深层神经网络与浅层神经网络的比较：

+ Ability to approximate different functions 
+ Number of linear regions per parameter 
+ Depth efficiency（深层优势，作者在末尾的note部分举了一些理论的例子）
+ Large, structured inputs (如图像，实际上只能深层)
+ Training and generalization 

> 在本章中，我们首先考虑了当我们组合两个浅层网络时会发生什么。我们论证了第一个网络"折叠"输入空间，而第二个网络则对其应用分段线性函数。当输入空间被折叠到自身时，第二个网络的效果会被复制。
>
> 我们接着证明了这种浅层网络的组合是深度网络中两层的特例。我们将每一层中的 ReLU 函数解释为在多个位置裁剪输入函数，并在输出函数中创造"新的点"。我们引入了超参数的概念，就目前所见的网络而言，超参数包括每层中隐藏单元的数量。
>
> 最后，我们比较了浅层和深度网络。我们发现：(i) 两种网络都能以足够的容量近似任何函数；(ii) 深度网络每个参数产生更多的线性区域，使得某些函数能够用深度网络更有效地近似；(iv) 大型、结构化的输入（如图像）最好在多个阶段处理；以及 (v) 在实践中，最好的结果往往是使用具有多层的深度网络获得的。
>
> 现在我们已经理解了深度和浅层网络模型，我们将注意力转向训练它们。在下一章中，我们将讨论损失函数。对于任何给定的参数值 $\varphi$，损失函数返回一个单一数值，表示模型输出与训练数据集的真实预测之间的不匹配程度。在第 6 章和第 7 章中，我们将讨论训练过程本身，即我们如何寻找能最小化这个损失的参数值。