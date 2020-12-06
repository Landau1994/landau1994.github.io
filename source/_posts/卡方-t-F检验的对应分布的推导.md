---
title: '卡方,t,F检验的对应分布的推导'
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - math
tags:
  - Probability
  - note
date: 2018-08-01 15:16:53
keywords: 数理统计 t检验 读书笔记
description: 数学分析和线性代数真的很有用
photos:
---
## 前言

假设检验是我们在日常研究中，经常碰到的统计问题。对于追求实用与效率的科研人员来说，各种不同的假设检验是可以用软件，点点鼠标，或者写写代码，就可以完成的。

不过，对于我们这些想要在生物信息领域深入和进阶，并且最终有所建树的学生来说，我们光会拧螺丝和用板子，用轮子是不够的，当有新的技术，新的需求出来之后，我们得要造新轮子，开发新方法。因此，我们还是得学学火箭是咋飞起来和板子以及轮子是咋造出来的知识。

我们以学习和介绍研究中比较基础的$\chi^2,t,F$三种检验的所对应的分布推导，开始我们的进阶之旅。

说明：本文的推导来自《概率统计讲义》第三版附录二，陈家鼎等编著，高等教育出版社出版。略微有所修改，阅读本文，只需修过本科阶段非数学专业的三门基础数学课：高等数学(不是很深，也不是很浅的数学分析)，线性代数，概率论与数理统计。
## 正交矩阵与正态分布

在线性代数课程中，我们知道，若$n$阶方阵$A=(a_{ij})_{n\times n}$满足$A^TA=I$，写成标量的形式就是：

$$
\sum_{k=1}^{n} a_{ki}a_{kj} = \begin{cases}  
1 & i=j\\
0 & i\ne j
\end{cases} \quad (1)
$$      

此时，我们称方阵$A$为正交矩阵。而且，通过线性代数的课程，我们知道，正交矩阵满足如下性质：
 + 1-1 设A是正交矩阵，则$AA^T=I$，并且结合（1）可得：

   $$
   \sum_{k=1}^{n} a_{kj} a_{lj} = \begin{cases} 1 & k=l\\
   0 & k\ne l
   \end{cases} \quad (2)
   $$
 
 + 1-2 设A是正交矩阵，则$A^T$也是正交矩阵，并且$|A|=1$或$|A|=-1$，其中$|\cdot|$表示行列式。

 + 1-3 若$A=(a_{ij})_{n\times n}$是正交矩阵，而$x_1,x_2,\dots,x_n$是任意n个实数，对于

   $$ y_i = \sum_{i=1}^{n} a_{ik}x_k  \quad i=1,2,\dots,n $$ 
   我们有
   $$ \sum_{i=1}^{n}y^2_i = \sum_{i=1}^{n} x_i^2 \quad (3)$$ 

很抱歉，开头罗列了这么多线性代数的事实，不过，也没办法，要做菜，我们得先备料不是吗。下面我们开始做菜了。

**定理1** 设$X_1,X_2,\dots,X_n$相互独立，且都服从$N(0,\sigma^2)$，又$A=(a_{ij})$是正交矩阵，构造随机变量
$$ Y_i = \sum_{i=1}^{n}a_{ij}X_j \quad (1\le i\le n) $$
**证明** 因$X_i$的分布密度是$\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{1}{2\sigma^2}x^2}$,且$X_i$是独立同分布样本（i.i.d.），故$X_1,X_2,\dots,X_n$联合密度为：$(\frac{1}{\sqrt{2\pi}\sigma})^ne^{-\frac{1}{2\sigma^2}\sum_{i=1}^{n}x^2}$

构造n维空间中的区域D:
$$ D=\{(x_1,x_2,\dots,x_n)|a_i < \sum_{j=1}^{n}a_{ij}x_j < b,i=1,2\dots,n\} $$
则有：

$$
\begin{aligned}
  & P\{a_1 < Y_1 < b_1,a_2 < Y_2 < b_2, \dots, a_n < Y_n < b_n \} \\
  & =P\{ (X_1,X_2,\dots,X_n)\in D \}\\
  & = \int_{D}(\frac{1}{\sqrt{2\pi}\sigma})^ne^{-\frac{1}{2\sigma^2}\sum_{i=1}^{n}x^2}dx_1dx_2\dots dx_n
\end{aligned}
$$

注意到
$$ y_i = \sum_{k=1}^{n}a_{ik}x_k \quad (i=1,2,\dots,n)$$
于是（利用正交矩阵的性质）
$$ x_i = \sum_{k=1}^{n}a_{ki}y_k \quad (i=1,2,\dots,n)$$
容易验证，变换的雅可比式为
$$ J(\frac{x_1,x_2,\dots,x_n}{y_1,y_2,\dots,y_n})= |A^T|=1,-1 $$

又$\sum_{k=1}^{n}x^2_{k}=\sum_{k=1}^n y_{k}^2$故

$$
\begin{aligned}
  &\int_D(\frac{1}{\sqrt{2\pi}\sigma})^ne^{-\frac{1}{2\sigma^2}\sum_{i=1}^{n}x^2}dx_1dx_2\dots dx_n \\
  &= \int_{a_1}^{b_1}\int_{a_2}^{b_2}\dots\int_{a_n}^{b_n}(\frac{1}{\sqrt{2\pi}\sigma})^ne^{-\frac{1}{2\sigma^2}\sum_{i=1}^{n}y^2} \cdot |J(\frac{x_1,x_2,\dots,x_n}{y_1,y_2,\dots,y_n})| dy_1dy_2\dots dy_n  \\
  &=\int_{a_1}^{b_1}\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{-y_1^2}{2\sigma^2}}dy_1\cdot\int_{a_2}^{b_2}\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{-y_2^2}{2\sigma^2}}dy_2 \cdots \int_{a_n}^{b_n}\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{-y_n^2}{2\sigma^2}}dy_n \\
  &=  P\{a_1 < Y_1 < b_1\} \cdot P\{a_2 < Y_2 < b_2\} \cdots P\{a_n < Y_n < b_n\}
\end{aligned}
$$
故$Y_1,Y_2,\dots,Y_n$相互独立，且不难看出，都服从$N(0,\sigma^2)$。定理1证毕。

**定理2**设$X_1,X_2,\dots,X_n$相互独立，且$X \sim N(\mu,\sigma^2)$。$A=(a_{ij})$是n阶正交矩阵,构造随机变量，
$$ Y_{i}=\sum_{k=1}^{n}a_{ik}X_k \quad (i=1,2,\dots,n) $$
则$Y_1,Y_2,\dots,Y_n$相互独立，且
$$ Y_{i} \sim N(\sum_{k=1}^{n}a_{ik}\mu_k,\sigma^2) $$
**证明**令$Z_i = X_i - \mu$，则$Z_1,Z_2,\dots,Z_n$相互独立，都服从$N(0,\sigma^2)$,根据定理1知，$\sum_{k=1}^{n}a_{1k}Z_k,\sum_{k=1}^{n}a_{2k}Z_k,\dots,\sum_{k=1}^{n}a_{nk}Z_k$相互独立。
且
$$ \sum_{k=1}^{n}a_{ik}Z_k \sim N(0,\sigma^2)$$
但是
$$ Y_i = \sum_{k=1}^{n}a_{ik}Z_k + \sum_{k=1}^{n}a_{ik}\mu_k $$
故$Y_1,Y_2,\dots,Y_n$相互独立，且$Y_i \sim N(\sum_{k=1}^{n}a_{ik}\mu_k,\sigma^2) \quad  (i=1,2,\dots,n)$

## 关于$\chi^2$分布

前面的的都是小菜，接下来上主菜。我们要开始证明一系列很fancy的定理

**定理3** 设$X_1,X_2,\dots,X_n$相互独立，并且都服从$N(0,1)$,则$\xi=\sum_{i=1}^n X^{2}_{i}$服从$n$个自由度的$\chi^2$分布，其PDF(probability density function)为
$$ p(u) = k_n(u)= \begin{cases}  
\frac{1}{2^{\frac{n}{2}}\Gamma(\frac{n}{2})}u^{\frac{n}{2}-1}e^{-\frac{u}{2}} & u>0\\
0 & u\le 0 
\end{cases} \quad (4)$$
**证明** 我们证明的策略是，先求出CDF(cumulative distribution function)$F(u)=P\{\xi \le u\}$，然后利用中值定理，证明$F^\prime(u)=k_n(u)$。

显然，当$u\le 0,F(u)=0,F^\prime(u)=k_n(u)$

当$u>0$时，由于$X_1,X_2,\dots,X_n$相互独立，故$X_1,X_2,\dots,X_n$联合密度为$(\frac{1}{\sqrt{2\pi}})^ne^{-\frac{1}{2}\sum_{i=1}^{n}x^2}$，
故
$$ 
\begin{aligned}
 F(u) &= P\{ \sum_{i=1}^{n}X_{i}^2 \le u \} \\
      &=  \int_{\sum_{i=1}^{n}x_{i}^2 \le u}(\frac{1}{\sqrt{2\pi}})^ne^{-\frac{1}{2}\sum_{i=1}^{n}x^2}dx_1dx_2\dots dx_n
\end{aligned}  
$$
故对于$h>0$,有
$$ 
\begin{aligned}
 F(u+h) - F(u)  =& \int_{u < \sum_{i=1}^{n}x_{i}^2 \le u+h}(\frac{1}{\sqrt{2\pi}})^ne^{-\frac{1}{2}\sum_{i=1}^{n}x^2}dx_1dx_2\dots dx_n \\
 \le& \int_{u < \sum_{i=1}^{n}x_{i}^2 \le u+h}(\frac{1}{\sqrt{2\pi}})^ne^{-\frac{u}{2}} dx_1dx_2\dots dx_n \\
 =& (\frac{1}{\sqrt{2\pi}})^n e^{-\frac{u}{2}} \int_{u < \sum_{i=1}^{n}x_{i}^2 \le u+h} dx_1dx_2\dots dx_n 
\end{aligned}  
$$
$$ 
\begin{aligned}
 F(u+h) - F(u)
\ge& \int_{u < \sum_{i=1}^{n}x_{i}^2 \le u+h}(\frac{1}{\sqrt{2\pi}})^ne^{-\frac{u+h}{2}} dx_1dx_2\dots dx_n \\
=& (\frac{1}{\sqrt{2\pi}})^n e^{-\frac{u+h}{2}} \int_{u < \sum_{i=1}^{n}x_{i}^2 \le u+h} dx_1dx_2\dots dx_n
\end{aligned}  
$$
令$S(x) = \int_{\sum_{i=1}^{n}x_{i}^2 \le x} dx_1dx_2\dots dx_n \quad (x>0)$
则
$$ 
(\frac{1}{\sqrt{2\pi}})^n e^{-\frac{u+h}{2}}\cdot\frac{S(u+h)-S(u)}{h} \\
\le \frac{F(u+h)-F(u)}{h} \le \\
(\frac{1}{\sqrt{2\pi}})^n e^{-\frac{u}{2}}\cdot\frac{S(u+h)-S(u)}{h}
 $$
问题现在变为如何求$S(x)$

做代换$x_i= y_i \sqrt{x}$,则
$$ dx_i = \sqrt{x} dy_i $$
由此 
$$ S(x) = \int_{\sum_{i=1}^{n}y_{i}^2 \le 1} (\sqrt{x})^n dy_1dy_2\dots dy_n = x^{\frac{n}{2}}\cdot C_n $$
有趣的是，我们可以看出
$$ C_n = \int_{\sum_{i=1}^{n}y_{i}^2 \le 1}  dy_1dy_2\dots dy_n $$ 
是$n$维单位球体的体积。不过在我们的问题中，我们可以看出它只和$n$有关的量。故$S(x)=\frac{n}{2}C_nx^{\frac{n}{2}-1}$

根据之前的不等式，结合中值定理：
$$ \lim_{h\rightarrow0^{+}} \frac{F(u+h)-F(u)}{h} = (\frac{1}{\sqrt{2\pi}})^nC_n\frac{n}{2}u^{\frac{n}{2}-1}e^{-\frac{n}{2}} $$
$$ \lim_{h\rightarrow0^{-}} \frac{F(u+h)-F(u)}{h} = (\frac{1}{\sqrt{2\pi}})^nC_n\frac{n}{2}u^{\frac{n}{2}-1}e^{-\frac{n}{2}} $$
所以
$$ F^{\prime}(u) = B_n u^{\frac{n}{2}-1}e^{-\frac{u}{2}} $$
综上
$$ p(u) = k_n(u)= \begin{cases}  
B_nu^{\frac{n}{2}-1}e^{-\frac{u}{2}} & u>0\\
0 & u\le 0 
\end{cases} \quad $$
由归一化条件$\int_{-\infty}^{+\infty}p(u)du=1$知$\int_{0}^{+\infty}B_nu^{\frac{n}{2}-1}e^{-\frac{u}{2}}du=1$
而在数学分析的知识告诉我们
$\int_{0}^{+\infty}u^{\frac{n}{2}-1}e^{-\frac{u}{2}}du=2^{\frac{n}{2}}\Gamma(\frac{n}{2})$

$$ B_n = \frac{1}{2^{\frac{n}{2}}\Gamma(\frac{n}{2})} $$
定理得证。

这个定理的一个副产物是，告诉了我们$n$维单位球体的体积$C_n=\frac{\pi^\frac{n}{2}}{\Gamma(\frac{n}{2}+1)}$

**推论** 若$\xi \sim \chi^2(n)$，则有$E(\xi)=n$

**证明** 由定理1，结合数学期望的性质，知$E(\xi)=E(\sum_{i=1}^{n}X_i^2)=\sum_{i=1}^{n}E(X_i^2)=\sum_{i=1}^{n}(D(X_i)+(E(X_i))^2)=\sum_{i=1}^{n}1=n$ 

$\quad \Box$

**定理4** 若$\xi$与$\eta$相互独立，且$\xi \sim\chi^2(n_1),\eta \sim\chi^2(n_2)$，则$\xi+\eta \sim \chi^2(n_1+n_2)$

**证明** 设$\xi,\eta，\xi+\eta$的分布函数分别为$p_1(x),p_2(x),p(x)$，我们先分别不加证明的引用概率论和Gamma函数的两个结论：

1).已知(X,Y)的联合密度是$p(x,y)$，$Z=Y+Y$的PDF为：
$$ 
p_z(z)=\int_{-\infty}^{\infty}p(x,z-x)dx
 $$
2).$\int_{0}^{1}v^{p-1}(1-v)^{q-1}dv=\frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}$(p,q为正整数)

下面开始证明：
当 $x \le 0$时，$P(\xi+\eta \le 0)=0,p(x)=0$,定理成立。

当 $x > 0$时，

$$ \begin{aligned}
p(x) &= \int_{0}^{x}\frac{1}{2^{\frac{n_1}{2}}\Gamma(\frac{n_1}{2})}u^{\frac{n_1}{2}-1}e^{-\frac{u}{2}}\frac{1}{2^{\frac{n_2}{2}}\Gamma(\frac{n_2}{2})}(x-u)^{\frac{n_2}{2}-1}e^{-\frac{x-u}{2}}du \\
&=\frac{e^{-\frac{x}{2}}x^{\frac{n_1+n_2}{2}-1}}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}\int_{0}^{1}v^{\frac{n_1}{2}-1}(1-v)^{\frac{n_2}{2}-1}dv （v=\frac{u}{x}) \\
&=\frac{e^{-\frac{x}{2}}x^{\frac{n_1+n_2}{2}-1}}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}\cdot\frac{\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}{\Gamma(\frac{n_1+n_2}{2})} \\
&=\frac{1}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}x^{\frac{n_1+n_2}{2}}e^{-\frac{x}{2}}
\end{aligned} 
$$

综上：

$$ 
p(x) = \begin{cases}
\frac{1}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}x^{\frac{n_1+n_2}{2}}e^{-\frac{x}{2}} & x>0\\
0 & x\le 0
\end{cases}
 $$ 
 $\Box$

 **定理5** 若$x_1,x_2,\dots,x_n$相互独立，且都服从分布$N(0,1)$,则有如下三条结论：
 
 1. $\bar{X}=\frac{X_1+X_2+\dots
 +X_n}{n}\sim N(0,\frac{1}{n})$ 
 2. $\sum_{i=1}^{n}(X_i-\bar{X}) \sim \chi^2(n-1)$
 3. $\bar{X}$与$\sum_{i=1}^{n}(X_i-\bar{X})$相互独立

**证明** 构造正交矩阵

$$ 
  \begin{bmatrix}
   \frac{1}{\sqrt{n}} & \frac{1}{\sqrt{n}} & \frac{1}{\sqrt{n}} & \dots & \frac{1}{\sqrt{n}}\\
   \frac{1}{\sqrt{1\cdot 2}} & \frac{-1}{\sqrt{1\cdot 2}} & 0 &\dots & 0  \\
   \frac{1}{\sqrt{2\cdot 3}} & \frac{1}{\sqrt{2\cdot 3}} & \frac{-2}{\sqrt{2\cdot 3}} &\dots & 0  \\
   \vdots & \vdots & \vdots & \ddots & 0 \\
   \frac{1}{\sqrt{(n-1)n}} & \frac{1}{\sqrt{(n-1)n}} & \frac{1}{\sqrt{(n-1)n}} &\dots & \frac{-(n-1)}{\sqrt{(n-1)n}} 
  \end{bmatrix} 
 $$
 由此正交矩阵，我们可以构造随机变量：

 $$ 
 \begin{aligned}
   Y_1 &= \frac{1}{\sqrt{n}}(X_1+X_2+\dots+X_n) \\
   Y_2 &= \frac{1}{\sqrt{1\cdot 2}}(X_1-X_2) \\
   Y_3 &= \frac{1}{\sqrt{2\cdot 3}}(X_1+X_2-2X_3) \\
   \dots \\
    Y_n &= \frac{1}{\sqrt{(n-1)n}}(X_1+X_2+\dots+X_{n}-(n-1)X_n)
 \end{aligned}
  $$
  有定理1可知，$Y_1,Y_2,\dots,Y_n$相互独立，且都服从$N(0,1)$，
  我们发现$Y_1 \sim N(0,1)$，因此$\bar{X}=\frac{1}{\sqrt{n}}Y_1\sim N(0,\frac{1}{n})$，第一条结论得证。

  由于$\sum_{i=1}^{n}X_i^{2}=\sum_{i=1}^{n}Y_i^{2}$故

  $$ 
 \begin{aligned}
   \sum_{i=1}^{n}(X_i-\bar{X}) &= \sum_{i=1}^{n}X_i^2-n\bar{X} \\
    &= \sum_{i=1}^{n}Y_i^2-n(\frac{1}{\sqrt{n}}Y_i)^2 \\
    &= \sum_{i=1}^{n}Y_i^2-Y^2_i \\
    &= \sum_{i=2}^{n}Y_i^2 \sim \chi^2(n-1)
 \end{aligned}
  $$
  第二条结论得证。

  由于$Y_1,Y_2,\dots,Y_n$相互独立，且
  $$ 
  \bar{X}=\frac{1}{\sqrt{n}}Y_1, \sum_{i=1}^{n}(X_i-\bar{X})^2=\sum_{i=2}^{n}Y^{2}_{i}
   $$
  故$\bar{X}$与$\sum_{i=1}^{n}(X_i-\bar{X})^2$独立，第三条结论得证
  $\Box$

  **推论** 若$x_1,x_2,\dots,x_n$相互独立，且都服从分布$N(\mu,\sigma^2)$,则有如下三条结论：
 
 1. $\bar{X}=\frac{X_1+X_2+\dots
 +X_n}{n}\sim N(\mu,\frac{\sigma^2}{n})$ 
 2. $\frac{1}{\sigma^2}\sum_{i=1}^{n}(X_i-\bar{X}) \sim \chi^2(n-1)$
 3. $\bar{X}$与$\sum_{i=1}^{n}(X_i-\bar{X})$相互独立

## 关于t分布
**定理6** 设$\xi,\eta$相互独立，且$\xi \sim N(0,1),\eta \sim \chi^2(n)$， 则$\zeta=\frac{\xi}{\sqrt{\frac{\eta}{n}}} \sim t(n)$，其PDF为：
$$ 
p(u)=t_n(u) = \frac{\Gamma(\frac{n+1}{2})}{\Gamma(\frac{n}{2})\sqrt{n\pi}}(1+\frac{u^2}{n})^{-\frac{n+1}{2}} \quad (5)
 $$
**证明** 与定理3证明的思路类似，设$F(u) = P\{\zeta \le u\}$证明$F^{\prime}(u)=t_n(u)$, 由已知：
$$ 
\begin{aligned}
  F(u) 
  &= P\{\frac{\xi}{\sqrt{\frac{\eta}{n}}} \le u\} \\
  &= P\{\frac{\xi}{\sqrt{\eta}} \le \frac{u}{\sqrt{n}}\} \\
  &= \iint_{\frac{x}{\sqrt{y}}\le\frac{u}{\sqrt{n}}} \frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}\cdot\frac{1}{2^{\frac{n}{2}}\Gamma(\frac{n}{2})}y^{\frac{n}{2}-1}e^{-\frac{y}{2}}dxdy \\
  &= \iint_{t \le \frac{u}{\sqrt{n}},y > 0} \frac{1}{\sqrt{2\pi}2^{\frac{n}{2}}\Gamma(\frac{n}{2})}s^{\frac{n}{2}-1}\cdot e^{-\frac{1}{2}(1+t^2)s}|J(\frac{x,y}{s,t})|dsdt \quad (t=\frac{x}{\sqrt{y}},y=s) \\
  &= \int_{-\infty}^{\frac{u}{\sqrt{n}}} \frac{dt}{\sqrt{\pi}\Gamma(\frac{n}{2})}\int_{0}^{+\infty}(\frac{s}{2})^{\frac{n-1}{2}}e^{-\frac{s}{2}(1+t^2)}d(\frac{s}{2}) \\
  &= \int_{-\infty}^{\frac{u}{\sqrt{n}}} \frac{1}{\sqrt{\pi}\Gamma(\frac{n}{2})} \cdot \frac{\Gamma(\frac{n+1}{2})}{(1+t^2)^{\frac{n+1}{2}}}dt \\
  &= \int_{-\infty}^{u} \frac{1}{\sqrt{n\pi}\Gamma(\frac{n}{2})} \cdot \frac{\Gamma(\frac{n+1}{2})}{(1+\frac{v^2}{n})^{\frac{n+1}{2}}}dv \quad (v=t\sqrt{n})
\end{aligned}
 $$
 故$F^{\prime}(u)= t_n(u) \quad \Box$ 
 
 定理5,6可以用来证明下面这个在统计学里很有作用的定理：

 **定理7** 设$X_1,X_2,\dots,X_n (n \ge 2)$相互独立，且都服从$N(\mu,\sigma^2)$,则$T=\frac{\bar{X}-\mu}{\sqrt{\frac{S^2}{n}}} \sim t(n-1)$其中
 $$ \bar{X}=\frac{X_1+X_2+\dots+X_n}{n}, S^2=\frac{1}{n-1}\sum_{i=1}^{n}(X_i-\bar{X})^2$$

 **证明** 构造随机变量
 $$ 
 \xi = \frac{\bar{X}-\mu}{\sqrt{\frac{\sigma^2}{n}}},\eta=\frac{1}{\sigma^2}\cdot\sum_{i=1}^{n}(X_i-\bar{X})^2
  $$
 根据定理5的推论，我们知道$\xi,\eta$相互独立，且$\xi\sim N(0,1),\eta\sim \chi^2(n-1)$
 故根据定理6，$\frac{\xi}{\sqrt{\frac{\eta}{n-1}}}\sim t(n-1)$
 故 
 $$ T=\frac{\bar{X}-\mu}{\sqrt{\frac{S^2}{n}}}=\frac{\frac{\bar{X}-\mu}{\sqrt{\frac{\sigma^2}{n}}}}{\sqrt{\frac{S^2}{\sigma^2}}}=\frac{\xi}{\sqrt{\frac{\eta}{n-1}}} \sim t(n-1) \quad \Box$$

 ## 关于F分布

 **定理8** 设$\xi,\eta$相互独立，且$\xi\sim \chi^2(n_1),\eta\sim \chi^2(n_2)$ 则 $\zeta=\frac{\frac{\xi}{n_1}}{\frac{\eta}{n_2}} \sim F(n_1,n_2)$ 其PDF为：
 $$ 
 p(u) = f_{n_1,n_2}(u) = \begin{cases}
   \frac{\Gamma(\frac{n_1+n_2}{2})}{\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}(\frac{n_1}{n_2})^{\frac{n_1}{2}}u^{\frac{n_1}{2}-1}(1+\frac{n_1}{n_2}u)^{-\frac{n_1+n_2}{2}} & u > 0 \\
   0 & u \le 0
 \end{cases}
  $$
**证明** 跟之前一样，令$F(u)=P\{\xi\le u\}$ ，证明$F^{\prime}(u)=f_{n_1,n_2}(u)$ 

当 $u \le 0, F(u)=0$,
$$ 
\begin{aligned}
  F(u) 
  &= P\{\xi\le u\} \\
  &= P\{\frac{\frac{\xi}{n_1}}{\frac{\eta}{n_2}}\le u\} \\
  &= P\{\frac{\xi}{\eta}\le \frac{n_1}{n_2}u\} \\
  &= \iint_{\frac{x}{y}\le\frac{n_1}{n_2}u,x>0,y>0}\frac{1}{2^{\frac{n_1}{2}}\Gamma(\frac{n_1}{2})}x^{\frac{n_1}{2}-1}e^{-\frac{x}{2}}\frac{1}{2^{\frac{n_2}{2}}\Gamma(\frac{n_2}{2})}y^{\frac{n_2}{2}-1}e^{-\frac{y}{2}}dxdy \\
  &= \iint_{0 < t \le \frac{n_1}{n_2}u,s>0} \frac{e^{-\frac{s}{2}(1+t)}s^{\frac{n_1+n_2}{2}-2}}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}t^{\frac{n_1}{2}-1}\cdot|J(\frac{x,y}{s,t})|dsdt \quad x=st,y=s \\
  &= \int_{0}^{\frac{n_1}{n_2}u}\frac{t^{\frac{n_1}{2}-1}}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}dt\int_{0}^{+\infty}s^{\frac{n_1+n_2}{2}-1}t^{\frac{n_1}{2}-1}dsdt \\
  &= \int_{0}^{\frac{n_1}{n_2}u}\frac{\Gamma(\frac{n_1+n_2}{2})}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}t^{\frac{n_1}{2}-1}(1+t)^{-\frac{n_1+n_2}{2}}dt \\
  &= \int_{0}^{u}\frac{\Gamma(\frac{n_1+n_2}{2})}{2^{\frac{n_1+n_2}{2}}\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})}(\frac{n_1}{n_2})^{\frac{n_1}{2}}v^{\frac{n_1}{2}-1}(1+\frac{n_1}{n_2}v)^{-\frac{n_1+n_2}{2}}dv \quad  (t=\frac{n_1}{n_2}v)
\end{aligned}
$$
故$F^{\prime}(u)=f_{n_1,n_2}(u)$  $\quad \Box$

**定理9** 设$X_1,X_2,\dots,X_{n_1},Y_1,Y_2,\dots,Y_n$, 这$n_1+n_2$个随机变量相互独立，且都服从$N(\mu,\sigma^2)$,则
$$ \zeta=\frac{\frac{1}{n_1-1}\sum_{i=1}^{n_1}(X_i-\bar{X})^2}{\frac{1}{n_2-1}\sum_{i=1}^{n_2}(Y_i-\bar{Y})^2} \sim F(n_1-1,n_2-1)$$ 

**证明** 构造随机变量

$$
\xi=\frac{1}{\sigma^2}\cdot\sum_{i=1}^{n}(X_i-\bar{Y})^2, \eta=\frac{1}{\sigma^2}\cdot\sum_{i=1}^{n}(Y_i-\bar{Y})^2 
$$
由之前的结论，我们知道$\xi\sim\chi^2(n_1-1),\eta\sim\chi^2(n_2-1)$, 接下来证明$\xi,\eta$的独立性，构造随机变量：
$$ U_i = \frac{X_i-\mu}{\sigma},V_i=\frac{Y_i-\mu}{\sigma} 
$$
则
$$ 
\xi=\sum_{i=1}^{n_1}(U_i-\bar{U}),\eta=\sum_{i=1}^{n_2}(V_i-\bar{V})
 $$
由已知$U_1,U_2,\dots,U_{n_1},V_1,V_2,\dots,V_{n_2}$ 相互独立，且都服从$N(0,1)$,于是其联合分布密度为
$$ 
p(u_1,u_2,\dots,u_{n_1},v_1,v_2,\dots,v_{n_2}) \\
=(\frac{1}{\sqrt{2\pi}})^{n_1+n_2}\cdot e^{-\frac{1}{2}(\sum_{i=1}^{n_1}u^2_{i}+\sum_{i=1}^{n_2}v^2_{i})}
$$
所以,对于任意的实数$a,b,c,d$
$$ 
\begin{aligned}
  P\{a< \xi < b, c < \eta < d \}
  =& \int\limits_{a<\sum_{i=1}^{n_1}(u_i-\bar{u})^2 < b,c < \sum_{i=1}^{n_1}(v_i-\bar{v})^2 < d  }p(u_1,u_2,\dots,u_{n_1}v_1,v_2,\dots,v_{n_2})du_1 \dots du_{n_1}dv_1\dots dv_{n_2} \\
  =& \int_{a<\sum_{i=1}^{n_1}(u_i-\bar{u})^2 < b}(\frac{1}{\sqrt{2\pi}})^{n_1}e^{-\frac{1}{2}\sum_{i=1}^{n_1}u_i^2}du_1\dots du_{n_1} \cdot \\
  &\int_{a<\sum_{i=1}^{n_1}(v_i-\bar{v})^2 < b}(\frac{1}{\sqrt{2\pi}})^{n_1}e^{-\frac{1}{2}\sum_{i=1}^{n_1}v_i^2}dv_1\dots dv_{n_1} \\
  =& P\{a < \xi < b\} \cdot P\{c < \eta < d\} 
\end{aligned}
 $$
 独立性得证。

 再结合定理8，
 
 $$\zeta=\frac{\frac{1}{n_1-1}\sum_{i=1}^{n_1}(X_i-\bar{X})^2}{\frac{1}{n_2-1}\sum_{i=1}^{n_2}(Y_i-\bar{Y})^2} = \frac{\frac{\xi}{n_1-1}}{\frac{\eta}{n_2-1}} \sim F(n_1-1,n_2-1) \quad \Box
 $$




