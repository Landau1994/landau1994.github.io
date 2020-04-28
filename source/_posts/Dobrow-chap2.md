---
title: Dobrow-chap2
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - math
tags:
  - stochastic Process
  - note
date: 2020-04-05 21:17:46
keywords:
description: 对于 `Introduction to stochastic processes with R`一书的笔记
photos:
---
对于 `Introduction to stochastic processes with R`一书的笔记

> Let us finish the article and the whole book with a good example of dependent trials, which approximately can be considered as a simple chain.
>                                            –Andrei Andreyevich Markov

Chap2: Markov Chains: First steps

本章讲马尔可夫链

### Introduction 

1. 引入的案例：
   这一节，用一个类似大富翁的游戏来引入马尔可夫夫性。
2. 马尔可夫链的形式化定义为
    > Markov Chain
    > Let $\mathcal{S}$ be a discrete set. A Markov chain is a sequence of random variables $X_0,X_1,\dots$ taking values in $\mathcal{S}$ with the property that
    > $$ P(X_{n+1}=j|X_0=x_0,\dots,X_{n-1}=x_{n-1},X_n=i) \\ = P(X_{n+1}=j|X_n=i) \tag{2.1} $$ 
    > for all $x_0,\dots,x_{n-1},i,j\in\mathcal{S},n\ge0$ The set $\mathcal{S}$ is the state space of the Markov chain. 、

3. $X_n=i$ 称为在时刻n到达状态i。
   
4. 时间齐性马尔可夫链：
   $$ P(X_{n+1}=j|X_{n}=i)=P(X_1=j|X_0=i) \tag{2.2}$$

5. transition matrix:
   	1. n步转移矩阵计算（矩阵乘法）
   	2. 若干例子：
   		1) 收敛于一个各行相等的矩阵；
   		2) 不收敛，进入跳跃的状态；
   		3) 收敛于一个各行不相等的矩阵

第五部分从直观上为下一章铺路。