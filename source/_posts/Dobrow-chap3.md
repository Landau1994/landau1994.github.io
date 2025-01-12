---
title: Dobrow-chap3
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - math
tags:
  - note
  - stochastic Process
description: Markov Chains for the Long term
date: 2020-05-06 22:15:53
keywords:
photos:
---


This is a note of the textbook `Introduction to stochastic processes with R`

> There exists everywhere a medium in things, determined by equilibrium.
>                                                  —Dmitri Mendeleev

承接上章最后的数值案例，本章主要讲转移步数趋于无穷时马尔可夫链的性质。

### 3.1 Limiting Distribution

#### 3.1.1 定义

作者给了有一个定义和三个等价定义。这个含义最清楚：

> A limiting distribution
for the Markov chain is a probability distribution 𝝀 with the property that, for any initial distribution $\boldsymbol{\alpha}$:
> $$\lim_{n\rightarrow\infty}\boldsymbol{\alpha}P^n=\boldsymbol{\lambda}$$

这个定义与原定义的等价性也很容易理解：

若对某一马尔可夫链的状态转移矩阵有：
$$\lim_{n\rightarrow\infty}P^n_{ij}=\lambda_j$$
并且设初始分布$\boldsymbol{\alpha}=(\alpha_1,\dots,\alpha_n)$,状态总数为$m$
则有：

$$
\begin{aligned}
  \lim_{n\rightarrow\infty}\boldsymbol{\alpha}\boldsymbol{P}^n &=(\alpha_1,\dots,\alpha_n)\begin{pmatrix}
   \lambda_1 & \lambda_2 & \cdots &  \lambda_m \\
   \lambda_1 & \lambda_2 & \cdots &  \lambda_m \\
   \vdots & \vdots & \ddots &  \lambda_m \\
   \lambda_1 & \lambda_2 & \cdots &  \lambda_m
  \end{pmatrix} \\ 
  & = \begin{pmatrix}
    \lambda_1\sum_{i=1}^n\alpha_i \\
    \lambda_2\sum_{i=1}^n\alpha_i \\
    \cdots \\
    \lambda_m\sum_{i=1}^n\alpha_i
  \end{pmatrix} \\
  & = \begin{pmatrix}
    \lambda_1 \\
    \lambda_2 \\
    \cdots \\
    \lambda_m
  \end{pmatrix} \\
  & = \boldsymbol{\lambda}
\end{aligned}
$$


Ex3.1 Two state Markov Chain
计算案例；

#### 3.1.2 Proportion of Time in Each State

利用计算条件期望的技术，来从状态在过程中所占的比例来理解Limit distribution



Ex3.2

```r
###### Simulate discrete-time Markov chain ########################
# Simulates n steps of a Markov chain 
# markov(init,mat,n,states)
# Generates X0, ..., Xn for a Markov chain with initiial
#  distribution init and transition matrix mat
# Labels can be a character vector of states; default is 1, .... k

markov <- function(init,mat,n,labels) { 
	if (missing(labels)) labels <- 1:length(init)
simlist <- numeric(n+1)
states <- 1:length(init)
simlist[1] <- sample(states,1,prob=init)
for (i in 2:(n+1)) 
	{ simlist[i] <- sample(states,1,prob=mat[simlist[i-1],]) }
labels[simlist]
}
####################################################
P <- matrix(c(0.1,0.2,0.4,0.3,0.4,0,0.4,0.2,0.3,0.3,0,0.4,0.2,0.1,0.4,0.3),
  nrow=4, byrow=TRUE)
lab <- c("Aerobics","Massage","Weights","Yoga")
rownames(P) <- lab
colnames(P) <- lab
P
init <- c(1/4,1/4,1/4,1/4) # initial distribution
states <- c("a","m","w","y")
# simulate chain for 100 steps
simlist <- markov(init,P,100,states)
simlist
table(simlist)/100
steps <- 1000000
simlist <- markov(init,P,steps,states)
table(simlist)/steps
```

### 3.2 Stationary Distribution

#### 3.2.1 定义
注意Stationary Distribution的定义没有出现极限。

> If the initial distributino is a stationary distribution, Then $X_0,X_1,\cdots,X_n$ is a sequence of identically distributed random variables. But it doesn't mean that the random variables are independent.

> Lemma 3.1: Limiting Distributions are stationary Distribution

The reverse is false, 反例：

$$
  P = \begin{pmatrix}
    0 & 1 \\
    1 & 0
  \end{pmatrix}
  ，\pi=(\frac{1}{2},\frac{1}{2}) 
$$

以及

$$
   P = \begin{pmatrix}
    1 & 0 \\
    0 & 1
  \end{pmatrix}
$$

#### 3.2.2 Regular Matrices

一个自然的想法是问，什么样的条件下，一个马尔可夫链有极限分布，而且极限分布就是平稳分布呢。

满足如下性质的马尔可夫链是符合这个要求的

> Regular Transition Matrix
> A transition matrix $\boldsymbol{P}$ is said to be regular if some power of $\boldsymbol{P}$ is positive. That is $\boldsymbol{P}^n > 0 $, for some $n\ge 1$

有定理：

> Theorem 3.2: A markov chain whose transition matrix $\boldsymbol{P}$ is regular has a limiting distribution, which is teh unique, positive, stationary distribution of the chanin.

Ex 3.3-3.4 具体算例，

#### 3.2.3 Finding the stationary distribution

本质上这是一个特征值问题。

```r
### Stationary distribution of discrete-time Markov chain
###  (uses eigenvectors)
###
stationary <- function(mat) {
x = eigen(t(mat))$vectors[,1]
as.double(x/sum(x))
}
```

Ex 3.5-3.6; 计算技巧，令$x_1=1$

Ex 3.7 The Ehrenfest dog-flea model

Ex 3.8 Random walk on a graph;

On weighted graph

> Stationiary Distribution for Random walk on a weighted graph
> 
> Let $G=(V,E)$ be  a weighted graph with edge weight function $w(i,j)$. For random walk on G, the stationary distribution $pi$ is proportion to the sum of teh edge weights incident to each vertex. That is.
> $$ \pi = \frac{w(v)}{\sum_{z}w(z)},\forall v\in V$$ 
> where
> $$ w(v) = \sum_{z \sim v }w(v,z) $$

On simple graph

> Stationary Distribution for simple Random Walk on a graph
> 
> For simple random walk on a weighted graph, set $w(i,j)=1,\forall i,j \in V $, then, $w(v)=deg(v)$,which gives
> $$ \pi_{v}=\frac{\deg(v)}{\sum_z \deg(z)}=\frac{\deg(v)}{2e} $$

Ex 3.9-3.10 如何计算的案例；

#### 3.2.4 The Eigenvalue Connection

转置，看出与特征值的关联。

Ex 3.10 理论案例： random walk in regular graph.。

### 3.3 Can you find the way to state $a$

#### 3.3.1 状态可到达与状态互通

> Say that state $j$ is accessible from state i, if $P_{ij}^n > 0$. That is,there is positive probability of reaching $j$ from $i$ in a finite number of steps. State $i$ and $j$ communicate if $i$ is accessible from $j$ and $j$ is accessible from $i$

eg 3.11 本例讲述了用 Transition graphs 展示 Communication classes

#### 3.3.2 不可约

> Irreducibility
> A Markov chain is called irreducible if it has exactly one cmmunication class. That is, all states communicate with each other

Ex 3.12 一个不可约链的例子；

#### 3.3.3 Recurrence and Transience

> Given a Markov chain $X_0,X_1,\dots$, let $T_j=\min\{n>0:X_n=j\}$be the first passage time to state $j$. If $X_n\ne j,\forall n>0$, see$T_j=\infty$. Let
> $$ f_j = P(T_j < \infty | X_0=j)$$
> be the probability started in $j$ eventually returns to $j$.
> 
> State $j$ is said to be recurrent if the Markov chain started in $j$ eventually revists $j$. That is $f_j=1$
> 
> State $j$ is said to be transient if there is positive probability that the Markov chain started in j never returns to $j$. That is $f_j < 1$

如何根据状态转移矩阵判定，某一个状态是Recurrent或Transient States。用示性函数，

$E(\sum_{n=0}^{\infty}I_n)=\sum_{n=0}^{\infty}E(I_n)=\sum_{n=0}^{\infty}P(X_n=j|X_0=i)=\sum_{n=0}^{\infty}P_{ij}^{n}$

由此可以推出另外一个判定条件；

> Recurrence, Transience
> 
> (i) State $j$ is recurrent if and only if
> $$ \sum_{n=0}^{\infty}P_{ij}^n=\infty $$
> 
> (ii) State j is transient if and only if
> 
> $$\sum_{n=0}^{\infty}P_{ij}^n<\infty$$

> Recuurence and Transience are Class Properties
> 
> Theorem 3.3 The states of a communication class are either all recurrent or all transient.
> Corollary 3.4 For a finite irreducible Markov chain, all states are recuurent.

Ex 3.13 接下来的例子是简单的一维随机游走；这个例子可以推广到高维。

#### 3.3.4 Canonical Decomposition

Closed Communication Class

> Lemma 3.5 A communication class is closed if it consists of all recurrent states. A finite communication class is closed only if it consits of all recurrent states.

反证法即可证得；最后便可以得到，我们想定义的；

> The state space S of a finite Markov chain can be partitioned into transient and reccurent states as $S=T \cup R_1 \cup \cdots R_m$, where T is the set of all transient states and $R_i$ are closed communiction classes of recurrent states. This is called the canonical decomposition.

注：由等价类的定义可以保障这么重排状态转移矩阵，是与原矩阵等价的。

> Given a canonical decomposition, the state space can be reordered so that the Markov transition matrix has the block matrix form

$$
	\boldsymbol{P}=
	\left(
	\begin{array}{c|c}
	\boldsymbol{Q} & \ast & \ast & \cdots & \ast \\ \hline 
	\boldsymbol{O} & \boldsymbol{P_1} & \boldsymbol{O} &\cdots & \boldsymbol{O}\\
  \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{P_2} &\cdots & \boldsymbol{O} \\
  \vdots & \vdots & \vdots &\ddots & \vdots \\
  \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} &\cdots & \boldsymbol{P_m}
	\end{array}
	\right)
$$

其中$\boldsymbol{O}=(p_{ij}=0),\boldsymbol{Q}=(p_{ij})_{i,j \in T},\boldsymbol{P_l}=(p_{ij})_{i,j \in R_l},l=1,2,\cdots,m$

Ex3.14 具体 case;

更进一步有：

$$
	\lim_{n\rightarrow\infty} \boldsymbol{P}^n=
	\left(
	\begin{array}{c|c}
	\boldsymbol{O} & \ast & \ast & \cdots & \ast \\ \hline 
	\boldsymbol{O} & \lim_{n\rightarrow\infty}\boldsymbol{P_1}^n & \boldsymbol{O} &\cdots & \boldsymbol{O}\\
  \boldsymbol{O} & \boldsymbol{O} & \lim_{n\rightarrow\infty}\boldsymbol{P_2}^n &\cdots & \boldsymbol{O} \\
  \vdots & \vdots & \vdots &\ddots & \vdots \\
  \boldsymbol{O} & \boldsymbol{O} & \boldsymbol{O} &\cdots & \lim_{n\rightarrow\infty}\boldsymbol{P_m}^n
	\end{array}
	\right)
$$

---
edit: 2020-05-05

### 3.4 Irreducible Markov Chains

#### 3.4.1 定理

假设地震台通过统计数据得到了一个不同级别的地震的状态概率矩阵，政府和大众关心的问题是，从上一次地震到即将发生的下一次地震，大概需要多少年？

一个显然的直觉上的感觉是，根据地震的强度而定。不同强度不同。

可以用马尔可夫链作为上述问题的数学模型。

回顾停时的定义$T_j = \min\{n>0:X_n=j\}$，有如下定理：

> **Theorem 3.6 Limit Theorem for Finite Irreducible Markov Chains**. Assume that $X_0,X_1,\dots$, is a finite irreducible Markov Chain. For each state $j$. let $\mu_j=E(T_j|X_0=j)$ be the expected return time to $j$. Then, $\mu_j$ is finite, and there exisits a unique, positve stationary distribution $\boldsymbol{\pi}$ such that 
> $$ \pi_j =\frac{1}{\mu_j} , \forall j$$
> Furthermore, $\forall i$
> $$ \pi_j = \lim_{n\rightarrow\infty}\frac{1}{n}\sum_{m=0}^{n-1}P^m_{ij} $$


Ex 3.15 地震的 recurrence的例子；

Ex 3.16 Frog jumping 算例；

#### 3.4.2 First-Step Analysis

Ex 3.17 一个计算Expected return time的例子；

### 3.5 Periodicity

如果该马尔可夫链没有平稳分布，该如何研究？

具体的案例可以看第一章最后的那几个算例。

答案是引入周期的概念。

> **Lemma 3.7** The states of a communication class all have the same period

根据这条引理，有

> A Markov chain is called perioidic if is irreducible and all states have period greater than 1. A Markov chain is called aperiodic if it is irreducible and all states have period equal 1.

### 3.6 Ergoidc Markov Chains 

A Markov chain is called ergodic if it is irreducible, aperiodic, and all states have finite expected return times.

> **Fundamental Limit Theorem for Ergodic Markov Chains**
> **Therem 3.8** Let $X_0,X_1,\cdots$ be an ergodic Markov chain. There exists a unique, positive, stationary distribution $\boldsymbol{\pi}$,which is the limiting distribution of the chain. That is, 
> $$\pi_j = \lim_{n\rightarrow\infty}P^n_{ij},\forall i,j$$

Ex 3.19（Modified Ehrenfrest Model）; 

Ex 3.20 算例;

Ex 3.21 PageRank; (damping factor,p=0.85)

---
edit: 2020-05-06

### 3.7 Time reversibility

#### 3.7.1 Definition

The property of time reversibility can be explained intuitively as follows. If you were to take a movie of Markov chain moving forward in time and then run the movie backwards, you could not tell the difference between the two.

换成数学上的语言就是，如果假设马尔可夫链处于稳态，这时存在：

$$ P(X_0=i,X_1=j)=P(X_0=j,X_1=i) $$

由全概率公式可知；

> **Time Reversibility**
> 
> An irreducible Markov chain with transition matrix P and stationary distribution $\boldsymbol{\pi}$, if
> $$ \pi_iP_{ij}=\pi_jP_{ji},\forall i,j $$


Ex 3.23; 

Ex 3.24；Simple random walk on a graph is time 



#### 3.7.2 Reversible Markov Chains and Radom walk

Every reversible Markov chain can be considered as a random walk on a weighted graph

Ex 3.25 算例；

#### 3.7.3 The key benifit of reversibility

> **Proposition 3.9** Let $\boldsymbol{P}$ be the transition matrix of a Markov chain. If $\boldsymbol{x}$ is a probability distribution which satisfies
> $$ x_iP_{ij}=x_jP_{ji},\forall i,j $$
> then $\boldsymbol{x}$ is the stationary distribution, and the markov chain is reversible.

Ex 3.26 **Birth-and-death chain**

Case: random walk with a partialy relecting boundary.

### 3.8 Absorbing Chains

#### 3.8.1 定义

> **Absorbing State, Absorbing Chain**
> 
> State $i$ is an absorbing state if $P_{ii}=1$. A Markov chain is called an absorbing chain if it has at leat one absorbing state.

根据这个定义，一个吸收的马尔可夫链的canoical decompostion可以写为：
$$
  \boldsymbol{P} = \lim_{n\rightarrow\infty} \boldsymbol{P}^n=
	\left(
	\begin{array}{c|c}
	\boldsymbol{Q} & \boldsymbol{R} \\
  \hline
  \boldsymbol{O} & \boldsymbol{I}
	\end{array}
	\right)
$$

Ex 3.31

#### 3.8.2 Expected Number of Visits to Transient States

> **Theorem 3.11** Consider an absorbing Markov chain with t transient states. Let $\boldsymbol{F}$ be a $t\times t$ matrix indexed by transient states. where $F_{ij}$ is the expected number of visits to $j$ given that the chain starts in $i$, Then, 
> $$F=(I-Q)^{-1}$$

#### 3.8.3 Expected Time to Absorption

> **Absorbing Markov Chains**
> For an absorbing Markov chain with all states either transient or absorbing. Let $F=(I-Q)^{-1}$
> 1. (Absorption probability) The probability that from transient state $i$ the chain is absorbed in state $j$ is $(FR)_{ij}$
> 2. (Absorption time) The expected number of steps from transient state $i$ until the chain is absorbed in some absorbing state is $(F1)_i$  

#### 3.8.4 Expected Hitting Time for Irreducible chain

#### 3.8.5 Patterns in Sequences

#### 3.9 Regenration and strong Markov property

#### 3.10 Proofs of limiting Theorem

剩下的都是常规内容。不再赘述。

总的来说，这一章的章节组织很有条理，对初学者友好。而且选的例子很有启发性。

---
题外话：

有个值得思考的问题，这种偏向应用的数学内容的教材，如何平衡论述的逻辑，理论以及思考的深度，以及应用的广度，以及对于读者的吸引性和实用性，都是很需要功底的。但是这方面做的好的教材真的太少了。

理想的大学教师，是学术成就和教学成就都很出色，但是这毕竟是少数。

有一对儿相互矛盾的命题：

+ 大学生应该提高自学能力，不要指望老师手把手的教？
+ 国家给了大学老师工资，大学生也付了学费，如果都靠自学的话，要大学干什么？

依笔者看来，大学提供的是一套适合学习和研究的硬件设施，一张平静的书桌，一群志同道合的良师益友。这些环境和人，是任何其它机构或者网课代替不了的。