---
title: Dobrow_chap1
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
date: 2020-04-04 19:41:57
mathjax: true
categories:
    - math
tags:
    - stochastic Process
    - note
keywords:
description: This is a note of the book `Introduction to stochastic processes with R`
photos:
---
This is a note of the book `Introduction to stochastic processes with R`

### Introduction and preview

> We demand rigidly defined areas of doubt and uncertatinty  –Douglas Adams, The Hitchhiker’s Guide to the Galaxy

#### 1.1 DETERMINISTIC AND STOCHASTIC MODELS

+ Consider a simple exponential growthmodel
where the random arises？
The deterministic model does not address the uncertainty present in the reproduction
rate of individual organisms.

+ In many biological processes, the exponential distribution is a common choice for modeling the times of births and deaths.

+ Ex1.1 PageRank: random walks on graphs

+ Ex1.2 Spread of infectious disease
SIR model: Susceptible-infected-removed
Reed-Frost model: Stochastic SIR model in discrete time.

#### 1.2 What is a stochastic process

+ The author said
    > A stochastic process, also called a random process, is simply one in which outcomes are uncertain. By contrast, in a deterministic system there is no randomness. In a deterministic system, the same output is always produced from a given input.

+ A stochastic process is specified by its index and state sapce, and by the dependency relationships among its random variables
    > Stochastic process: A stochastic process is a collection of random variables $\{X_t,t \in I\}$.The set I is the index set of the process. The random variables are defined on a commmon state space S.

+ Ex 1.3 Monopoly
+ EX 1.4 Discrete time, continous state space
+ Ex 1.5 Continuous time, discrete state space
: arrival process, Poisson process
+ Ex 1.6 Random walk and gambler's ruin: Random walk, discrete-time stochastic process whose state space is $\mathbb{Z}$
+ Ex 1.7 Brownian motion: Brownian motion is a continuous-time, contiuous state space stochastic process

#### 1.3 Monte Carlo Simulation

+ Given a random experiment and event A, a Monte Carlo estimate of $P(A)$ is obtained by repeating random experiment many times and taking the proportion of trials in which A occurs as an approximation for $P(A)$

+ Strong law of large numbers
     
     $$ \lim_{n\rightarrow\infty}\frac{X_1+\dots+X_n}{n}=P(A), a.s. \tag{1.1}$$

#### 1.4 Conditional Probability

> The simplest stochastic process is a sequence of i.i.d. random variables. Such
sequences are often used to model random samples in statistics. However, most
real-world systems exhibit some type of dependency between variables, and an
independent sequence is often an unrealistic model.

> Thus, the study of stochastic processes really begins with conditional
probability—conditional distributions and conditional expectation. These will
become essential tools for all that follows.

+ Conditional Probability: 
  
  $P(A|B) = \frac{P(A\bigcap B)}{P(B)}$

+ Law of Total probability:
  
  Let $B_1,\dots,B_k$ be a sequence of events that partition the sample space. That is, the $B_i$ are mutually exclusive(disjoint) and their union is equal to $\Omega$. Then, for many event A,
  $$ P(A)=\sum_{i=1}^{k}P(A\bigcap B_i)=\sum_{i=1}^{k}P(A|B_i)P(B_i) $$

+ Ex1.8 Disease tests

+ Ex1.9 Find the probability that it is a heart

+ Ex1.10 Gambler's ruin
  let $p_k$ denote the probability of reaching n when the gambler's fortune is k.
  $$ p_k = p_{k+1}(\frac{1}{2})+p_{k-1}(\frac{1}{2}) $$
  or 
  $$ p_{k+1}-p_k = p_k - p_{k-1},  k = 1,\dots, n-1 \tag{1.2}$$
  
  using $p_0 =0, p_n = 1$

  we have $p_k=kp_1=\frac{k}{n},k=0,\dots,n$
  
  The gambler's ruin is \frac{n-k}{n} 

+ Bayes Rule

  Given a countable sequence of events $B_1,B_2,\dots$ which partition the sample space, a more general form of Bayes' rule is

  $$P(B_i|A)=\frac{P(A|B_i)P(B_i)}{\sum_jP(A|B_j)P(B_j)}$$ 

+ Ex 1.11 The probability that teh employee is in fact lying.

+ Conditional Distribution
  
  joint density function:

  $P(X \le x, Y \le y) = \int_{-\infty}^{x}\int_{-\infty}^{y}f(x,t)dtds$
  
  + Discrete Case:
    
    P(Y=y|X=x)=\frac{X=x,Y=y}{P(X=x)}
    
  + Continuous Case:
  
    For continuous randomo variables X and Y, the conditional density function of Y given X=x is 

    $f_{Y|X}(y|x)=\frac{f(x,y)}{f_X(x)}$

    $P(Y \in R | X = x)=\int_Rf_{Y|X}(y|x)dy$

#### 1.5 Conditional Expectation of Y given X=x

+ Conditional Expectation of Y given X=x

    $$E(Y|X=x)= \begin{cases}\sum_{y}yP(Y=y|X=x),&\text{discrete} \\
    \int_{-\infty}^{\infty}yf_{Y|X}(y|x)dy, & \text{continuous}   
    \end{cases} $$