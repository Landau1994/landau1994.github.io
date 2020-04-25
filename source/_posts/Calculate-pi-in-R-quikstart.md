---
title: Calculate pi in R quikstart
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - R
keywords:
  description: null
date: 2020-04-24 22:25:00
description:
photos:
---


R中也可以用`Rmpfr`包实现多精度的计算。例如，我们可以用如下代码实现AGM算法计算Pi到小数点后256位。

``` r
library(Rmpfr)
piMpfr <- function(prec=256, itermax = 100, verbose=TRUE) {
  m2 <- mpfr(2, prec) # '2' as mpfr number
  ## -> all derived numbers are mpfr (with precision 'prec')
  p <- m2 + sqrt(m2) # 2 + sqrt(2) = 3.414..
  y <- sqrt(sqrt(m2)) # 2^ {1/4}
  x <- (y+1/y) / m2
  it <- 0L
  repeat {
    p.old <- p
    it <- it+1L
    p <- p * (1+x) / (1+y)
    if(verbose) cat(sprintf("it=%2d, pi^ = %s, |.-.|/|.|=%e\n",
                            it, formatMpfr(p, min(50, prec/log2(10))), 1-p.old/p))
    if (abs(p-p.old) <= m2^(-prec))
      break
    if(it > itermax) {
      warning("not converged in", it, "iterations") ; break
    }
    ## else
    s <- sqrt(x)
    y <- (y*s + 1/s) / (1+y)
    x <- (s+1/s)/2
  }
  p
}
piMpfr(prec = 256)
```

    ## it= 1, pi^ = 3.1426067539416226007907198236183018919713562462772, |.-.|/|.|=-8.642723e-02
    ## it= 2, pi^ = 3.1415926609660442304977522351203396906792842568645, |.-.|/|.|=-3.227958e-04
    ## it= 3, pi^ = 3.1415926535897932386457739917571417940347896238675, |.-.|/|.|=-2.347934e-09
    ## it= 4, pi^ = 3.1415926535897932384626433832795028841972241204666, |.-.|/|.|=-5.829228e-20
    ## it= 5, pi^ = 3.1415926535897932384626433832795028841971693993751, |.-.|/|.|=-1.741826e-41
    ## it= 6, pi^ = 3.1415926535897932384626433832795028841971693993751, |.-.|/|.|=0.000000e+00

    ## 1 'mpfr' number of precision  256   bits 
    ## [1] 3.141592653589793238462643383279502884197169399375105820974944592307816406286163
