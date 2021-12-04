---
title: learn_python_001
author: 夏目沉吟
avatar: /images/faceicon.png
authorLink: 'https://github.com/Landau1994'
authorAbout: 'https://github.com/Landau1994'
authorDesc: A PhD student in bioinformatics
mathjax: true
categories:
  - implementation
tags:
  - note
  - Python
date: 2021-12-04 21:19:59
keywords:
description:
photos:
---


# Learn_python_001 Top K problem

### &#x1F7E7;&#x2753; The problem

Top K question:

输入整数数组 arr ，找出其中最小的 k 个数。例如，输入4、5、1、6、2、7、3、8这8个数字，则最小的4个数字是1、2、3、4。

示例 1：

输入：arr = [3,2,1], k = 2
输出：[1,2] 或者 [2,1]
示例 2：

输入：arr = [0,1,2,1], k = 1
输出：[0]
 
限制：

0 <= k <= arr.length <= 10000
0 <= arr[i] <= 10000

来源：力扣（LeetCode）
链接：https://leetcode-cn.com/problems/zui-xiao-de-kge-shu-lcof
著作权归领扣网络所有。商业转载请联系官方授权，非商业转载请注明

### &#x1F4C4; code

有三种解法，解法二，三，是符合题目要求的两种（因为题目也考察了排序算法）。详解见https://leetcode-cn.com/problems/zui-xiao-de-kge-shu-lcof/solution/jian-zhi-offer-40-zui-xiao-de-k-ge-shu-j-9yze/

```python
class Solution1:
    def getLeastNumbers(self, arr: List[int], k: int) -> List[int]:
        lstStd = arr
        lstStd.sort()
        res = lstStd[:k]
        return res
### solution2:
### wirte your own quick sort
### This was taken from https://leetcode-cn.com/problems/zui-xiao-de-kge-shu-lcof/solution/jian-zhi-offer-40-zui-xiao-de-k-ge-shu-j-9yze/
class Solution2:
    def getLeastNumbers(self, arr: List[int], k: int) -> List[int]:
        def quick_sort(arr, l, r):
            if l >= r: return
            i, j = l, r
            while i < j:
                while i < j and arr[j] >= arr[l]: j -= 1
                while i < j and arr[i] <= arr[l]: i += 1
                arr[i], arr[j] = arr[j], arr[i]
            arr[l], arr[i] = arr[i], arr[l]
            quick_sort(arr, l, i - 1)
            quick_sort(arr, i + 1, r)
        quick_sort(arr, 0, len(arr) - 1)
        return arr[:k]
class Solution3:
    def getLeastNumbers(self, arr: List[int], k: int) -> List[int]:
        if k >= len(arr): return arr
        def quick_sort(l, r):
            i, j = l, r
            while i < j:
                while i < j and arr[j] >= arr[l]: j -= 1
                while i < j and arr[i] <= arr[l]: i += 1
                arr[i], arr[j] = arr[j], arr[i]
            arr[l], arr[i] = arr[i], arr[l]
            if k < i: return quick_sort(l, i - 1) 
            if k > i: return quick_sort(i + 1, r)
            return arr[:k]     
        return quick_sort(0, len(arr) - 1)
```
### &#x1F4cF; 测试

我们用如下代码测试：

```python
### import require package
import numpy as np
import random 
import matplotlib.pyplot as plt
import time
import seaborn as sns
from typing import List, Dict, Tuple, Sequence
def ProgramTime(N,func):
    lst = [random.randrange(10**7) for n in range(N)]
    start = time.perf_counter()
    func(lst,10)
    runtime = (time.perf_counter() - start)
    return runtime
ProgramTimeVec = np.vectorize(ProgramTime)

### define theory function
def f1(n, k):
    return k*n
def f2(n, k):
    return k*n*np.log(n)

### plot test curve
n = np.arange(1, 2000)
colors = sns.color_palette("Set1")

plt.plot(n, f1(n, 1e-7), c=colors[0])
plt.plot(n, f2(n, 1e-7), c=colors[1])
plt.plot(n, ProgramTimeVec(n,sol1.getLeastNumbers),c=colors[2])
plt.plot(n, ProgramTimeVec(n,sol2.getLeastNumbers),c=colors[3])
plt.plot(n, ProgramTimeVec(n,sol3.getLeastNumbers),c=colors[4])
plt.xlabel('Size of input (n)', fontsize=16)
plt.ylabel('Time', fontsize=16)
#plt.legend(['$\mathcal{O}(n^2)$', '$\mathcal{O}(n \log n)$'], loc='best', fontsize=20)
plt.legend(['$\mathcal{O}(n)$', '$\mathcal{O}(n \log n)$','sol1', 'sol2','sol3'], loc='best', fontsize=20)
fig = plt.gcf()
fig.set_size_inches(8, 6)
plt.savefig("../fig/test.png",dpi=300)
```
结果如下, 可以看出，使用解法三，也就是基于快速排序的数组划分，可以实现线性时间：

![figure](https://s3.bmp.ovh/imgs/2021/12/90470f3f3dcdc66c.png)
