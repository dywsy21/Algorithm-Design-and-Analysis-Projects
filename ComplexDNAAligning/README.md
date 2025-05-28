---
title: "复杂DNA序列比对算法 实验报告"
author: "王思宇"
date: \today
using_title: true
using_table_of_content: true
---

# 复杂DNA序列比对算法 实验报告

## 项目地址

本项目托管在GitHub上：[dywsy21/Algorithm-Design-and-Analysis-Projects](https://github.com/dywsy21/Algorithm-Design-and-Analysis-Projects)

## 项目概述

本项目实现了一个高效的DNA序列比对算法，能够处理长序列和短序列的比对需求。算法基于k-mer锚点检测、动态规划延展和智能间隙填充策略，在保证比对精度的同时优化了覆盖率。

## 算法整体思路

使用的算法本质上是图匹配算法，将DNA序列比对问题转化为在序列图中寻找最优路径的问题：

1. 将参考序列和查询序列的k-mer作为图中的节点
2. 相同的k-mer之间建立连接边，形成潜在的比对锚点
3. 通过延展锚点找到连续的比对路径
4. 选择非重叠的最优路径组合，最大化覆盖率

这种图匹配方法避免了传统动态规划的二次时间复杂度，同时保持了高精度的比对质量。

## 算法核心思想

### 1. 分层策略设计

算法根据查询序列长度自适应选择不同的比对策略：

1. 长序列策略（>3000bp）：注重稳定性和全局覆盖
2. 短序列策略（≤3000bp）：注重敏感性和局部精度

### 2. 多k-mer锚点检测

使用多种k-mer大小进行锚点检测，兼顾敏感性和特异性：

- 长序列：k=[8,15]，适中采样（shift=1）
- 短序列：k=[6,8,10,12,15]，密集采样（shift=1）

#### 锚点检测伪代码

伪代码如下：

```
GENERATE-ANCHORS(ref, query, kmer_sizes)
  anchors ← ∅
  for each k ∈ kmer_sizes
      do anchor_set ← HASHTABLE-KMER-SEARCH(ref, query, k, 1)
         anchors ← anchors ∪ anchor_set
  return anchors
```

```
HASHTABLE-KMER-SEARCH(ref, query, k, shift)
  hash_table ← BUILD-KMER-HASHTABLE(ref, k)
  rc_query ← REVERSE-COMPLEMENT(query)
  anchors ← ∅
  for i ← 0 to length[query] - k by shift
      do kmer ← query[i..i+k-1]
         if HASH(kmer) ∈ hash_table
             then for each pos ∈ hash_table[HASH(kmer)]
                      do anchors ← anchors ∪ {(i, pos, +1, k)}
         rc_kmer ← rc_query[length[query]-i-k..length[query]-i-1]
         if HASH(rc_kmer) ∈ hash_table
             then for each pos ∈ hash_table[HASH(rc_kmer)]
                      do anchors ← anchors ∪ {(i, pos, -1, k)}
  return anchors
```

### 3. 保守延展策略

对每个锚点进行双向延展，严格控制错误匹配率：

- 实时监控错误匹配率，超过阈值立即停止延展
- 正向链和反向互补链分别处理

#### 锚点延展伪代码

伪代码如下：

```
EXTEND-ANCHOR(query, ref, qpos, rpos, strand, init_len, params)
  if strand = +1
      then return EXTEND-FORWARD-STRAND(query, ref, qpos, rpos, init_len, params)
      else return EXTEND-REVERSE-STRAND(query, ref, qpos, rpos, init_len, params)
```

```
EXTEND-FORWARD-STRAND(query, ref, qpos, rpos, init_len, params)
   left_ext ← 0
   mismatches ← 0
   while qpos - left_ext > 0 and rpos - left_ext > 0 and left_ext < params.extension_limit
       do if query[qpos - left_ext - 1] = ref[rpos - left_ext - 1]
              then left_ext ← left_ext + 1
              else mismatches ← mismatches + 1
                   if left_ext > params.mismatch_threshold and 
                      mismatches / left_ext > params.mismatch_rate
                       then break
                   left_ext ← left_ext + 1
   right_ext ← EXTEND-RIGHT-DIRECTION(query, ref, qpos, rpos, init_len, params)
   return (qpos - left_ext, qpos + init_len + right_ext, 
           rpos - left_ext, rpos + init_len + right_ext)
```

## 详细算法实现

### 阶段一：锚点生成

符号说明:
- n：参考序列长度
- m：查询序列长度  
- k：k-mer数量

**时间复杂度**：O(n + m·k)

### 阶段二：锚点延展

参数设置为: 
- 短序列：延展限制80bp，错误率阈值8%
- 长序列：延展限制200bp，错误率阈值5%

符号说明: L：平均延展长度
**时间复杂度**：O(L)

### 阶段三：片段选择与优化

#### 片段选择伪代码

伪代码如下：

```
SELECT-NON-OVERLAPPING-SEGMENTS(segments)
   valid_segments ← ∅
   for each segment ∈ segments
       do if VALIDATE-SEGMENT(segment)
              then valid_segments ← valid_segments ∪ {segment}
   SORT(valid_segments, key = query_start)
   selected ← ∅
   last_end ← 0
   for each segment ∈ valid_segments
       do if segment.query_start ≥ last_end
              then selected ← selected ∪ {segment}
                   last_end ← segment.query_end
   return selected
```

```
VALIDATE-SEGMENT(query, ref, qstart, qend, rstart, rend)
   if qend - qstart < MIN_SEGMENT_LENGTH
       then return FALSE
   edit_distance ← CALCULATE-EDIT-DISTANCE(ref, query, rstart, rend, qstart, qend)
   edit_rate ← edit_distance / (qend - qstart)
   return edit_rate ≤ MAX_EDIT_RATE
```

验证条件的参数设置为：
- 片段长度 ≥ 30bp
- 编辑距离率 ≤ 10%
- 无查询位置重叠

### 阶段四：间隙填充

保守的间隙填充策略，仅在精确匹配时合并片段：

#### 间隙填充伪代码

伪代码如下：

```
CONSERVATIVE-GAP-FILLING(segments, query, ref, params)
   filled ← ∅
   for i ← 1 to length[segments]
       do segment ← segments[i]
          if i = 1 and segment.query_start > 0
              then segment ← TRY-EXTEND-TO-START(segment, query, ref, params)
          if filled ≠ ∅
              then gap_size ← segment.query_start - filled[last].query_end
                   if 0 < gap_size ≤ params.gap_fill_max and EXACT-MATCH-CHECK(...)
                       then filled[last] ← MERGE-SEGMENTS(filled[last], segment)
                             continue
          filled ← filled ∪ {segment}
   if filled ≠ ∅ and filled[last].query_end < length[query]
       then filled[last] ← TRY-EXTEND-TO-END(filled[last], query, ref, params)
   return VALIDATE-NO-OVERLAPS(filled)
```

## 复杂度分析

### 总的符号说明

- n：参考序列长度
- m：查询序列长度
- k：k-mer总数量
- A：锚点数量
- L：平均延展长度
- S：有效片段数量
- V：单个片段验证时间
- G：平均间隙大小

### 时间复杂度

1. 锚点生成：O(n + m·k)
   - 构建哈希表：O(n)
   - 查询匹配：O(m·k)

2. 锚点延展：O(A·L)
   - A为锚点数量，L为平均延展长度
   - 实际中A << m，L << n

3. 片段选择：O(S·log S + S·V)
   - S为有效片段数，V为验证时间
   - 排序：O(S·log S)
   - 验证：O(S·V)

4. 间隙填充：O(S·G)
   - G为平均间隙大小，通常很小

**总时间复杂度**：O(n + m·k + A·L + S·V)

在实际应用中，由于A << m，S << A，算法表现接近线性时间。

## 性能表现

在标准测试数据集上的表现：

- **长序列（T1）**：得分29710/29820（99.6%）
- **短序列（T2）**：得分1983/2090 (94.9%)

两序列均仅差100分达到基线，实现了优秀的覆盖率和运行效率。

## 参数配置

所有关键参数已提取到`config.py`中，便于针对不同应用场景进行调优.

## 使用方法

直接跑`python main.py`即可.

## 总结

本算法通过图匹配的思想，结合分层策略设计、多k-mer锚点检测、保守延展和智能间隙填充，实现了高效准确的DNA序列比对。算法具有良好的时间复杂度、空间效率和实际性能，适用于各种规模的序列比对任务。
