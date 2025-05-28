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

我们实现的算法的时间复杂度在近似线性的范围内。

## 算法概述

使用的算法本质上是图匹配算法，将DNA序列比对问题转化为在序列图中寻找最优路径的问题：

1. 将参考序列和查询序列的k-mer作为图中的节点
2. 相同的k-mer之间建立连接边，形成潜在的比对锚点
3. 通过延展锚点找到连续的比对路径
4. 选择非重叠的最优路径组合，最大化覆盖率

这种图匹配方法避免了传统动态规划的二次时间复杂度，同时保持了高精度的比对质量。

## 详细算法实现

### 阶段一：多k-mer锚点检测

使用多种k-mer大小进行锚点检测，兼顾敏感性和特异性。小k-mer可提高敏感性，捕获短匹配，大k-mer可提高特异性，减少假阳性。这样的差异化采样可以平衡计算效率和检测精度。

#### 锚点检测伪代码

伪代码如下：

```
GENERATE-ANCHORS(ref, query, kmer_sizes)
  anchors ← empty set
  for each k in kmer_sizes
      do if k <= 5 then shift ← 3
         else if k <= 8 then shift ← 1  
         else shift ← 1
         anchor_set ← HASHTABLE-KMER-SEARCH(ref, query, k, shift)
         anchors ← anchors union anchor_set
  return anchors
```

```
HASHTABLE-KMER-SEARCH(ref, query, k, shift)
  hash_table ← BUILD-KMER-HASHTABLE(ref, k)
  rc_query ← REVERSE-COMPLEMENT(query)
  anchors ← empty set
  for i ← 0 to length[query] - k by shift
      do kmer ← query[i..i+k-1]
         if HASH(kmer) in hash_table
             then for each pos in hash_table[HASH(kmer)]
                      do anchors ← anchors union {(i, pos, +1, k)}
         rc_kmer ← rc_query[length[query]-i-k..length[query]-i-1]
         if HASH(rc_kmer) in hash_table
             then for each pos in hash_table[HASH(rc_kmer)]
                      do anchors ← anchors union {(i, pos, -1, k)}
  return anchors
```

### 阶段二：保守延展策略

对每个锚点进行双向延展，通过以下实时监控错误匹配率，超过阈值立即停止延展的方法严格控制错误匹配率，并且将正向链和反向互补链分别处理。短序列支持边界延展优化。

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
   
   // 短序列边界优化
   if length[query] < SHORT_THRESHOLD
       then if qpos - left_ext <= BOUNDARY_TOLERANCE
                then left_ext ← max(left_ext, min(qpos, rpos))
            if qpos + init_len + right_ext ≥ length[query] - BOUNDARY_TOLERANCE
                then right_ext ← max(right_ext, length[query] - qpos - init_len)
   
   return (qpos - left_ext, qpos + init_len + right_ext, 
           rpos - left_ext, rpos + init_len + right_ext)
```

### 阶段三：片段选择与优化

#### 片段选择伪代码

伪代码如下：

```
SELECT-NON-OVERLAPPING-SEGMENTS(segments)
   valid_segments ← empty set
   for each segment in segments
       do if VALIDATE-SEGMENT(segment)
              then valid_segments ← valid_segments union {segment}
   SORT(valid_segments, key = query_start)
   selected ← empty set
   last_end ← 0
   for each segment in valid_segments
       do if segment.query_start ≥ last_end
              then selected ← selected union {segment}
                   last_end ← segment.query_end
   return selected
```

```
VALIDATE-SEGMENT(query, ref, qstart, qend, rstart, rend)
   if qend - qstart < MIN_SEGMENT_LENGTH
       then return FALSE
   edit_distance ← CALCULATE-EDIT-DISTANCE(ref, query, rstart, rend, qstart, qend)
   edit_rate ← edit_distance / (qend - qstart)
   return edit_rate <= MAX_EDIT_RATE
```

### 阶段四：智能间隙填充

我们采用了基于序列类型的差异化间隙填充策略：

1. 短序列间隙填充：

- 小间隙容忍：<=20bp，允许15%错误率

2. 长序列间隙填充：

- 大间隙容忍：<=60bp，要求85%匹配度
- 合并片段验证：允许12%编辑距离率

#### 间隙填充伪代码

伪代码如下：

```
ADAPTIVE-GAP-FILLING(segments, query, ref, params)
   filled ← BASIC-GAP-FILLING(segments, query, ref, params)
   
   if length[query] < SHORT_THRESHOLD
       then return SHORT-SEQUENCE-GAP-FILLING(filled, query, ref)
       else return LONG-SEQUENCE-GAP-FILLING(filled, query, ref)
```

```
LONG-SEQUENCE-GAP-FILLING(segments, query, ref)
   final_filled ← empty set
   for each segment in segments
       do if final_filled ≠ empty set
              then gap_size ← segment.query_start - final_filled[last].query_end
                   if 0 < gap_size <= LARGE_GAP_MAX
                       then match_rate ← CALCULATE-MATCH-RATE(gap_region)
                            if match_rate ≥ LARGE_GAP_MATCH_RATE
                                then merged ← MERGE-WITH-VALIDATION(final_filled[last], segment)
                                     if VALID-MERGED-SEGMENT(merged)
                                         then final_filled[last] ← merged
                                              continue
          final_filled ← final_filled union {segment}
   return final_filled
```

### 自适应多策略设计

这是我们的参数部分。算法根据序列长度自动选择最优策略：

1. 长序列策略（>=3000bp）：

- k-mer大小：[7, 8, 9, 10, 14, 15, 18]
- 采样密度：shift=1（密集采样）
- 延展限制：198bp，错误率<=9%
- 间隙填充：最大60bp，要求85%匹配度

2. 短序列策略（<3000bp）：

- k-mer大小：[4, 5, 6, 7, 8, 9, 11, 13, 15]
- 采样策略：k<=5用shift=3，k<=8用shift=1，其他用shift=1
- 延展限制：74bp，错误率<=8%
- 间隙填充：最大40bp，支持边界延展

## 复杂度分析

### 符号说明

- n：参考序列长度
- m：查询序列长度
- k：k-mer总数量
- A：锚点数量
- L：平均延展长度
- S：有效片段数量
- V：单个片段验证时间
- G：平均间隙大小

### 各阶段时间复杂度

1. 锚点生成：$O(n + m·k)$
   - 构建哈希表：$O(n)$
   - 多k-mer查询匹配：$O(m·k)$

2. 锚点延展：$O(A·L)$
   - A为锚点数量，L为平均延展长度
   - 实际中A << m，L << n

3. 片段选择：$O(S·log S + S·V)$
   - S为有效片段数，V为验证时间
   - 排序：$O(S·log S)$
   - 验证：$O(S·V)$

4. 间隙填充：$O(S·G)$
   - G为平均间隙大小，通常很小

### 总时间复杂度

**总时间复杂度**：$O(n + m·k + A·L + S·V)$

在实际应用中，由于A << m，S << A，算法表现接近线性时间。

## 性能表现

算法输出结果如下：

测试1：

```
[(0, 380, 0, 380), (380, 783, 380, 783), (783, 1197, 783, 1197), (1197, 1608, 1197, 1608), (1608, 2014, 1608, 2014), (2014, 2419, 2014, 2419), (2419, 2830, 2419, 2830), (2830, 3240, 2830, 3240), (3240, 3645, 3240, 3645), (3645, 4048, 3645, 4048), (4048, 4453, 4048, 4453), (4453, 4858, 4453, 4858), (4858, 5272, 4858, 5272), (5272, 5675, 5272, 5675), (5675, 6078, 5675, 6078), (6078, 6489, 6078, 6489), (6489, 6903, 22926, 23340), (6903, 7286, 22543, 22926), (7287, 7511, 22318, 22542), (7511, 7851, 21978, 22318), (7852, 8093, 21736, 21977), (8093, 8497, 21332, 21736), (8497, 8902, 20927, 21332), (8902, 9307, 20522, 20927), (9307, 9717, 20112, 20522), (9717, 9931, 19898, 20112), (9932, 10225, 19604, 19897), (10225, 10630, 19199, 19604), (10630, 11034, 18795, 19199), (11034, 11445, 18384, 18795), (11445, 11855, 17974, 18384), (11855, 12194, 17635, 17974), (12194, 12600, 17228, 17634), (12600, 13003, 16825, 17228), (13003, 13413, 16415, 16825), (13413, 13817, 16011, 16415), (13817, 14223, 15605, 16011), (14223, 14629, 15199, 15605), (14629, 15040, 14788, 15199), (15040, 15371, 14457, 14788), (15372, 15695, 14133, 14456), (15695, 16100, 13728, 14133), (16100, 16510, 13318, 13728), (16510, 16920, 12908, 13318), (16920, 17323, 12505, 12908), (17323, 17726, 12102, 12505), (17726, 18129, 11699, 12102), (18129, 18540, 11288, 11699), (18540, 18951, 10877, 11288), (18951, 19362, 10466, 10877), (19362, 19678, 10150, 10466), (19679, 20037, 9791, 10149), (20037, 20451, 9377, 9791), (20451, 20857, 8971, 9377), (20857, 21262, 8566, 8971), (21262, 21665, 8163, 8566), (21665, 22075, 7753, 8163), (22075, 22486, 7342, 7753), (22486, 22900, 6928, 7342), (22900, 23303, 6525, 6928), (23303, 23717, 6111, 6525), (23717, 24127, 23718, 24128), (24127, 24532, 24128, 24533), (24532, 24938, 24533, 24939), (24938, 25344, 24939, 25345), (25344, 25755, 25345, 25756), (25755, 26165, 25756, 26166), (26165, 26575, 26166, 26576), (26575, 26985, 26576, 26986), (26985, 27390, 26986, 27391), (27390, 27801, 27391, 27802), (27801, 28211, 27802, 28212), (28211, 28621, 28212, 28622), (28621, 29031, 28622, 29032), (29031, 29441, 29032, 29442), (29441, 29829, 29442, 29830)]
```

测试2：

```
[(0, 195, 0, 195), (196, 236, 196, 236), (237, 282, 237, 282), (282, 334, 382, 434), (335, 401, 435, 501), (401, 503, 497, 599), (505, 558, 705, 758), (558, 719, 658, 819), (719, 810, 719, 810), (810, 903, 710, 803), (903, 1003, 697, 797), (1003, 1102, 703, 802), (1103, 1201, 803, 901), (1201, 1301, 899, 999), (1301, 1394, 901, 994), (1396, 1457, 396, 457), (1458, 1502, 458, 502), (1513, 1618, 1013, 1118), (1619, 1657, 1319, 1357), (1658, 1701, 1358, 1401), (1710, 1743, 1210, 1243), (1744, 1806, 1244, 1306), (1807, 1904, 1107, 1204), (1904, 1998, 1404, 1498), (2288, 2388, 1488, 1588), (2388, 2500, 1588, 1700)]
```

其中测试1得分29821，超过基线1分；测试2得分2065，距离基线仅25分，获得了98.8%的分数。

测试1大约需半分钟跑完，测试2半秒左右即可跑完。这也侧面映证了确实是接近线性的时间复杂度。

## 使用方法

直接运行`python main.py`即可开始序列比对。算法会自动根据序列长度选择最优策略。

## 总结

本算法通过图匹配的思想，结合自适应策略设计、多k-mer锚点检测、保守延展和智能间隙填充，实现了高效准确的DNA序列比对。算法具有良好的时间复杂度、空间效率和实际性能，适用于各种规模的序列比对任务。创新的自适应策略使算法能够针对不同长度的序列优化性能，达到了接近最优的比对效果。
