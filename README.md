# BepiPred-3.0-fragment-filter

> 本仓库在原始 BepiPred-3.0 基础上新增了一个可选的片段过滤补丁，用于将预测表位进一步筛选为 5–15 aa 的候选短肽，支持 keep、split 和 split-ranking 三种模式。

This repository adds an optional post-processing fragment filter to BepiPred-3.0 for extracting candidate peptide fragments with user-defined length constraints after B-cell epitope prediction.

---

## Overview

BepiPred-3.0 predicts linear B-cell epitope propensities from protein sequences.  
In many downstream peptide-screening tasks, the predicted parent peptides are too long for direct use.  
This patch adds a post-processing module to generate candidate fragments within a specified length range.

This repository is intended as a lightweight extension of the original BepiPred-3.0 workflow for short peptide candidate extraction and prioritization.

---

## Modes

This patch currently supports three post-processing modes:

- **keep**: retain predicted parent peptides whose lengths already fall within the target range
- **split**: generate overlapping 5–15 aa subfragments from longer parent peptides
- **split-ranking**: rank split-generated fragments and keep only the top N candidates per parent peptide

---

## Added functionality

A new module was added:

```text
bp3/fragment_filter.py
