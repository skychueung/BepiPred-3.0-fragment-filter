# BepiPred-3.0-fragment-filter

> 本仓库在原始 BepiPred-3.0 基础上新增了两个可选的后处理补丁模块：
> 1. **片段过滤 (Fragment Filter)** — 将预测表位筛选为 5–15 aa 的候选短肽，支持 keep、split 和 split-ranking 三种模式
> 2. **理化性质过滤 (Physicochemical Filter)** — 在片段过滤基础上进一步按电荷、疏水性、pI、半胱氨酸等理化指标筛选

This repository adds two optional post-processing patches to BepiPred-3.0 for extracting and filtering candidate peptide fragments after B-cell epitope prediction.

---

## Overview

BepiPred-3.0 predicts linear B-cell epitope propensities from protein sequences.  
In many downstream peptide-screening tasks, the predicted parent peptides are too long for direct use.  
This patch adds two post-processing modules:

1. **Fragment Filter** — generate candidate fragments within a specified length range
2. **Physicochemical Filter** — further filter fragments by physicochemical properties (net charge, GRAVY, pI, cysteine count, disulfide risk)

---

## Modules

### 1. Fragment Filter (`bp3/fragment_filter.py`)

Post-processing module for fragment length filtering and ranking.

#### Modes

| Mode | Description |
|------|-------------|
| `keep` | Retain predicted parent peptides whose lengths already fall within the target range |
| `split` | Generate overlapping subfragments from longer parent peptides using a sliding window |
| `split + top_n` | Split mode + keep only the top N ranked fragments per parent peptide |

#### Ranking Heuristic

Fragments are ranked by a heuristic score based on:
- Aromatic residues (F, W, Y) — positive contribution
- Positively charged residues (K, R, H) — positive contribution
- Glycine/Serine/Proline/Asn/Gln — negative penalty
- Offset within parent peptide — slight penalty

#### Output

| Mode | Output Filename |
|------|----------------|
| keep | `fragment_candidates_keep_{min}_{max}.tsv` |
| split | `fragment_candidates_split_{min}_{max}.tsv` |
| split + top_n | `fragment_candidates_split_{min}_{max}_top{N}.tsv` |

---

### 2. Physicochemical Filter (`bp3/physicochemical_filter.py`)

Chains after fragment filter output to apply physicochemical property filtering.

#### Calculated Properties

Each fragment is analyzed for:

| Property | Description |
|----------|-------------|
| `Net_Charge_pH7.4` | Net charge at pH 7.4 (calculated from Henderson-Hasselbalch pKa values) |
| `GRAVY` | Grand Average of Hydropathicity (Kyte-Doolittle scale via BioPython) |
| `pI` | Isoelectric point (via BioPython ProtParam) |
| `Cys_Count` | Number of cysteine residues |
| `Disulfide_Risk` | Disulfide bond risk prediction: none / low / medium / higher / high |

#### Filter Criteria

| Parameter | CLI Flag | Default | Description |
|-----------|----------|---------|-------------|
| Min charge | `--phys_min_charge` | 0.0 | Minimum net charge at pH 7.4 |
| Max GRAVY | `--phys_max_gravy` | 0.0 | Maximum GRAVY (lower = more hydrophilic) |
| Min pI | `--phys_min_pi` | None | Minimum isoelectric point |
| Max pI | `--phys_max_pi` | None | Maximum isoelectric point |
| Max Cys | `--phys_max_cys` | None | Maximum cysteine count |
| Disulfide risk | `--phys_exclude_high_disulfide_risk` | False | Exclude "higher" and "high" risk peptides |

#### Output

Output filename: `{fragment_output_name}_physicochemical.tsv`

Example: `fragment_candidates_keep_5_15_physicochemical.tsv`

---

## CLI Usage

### Prerequisites

```bash
pip install torch fair-esm biopython pandas plotly
```

### Basic BepiPred-3.0 prediction

```bash
python bepipred3_CLI.py -i input.fasta -o output_dir -pred vt_pred
```

### A. Fragment filter — keep mode

Keep only parent peptides within length range 5–15:

```bash
python bepipred3_CLI.py \
  -i P11311_full_sequence.fasta \
  -o output_test_keep \
  -pred vt_pred \
  --fragment_filter \
  --frag_input_table iedb_p1_results.tsv \
  --frag_min_len 5 \
  --frag_max_len 15 \
  --frag_mode keep
```

### B. Fragment filter — split + top N mode

Split long peptides into fragments, keep top 3 per parent:

```bash
python bepipred3_CLI.py \
  -i P11311_full_sequence.fasta \
  -o output_test_split \
  -pred vt_pred \
  --fragment_filter \
  --frag_input_table iedb_p1_results.tsv \
  --frag_min_len 5 \
  --frag_max_len 15 \
  --frag_mode split \
  --frag_step 1 \
  --frag_top_n 3
```

### C. Fragment filter + physicochemical filter

Chain both filters — keep mode then physicochemical screening:

```bash
python bepipred3_CLI.py \
  -i P11311_full_sequence.fasta \
  -o output_test_phys \
  -pred vt_pred \
  --fragment_filter \
  --frag_input_table iedb_p1_results.tsv \
  --frag_min_len 5 \
  --frag_max_len 15 \
  --frag_mode keep \
  --physicochemical_filter \
  --phys_min_charge 0 \
  --phys_max_gravy 0 \
  --phys_min_pi 8.0 \
  --phys_max_pi 11.5 \
  --phys_max_cys 1 \
  --phys_exclude_high_disulfide_risk
```

---

## CLI Arguments Reference

### Fragment Filter Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--fragment_filter` | flag | — | Enable fragment filter post-processing |
| `--frag_input_table` | Path | auto-detect | Explicit input peptide table (.tsv/.csv) |
| `--frag_min_len` | int | 5 | Minimum fragment length |
| `--frag_max_len` | int | 15 | Maximum fragment length |
| `--frag_mode` | str | keep | `keep` or `split` |
| `--frag_step` | int | 1 | Sliding window step for split mode |
| `--frag_top_n` | int | None | Keep top N fragments per parent (split mode) |
| `--frag_best_per_parent` | flag | False | Keep only the best fragment per parent |
| `--frag_output_name` | str | auto | Custom output filename |

### Physicochemical Filter Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--physicochemical_filter` | flag | — | Enable physicochemical filter (requires `--fragment_filter`) |
| `--phys_min_charge` | float | 0.0 | Minimum net charge at pH 7.4 |
| `--phys_max_gravy` | float | 0.0 | Maximum GRAVY score |
| `--phys_min_pi` | float | None | Minimum isoelectric point |
| `--phys_max_pi` | float | None | Maximum isoelectric point |
| `--phys_max_cys` | int | None | Maximum cysteine count |
| `--phys_exclude_high_disulfide_risk` | flag | False | Exclude high/higher disulfide risk peptides |

---

## File Structure

```
BepiPred-3.0-fragment-filter/
├── bepipred3_CLI.py              # Main CLI entry point
├── bp3/
│   ├── __init__.py
│   ├── bepipred3.py              # Core BepiPred-3.0 prediction engine
│   ├── fragment_filter.py        # Fragment filter module (patched)
│   ├── physicochemical_filter.py # Physicochemical filter module (patched)
│   └── BP3Models/                # Pre-trained model weights
├── esm_encodings/                # ESM-2 cached encodings
├── example_antigens/             # Example input files
├── requirements.txt
└── README.md
```

---

## Bug Fixes

### `fragment_filter.py` — missing `top_n` parameter (fixed)

The CLI passes `top_n=frag_top_n` to `run_fragment_filter()`, but the original function signature did not accept this parameter, causing a `TypeError` when using `--frag_top_n`.

**Fix applied:**
- Added `top_n: Optional[int] = None` parameter to `run_fragment_filter()`
- Added `top_n_per_parent()` function for top-N-per-parent filtering
- Updated output filename generation to include `_top{N}` suffix when `top_n` is set

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `torch` | PyTorch backend for ESM-2 and BP3 models |
| `fair-esm` | ESM-2 protein language model for sequence encoding |
| `biopython` | ProteinAnalysis for GRAVY and pI calculation |
| `pandas` | TSV/CSV table reading and processing |
| `plotly` | Interactive prediction plots |

---

## Tested Scenarios

| Test | Mode | Output | Status |
|------|------|--------|--------|
| A | keep | `fragment_candidates_keep_5_15.tsv` | ✅ Pass |
| B | split + top_n=3 | `fragment_candidates_split_5_15_top3.tsv` | ✅ Pass |
| C | keep + physicochemical | `fragment_candidates_keep_5_15_physicochemical.tsv` | ✅ Pass |

All physicochemical columns verified: `Net_Charge_pH7.4`, `pI`, `GRAVY`, `Cys_Count`, `Disulfide_Risk`.
