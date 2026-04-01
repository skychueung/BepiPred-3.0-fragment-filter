from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict
import csv
import re

try:
    import pandas as pd
except ImportError as e:
    raise ImportError("fragment_filter.py requires pandas. Please install pandas first.") from e


COLUMN_ALIASES: Dict[str, List[str]] = {
    "parent_rank": ["parent_rank", "rank", "no", "index", "#"],
    "parent_start": ["parent_start", "start", "from", "begin"],
    "parent_end": ["parent_end", "end", "to", "stop"],
    "parent_peptide": ["parent_peptide", "peptide", "sequence", "epitope", "seq"],
    "parent_length": ["parent_length", "length", "len"],
    "parent_score": ["parent_score", "score", "avg_score", "mean_score", "probability", "confidence"],
}


@dataclass
class FragmentFilterConfig:
    min_len: int = 5
    max_len: int = 15
    mode: str = "keep"  # keep | split
    step: int = 1
    best_per_parent: bool = False
    input_table: Optional[Path] = None
    output_name: Optional[str] = None


def _normalize_colname(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def _guess_sep(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    if path.suffix.lower() == ".csv":
        return ","
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except Exception:
        return "\t" if "\t" in sample else ","


def _resolve_columns(df: pd.DataFrame) -> pd.DataFrame:
    raw_to_norm = {c: _normalize_colname(c) for c in df.columns}
    norm_to_raw = {v: k for k, v in raw_to_norm.items()}

    resolved = {}
    for target, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            alias_norm = _normalize_colname(alias)
            if alias_norm in norm_to_raw:
                resolved[target] = norm_to_raw[alias_norm]
                break

    required = ["parent_start", "parent_end", "parent_peptide"]
    missing = [c for c in required if c not in resolved]
    if missing:
        raise ValueError(
            f"Could not resolve required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["parent_rank"] = df[resolved["parent_rank"]] if "parent_rank" in resolved else range(1, len(df) + 1)
    out["parent_start"] = df[resolved["parent_start"]]
    out["parent_end"] = df[resolved["parent_end"]]
    out["parent_peptide"] = df[resolved["parent_peptide"]]
    out["parent_length"] = df[resolved["parent_length"]] if "parent_length" in resolved else df[resolved["parent_peptide"]].astype(str).str.len()
    out["parent_score"] = df[resolved["parent_score"]] if "parent_score" in resolved else None

    out["parent_start"] = pd.to_numeric(out["parent_start"], errors="coerce")
    out["parent_end"] = pd.to_numeric(out["parent_end"], errors="coerce")
    out["parent_length"] = pd.to_numeric(out["parent_length"], errors="coerce")
    out["parent_peptide"] = out["parent_peptide"].astype(str).str.strip()

    out = out.dropna(subset=["parent_start", "parent_end", "parent_length"])
    out = out[out["parent_peptide"].str.len() > 0].copy()

    out["parent_start"] = out["parent_start"].astype(int)
    out["parent_end"] = out["parent_end"].astype(int)
    out["parent_length"] = out["parent_length"].astype(int)
    return out


def read_prediction_table(path: Path) -> pd.DataFrame:
    sep = _guess_sep(path)
    df = pd.read_csv(path, sep=sep)
    return _resolve_columns(df)


def find_candidate_table(out_dir: Path) -> Optional[Path]:
    preferred_names = [
        "predicted_peptides.tsv",
        "predicted_peptides.csv",
        "linear_epitopes.tsv",
        "linear_epitopes.csv",
        "epitopes.tsv",
        "epitopes.csv",
        "raw_output.tsv",
        "raw_output.csv",
    ]

    for name in preferred_names:
        p = out_dir / name
        if p.exists():
            return p

    all_tables = list(out_dir.glob("*.tsv")) + list(out_dir.glob("*.csv"))
    if not all_tables:
        return None

    # Prefer files with likely names
    keywords = ("peptide", "epitope", "output", "raw")
    scored = []
    for p in all_tables:
        score = sum(k in p.name.lower() for k in keywords)
        scored.append((score, p))
    scored.sort(key=lambda x: (-x[0], x[1].name))
    return scored[0][1]


def heuristic_rank_score(fragment: str, offset_in_parent: int = 0) -> float:
    seq = fragment.upper()
    length = len(seq)
    if length == 0:
        return 0.0

    aromatic = sum(aa in "FWY" for aa in seq)
    positive = sum(aa in "KRH" for aa in seq)
    gly_ser_penalty = sum(aa in "GSPNQ" for aa in seq)

    base = 1.0
    base += 0.10 * aromatic
    base += 0.04 * positive
    base -= 0.02 * gly_ser_penalty
    base -= 0.005 * offset_in_parent
    return round(base / length * 10, 4)


def keep_mode(df: pd.DataFrame, min_len: int, max_len: int) -> pd.DataFrame:
    out = df[(df["parent_length"] >= min_len) & (df["parent_length"] <= max_len)].copy()
    out["fragment_start"] = out["parent_start"]
    out["fragment_end"] = out["parent_end"]
    out["fragment"] = out["parent_peptide"]
    out["fragment_length"] = out["parent_length"]
    out["offset_in_parent"] = 0
    out["heuristic_rank_score"] = [
        heuristic_rank_score(seq, 0) for seq in out["fragment"].tolist()
    ]
    return out


def split_mode(df: pd.DataFrame, min_len: int, max_len: int, step: int) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        parent_rank = row["parent_rank"]
        parent_start = int(row["parent_start"])
        parent_end = int(row["parent_end"])
        parent_peptide = str(row["parent_peptide"])
        parent_length = int(row["parent_length"])
        parent_score = row["parent_score"]

        seq = parent_peptide
        seq_len = len(seq)

        if seq_len < min_len:
            continue

        upper = min(max_len, seq_len)

        for frag_len in range(upper, min_len - 1, -1):
            for offset in range(0, seq_len - frag_len + 1, step):
                frag = seq[offset: offset + frag_len]
                frag_start = parent_start + offset
                frag_end = frag_start + frag_len - 1

                rows.append({
                    "parent_rank": parent_rank,
                    "parent_start": parent_start,
                    "parent_end": parent_end,
                    "parent_peptide": parent_peptide,
                    "parent_length": parent_length,
                    "parent_score": parent_score,
                    "fragment_start": frag_start,
                    "fragment_end": frag_end,
                    "fragment": frag,
                    "fragment_length": frag_len,
                    "offset_in_parent": offset,
                    "heuristic_rank_score": heuristic_rank_score(frag, offset),
                })

    return pd.DataFrame(rows)


def best_per_parent(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    sort_cols = ["parent_rank", "heuristic_rank_score", "fragment_length", "offset_in_parent"]
    ascending = [True, False, False, True]
    ranked = df.sort_values(sort_cols, ascending=ascending).copy()
    best = ranked.groupby("parent_rank", as_index=False).head(1).copy()
    best = best.sort_values(["parent_rank", "fragment_start", "fragment_end"]).reset_index(drop=True)
    return best


def write_fragment_table(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if out_path.suffix.lower() == ".tsv" else ","
    df.to_csv(out_path, sep=sep, index=False)


def run_fragment_filter(
    out_dir: Path,
    min_len: int = 5,
    max_len: int = 15,
    mode: str = "keep",
    step: int = 1,
    best_only: bool = False,
    input_table: Optional[Path] = None,
    output_name: Optional[str] = None,
) -> Path:
    if min_len < 1:
        raise ValueError("min_len must be >= 1")
    if max_len < min_len:
        raise ValueError("max_len must be >= min_len")
    if step < 1:
        raise ValueError("step must be >= 1")
    if mode not in {"keep", "split"}:
        raise ValueError("mode must be 'keep' or 'split'")

    source_table = input_table if input_table else find_candidate_table(out_dir)
    if source_table is None or not Path(source_table).exists():
        raise FileNotFoundError(
            "No input peptide table found. Please provide --frag_input_table explicitly."
        )

    source_table = Path(source_table)
    parents = read_prediction_table(source_table)

    if mode == "keep":
        result = keep_mode(parents, min_len=min_len, max_len=max_len)
    else:
        result = split_mode(parents, min_len=min_len, max_len=max_len, step=step)

    if best_only:
        result = best_per_parent(result)

    if output_name is None:
        suffix = "tsv"
        output_name = f"fragment_candidates_{mode}_{min_len}_{max_len}.{suffix}"

    out_path = out_dir / output_name
    write_fragment_table(result, out_path)

    print(f"[fragment_filter] input table: {source_table}")
    print(f"[fragment_filter] input rows: {len(parents)}")
    print(f"[fragment_filter] output rows: {len(result)}")
    print(f"[fragment_filter] saved to: {out_path}")
    return out_path