from __future__ import annotations

from pathlib import Path
from typing import Optional, Dict
import pandas as pd

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as e:
    raise ImportError("physicochemical_filter.py requires biopython. Please install biopython first.") from e


def calculate_net_charge_at_ph(sequence: str, ph: float = 7.4) -> float:
    sequence = sequence.upper()

    pka_asp = 3.65
    pka_glu = 4.25
    pka_his = 6.00
    pka_cys = 8.18
    pka_tyr = 10.07
    pka_lys = 10.53
    pka_arg = 12.48
    pka_n_term = 8.0
    pka_c_term = 3.1

    charge = 0.0

    charge += 1.0 / (1.0 + 10 ** (ph - pka_n_term))
    charge -= 1.0 / (1.0 + 10 ** (pka_c_term - ph))

    for aa, pka in [("D", pka_asp), ("E", pka_glu), ("C", pka_cys), ("Y", pka_tyr)]:
        count = sequence.count(aa)
        charge -= count * (1.0 / (1.0 + 10 ** (pka - ph)))

    for aa, pka in [("H", pka_his), ("K", pka_lys), ("R", pka_arg)]:
        count = sequence.count(aa)
        charge += count * (1.0 / (1.0 + 10 ** (ph - pka)))

    return round(charge, 4)


def predict_disulfide_risk(sequence: str) -> str:
    seq = sequence.upper()
    cys_positions = [i for i, aa in enumerate(seq) if aa == "C"]
    cys_count = len(cys_positions)

    if cys_count == 0:
        return "none"
    if cys_count == 1:
        return "low"
    if cys_count >= 3:
        return "high"

    gap = abs(cys_positions[1] - cys_positions[0]) - 1
    if 2 <= gap <= 8 and len(seq) <= 15:
        return "higher"
    return "medium"


def calculate_properties(sequence: str) -> dict:
    analyzer = ProteinAnalysis(sequence)
    return {
        "Net_Charge_pH7.4": calculate_net_charge_at_ph(sequence, ph=7.4),
        "GRAVY": round(analyzer.gravy(), 4),
        "pI": round(analyzer.isoelectric_point(), 4),
        "Cys_Count": sequence.upper().count("C"),
        "Disulfide_Risk": predict_disulfide_risk(sequence),
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _enrich_with_properties(df: pd.DataFrame, sequence_col: str) -> pd.DataFrame:
    """Add physicochemical property columns to every row."""
    rows = []
    for _, row in df.iterrows():
        seq = str(row[sequence_col]).strip()
        if not seq:
            continue
        props = calculate_properties(seq)
        out = dict(row)
        out.update(props)
        rows.append(out)
    return pd.DataFrame(rows)


def _save_layer(df: pd.DataFrame, output_table: Path, label: str) -> None:
    """Save a single layer result and print stats."""
    output_table.parent.mkdir(parents=True, exist_ok=True)
    out_sep = "\t" if output_table.suffix.lower() == ".tsv" else ","
    df.to_csv(output_table, sep=out_sep, index=False)
    print(f"[physicochemical_filter][{label}] output rows: {len(df)}")
    print(f"[physicochemical_filter][{label}] saved to: {output_table}")


# ---------------------------------------------------------------------------
# Single-file mode (original)
# ---------------------------------------------------------------------------

def run_physicochemical_filter(
    input_table: Path,
    output_table: Path,
    sequence_col: str = "fragment",
    length_col: str = "fragment_length",
    start_col: str = "fragment_start",
    end_col: str = "fragment_end",
    min_charge: float = 0.0,
    max_gravy: float = 0.0,
    min_pi: Optional[float] = None,
    max_pi: Optional[float] = None,
    max_cys: Optional[int] = None,
    exclude_high_disulfide_risk: bool = False,
) -> Path:
    input_table = Path(input_table)
    if not input_table.exists():
        raise FileNotFoundError(f"Input table not found: {input_table}")

    sep = "\t" if input_table.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(input_table, sep=sep)

    required = [sequence_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Available columns: {list(df.columns)}")

    rows = []
    for _, row in df.iterrows():
        seq = str(row[sequence_col]).strip()
        if not seq:
            continue
        props = calculate_properties(seq)
        out = dict(row)
        out.update(props)
        rows.append(out)

    prop_df = pd.DataFrame(rows)
    if prop_df.empty:
        output_table.parent.mkdir(parents=True, exist_ok=True)
        prop_df.to_csv(output_table, sep="\t" if output_table.suffix.lower() == ".tsv" else ",", index=False)
        return output_table

    mask = (
        (prop_df["Net_Charge_pH7.4"] > min_charge) &
        (prop_df["GRAVY"] < max_gravy)
    )
    if min_pi is not None:
        mask &= prop_df["pI"] >= min_pi
    if max_pi is not None:
        mask &= prop_df["pI"] <= max_pi
    if max_cys is not None:
        mask &= prop_df["Cys_Count"] <= max_cys
    if exclude_high_disulfide_risk:
        mask &= ~prop_df["Disulfide_Risk"].isin(["higher", "high"])

    filtered_df = prop_df[mask].copy()

    sort_cols = ["Net_Charge_pH7.4", "GRAVY", "pI"]
    ascending = [False, True, False]
    existing_sort_cols = [c for c in sort_cols if c in filtered_df.columns]
    if existing_sort_cols:
        filtered_df = filtered_df.sort_values(
            by=existing_sort_cols,
            ascending=ascending[:len(existing_sort_cols)]
        ).reset_index(drop=True)

    output_table.parent.mkdir(parents=True, exist_ok=True)
    out_sep = "\t" if output_table.suffix.lower() == ".tsv" else ","
    filtered_df.to_csv(output_table, sep=out_sep, index=False)

    print(f"[physicochemical_filter] input table: {input_table}")
    print(f"[physicochemical_filter] input rows: {len(prop_df)}")
    print(f"[physicochemical_filter] output rows: {len(filtered_df)}")
    print(f"[physicochemical_filter] saved to: {output_table}")
    return output_table


# ---------------------------------------------------------------------------
# Layered mode (new)
# ---------------------------------------------------------------------------

def run_physicochemical_filter_layered(
    input_table: Path,
    output_dir: Path,
    stem: str,
    sequence_col: str = "fragment",
    min_charge: float = 0.0,
    max_gravy: float = 0.0,
    min_pi: Optional[float] = None,
    max_pi: Optional[float] = None,
    max_cys: Optional[int] = None,
    exclude_high_disulfide_risk: bool = False,
) -> Dict[str, Path]:
    """Run physicochemical filtering with independent per-layer outputs.

    Each layer applies ONLY its own criterion on the full property-enriched
    input, so the user can compare the effect of each filter independently.
    A final layer applies all criteria combined.

    Output files:
        {stem}_charge.tsv   - charge filter only
        {stem}_pi.tsv       - pI filter only
        {stem}_gravy.tsv    - GRAVY filter only
        {stem}_cys.tsv      - Cys count / disulfide risk filter only
        {stem}_final.tsv    - all criteria combined

    Returns a dict mapping layer name -> output Path.
    """
    input_table = Path(input_table)
    if not input_table.exists():
        raise FileNotFoundError(f"Input table not found: {input_table}")

    sep = "\t" if input_table.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(input_table, sep=sep)

    required = [sequence_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Available: {list(df.columns)}")

    print(f"[physicochemical_filter][layered] input table: {input_table}")
    print(f"[physicochemical_filter][layered] input rows: {len(df)}")

    # Enrich every row with all properties (no filtering yet)
    prop_df = _enrich_with_properties(df, sequence_col)
    print(f"[physicochemical_filter][layered] enriched rows: {len(prop_df)}")

    output_dir = Path(output_dir)
    results = {}

    # --- Layer 1: charge only ---
    mask_charge = prop_df["Net_Charge_pH7.4"] > min_charge
    df_charge = prop_df[mask_charge].copy()
    path_charge = output_dir / f"{stem}_charge.tsv"
    _save_layer(df_charge, path_charge, "charge")
    results["charge"] = path_charge

    # --- Layer 2: pI only ---
    mask_pi = pd.Series(True, index=prop_df.index)
    if min_pi is not None:
        mask_pi &= prop_df["pI"] >= min_pi
    if max_pi is not None:
        mask_pi &= prop_df["pI"] <= max_pi
    df_pi = prop_df[mask_pi].copy()
    path_pi = output_dir / f"{stem}_pi.tsv"
    _save_layer(df_pi, path_pi, "pi")
    results["pi"] = path_pi

    # --- Layer 3: GRAVY only ---
    mask_gravy = prop_df["GRAVY"] < max_gravy
    df_gravy = prop_df[mask_gravy].copy()
    path_gravy = output_dir / f"{stem}_gravy.tsv"
    _save_layer(df_gravy, path_gravy, "gravy")
    results["gravy"] = path_gravy

    # --- Layer 4: Cys / disulfide risk only ---
    mask_cys = pd.Series(True, index=prop_df.index)
    if max_cys is not None:
        mask_cys &= prop_df["Cys_Count"] <= max_cys
    if exclude_high_disulfide_risk:
        mask_cys &= ~prop_df["Disulfide_Risk"].isin(["higher", "high"])
    df_cys = prop_df[mask_cys].copy()
    path_cys = output_dir / f"{stem}_cys.tsv"
    _save_layer(df_cys, path_cys, "cys")
    results["cys"] = path_cys

    # --- Layer 5: final (all criteria combined) ---
    mask_final = (
        (prop_df["Net_Charge_pH7.4"] > min_charge) &
        (prop_df["GRAVY"] < max_gravy)
    )
    if min_pi is not None:
        mask_final &= prop_df["pI"] >= min_pi
    if max_pi is not None:
        mask_final &= prop_df["pI"] <= max_pi
    if max_cys is not None:
        mask_final &= prop_df["Cys_Count"] <= max_cys
    if exclude_high_disulfide_risk:
        mask_final &= ~prop_df["Disulfide_Risk"].isin(["higher", "high"])

    df_final = prop_df[mask_final].copy()
    sort_cols = ["Net_Charge_pH7.4", "GRAVY", "pI"]
    ascending = [False, True, False]
    existing = [c for c in sort_cols if c in df_final.columns]
    if existing:
        df_final = df_final.sort_values(
            by=existing, ascending=ascending[:len(existing)]
        ).reset_index(drop=True)

    path_final = output_dir / f"{stem}_final.tsv"
    _save_layer(df_final, path_final, "final")
    results["final"] = path_final

    return results
