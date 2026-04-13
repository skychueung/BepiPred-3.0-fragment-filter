from __future__ import annotations

from pathlib import Path
from typing import Optional
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