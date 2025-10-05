#!/usr/bin/env python3
"""Create a Manhattan plot from a CAPER .simple file."""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_line,
    element_text,
    geom_hline,
    geom_point,
    ggplot,
    ggtitle,
    scale_color_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_bw,
)


SIMPLE_METRICS = {
    "empirical-p": ("Empirical_P", "Empirical P-value"),
    "empirical-midp": ("Empirical_MidP", "Empirical MidP"),
    "mgit-p": ("MGIT_P", "MGIT P-value"),
    "mgit-midp": ("MGIT_MIDP", "MGIT MidP"),
    "analytic-p": ("Analytic_P", "Analytic P-value"),
}


@dataclass(frozen=True)
class GeneCoordinate:
    chromosome: str
    position: int


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a publication-quality Manhattan plot from a CAPER .simple "
            "file. The script can infer genomic coordinates from a GFF3 "
            "annotation file or from columns already present in the simple file."
        )
    )
    parser.add_argument(
        "simple",
        type=Path,
        help="Path to the CAPER .simple file",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Destination for the generated plot (extension determines format)",
    )
    parser.add_argument(
        "--gff3",
        type=Path,
        help="Optional path to a GFF3 file for resolving gene coordinates",
    )
    parser.add_argument(
        "--title",
        help="Optional plot title",
    )
    parser.add_argument(
        "--metric",
        choices=sorted(SIMPLE_METRICS.keys()),
        default="empirical-p",
        help=(
            "Which statistic from the simple file to visualise. Choices mirror the "
            "columns written by CAPER (default: empirical-p)."
        ),
    )
    parser.add_argument(
        "--value-column",
        help=(
            "Explicit simple file column to plot. Overrides --metric. Useful for "
            "augmented tables that include additional association statistics."
        ),
    )
    parser.add_argument(
        "--gene-column",
        default="Gene",
        help="Column containing gene identifiers (default: Gene)",
    )
    parser.add_argument(
        "--transcript-column",
        default="Transcript",
        help="Column containing transcript identifiers (default: Transcript)",
    )
    parser.add_argument(
        "--chrom-column",
        help="Name of a column that already stores chromosome labels",
    )
    parser.add_argument(
        "--position-column",
        help="Name of a column that already stores genomic positions",
    )
    parser.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Keep multiple rows per gene instead of selecting the best p-value",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=6.5,
        help="Width of the output figure in inches (default: 6.5)",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=6.5,
        help="Height of the output figure in inches (default: 6.5)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution of raster outputs in DPI (default: 300)",
    )
    parser.add_argument(
        "--suggestive-line",
        type=float,
        metavar="P",
        help="Draw a suggestive significance line at the given p-value",
    )
    parser.add_argument(
        "--significance-line",
        type=float,
        metavar="P",
        help="Draw a genome-wide significance line at the given p-value",
    )

    args = parser.parse_args(argv)

    if args.chrom_column and not args.position_column:
        parser.error("--chrom-column requires --position-column to be specified as well")
    if args.position_column and not args.chrom_column:
        parser.error("--position-column requires --chrom-column to be specified as well")

    return args


def read_simple(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Simple file '{path}' does not exist")

    try:
        df = pd.read_csv(
            path,
            comment="#",
            delim_whitespace=True,
            dtype=str,
        )
    except Exception as exc:  # pragma: no cover - defensive programming
        raise RuntimeError(f"Unable to read simple file '{path}': {exc}")

    if df.empty:
        raise ValueError("The provided simple file does not contain any records")

    return df


def parse_attributes(field: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for item in field.split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key.strip()] = value.strip()
    return attributes


def load_gene_coordinates(path: Path) -> Tuple[Dict[str, GeneCoordinate], Dict[str, str]]:
    gene_coords: Dict[str, GeneCoordinate] = {}
    transcript_to_gene: Dict[str, str] = {}

    if path is None:
        return gene_coords, transcript_to_gene

    if not path.exists():
        raise FileNotFoundError(f"GFF3 file '{path}' does not exist")

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 9:
                continue
            chrom, _source, feature_type, start, _end, _score, _strand, _phase, attributes = parts
            attrs = parse_attributes(attributes)

            if feature_type == "gene":
                identifiers = {
                    attrs.get("ID"),
                    attrs.get("Name"),
                    attrs.get("gene_id"),
                    attrs.get("gene_name"),
                }
                identifiers = {value for value in identifiers if value}
                if not identifiers:
                    continue
                try:
                    position = int(start)
                except ValueError:
                    continue
                coordinate = GeneCoordinate(chromosome=chrom, position=position)
                for identifier in identifiers:
                    gene_coords.setdefault(identifier, coordinate)
                    # Also store versionless identifiers for convenience.
                    if "." in identifier:
                        gene_coords.setdefault(identifier.split(".")[0], coordinate)

            elif feature_type in {"mRNA", "transcript"}:
                transcript_id = attrs.get("ID")
                parents = attrs.get("Parent", "")
                if not transcript_id or not parents:
                    continue
                parent_gene = parents.split(",")[0]
                transcript_to_gene.setdefault(transcript_id, parent_gene)
                if "." in transcript_id:
                    transcript_to_gene.setdefault(transcript_id.split(".")[0], parent_gene)

    return gene_coords, transcript_to_gene


def resolve_coordinate(
    row: pd.Series,
    identifier_columns: Sequence[str],
    gene_coords: Dict[str, GeneCoordinate],
    transcript_to_gene: Dict[str, str],
) -> Optional[GeneCoordinate]:
    identifiers: List[str] = []
    for column in identifier_columns:
        value = row.get(column)
        if isinstance(value, str) and value:
            identifiers.append(value)
            if "." in value:
                identifiers.append(value.split(".")[0])

    for identifier in identifiers:
        if identifier in gene_coords:
            return gene_coords[identifier]
        if identifier in transcript_to_gene:
            parent = transcript_to_gene[identifier]
            if parent in gene_coords:
                return gene_coords[parent]
            if "." in parent and parent.split(".")[0] in gene_coords:
                return gene_coords[parent.split(".")[0]]

    return None


def coerce_numeric(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return values


def clean_chromosome(value: str) -> str:
    if not isinstance(value, str):
        return str(value)
    val = value.strip()
    if val.lower().startswith("chr"):
        val = val[3:]
    if val in {"MT", "Mt", "mt"}:
        return "MT"
    return val


def chromosome_sort_key(chrom: str) -> Tuple[int, str]:
    clean = clean_chromosome(chrom)
    upper = clean.upper()
    special = {"X": 1000, "Y": 1001, "MT": 1002, "M": 1002}
    if upper in special:
        return (special[upper], "")
    try:
        return (int(clean), "")
    except ValueError:
        return (2000, upper)


def prepare_dataframe(
    df: pd.DataFrame,
    value_column: str,
    identifier_columns: Sequence[str],
    chrom_column: Optional[str],
    position_column: Optional[str],
    keep_duplicates: bool,
    gene_coords: Dict[str, GeneCoordinate],
    transcript_to_gene: Dict[str, str],
) -> Tuple[pd.DataFrame, List[float], List[str]]:
    if value_column not in df.columns:
        raise ValueError(f"Column '{value_column}' was not found in the simple file")

    df = df.copy()

    df[value_column] = coerce_numeric(df[value_column])
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=[value_column])

    if not keep_duplicates:
        for column in identifier_columns:
            if column in df.columns:
                df = df.sort_values(value_column, ascending=True, kind="mergesort")
                df = df.drop_duplicates(subset=column, keep="first")
                break

    if chrom_column and position_column:
        if chrom_column not in df.columns or position_column not in df.columns:
            raise ValueError(
                "Both --chrom-column and --position-column must exist in the simple file"
            )
        df["chromosome"] = df[chrom_column].astype(str)
        df["position"] = coerce_numeric(df[position_column])
    else:
        if not identifier_columns:
            raise ValueError(
                "Unable to resolve genomic coordinates. Provide --gene-column/"
                "--transcript-column names that exist in the table, or supply "
                "--chrom-column/--position-column."
            )
        coordinates = df.apply(
            resolve_coordinate,
            axis=1,
            identifier_columns=identifier_columns,
            gene_coords=gene_coords,
            transcript_to_gene=transcript_to_gene,
        )
        df["chromosome"] = [coord.chromosome if coord else np.nan for coord in coordinates]
        df["position"] = [coord.position if coord else np.nan for coord in coordinates]

    df = df.dropna(subset=["chromosome", "position"])
    df["position"] = coerce_numeric(df["position"])
    df = df.dropna(subset=["position"])

    if df.empty:
        raise ValueError(
            "No rows with genomic coordinates were available. "
            "Provide --gff3 or specify --chrom-column/--position-column."
        )

    df["chromosome"] = df["chromosome"].map(clean_chromosome)

    min_positive = df.loc[df[value_column] > 0, value_column].min()
    if pd.isna(min_positive):
        min_positive = 1e-300
    replacement = min_positive / 10 if min_positive > 0 else 1e-300
    df.loc[df[value_column] <= 0, value_column] = replacement

    df["neg_log_value"] = -np.log10(df[value_column])

    chrom_order = sorted(df["chromosome"].unique(), key=chromosome_sort_key)
    df["chromosome"] = pd.Categorical(df["chromosome"], categories=chrom_order, ordered=True)
    df = df.sort_values(["chromosome", "position"], ascending=[True, True])

    cumulative_positions: List[float] = []
    tick_positions: List[float] = []
    tick_labels: List[str] = []
    current_offset = 0.0

    for index, chrom in enumerate(chrom_order):
        chrom_mask = df["chromosome"] == chrom
        chrom_positions = df.loc[chrom_mask, "position"]
        if chrom_positions.empty:
            continue
        start = float(chrom_positions.min())
        offset = current_offset - start
        adjusted = chrom_positions.astype(float) + offset
        df.loc[chrom_mask, "cumulative_position"] = adjusted

        min_pos = float(adjusted.min())
        max_pos = float(adjusted.max())
        tick_positions.append(min_pos + (max_pos - min_pos) / 2.0)
        tick_labels.append(str(chrom))

        current_offset = max_pos + 1_000_000.0

    df["color_group"] = df["chromosome"].cat.codes % 2

    return df, tick_positions, tick_labels


def build_plot(
    df: pd.DataFrame,
    tick_positions: Sequence[float],
    tick_labels: Sequence[str],
    title: Optional[str],
    suggestive_line: Optional[float],
    significance_line: Optional[float],
    y_label: str,
) -> ggplot:
    palette = ["#1f78b4", "#33a02c"]

    plot = (
        ggplot(
            df,
            aes(x="cumulative_position", y="neg_log_value", color="factor(color_group)"),
        )
        + geom_point(size=1.1, alpha=0.8)
        + scale_color_manual(values=palette, guide=False)
        + scale_x_continuous(name="Chromosome", breaks=tick_positions, labels=tick_labels)
        + scale_y_continuous(name=y_label, expand=(0.02, 0))
        + theme_bw()
        + theme(
            panel_grid_major_x=element_blank(),
            panel_grid_minor_x=element_blank(),
            panel_grid_minor_y=element_line(color="#e0e0e0", size=0.3),
            panel_grid_major_y=element_line(color="#cccccc", size=0.4),
            axis_text_x=element_text(rotation=45, ha="right"),
            axis_title_x=element_text(margin={"t": 10}),
            axis_title_y=element_text(margin={"r": 10}),
            figure_size=(8, 4.5),
        )
    )

    if title:
        plot += ggtitle(title)

    for line_value, style in (
        (suggestive_line, {"linetype": "dashed", "color": "#6a3d9a"}),
        (significance_line, {"linetype": "solid", "color": "#e31a1c"}),
    ):
        if line_value is None:
            continue
        if line_value <= 0 or math.isnan(line_value):
            continue
        neg_log = -math.log10(line_value)
        plot += geom_hline(yintercept=neg_log, **style)

    return plot


def save_plot(plot: ggplot, output: Path, width: float, height: float, dpi: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        plot.save(filename=str(output), width=width, height=height, dpi=dpi)
    except Exception as exc:  # pragma: no cover - defensive programming
        raise RuntimeError(f"Failed to save plot to '{output}': {exc}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    df = read_simple(args.simple)
    gene_coords, transcript_to_gene = load_gene_coordinates(args.gff3)

    identifier_columns: List[str] = []
    for column in (args.gene_column, args.transcript_column):
        if column:
            identifier_columns.append(column)

    if args.value_column:
        value_column = args.value_column
        y_label = f"-log10({args.value_column})"
    else:
        value_column, friendly_name = SIMPLE_METRICS[args.metric]
        y_label = f"-log10({friendly_name})"

    processed, tick_positions, tick_labels = prepare_dataframe(
        df,
        value_column=value_column,
        identifier_columns=identifier_columns,
        chrom_column=args.chrom_column,
        position_column=args.position_column,
        keep_duplicates=args.keep_duplicates,
        gene_coords=gene_coords,
        transcript_to_gene=transcript_to_gene,
    )

    plot = build_plot(
        processed,
        tick_positions=tick_positions,
        tick_labels=tick_labels,
        title=args.title,
        suggestive_line=args.suggestive_line,
        significance_line=args.significance_line,
        y_label=y_label,
    )

    save_plot(plot, args.output, width=args.width, height=args.height, dpi=args.dpi)

    return 0


if __name__ == "__main__":  # pragma: no cover - entry point
    try:
        sys.exit(main())
    except Exception as exc:  # pragma: no cover - provide friendly error messages
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
