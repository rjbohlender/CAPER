#!/usr/bin/env python3
"""Merge CAPER .simple files from sharded runs and re-rank the results."""
from __future__ import annotations

import argparse
import io
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence


@dataclass
class SimpleFile:
    comments: List[str]
    header: List[str]
    rows: List[List[str]]
    base_length: int


def parse_simple(path: Path) -> SimpleFile:
    """Parse a .simple file."""
    comments: List[str] = []
    header: Optional[List[str]] = None
    rows: List[List[str]] = []
    base_length: Optional[int] = None

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                continue
            if header is None:
                header = line.split()
                base_length = len(header)
                continue
            fields = line.split()
            if base_length is None:
                raise ValueError("Internal error: base_length should be initialized")

            expected_columns = len(header)
            if len(fields) != expected_columns:
                raise ValueError(
                    "Encountered a ragged row while parsing '{}'\n"
                    "Expected {} columns based on the header, but found {}.\n"
                    "Offending line: {}".format(
                        path,
                        expected_columns,
                        len(fields),
                        raw_line.rstrip(),
                    )
                )

            rows.append(fields)

    if header is None or base_length is None:
        raise ValueError(f"File '{path}' does not contain a header row")

    return SimpleFile(comments=comments, header=header, rows=rows, base_length=base_length)


def determine_sort_column(columns: Sequence[str], explicit: Optional[str] = None) -> str:
    if explicit:
        if explicit not in columns:
            raise ValueError(f"Column '{explicit}' is not present in the merged table")
        return explicit

    for candidate in ("Empirical_P", "Analytic_P", "Test_Statistic"):
        if candidate in columns:
            return candidate

    raise ValueError(
        "Unable to determine a sort column automatically. "
        "Use --sort-column to choose one."
    )


def parse_numeric(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float("inf")


def merge_simples(
    files: List[SimpleFile], sort_column: str, columns: Sequence[str]
) -> SimpleFile:
    if not files:
        raise ValueError("No input files were provided")

    header = list(columns)
    if not header or header[0] != "Rank":
        raise ValueError("The merged table must begin with a 'Rank' column")

    header_index = {name: idx for idx, name in enumerate(header)}
    sort_index = header_index[sort_column]

    merged_rows: List[List[str]] = []
    for simple in files:
        mapping: List[int] = []
        for name in simple.header:
            if name not in header_index:
                raise ValueError(
                    f"Column '{name}' from one of the inputs is missing in the merged header"
                )
            mapping.append(header_index[name])

        for row in simple.rows:
            if len(row) != len(mapping):
                raise ValueError("Encountered a row with an unexpected number of columns")

            aligned = ["" for _ in header]
            for src_idx, dest_idx in enumerate(mapping):
                aligned[dest_idx] = row[src_idx]
            merged_rows.append(aligned)

    merged_rows.sort(key=lambda row: parse_numeric(row[sort_index]))

    for idx, row in enumerate(merged_rows, start=1):
        row[0] = str(idx)

    return SimpleFile(
        comments=[],
        header=header,
        rows=merged_rows,
        base_length=files[0].base_length,
    )


def write_simple(simple: SimpleFile, destination: Optional[Path]) -> None:
    output: io.TextIOBase
    close_after = False

    if destination is None:
        output = sys.stdout
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        output = destination.open("w", encoding="utf-8", newline="\n")
        close_after = True

    try:
        for comment in simple.comments:
            output.write(comment + "\n")
        output.write("\t".join(simple.header) + "\n")
        for row in simple.rows:
            output.write("\t".join(row) + "\n")
    finally:
        if close_after:
            output.close()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Merge CAPER .simple files from sharded runs, re-rank the results, and "
            "output a single aggregated .simple file."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Paths to the input .simple files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="File to write merged results to (defaults to stdout)",
    )
    parser.add_argument(
        "--sort-column",
        help="Column to use when ordering the merged results",
    )

    args = parser.parse_args(argv)

    simples: List[SimpleFile] = []
    for path in args.inputs:
        if not path.exists():
            parser.error(f"Input file '{path}' does not exist")
        simples.append(parse_simple(path))

    column_order: List[str] = []
    for simple in simples:
        for column in simple.header:
            if column not in column_order:
                column_order.append(column)

    sort_column = determine_sort_column(column_order, args.sort_column)
    merged = merge_simples(simples, sort_column, column_order)
    write_simple(merged, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
