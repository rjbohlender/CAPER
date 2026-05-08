#!/usr/bin/env python3
"""Create regenie gene-based supplemental files from a CAPER matrix."""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


FIELD_NAMES = (
    "chrom",
    "start",
    "end",
    "ref",
    "alt",
    "type",
    "gene",
    "transcript",
    "region",
    "function",
    "annotation",
)

ANNOTATION_FIELDS = ("type", "region", "function", "annotation")


@dataclass
class MatrixVariant:
    chrom: str
    start: str
    end: str
    ref: str
    alt: str
    type: str
    gene: str
    transcript: str
    region: str
    function: str
    annotation: str


@dataclass
class SetRecord:
    chrom: str
    position: int
    variants: OrderedDict[str, None] = field(default_factory=OrderedDict)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the metadata columns from a CAPER matrix into regenie "
            "--anno-file, --set-list, and --mask-def inputs. The matrix is "
            "parsed strictly with CAPER's parser layout: chrom, start, end, "
            "ref, alt, type, gene, transcript, region, function, annotation, "
            "then sample genotype columns."
        )
    )
    parser.add_argument(
        "matrix",
        type=Path,
        help="CAPER matrix file, optionally gzip-compressed",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        required=True,
        help=(
            "Output prefix. Defaults are '<prefix>.anno.txt', "
            "'<prefix>.set-list.txt', and '<prefix>.mask-def.txt'."
        ),
    )
    parser.add_argument(
        "--anno-file",
        type=Path,
        help="Explicit destination for the regenie annotation file",
    )
    parser.add_argument(
        "--set-list",
        type=Path,
        help="Explicit destination for the regenie set-list file",
    )
    parser.add_argument(
        "--mask-def",
        type=Path,
        help="Explicit destination for the regenie mask-definition file",
    )
    parser.add_argument(
        "--annotation-column",
        choices=ANNOTATION_FIELDS,
        default="function",
        help=(
            "CAPER metadata column to use as regenie annotation category "
            "(default: function)"
        ),
    )
    parser.add_argument(
        "--missing-annotation",
        default="unannotated",
        help="Category to use when the selected annotation column is empty",
    )
    parser.add_argument(
        "--mask-column",
        choices=ANNOTATION_FIELDS,
        default="function",
        help=(
            "CAPER metadata column used to group annotation categories into "
            "regenie masks (default: function)"
        ),
    )
    parser.add_argument(
        "--missing-mask",
        default="unmasked",
        help="Mask name to use when the selected mask column is empty",
    )
    parser.add_argument(
        "--split-annotation-values",
        default="",
        metavar="CHARS",
        help=(
            "Optional characters that split one CAPER annotation field into "
            "multiple regenie categories. For example, ';' splits A;B."
        ),
    )
    parser.add_argument(
        "--mask-mode",
        choices=("all", "per-mask", "both"),
        default="per-mask",
        help=(
            "Which mask definitions to write. 'per-mask' writes masks grouped "
            "by --mask-column (default)."
        ),
    )
    parser.add_argument(
        "--all-mask-name",
        default="all",
        help="Mask name used when --mask-mode is all or both",
    )
    parser.add_argument(
        "--variant-id-format",
        choices=("chr-pos-ref-alt", "physical", "caper-fields"),
        default="chr-pos-ref-alt",
        help=(
            "How to derive variant IDs. Use the format that matches the .bim "
            "IDs used by regenie. Defaults to chrom_start_ref_alt."
        ),
    )
    parser.add_argument(
        "--variant-id-template",
        help=(
            "Custom Python format string for variant IDs, using parser-layout "
            "fields like {chrom}, {start}, {ref}, {alt}. Overrides "
            "--variant-id-format."
        ),
    )
    parser.add_argument(
        "--set-id-format",
        choices=("gene", "transcript", "gene-transcript"),
        default="gene",
        help="How to name regenie sets (default: gene)",
    )
    parser.add_argument(
        "--set-id-template",
        help=(
            "Custom Python format string for set IDs, using parser-layout "
            "fields like {gene} and {transcript}. Overrides --set-id-format."
        ),
    )
    return parser.parse_args(argv)


def open_text(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(2)

    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")

    return path.open("r", encoding="utf-8")


def output_paths(args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    prefix = args.out_prefix
    anno_file = args.anno_file or prefix.with_suffix(prefix.suffix + ".anno.txt")
    set_list = args.set_list or prefix.with_suffix(prefix.suffix + ".set-list.txt")
    mask_def = args.mask_def or prefix.with_suffix(prefix.suffix + ".mask-def.txt")
    return anno_file, set_list, mask_def


def parse_matrix_variant(raw_line: str, line_no: int) -> MatrixVariant:
    fields = raw_line.rstrip("\n").split("\t", 11)
    if len(fields) < len(FIELD_NAMES):
        raise ValueError(
            f"Line {line_no} has {len(fields)} tab-delimited fields; "
            f"expected at least {len(FIELD_NAMES)} CAPER metadata fields"
        )

    values = fields[: len(FIELD_NAMES)]
    return MatrixVariant(**dict(zip(FIELD_NAMES, values)))


def field_map(variant: MatrixVariant) -> Dict[str, str]:
    return {name: getattr(variant, name) for name in FIELD_NAMES}


def variant_id(variant: MatrixVariant, args: argparse.Namespace) -> str:
    values = field_map(variant)
    if args.variant_id_template:
        return args.variant_id_template.format(**values)

    if args.variant_id_format == "chr-pos-ref-alt":
        return f"{variant.chrom}_{variant.start}_{variant.ref}_{variant.alt}"
    if args.variant_id_format == "physical":
        return "_".join(
            (
                variant.chrom,
                variant.start,
                variant.end,
                variant.ref,
                variant.alt,
                variant.type,
            )
        )
    if args.variant_id_format == "caper-fields":
        return "_".join(
            (
                variant.chrom,
                variant.start,
                variant.end,
                variant.ref,
                variant.alt,
                variant.type,
                variant.gene,
                variant.transcript,
            )
        )

    raise ValueError(f"Unsupported variant ID format: {args.variant_id_format}")


def set_id(variant: MatrixVariant, args: argparse.Namespace) -> str:
    values = field_map(variant)
    if args.set_id_template:
        return args.set_id_template.format(**values)

    if args.set_id_format == "gene":
        return variant.gene
    if args.set_id_format == "transcript":
        return variant.transcript
    if args.set_id_format == "gene-transcript":
        return f"{variant.gene}_{variant.transcript}"

    raise ValueError(f"Unsupported set ID format: {args.set_id_format}")


def sanitize_category(value: str, missing: str) -> str:
    value = value.strip() or missing
    value = re.sub(r"\s+", "_", value)
    value = value.replace(",", "_")
    return value


def annotation_categories(
    variant: MatrixVariant, args: argparse.Namespace
) -> List[str]:
    raw_value = getattr(variant, args.annotation_column)
    if args.split_annotation_values and raw_value.strip():
        splitter = "[" + re.escape(args.split_annotation_values) + "]"
        raw_values = [
            value for value in re.split(splitter, raw_value) if value.strip()
        ]
    else:
        raw_values = [raw_value]

    categories: List[str] = []
    seen: set[str] = set()
    for value in raw_values:
        category = sanitize_category(value, args.missing_annotation)
        if category not in seen:
            categories.append(category)
            seen.add(category)
    return categories


def mask_name(variant: MatrixVariant, args: argparse.Namespace) -> str:
    return sanitize_category(getattr(variant, args.mask_column), args.missing_mask)


def parse_start(variant: MatrixVariant, line_no: int) -> int:
    try:
        return int(variant.start)
    except ValueError as exc:
        raise ValueError(
            f"Line {line_no} has a non-integer start position: {variant.start!r}"
        ) from exc


def ensure_no_whitespace(value: str, kind: str, line_no: int) -> None:
    if re.search(r"\s", value):
        raise ValueError(
            f"Line {line_no} produced {kind} with whitespace: {value!r}. "
            "Use a custom template or clean the matrix field before running regenie."
        )


def ensure_valid_variant_id(value: str, line_no: int) -> None:
    ensure_no_whitespace(value, "variant ID", line_no)
    if "," in value:
        raise ValueError(
            f"Line {line_no} produced a variant ID containing ',': {value!r}. "
            "regenie's set-list uses commas to separate variant IDs, so choose "
            "a comma-free --variant-id-template."
        )


def convert(args: argparse.Namespace) -> Tuple[int, int, int, int]:
    sets: OrderedDict[str, SetRecord] = OrderedDict()
    anno_rows: List[Tuple[str, str, str]] = []
    anno_seen: set[Tuple[str, str, str]] = set()
    variant_seen: set[str] = set()
    observed_categories: OrderedDict[str, None] = OrderedDict()
    masks: OrderedDict[str, OrderedDict[str, None]] = OrderedDict()

    n_rows = 0

    with open_text(args.matrix) as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Matrix file is empty: {args.matrix}")

        header_fields = header.rstrip("\n").split("\t", 11)
        if len(header_fields) < len(FIELD_NAMES) + 1:
            raise ValueError(
                "Matrix header does not look like CAPER parser layout plus "
                "sample columns"
            )

        for line_no, raw_line in enumerate(handle, start=2):
            if not raw_line.strip():
                continue

            variant = parse_matrix_variant(raw_line, line_no)
            n_rows += 1

            vid = variant_id(variant, args)
            sid = set_id(variant, args)
            ensure_valid_variant_id(vid, line_no)
            ensure_no_whitespace(sid, "set ID", line_no)
            if vid in variant_seen:
                continue
            variant_seen.add(vid)

            start = parse_start(variant, line_no)
            if sid not in sets:
                sets[sid] = SetRecord(chrom=variant.chrom, position=start)
            else:
                if sets[sid].chrom != variant.chrom:
                    raise ValueError(
                        f"Line {line_no} assigns set {sid!r} to chromosome "
                        f"{variant.chrom!r}, but it was already seen on "
                        f"chromosome {sets[sid].chrom!r}"
                    )
                sets[sid].position = min(sets[sid].position, start)

            sets[sid].variants.setdefault(vid, None)

            mask = mask_name(variant, args)
            ensure_no_whitespace(mask, "mask name", line_no)
            masks.setdefault(mask, OrderedDict())

            for category in annotation_categories(variant, args):
                ensure_no_whitespace(category, "annotation category", line_no)
                observed_categories.setdefault(category, None)
                masks[mask].setdefault(category, None)
                anno_row = (vid, sid, category)
                if anno_row not in anno_seen:
                    anno_rows.append(anno_row)
                    anno_seen.add(anno_row)

    if n_rows == 0:
        raise ValueError(f"Matrix file contains no variant rows: {args.matrix}")

    anno_file, set_list, mask_def = output_paths(args)
    for path in (anno_file, set_list, mask_def):
        path.parent.mkdir(parents=True, exist_ok=True)

    with anno_file.open("w", encoding="utf-8") as handle:
        for vid, sid, category in anno_rows:
            handle.write(f"{vid}\t{sid}\t{category}\n")

    with set_list.open("w", encoding="utf-8") as handle:
        for sid, record in sets.items():
            variants = ",".join(record.variants.keys())
            handle.write(f"{sid}\t{record.chrom}\t{record.position}\t{variants}\n")

    categories = list(observed_categories.keys())
    with mask_def.open("w", encoding="utf-8") as handle:
        if args.mask_mode in ("all", "both"):
            handle.write(f"{args.all_mask_name}\t{','.join(categories)}\n")
        if args.mask_mode in ("per-mask", "both"):
            for mask, mask_categories in masks.items():
                handle.write(f"{mask}\t{','.join(mask_categories.keys())}\n")

    n_set_variants = sum(len(record.variants) for record in sets.values())
    return n_rows, len(sets), n_set_variants, len(categories)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        n_rows, n_sets, n_set_variants, n_categories = convert(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    anno_file, set_list, mask_def = output_paths(args)
    print(
        "Wrote regenie supplemental files:\n"
        f"  anno-file: {anno_file}\n"
        f"  set-list:  {set_list}\n"
        f"  mask-def:  {mask_def}\n"
        f"Parsed {n_rows} matrix rows into {n_sets} sets, "
        f"{n_set_variants} set-variant entries, and {n_categories} categories.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
