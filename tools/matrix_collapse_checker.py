#!/usr/bin/env python3
"""Check CAPER matrix variant-collapsing counts.

This mirrors the current C++ Gene path for no-call handling, minor-allele
conversion, variant collapsing, and post-collapse MAF/AAF/MAC filtering closely
enough to diagnose whether zero-count/no-call variants change a collapsed
pseudo-variant.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


CHROM = 0
START = 1
END = 2
REF = 3
ALT = 4
TYPE = 5
GENE = 6
TRANSCRIPT = 7
FUNCTION = 9
FIRST_SAMPLE = 11


@dataclass
class Variant:
    line: str
    fields: list[str]
    genotypes: list[float]
    missing: list[bool]
    removed_by_function_filter: bool = False

    @property
    def variant_id(self) -> str:
        return ",".join(
            [
                self.fields[CHROM],
                self.fields[START],
                self.fields[END],
                self.fields[REF],
                self.fields[ALT],
                self.fields[TYPE],
                self.fields[GENE],
                self.fields[TRANSCRIPT],
            ]
        )

    @property
    def kind(self) -> str:
        return self.fields[TYPE]


@dataclass
class Transcript:
    gene: str
    transcript: str
    variants: list[Variant] = field(default_factory=list)


@dataclass
class CollapseReport:
    gene: str
    transcript: str
    input_variants: int
    collapse_candidate_variants: int
    no_call_cells_before: int
    no_call_samples_before: int
    non_no_call_carriers_before: int
    carriers_after_collapse: int
    collapsed_removed_by_mac: bool
    carriers_after_filtering: int
    remaining_variants_after_filtering: int
    case_ref: int | None = None
    case_alt: int | None = None
    cont_ref: int | None = None
    cont_alt: int | None = None
    collapse_variant_patterns: list[str] = field(default_factory=list)
    collapse_variant_lines: list[str] = field(default_factory=list)


def open_text(path: str):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def open_output(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "wt", newline="")


def read_sample_filter(path: str | None) -> set[str] | None:
    if path is None:
        return None
    samples: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                samples.add(line.split()[0])
    return samples


def read_ped(path: str | None, quantitative: bool = False) -> dict[str, float] | None:
    if path is None:
        return None
    phenotypes: dict[str, float] = {}
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                raise SystemExit(f"{path}:{line_no}: expected at least 6 PED fields")
            sample = fields[1]
            try:
                phenotype = float(fields[5]) if quantitative else int(fields[5]) - 1
            except ValueError as exc:
                raise SystemExit(
                    f"{path}:{line_no}: could not parse phenotype {fields[5]!r}"
                ) from exc
            phenotypes[sample] = phenotype
    return phenotypes


def read_filter_whitelist(path: str | None, method: str) -> dict[str, set[str]]:
    allowed = {"TYPE": set(), "FUNCTION": set()}
    if path is None:
        return allowed

    with open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        try:
            method_col = header.index(method)
        except ValueError as exc:
            raise SystemExit(f"method {method!r} not found in {path}") from exc

        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            key = row[0].split(":", 1)
            if len(key) != 2:
                continue
            category, value = key
            if method_col < len(row) and row[method_col] == "1":
                allowed.setdefault(category, set()).add(value)
    return allowed


def allow_variant(variant: Variant, allowed: dict[str, set[str]]) -> bool:
    if not allowed["TYPE"] and not allowed["FUNCTION"]:
        return True
    return (
        variant.fields[TYPE] in allowed["TYPE"]
        and variant.fields[FUNCTION] in allowed["FUNCTION"]
    )


def parse_genotypes(
    fields: list[str], sample_columns: list[int], header_sample_count: int
) -> tuple[list[float], list[bool]]:
    if len(fields) == FIRST_SAMPLE + header_sample_count:
        raw = [fields[col] for col in sample_columns]
    else:
        # CAPER's C++ parser keeps all sample genotype characters in the final
        # split field and indexes it by header column offset.
        genotype_blob = "".join(fields[FIRST_SAMPLE:])
        raw = [genotype_blob[col - FIRST_SAMPLE] for col in sample_columns]

    genotypes: list[float] = []
    missing: list[bool] = []
    for token in raw:
        if not token:
            value = 0
        else:
            value = ord(token[0]) - ord("0")
        if value > 2:
            genotypes.append(0.0)
            missing.append(True)
        else:
            genotypes.append(float(value))
            missing.append(False)
    return genotypes, missing


def parse_matrix(
    path: str,
    sample_filter: set[str] | None,
    allowed: dict[str, set[str]],
    whole_gene: bool,
) -> tuple[str, list[str], dict[tuple[str, str], Transcript]]:
    transcripts: dict[tuple[str, str], Transcript] = {}
    samples: list[str] = []
    sample_columns: list[int] = []
    header_sample_count = 0
    header_line = ""

    with open_text(path) as handle:
        for line_no, line in enumerate(handle):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if line_no == 0:
                header_line = line
                header_samples = fields[FIRST_SAMPLE:]
                header_sample_count = len(header_samples)
                for offset, sample in enumerate(header_samples, start=FIRST_SAMPLE):
                    if sample_filter is None or sample in sample_filter:
                        samples.append(sample)
                        sample_columns.append(offset)
                continue

            if len(fields) < FIRST_SAMPLE + 1:
                raise SystemExit(f"line {line_no + 1}: expected matrix fields")

            gene = fields[GENE]
            transcript = gene if whole_gene else fields[TRANSCRIPT]
            key = (gene, transcript)
            if key not in transcripts:
                transcripts[key] = Transcript(gene=gene, transcript=transcript)

            genotypes, missing = parse_genotypes(fields, sample_columns, header_sample_count)
            variant = Variant(line=line, fields=fields, genotypes=genotypes, missing=missing)
            variant.removed_by_function_filter = not allow_variant(variant, allowed)
            transcripts[key].variants.append(variant)

    return header_line, samples, transcripts


def convert_to_minor_allele_counts(
    genotypes: list[list[float]], missing: list[list[bool]], types: list[str]
) -> list[list[float]]:
    converted = [row[:] for row in genotypes]
    if not converted:
        return converted

    nsamples = len(converted[0])
    for idx, row in enumerate(converted):
        if idx < len(types) and types[idx] == "Collapsed":
            continue
        maf = sum(row) / (2.0 * nsamples)
        if maf > 0.5:
            flipped = []
            for sample_idx, value in enumerate(row):
                if missing[idx][sample_idx]:
                    flipped.append(0.0)
                elif value == 0:
                    flipped.append(2.0)
                elif value == 2:
                    flipped.append(0.0)
                else:
                    flipped.append(value)
            converted[idx] = flipped
    return converted


def filter_after_collapse(
    rows: list[list[float]],
    types: list[str],
    function_removed: list[bool],
    maf_threshold: float,
    mac_threshold: float,
    aaf_filter: bool,
) -> tuple[int, bool]:
    kept = 0
    collapsed_removed_by_mac = False
    nsamples = len(rows[0]) if rows else 0

    filtered_rows = rows
    if not aaf_filter:
        filtered_rows = convert_to_minor_allele_counts(
            rows, [[False] * nsamples for _ in rows], types
        )

    for idx, row in enumerate(filtered_rows):
        is_collapsed = types[idx] == "Collapsed"
        allele_sum = sum(row)
        frequency = allele_sum / (2.0 * nsamples) if nsamples else 0.0
        remove_by_mac = allele_sum > mac_threshold
        if aaf_filter:
            remove_by_frequency = (not is_collapsed) and frequency >= maf_threshold
        else:
            remove_by_frequency = (not is_collapsed) and frequency > maf_threshold
        remove = remove_by_mac or remove_by_frequency or function_removed[idx]
        if remove and is_collapsed and remove_by_mac:
            collapsed_removed_by_mac = True
        if not remove:
            kept += 1
    return kept, collapsed_removed_by_mac


def analyze_transcript(
    transcript: Transcript,
    nsamples: int,
    phenotypes: list[float] | None,
    collapse_threshold: int,
    maf_threshold: float,
    mac_threshold: float,
    aaf_filter: bool,
) -> CollapseReport:
    genotypes = [variant.genotypes for variant in transcript.variants]
    missing = [variant.missing for variant in transcript.variants]
    types = [variant.kind for variant in transcript.variants]
    function_removed = [variant.removed_by_function_filter for variant in transcript.variants]

    collapse_genotypes = convert_to_minor_allele_counts(genotypes, missing, types)
    allele_sums = [sum(row) for row in collapse_genotypes]
    rare_indices = [
        idx
        for idx, allele_sum in enumerate(allele_sums)
        if 0 < allele_sum <= collapse_threshold
    ]

    no_call_cells = sum(
        1 for idx in rare_indices for value in missing[idx] if value
    )
    no_call_samples = sum(
        1 for sample_idx in range(nsamples) if any(missing[idx][sample_idx] for idx in rare_indices)
    )
    carriers_before = sum(
        1
        for sample_idx in range(nsamples)
        if any(
            collapse_genotypes[idx][sample_idx] > 0 and not missing[idx][sample_idx]
            for idx in rare_indices
        )
    )

    pseudo = [
        1.0
        if any(collapse_genotypes[idx][sample_idx] > 0 for idx in rare_indices)
        else 0.0
        for sample_idx in range(nsamples)
    ]
    carriers_after_collapse = int(sum(pseudo))

    common_indices = [idx for idx in range(len(genotypes)) if idx not in set(rare_indices)]
    post_rows = [genotypes[idx][:] for idx in common_indices]
    post_types = [types[idx] for idx in common_indices]
    post_function_removed = [function_removed[idx] for idx in common_indices]
    if rare_indices:
        post_rows.append(pseudo)
        post_types.append("Collapsed")
        post_function_removed.append(False)

    remaining, collapsed_removed_by_mac = filter_after_collapse(
        post_rows,
        post_types,
        post_function_removed,
        maf_threshold,
        mac_threshold,
        aaf_filter,
    )
    carriers_after_filtering = 0 if collapsed_removed_by_mac else carriers_after_collapse
    case_ref = case_alt = cont_ref = cont_alt = None
    if phenotypes is not None:
        case_alt = int(
            sum(value for value, phenotype in zip(pseudo, phenotypes) if phenotype == 1)
        )
        cont_alt = int(
            sum(value for value, phenotype in zip(pseudo, phenotypes) if phenotype == 0)
        )
        ncases = sum(1 for phenotype in phenotypes if phenotype == 1)
        ncontrols = sum(1 for phenotype in phenotypes if phenotype == 0)
        case_ref = ncases - case_alt
        cont_ref = ncontrols - cont_alt

    return CollapseReport(
        gene=transcript.gene,
        transcript=transcript.transcript,
        input_variants=len(transcript.variants),
        collapse_candidate_variants=len(rare_indices),
        no_call_cells_before=no_call_cells,
        no_call_samples_before=no_call_samples,
        non_no_call_carriers_before=carriers_before,
        carriers_after_collapse=carriers_after_collapse,
        collapsed_removed_by_mac=collapsed_removed_by_mac,
        carriers_after_filtering=carriers_after_filtering,
        remaining_variants_after_filtering=remaining,
        case_ref=case_ref,
        case_alt=case_alt,
        cont_ref=cont_ref,
        cont_alt=cont_alt,
        collapse_variant_patterns=[
            "\t".join(transcript.variants[idx].fields[:FIRST_SAMPLE])
            for idx in rare_indices
        ],
        collapse_variant_lines=[transcript.variants[idx].line for idx in rare_indices],
    )


def write_reports(reports: Iterable[CollapseReport]) -> None:
    reports = list(reports)
    writer = csv.writer(sys.stdout, dialect="excel-tab")
    writer.writerow(
        [
            "gene",
            "transcript",
            "input_variants",
            "collapse_candidate_variants",
            "no_call_cells_before",
            "no_call_samples_before",
            "non_no_call_carriers_before",
            "carriers_after_collapse",
            "collapsed_removed_by_mac",
            "carriers_after_filtering",
            "remaining_variants_after_filtering",
            "case_ref",
            "case_alt",
            "cont_ref",
            "cont_alt",
        ]
    )
    for report in reports:
        writer.writerow(
            [
                report.gene,
                report.transcript,
                report.input_variants,
                report.collapse_candidate_variants,
                report.no_call_cells_before,
                report.no_call_samples_before,
                report.non_no_call_carriers_before,
                report.carriers_after_collapse,
                int(report.collapsed_removed_by_mac),
                report.carriers_after_filtering,
                report.remaining_variants_after_filtering,
                "" if report.case_ref is None else report.case_ref,
                "" if report.case_alt is None else report.case_alt,
                "" if report.cont_ref is None else report.cont_ref,
                "" if report.cont_alt is None else report.cont_alt,
            ]
        )

    print()
    print("## Collapse candidate grep patterns")
    print("# Use the non-comment lines below with: grep -F -f patterns.txt matrix.tsv")
    for report in reports:
        for pattern in report.collapse_variant_patterns:
            print(pattern)


def write_collapsed_matrix(path: str, header_line: str, reports: Iterable[CollapseReport]) -> None:
    with open_output(path) as handle:
        handle.write(header_line)
        handle.write("\n")
        for report in reports:
            for line in report.collapse_variant_lines:
                handle.write(line)
                handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report CAPER-style collapsed variant no-call and carrier counts."
    )
    parser.add_argument("matrix", help="CAPER matrix file, optionally .gz, or '-'")
    parser.add_argument(
        "--samples",
        help="Optional sample include file. First whitespace-delimited field is used.",
    )
    parser.add_argument(
        "--ped",
        help="Optional PED file. Uses column 2 as sample ID and column 6 as binary phenotype.",
    )
    parser.add_argument(
        "--collapse-threshold",
        type=int,
        default=10,
        help="Allele-count threshold used by --var_collapsing. Default: 10.",
    )
    parser.add_argument(
        "--maf",
        type=float,
        default=0.5,
        help="Post-collapse MAF/AAF frequency threshold. Default: 0.5.",
    )
    parser.add_argument(
        "--mac",
        type=float,
        default=math.inf,
        help="Post-collapse allele-count filter threshold. Default: no MAC filter.",
    )
    parser.add_argument(
        "--ma-count",
        action="store_true",
        help="Mirror CAPER --ma_count: use MAF filter instead of AAF filter after collapsing.",
    )
    parser.add_argument(
        "--whole-gene",
        action="store_true",
        help="Group all rows by gene instead of transcript.",
    )
    parser.add_argument(
        "--method",
        default="SKATO",
        help="Method column for the function/type whitelist. Default: SKATO.",
    )
    parser.add_argument(
        "--filter-whitelist",
        default=str(Path(__file__).resolve().parents[1] / "filter" / "filter_whitelist.csv"),
        help="CAPER filter whitelist CSV. Use 'none' to disable. Default: repo filter whitelist.",
    )
    parser.add_argument(
        "--collapsed-matrix-out",
        help="Write a matrix row subset containing only collapse-candidate variants.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sample_filter = read_sample_filter(args.samples)
    ped_phenotypes = read_ped(args.ped)
    if ped_phenotypes is not None:
        ped_samples = set(ped_phenotypes)
        sample_filter = (
            ped_samples if sample_filter is None else sample_filter.intersection(ped_samples)
        )
    whitelist = None if args.filter_whitelist == "none" else args.filter_whitelist
    allowed = read_filter_whitelist(whitelist, args.method)
    header_line, samples, transcripts = parse_matrix(
        args.matrix, sample_filter, allowed, whole_gene=args.whole_gene
    )
    phenotypes = None
    if ped_phenotypes is not None:
        phenotypes = [ped_phenotypes[sample] for sample in samples]
    reports = [
        analyze_transcript(
            transcript,
            nsamples=len(samples),
            phenotypes=phenotypes,
            collapse_threshold=args.collapse_threshold,
            maf_threshold=args.maf,
            mac_threshold=args.mac,
            aaf_filter=not args.ma_count,
        )
        for transcript in transcripts.values()
    ]
    write_reports(reports)
    if args.collapsed_matrix_out:
        write_collapsed_matrix(args.collapsed_matrix_out, header_line, reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
