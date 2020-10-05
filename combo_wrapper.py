"""
Ryan Bohlender 2020

Wrapper for CARVA that runs a low pass whole-dataset run, and a high pass run on the remaining (more significant) genes.
This automates the process of running CARVA multiple times, and automatically combines the output. This should result
in performance benefits when the number of remaining genes to be processed at a high number of permutations is less
than the number of cores available for processing.

# For initial versions
 - MGIT will be excluded. We don't currently have a good way to combine the output across transcripts.
"""

import argparse as ap
import subprocess as sp
from pathlib import Path
import sys
from scipy.stats import chi2
import numpy as np


SIMPLE_MAP = {
    "Rank": 0,
    "Gene": 1,
    "Transcript": 2,
    "Test_Statistic": 3,
    "Exact_P": 4,
    "Empirical_P": 5,
    "Empirical_P_CI": 6,
    "Empirical_MidP": 7,
    "Empirical_MidP_CI": 8,
    "MGIT_P": 9,
    "Successes": 10,
    "MGIT_Successes": 11,
    "Permutations": 12
}


def run_carva(cmd):
    """

    :param cmd: CARVA run command
    :return:
    """
    proc = sp.Popen(cmd, stderr=sp.PIPE)
    found_cmd = False
    next_cmd = None
    for l in proc.stderr:
        l = l.decode()
        print(l, file=sys.stderr, end='')
        if found_cmd:
            next_cmd = l
            break
        if l.startswith("Some genes did not reach the success threshold. Run the following command to check those genes."):
            found_cmd = True
    proc.wait()
    return next_cmd


def build_cmd(args, max_perms=None):
    """

    :type args: str
    :param args:
    :param max_perms: Optional number of maximum permutations. Used to set the total in the second pass.
    :return:
    """
    cmd = args.strip().split()
    if max_perms:
        try:
            max_perms_idx = cmd.index('--max_perms') + 1
            cmd[max_perms_idx] = f'{max_perms}'
        except ValueError:
            cmd.append('--max_perms')
            cmd.append(str(max_perms))
    return cmd


def handle_pass_one_output(args):
    """Rename files to preserve them for the second run.

    :type args: str
    :param args: Argument string for CARVA call.
    :return:
    """
    fpath = extract_output(args)
    mod_path = extract_output(args).with_suffix('')
    fpath = fpath.rename(mod_path)

    return fpath


def extract_output(args):
    """Construct output path from CARVA arguments

    :type args: str
    :param args: Argument string for CARVA call. Used to build output.
    :return: Path to file
    """
    split_args = args.strip().split()
    try:
        out_idx = split_args.index('-o') + 1
    except ValueError:
        out_idx = split_args.index('--output') + 1
    # If the ValueError branch fails then we want to raise an error because the command is incomplete
    out_path = split_args[out_idx]

    try:
        method_idx = split_args.index('-m') + 1
    except ValueError:
        try:
            method_idx = split_args.index('--method') + 1
        except ValueError:  # Default branch
            method_idx = None
    if method_idx:
        method = split_args[method_idx]
    else:
        method = "VAAST"

    try:
        grouping_idx = split_args.index('-g') + 1
    except ValueError:
        try:
            grouping_idx = split_args.index('--grouping') + 1
        except ValueError:
            grouping_idx = None
    if grouping_idx:
        grouping = split_args[grouping_idx]
    else:
        grouping = 0

    try:
        testable_idx = split_args.index('--testable')
    except ValueError:
        testable_idx = None

    try:
        biallelic_idx = split_args.index('--biallelic')
    except ValueError:
        biallelic_idx = None

    try:
        range_idx = (split_args.index('--range') + 1, split_args.index('--range') + 2)
    except ValueError:
        range_idx = None

    constructed_path = out_path
    constructed_path += '/' + method
    if grouping != 0:
        constructed_path += '.' + f'g{grouping}'
    if testable_idx:
        constructed_path += '.' + 'testable'
    if biallelic_idx:
        constructed_path += '.' + 'biallelic'
    if range_idx:
        constructed_path += '.' + f'{range_idx[0]}' + '.' + f'{range_idx[1]}'
    constructed_path += '.simple'
    return Path(constructed_path)


def combine_output(pass_one, pass_two):
    """

    :type pass_one: Path
    :type pass_two: Path
    :param pass_one: Output of the initial pass through the whole dataset
    :param pass_two: Output of the second pass through the remaining genes
    :return:
    """
    data_one = parse_simple(pass_one)
    data_two = parse_simple(pass_two)

    combined = []
    for k, v in data_one.items():
        if k in data_two:
            new_line = data_two[k].split()
            succ_1 = int(v.split()[SIMPLE_MAP["Successes"]])
            succ_2 = int(data_two[k].split()[SIMPLE_MAP["Successes"]])
            perm_1 = int(v.split()[SIMPLE_MAP["Permutations"]])
            perm_2 = int(data_two[k].split()[SIMPLE_MAP["Permutations"]])
            new_line[SIMPLE_MAP["Successes"]] = f'{succ_1 + succ_2}'
            new_line[SIMPLE_MAP["Permutations"]] = f'{perm_1 + perm_2}'
            new_line[SIMPLE_MAP["Empirical_P"]] = f'{calc_empirical_p(new_line)}'
            new_line[SIMPLE_MAP["Empirical_P_CI"]] = f'{calc_empirical_ci(new_line)}'
            combined.append(new_line)
        else:
            combined.append(v.strip().split())

    pvals = np.zeros(len(combined))
    for i, v in enumerate(combined):
        pvals[i] = float(v[SIMPLE_MAP["Empirical_P"]])

    out_path = Path(str(pass_two) + '.combine')

    order = np.argsort(pvals)
    rank = 0
    with out_path.open('w') as f:
        for i in order:
            rank += 1
            print(format_output(rank, combined[i]), file=f)


def format_output(rank, line):
    """

    :type rank: int
    :type line: List[str]
    :param rank: The rank of the transcript being formatted
    :param line: The line to be formatted
    :return: Formatted line
    """
    return f'{rank}\t{line[SIMPLE_MAP["Gene"]]}\t{line[SIMPLE_MAP["Transcript"]]}\t{line[SIMPLE_MAP["Empirical_P"]]}\t{line[SIMPLE_MAP["Empirical_P_CI"]]}\t{line[SIMPLE_MAP["Successes"]]}\t{line[SIMPLE_MAP["Permutations"]]}'


def parse_simple(fpath):
    """

    Transcript id in column 2, zero-indexed.
    :type fpath: Path
    :param fpath:
    :return: dict containing values identified by transcript
    """
    data = {}
    with fpath.open('r') as f:
        for l in f:
            if l.startswith('#') or l.startswith('Rank'):
                continue
            sp_l = l.strip().split()
            data[sp_l[SIMPLE_MAP["Transcript"]]] = l.strip()
    return data


def calc_empirical_p(line):
    """

    :type line: List[str]
    :param line: Entry from the .simple output
    :return: Calculated empirical p-value
    """
    succ = int(line[SIMPLE_MAP["Successes"]])
    perm = int(line[SIMPLE_MAP["Permutations"]])

    return (succ + 1) / (perm + 1)


def calc_empirical_ci(line):
    """

    :type line: List[str]
    :param line: Entry from the .simple output
    :return: Calculated empirical p-value
    """
    succ = float(line[SIMPLE_MAP["Successes"]])
    perm = int(line[SIMPLE_MAP["Permutations"]])

    lo_ci = 0
    df = 2 * succ

    if succ > 0:
        lo_ci = chi2.ppf(0.025, df) / 2.
    hi_ci = chi2.ppf(0.975, df + 2) / 2.

    return f'{lo_ci / perm},{hi_ci / perm}'


def main():
    """Entrypoint.
    Handles program arguments and program run logic.
    """
    parser = ap.ArgumentParser()
    parser.add_argument('--carva_call', required=True, help="Command used to execute inital CARVA run.")
    parser.add_argument('--maximum_permutations', required=True,
                        help="Total number of permutations to run in the second pass.")
    args = parser.parse_args()

    call_parts = args.carva_call.strip().split()
    try:
        nperm = int(call_parts[call_parts.index('--nperm') + 1])
    except ValueError:
        nperm = 10000
    total_perms = int(args.maximum_permutations) - nperm

    cmd = build_cmd(args.carva_call)
    next_cmd = run_carva(cmd)
    pass_one_output = handle_pass_one_output(args.carva_call)
    cmd = build_cmd(next_cmd, str(total_perms))
    run_carva(cmd)
    pass_two_output = extract_output(" ".join(cmd))

    combine_output(pass_one_output, pass_two_output)


if __name__ == "__main__":
    main()
