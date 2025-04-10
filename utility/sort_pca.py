import argparse as ap
from pathlib import Path

if __name__ == "__main__":
    PARSER = ap.ArgumentParser()
    PARSER.add_argument("genotypes", type=Path)
    PARSER.add_argument("PCA", type=Path)
    ARGS = PARSER.parse_args()

    with ARGS.genotypes.open('r') as f:
        header = f.readline()

    header = header.strip().split()[3:]

    data = []
    with ARGS.PCA.open('r') as f:
        for l in f:
            data.append(l.strip().split())

    for sample in header:
        for row in data:
            if row[0] == sample:
                if "TCGA" in sample:
                    print("1\t" + row[1] + "\t" + row[2])
                else:
                    print("0\t" + row[1] + "\t" + row[2])
