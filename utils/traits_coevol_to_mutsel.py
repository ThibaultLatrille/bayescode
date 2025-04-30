#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np


def transform(s):
    try:
        x = float(s)
        if x < 0:
            return np.nan
        else:
            return np.log(x)
    except ValueError:
        return s


def clean_line(line):
    while "  " in line:
        line = line.replace("  ", " ")
    while "\t\t" in line:
        line = line.replace("\t\t", "\t")
    return line.replace(" ", "\t").strip()


def main(args):
    r = [clean_line(line) for line in open(args.input, 'r').readlines()]
    assert r[0].strip() == '#TRAITS'
    n_rows = int(r[1].split('\t')[0])
    n_cols = int(r[1].split('\t')[1])
    print(f"\nExpect {n_rows} taxa and {n_cols} traits from the header.\n")
    columns = r[1].split('\t')[2:]
    if len(columns) == n_cols:
        print(f"✅ Found {len(columns)} traits (columns).")
    else:
        print(f"❌ Found {len(columns)} traits (columns) instead of {n_cols} expected from the header. "
              f"You might want to check your input and output files.")
    table = [[transform(s) for s in line.split('\t')] for line in r[2:]]
    if len(table) == n_rows:
        print(f"✅ Found {len(table)} taxa (rows).")
    else:
        print(f"❌ Found {len(table)} taxa (rows) instead of {n_rows} expected from the header. "
              f"You might want to check your input and output files.")

    df = pd.DataFrame(table, columns=["TaxonName"] + columns)
    assert len(df) == len(table)
    assert len(df.columns) == len(columns) + 1
    df.to_csv(args.output, index=False, na_rep="NaN", sep="\t")
    print(f"\nData have been transformed to match the MutSel format with the function y=log(x).")
    print(f"Converted file '{args.input}' from CoEvol format (in natural-space) to MutSel format (in log-space).")
    print(f"\n'{args.output}' written, containing {len(table)} taxa and {len(columns)} traits.")


if __name__ == '__main__':
    descr = ("A script to convert a CoEvol trait file (.nexus) to a BayesCode trait file (.tsv).\n"
             "The script will read the CoEvol file and write it in a new file in MutSel format.\n"
             "The CoEvol file must contain a header with the number of taxa and the number of traits.\n"
             "The script will convert the data to log-space using the function y=log(x).\n"
             "The output file will be in tab-separated format with the first column containing the taxon names.\n"
             "The script will also check that the number of taxa and traits in the header matches the number of taxa and traits in the data.\n")

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument('--input', required=True, type=str, dest="input", help="Input trait file (CoEvol format).")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output trait file (BayesCode format).")
    main(parser.parse_args())
