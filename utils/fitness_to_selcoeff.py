#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

# Define the codon table
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})


def sel_coeff(profiles, anc_aa, der_aa, codon_site):
    # The selection coefficient is the difference between the fitness of the
    # derived amino acid and the fitness of the ancestral amino acid

    # The index of the site is 1-based in the input file, so we need to
    # subtract 1 to get the correct index in 0-based Python
    site = codon_site - 1

    if der_aa in ["-", "X"] or anc_aa in ["-", "X"]:
        # If one of the codon is a stop codon or an indeterminate codon, return "NaN"
        return "NaN"
    elif der_aa == anc_aa:
        # If the transition is synonymous, return 0
        return 0
    else:
        # If the transition is non-synonymous, return the selection coefficient as the log of the ratio of the
        # Wrightian fitness of the derived allele to the Wrightian fitness of the ancestral allele
        return round(np.log(profiles[der_aa][site] / profiles[anc_aa][site]), 6)


def main(args):
    assert args.sep in ["\t", ","]

    # Read the fitness profiles (Wrightian fitness) from the tab-delimited input file
    profiles = pd.read_csv(args.input_profiles, sep="\t")
    # Read the transitions from the tab-delimited input file
    df = pd.read_csv(args.input_transitions, sep=args.sep, dtype={"site": int, "anc": "string", "der": 'string'})

    # Assert that the indices in the transitions match the number of sites in the profiles
    assert min(df["site"]) >= 1, "The first index in the transitions file must be at least 1 (1-based counting)."
    assert max(df["site"]) <= len(profiles["site"]), f"The codon indices in the transitions ({max(df['site'])}) " \
                                                     f"is greater than the number of codon sites in the " \
                                                     f"profiles (n={len(profiles['site'])})."

    # Create new columns for the ancestral and derived amino acids
    df["anc_aa"] = df.apply(lambda row: codontable[row["anc"]], axis=1)
    df["der_aa"] = df.apply(lambda row: codontable[row["der"]], axis=1)
    # Create a new column for the selection coefficient
    df["S"] = df.apply(lambda row: sel_coeff(profiles, row["anc_aa"], row["der_aa"], row["site"]), axis=1)

    # Write the output as a tab-delimited file
    df.to_csv(args.output, sep=args.sep, index=False)


if __name__ == '__main__':
    # parse the arguments
    descr = (
        "A script to convert a fitness profile to selection coefficients.\n"
        "The script will read the fitness profiles and the transitions file and write it in a new file with the selection coefficients.\n"
        "--input_profiles: the fitness profiles obtained by running readnodemutsel or readaamutsel with the option --ss.\n"
        "--input_transitions: the transitions file, which must contain the columns 'site', 'anc' and 'der'.\n"
        "Selection coefficients are calculated as the log of the ratio of the fitness of the derived allele to the ancestral allele:\n"
        "S = log(W(der)/W(anc))\n"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument('--input_transitions', required=True, type=str, dest="input_transitions",
                        help="Input transitions file, must contain the columns 'site', 'anc' , and 'der'.\n"
                             "'site' is the 1-based index of the codon site, between 1 and n,"
                             " where n is the number of codons in the alignment.\n"
                             "'anc' is the 3-letter ancestral codon.\n"
                             "'der' is the 3-letter derived codon.\n")
    parser.add_argument('--input_profiles', required=True, type=str, dest="input_profiles",
                        help="Input fitness profiles file (obtained by running readnodemutsel or readaamutsel with the option --ss).\n"
                             "Contains n fitness profiles where n is the number of codons in the alignment.")
    parser.add_argument('--output', required=True, type=str, dest="output",
                        help="Output transition file also containing the columns 'anc_aa' (ancestral amino acid), "
                             "'der_aa' (derived amino acid), and 'S' (scaled selection coefficient).")
    parser.add_argument('--sep', required=False, type=str, dest="sep", default="\t",
                        help="Separator (either '\\t' or ',') for the input and output transition files.")

    # run the main function with the arguments from the command line
    main(parser.parse_args())
