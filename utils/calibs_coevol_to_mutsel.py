#!/usr/bin/env python3
import argparse
import pandas as pd
from ete3 import Tree


def is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def name_internal_node(tree):
    names = set()
    changes = list()
    for node in tree.traverse():
        if not node.name or is_float(node.name):
            if node.is_root():
                name = "Root"
            else:
                leaves = node.get_leaf_names()
                name = "".join([i[0:(int(12 / len(leaves)) + 1)] for i in leaves])
                while name in names:
                    name += "B"
                changes.append(["None" if not node.name else node.name, name])
            names.add(name)
            node.name = name
    return tree, changes


def main(args):
    # Open the tree
    tree = Tree(args.tree, format=1)
    tree, changes = name_internal_node(tree)
    if len(changes) > 0:
        tree.write(format=1, outfile=args.tree)
        print(f"The tree ({args.tree}) has been re-written because some internal nodes names have been changed:")
        print("\n".join([f"\t{c[0]} → {c[1]}." for c in changes]))

    df = []
    # Open the calibration file
    csv = pd.read_csv(args.input, sep=' ', skiprows=[0], header=None)
    # For each calibration, find the node name and add the calibration to the dataframe
    for index, (leaf_1, leaf_2, upper, lower) in csv.iterrows():
        # Get the common ancestor of the two leaves
        name = tree.get_common_ancestor([leaf_1, leaf_2]).name
        age = (lower + upper) / 2
        if lower == -1:
            age = upper
            lower = 0
        if upper == -1:
            age = lower
            upper = float("inf")
        assert lower < upper, f"Lower bound ({lower}) is greater than upper bound ({upper})"
        df += [[name, age, lower, upper]]

    # Write the dataframe to a TSV file
    header = ["NodeName", "Age", "LowerBound", "UpperBound"]
    pd.DataFrame(df).to_csv(args.output, index=False, header=header, sep="\t")
    print(f"'{args.output}' written, containing {len(df)} calibration points.")


if __name__ == '__main__':
    descr = ("A script to convert a CoEvol calibration file to a BayesCode calibration file (.tsv).\n"
             "The script needs also a newick tree file to find the most recent common ancestors.\n"
             "If internal node names are not present in the tree, they will be generated and the tree will be rewritten.\n"
             )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument('--tree', required=True, type=str, dest="tree", help="Input newick tree")
    parser.add_argument('--input', required=True, type=str, dest="input", help="Input calibration file (CoEvol format)")
    parser.add_argument('--output', required=True, type=str, dest="output",
                   help="Output calibration file (BayesCode format)")
    main(parser.parse_args())
