#!/usr/bin/env python3
import argparse
from collections import defaultdict
from ete3 import Tree
import pandas as pd


def name_internal_nodes(tree: Tree) -> Tree:
    node_i = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    return tree


def main(input_tree):
    # Convert newick to nexus including root
    tree = name_internal_nodes(Tree(input_tree, format=1))
    # Write as a .tsv file, one line per node
    dico_output = defaultdict(list)
    feature_list = sorted(set([f for node in tree.traverse() for f in node.features]) - {"support", "name"})
    for node in tree.traverse():
        dico_output["node_name"].append(node.name)
        for feature in feature_list:
            dico_output[feature].append(getattr(node, feature))

    df = pd.DataFrame(dico_output)
    df.to_csv(input_tree.replace(".nhx", ".tsv"), sep="\t", index=False)
    tree.write(format=1, outfile=input_tree.replace(".nhx", ".tree"))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, default="nodeomega_gal4.Omega.nhx", type=str, dest="input")
    args = parser.parse_args()
    main(args.input)