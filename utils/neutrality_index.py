#!/usr/bin/env python3
import os
import argparse
import itertools
from gzip import open as gzopen
from collections import defaultdict
import pandas as pd
import numpy as np
from ete3 import Tree


def open_tree(tree_path: str, format_ete3: int = 1) -> Tree:
    if tree_path.endswith(".gz"):
        newick = gzopen(tree_path).read().decode()
        return Tree(newick, format=format_ete3)
    else:
        return Tree(tree_path, format=format_ete3)


def prune_tree(input_tree: Tree, list_taxa: list = None) -> Tree:
    tree = input_tree.copy()

    # Prune tree
    if list_taxa is not None:
        tree.prune(list_taxa, preserve_branch_length=True)
        assert len(tree.get_leaves()) == len(list_taxa), f"Pruning failed: {len(tree.get_leaves())} != {len(list_taxa)}"

    # Add polytomies if branch length are 0
    remove_nodes = set([n for n in tree.traverse() if (n.dist == 0.0 and not n.is_root())])
    for n in remove_nodes:
        n.delete()
    assert len(
        set([n for n in tree.traverse()]).intersection(remove_nodes)) == 0, "Some polytomies could not be removed"
    for n in tree.traverse():
        if not n.is_root():
            assert n.dist > 0.0, f"Branch length is 0.0 for node {n.name}"
    return tree


def brownian_fitting(tree: Tree, trait: str) -> (float, float):
    leaves = tree.get_leaves()
    n = len(leaves)
    C = np.zeros((n, n))
    x = np.zeros(n)
    for i, leaf in enumerate(leaves):
        x[i] = float(getattr(leaf, trait))

    for node in tree.traverse("levelorder"):
        if node.is_root():
            node.t = 0
        else:
            # Don't forget the 4 factor (sigma^2_phy = Var(x)/4t)
            node.t = 4 * node.dist + node.up.t

    for i, j in itertools.product(range(n), range(n)):
        if i == j:
            C[i, j] = leaves[i].t
        elif i < j:
            ancestor = leaves[i].get_common_ancestor(leaves[j])
            C[i, j] = ancestor.t
            C[j, i] = ancestor.t

    invC = np.linalg.inv(C)
    ones = np.ones(n)
    v = np.dot(ones.T, invC)
    anc_z = float(np.dot(v, x)) / float(np.dot(v, ones))
    assert np.isfinite(anc_z)

    d = (x - anc_z * ones)
    var = np.dot(d.T, np.dot(invC, d)) / (n - 1)
    return anc_z, var


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(input_traits: str, input_tree: str, input_var_within: str, output_tsv: str):
    for path in [input_traits, input_tree, input_var_within]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_tsv]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Open tree, branch length should be in substitution per site
    tree = open_tree(input_tree)
    total_tree_length = sum([n.dist for n in tree.traverse() if not n.is_root()])
    print(f"The total tree length is {total_tree_length} substitutions per site")
    set_taxa = set(tree.get_leaf_names())

    # Open mean trait file
    # The column names should be in the format "*_mean" (where * is the trait name)
    # The first column should be "TaxonName"
    # "*_mean" is mandatory, but can contain missing values (ideally not)
    trait_df = pd.read_csv(input_traits, sep="\t")
    traits = [replace_last(col, "_mean", "") for col in trait_df.columns if col.endswith("_mean")]
    assert len(traits) > 0, f"No column ending with '_mean' found in {input_traits}."
    assert "TaxonName" == trait_df.columns[0], f"Column 'TaxonName' is not the first column in {input_traits}."
    trait_set_taxa = set(trait_df["TaxonName"])
    for taxa in trait_set_taxa:
        assert taxa in set_taxa, f"Taxon {taxa} in file {input_traits} not found in {input_tree}."

    # Open variance file containing the population variance and heritability
    # The column names should be in the format "*_variance" and "*_heritability" (where * is the trait name)
    # The trait names should be the same as in the mean trait file
    # The first column should be "TaxonName"
    # The second column should be "Nucleotide_diversity", and not contain any missing value
    # "*_variance" is mandatory, but can contain missing values (ideally not)
    # "*_heritability" is optional, if not found, it will be set to 1.0
    var_pop_df = pd.read_csv(input_var_within, sep="\t")
    assert "TaxonName" == var_pop_df.columns[0], f"Column 'TaxonName' is not the first column in {input_var_within}."
    assert "Nucleotide_diversity" == var_pop_df.columns[
        1], f"Column 'Nucleotide_diversity' is not the second column in {input_var_within}."
    for trait in traits:
        assert f"{trait}_variance" in var_pop_df.columns, f"Column {trait}_variance not found in {input_var_within}."
    for taxa in var_pop_df["TaxonName"]:
        assert taxa in set_taxa, f"Taxon {taxa} in file {input_var_within} not found in {input_tree}."
        assert taxa in trait_set_taxa, f"Taxon {taxa} in file {input_var_within} not found in {input_traits}."
    print("The traits found are:")
    print("\t" + "\n\t".join(traits))
    output_dict = defaultdict(list)
    for trait in traits:
        print(f"\nProcessing phenotype {trait}.")
        h, h_low, h_up = f"{trait}_heritability", f"{trait}_heritability_lower", f"{trait}_heritability_upper"
        if h_low in var_pop_df.columns and h_up in var_pop_df.columns:
            print(f"Found upper and lower bound for heritability of {trait}.")
            print(f"Using the mean of the confidence interval.")
            var_pop_df[h] = var_pop_df[[h_low, h_up]].mean(axis=1)
        elif f"{trait}_heritability" in var_pop_df.columns:
            print(f"Found heritability for {trait}.")
        else:
            print(f"Warning: column {trait}_heritability not found in {input_var_within}.")
            print("Assuming heritability = 1.0.")
            var_pop_df[h] = 1.0
        notna = (np.isfinite(var_pop_df[f"{trait}_variance"]) & (var_pop_df[f"{trait}_variance"] > 0.0))
        notna = (notna & (np.isfinite(var_pop_df[h]) & (var_pop_df[h] >= 0.0) & (var_pop_df[h] <= 1.0)))
        # Computing the genetic variance (geno = h² * pheno)
        genetic_variance_array = var_pop_df[f"{trait}_heritability"][notna] * var_pop_df[f"{trait}_variance"][notna]
        var_within_array = genetic_variance_array / (var_pop_df["Nucleotide_diversity"][notna])
        print(f'Found {len(genetic_variance_array)} species with a variance for {trait}.')
        var_within = np.mean(var_within_array)

        sp_mean_pheno = {sp: v for sp, v in zip(trait_df["TaxonName"], trait_df[f"{trait}_mean"]) if np.isfinite(v)}
        keep_leaf = [leaf.name for leaf in tree.get_leaves() if leaf.name in sp_mean_pheno]
        print(f'Found {len(keep_leaf)} species with a mean {trait}.')
        pruned_tree = prune_tree(tree, keep_leaf)
        assert len(pruned_tree.get_leaves()) == len(keep_leaf), "Error in the pruning of the tree."
        for leaf in pruned_tree.get_leaves():
            assert leaf.name in sp_mean_pheno, f"Leaf {leaf.name} not found in {input_traits}."
            setattr(leaf, trait, sp_mean_pheno[leaf.name])

        # Fitting the brownian model to obtain the phylogenetic variance
        anc_phenotype_mean, var_between = brownian_fitting(pruned_tree, trait=trait)
        ratio = var_between / var_within
        print(f"var_between = {var_between}")
        print(f"var_within = {var_within}.")
        print(f"ratio = {ratio}.")
        output_dict["trait"].append(trait)
        output_dict["nbr_taxa_between"].append(len(keep_leaf))
        output_dict["var_between"].append(var_between)
        output_dict["nbr_taxa_within"].append(len(genetic_variance_array))
        output_dict["var_within"].append(var_within)
        output_dict["ratio"].append(ratio)
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(output_tsv, sep="\t", index=False)


if __name__ == '__main__':
    descr = ("A script to obtain the neutrality index using maximum likelihood computation.\n"
             "It requires three input files:\n"
             "--tree: a phylogenetic tree in newick format, with branch lengths in number of substitutions per site (neutral markers).\n"
             "--traitsfile: a file containing the mean trait values for each species.\n"
             "--var_within: a file containing the variation within species for each trait and the genetic variation within species (neutral markers).\n"
             "The script will output a file containing the neutrality index for each trait (--output).\n"
             )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument("--traitsfile", help="Input trait file", required=True)
    parser.add_argument("--var_within", help="Input var_within file", required=True)
    parser.add_argument("--tree", help="Input tree file", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()
    main(args.traitsfile, args.tree, args.var_within, args.output)
