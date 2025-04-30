#!/usr/bin/env python3
import argparse
from glob import glob
from itertools import chain

import numpy as np
import pandas as pd
from ete3 import Tree
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
import matplotlib.colors as colors
from matplotlib import colormaps

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def format_float(x):
    if 0.001 < x < 10:
        return "{:6.3f}".format(x)
    elif x < 10000:
        return "{:6.1f}".format(x)
    else:
        s = "{:6.2g}".format(x)
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + 'e$^{' + str(int(exp)) + '}$'
            s = " " * (5 + 6 - len(s)) + s
        return s


def label_transform(s):
    if s == "PopulationSize" or s == "LogNe":
        return 'Effective population size ($N_{\\mathrm{e}} $)'
    elif s == "MutationRate" or s == "LogMutationRatePerTime":
        return 'Mutation rate per unit of time ($\\mu$)'
    else:
        return s


def get_annot(n: Tree, f: str) -> float:
    return float(getattr(n, f))


def plot_tree(path, feature, font_size=14, line_type="-", vt_line_width=0.5, hz_line_width=0.2):
    tree = Tree(path, format=1)
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])
    vlinec, vlines, vblankline, hblankline, nodes, nodex, nodey = [], [], [], [], [], [], []

    if len(tree) < 50:
        fig = plt.figure(figsize=(16, 9))
    elif len(tree) < 70:
        fig = plt.figure(figsize=(16, 12))
    else:
        fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111)

    min_annot = min(get_annot(n, feature) for n in tree.iter_leaves() if feature in n.features)
    max_annot = max(get_annot(n, feature) for n in tree.iter_leaves() if feature in n.features)

    cmap = colormaps["inferno"]
    norm = colors.LogNorm(vmin=min_annot, vmax=max_annot)
    color_map = ScalarMappable(norm=norm, cmap=cmap)

    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))

    max_name_size = max(len(n.name) for n in tree)
    # draw tree
    rows = []
    for n in node_list:
        x = sum(n2.dist for n2 in n.iter_ancestors()) + n.dist

        min_node_annot, max_node_annot = False, False
        if (feature + "_min" in n.features) and (feature + "_max" in n.features):
            min_node_annot = get_annot(n, feature + "_min")
            max_node_annot = get_annot(n, feature + "_max")

        if n.is_leaf():
            y = node_pos[n]

            node_name = " " + n.name
            row = {"Taxon": n.name}
            if len(n.name) != max_name_size:
                node_name += " " * (max_name_size - len(n.name))
            if feature in n.features:
                node_name += " " + format_float(get_annot(n, feature))
                row[feature] = get_annot(n, feature)
            if min_node_annot and max_node_annot:
                row[feature + "Lower"] = min_node_annot
                row[feature + "Upper"] = max_node_annot
            ax.text(x, y, node_name, va='center', size=font_size)
            rows.append(row)
        else:
            y = sum([node_pos[n2] for n2 in n.children]) / len(n.children)
            node_pos[n] = y

            if feature in n.features:
                node_annot = get_annot(n, feature)
                # draw vertical line
                vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                vlines.append(node_annot)

                # draw horizontal lines
                for child in n.children:
                    child_annot = get_annot(child, feature)
                    h = node_pos[child]
                    xs = [[x, x], [x + child.dist, x + child.dist]]
                    ys = [[h - hz_line_width, h + hz_line_width], [h - hz_line_width, h + hz_line_width]]
                    zs = [[node_annot, node_annot], [child_annot, child_annot]]
                    ax.pcolormesh(xs, ys, zs, cmap=cmap, norm=norm, shading="gouraud")
            else:
                vblankline.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                for child in n.children:
                    h = node_pos[child]
                    hblankline.append(((x, h), (x + child.dist, h)))

        nodex.append(x)
        nodey.append(y)

    pd.DataFrame(rows).to_csv(path.replace('.nhx', '.tsv'), index=False, sep='\t', header=list(rows[0].keys()))
    vline_col = LineCollection(vlinec, colors=[color_map.to_rgba(l) for l in vlines],
                               linestyle=line_type,
                               linewidth=vt_line_width * 2)
    ax.add_collection(LineCollection(hblankline, colors='black', linestyle=line_type, linewidth=hz_line_width * 2))
    ax.add_collection(LineCollection(vblankline, colors='black', linestyle=line_type, linewidth=vt_line_width * 2))
    ax.add_collection(vline_col)

    # scale line
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    diffy = ymax - ymin
    dist = round((xmax - xmin) / 4, 1)
    padding = 200.
    ymin -= diffy / padding
    ax.plot([xmin, xmin + dist], [ymin, ymin], color='k')
    ax.plot([xmin, xmin], [ymin - diffy / padding, ymin + diffy / padding], color='k')
    ax.plot([xmin + dist, xmin + dist], [ymin - diffy / padding, ymin + diffy / padding],
            color='k')
    ax.text((xmin + xmin + dist) / 2, ymin - diffy / padding, dist, va='top',
            ha='center', size=font_size)

    ax.set_axis_off()
    # plt.tight_layout()
    cbar = fig.colorbar(color_map, orientation='horizontal', pad=0, shrink=0.6, ax=ax)
    cbar.ax.xaxis.set_tick_params('major', labelsize=font_size * 1.8)
    cbar.ax.xaxis.set_tick_params('minor', labelsize=font_size)
    cbar.ax.set_xlabel(label_transform(feature), labelpad=5, size=font_size * 1.8)
    plt.tight_layout()
    for o in fig.findobj():
        o.set_clip_on(False)
    plt.savefig(path.replace('.nhx', '.pdf'), format="pdf", bbox_inches='tight')
    plt.close("all")


def main(args):
    for path in glob(f"{args.input}.*.nhx"):
        feature = path.split(".")[-2]
        if feature.startswith("Branch") or feature == "ContrastPopulationSize":
            continue
        plot_tree(path, feature)
        print("Done: " + path)


if __name__ == '__main__':
    descr = (
        "This script processes and visualizes phylogenetic trees in the New Hampshire X (.nhx) format from the input directory (--input).\n"
        "It reads .nhx files from the input directory, extracts features (e.g., annotations like population size or mutation rate) from tree nodes.\n"
        "Nodes are colored based on the values of these features and branch are colored as a gradient.\n"
        "The script outputs a .pdf file for each tree visualization and a .tsv file containing extracted feature values for leaves.\n"
        "It skips certain features (e.g., Branch or ContrastPopulationSize) and adjusts scaling and formatting for the visualizations.\n"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    args = parser.parse_args()
    main(args)
