#!python3
import os
import argparse
from glob import glob
from ete3 import Tree
from itertools import chain
import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
import matplotlib.colors as colors

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.font_manager as font_manager

mono_font = font_manager.FontProperties(family='Latin Modern Mono', style='normal')


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
    if s == "LogPopulationSize" or s == "LogNe":
        return 'Effective population size ($N_{\\mathrm{e}} $)'
    elif s == "LogMutationRate" or s == "LogMutationRatePerTime":
        return 'Mutation rate per unit of time ($\\mu$)'
    elif s == "TraitsLogGenomeSize":
        return "Genome size (Mb)"
    elif s == "LogOmega":
        return 'Non-synonymous relative substitution rate ($\\omega$)'
    elif s == "TraitsAdult_weight_":
        return 'Adult weight (g)'
    elif s == "TraitsFemale_maturity_":
        return 'Female maturity (days)'
    elif s == "TraitsMaximum_longevity_":
        return 'Maximum longevity (years)'
    elif s == "Traitsgeneration_time":
        return 'Generation time (days)'
    elif s == "Traitslongevity":
        return 'Longevity (years)'
    elif s == "Traitsmass":
        return 'Mass (kg)'
    elif s == "Traitsmaturity":
        return 'Maturity (days)'
    elif s == "TraitspiNpiS":
        return '$\\pi_{N} / \\pi_{S}$'
    elif s == "TraitspiS":
        return '$\\pi_{S}$'
    elif s == "aa-preferences":
        return 'Amino-acid preferences'
    elif s == "ContrastPopulationSize":
        return 'Contrast Population Size'
    elif s == "Log10BranchLength":
        return 'Branch length'
    elif s == "BranchTime":
        return 'Branch time'
    else:
        return s


def get_annot(n, f):
    return np.exp(float(getattr(n, f)))


def plot_tree(tree, feature, outputpath, font_size=14, line_type="-", vt_line_width=0.5, hz_line_width=0.2,
              max_circle_size=20, min_circle_size=4):
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])
    print(feature + "\n" + outputpath)
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

    cmap = plt.get_cmap("inferno")
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
                node_name += " [{0},{1}]".format(format_float(min_node_annot), format_float(max_node_annot))
                row[feature + "Lower"] = min_node_annot
                row[feature + "Upper"] = max_node_annot
            ax.text(x, y, node_name, va='center', size=font_size, name="Latin Modern Mono")
            rows.append(row)
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
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

    pd.DataFrame(rows).to_csv(outputpath[:-4] + ".tsv", index=None, header=rows[0].keys(), sep='\t')
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
    color_map._A = []
    ticks = "PopulationSize" in feature or "MutationRatePerTime" in feature or "Omega" in feature
    cbar = fig.colorbar(color_map, ticks=(
        [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50] if ticks else None),
                        orientation='horizontal', pad=0, shrink=0.6)
    cbar.ax.xaxis.set_tick_params('major', labelsize=font_size * 1.8)
    cbar.ax.xaxis.set_tick_params('minor', labelsize=font_size)
    if ticks:
        cbar.ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    cbar.ax.set_xlabel(label_transform(feature), labelpad=5, size=font_size * 1.8)
    plt.tight_layout()
    for o in fig.findobj():
        o.set_clip_on(False)
    plt.savefig(outputpath, format=outputpath[outputpath.rfind('.') + 1:], bbox_inches='tight')
    plt.close("all")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, default="", type=str, dest="tree")
    parser.add_argument('-f', '--feature', required=False, default="Omega", type=str, dest="feature")
    parser.add_argument('-o', '--output', required=False, default="tree.pdf", type=str, dest="output")
    args = parser.parse_args()
    t = Tree(args.tree, format=1)
    plot_tree(t, args.feature, args.output)
