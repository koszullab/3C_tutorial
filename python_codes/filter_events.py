# -*- coding: utf-8 -*-
"""
Script to analyse the contents of a 3C library in terms of loops, uncuts,
weirds events, and filter those events.
@author: cmdoret (reimplementation of Axel KournaK's code)
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import sys
import argparse
import pysam as ps
import os


def parse_args():
    """ Gets the arguments from the command line."""
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "input_file",
        type=str,
        help="The file containing the coordinates of Hi-C interacting pairs "
        "and, the indices of their restriction fragments (.dat.indices format). "
        "Use '-' to read from stdin.",
    )
    parser.add_argument(
        "output_file",
        type=str,
        nargs="?",
        default=sys.stdout,
        help="Path to the output file (filtered dat.indices). Defaults to stdout.",
    )
    group.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        required=False,
        help="Interactively shows plots and asks the user for thresholds. "
        "Overrides predefined thresholds. Disabled by default.",
    )
    group.add_argument(
        "-t",
        "--thresholds",
        nargs=2,
        type=int,
        required=False,
        help="The minimum number of restriction fragments between reads to "
        "consider loops and uncut events. Estimated automatically by default. "
        "Must be two integers separated by a space (-t <loop> <uncut>).",
    )
    parser.add_argument(
        "-p",
        "--plot_summary",
        required=False,
        action="store_true",
        default=False,
        help="Output a piechart summarizing library composition.",
    )
    return parser.parse_args()


def process_read_pair(line, sens2type):
    """
    Takes a record (line) from a dat.indices file as input, reorders the pair so that
    read 1 in intrachromosomal pairs always has the smallest genomic coordinate.
    Parameters
    ----------
    line : str
        Read pair from a dat.indices file.
    sens2type : dict
        Dictionary with entries mapping a tuple of int (sens1, sens2) to the
        name of the event type. e.g (16, 16): '-- (weirds)'.
    Returns
    -------
    dict
        Dictionary with reordered pair where each column from the line as an entry.

    >>> process_read_pair("a 1 16 3 b 2 16 4", {(16,16):'--'})
    {'chr1': 'a',
     'pos1': 1,
     'sens1': 16,
     'indice1': 3,
     'chr2': 'b',
     'pos2': 2,
     'sens2': 16,
     'indice2': 4,
     'n_sites': 1,
     'type': 'inter'}
    """
    # Split line by whitespace
    p = line.split()
    # Saving each column as a dictionary entry with meaningful key
    cols = ["chr1", "pos1", "sens1", "indice1", "chr2", "pos2", "sens2", "indice2"]
    p = {cols[i]: p[i] for i in range(len(cols))}
    # Transforming numeric columns to int
    for col in ["pos1", "sens1", "indice1", "pos2", "sens2", "indice2"]:
        p[col] = int(p[col])

    # invert records for intrachromosomal pairs where rec2 comes before
    # rec1 in genomic coordinates
    if p["chr1"] == p["chr2"] and p["pos2"] < p["pos1"]:
        p["sens1"], p["sens2"] = p["sens2"], p["sens1"]
        p["pos1"], p["pos2"] = p["pos2"], p["pos1"]
        p["indice1"], p["indice2"] = p["indice2"], p["indice1"]

    # Number of restriction sites separating reads in the pair
    p["n_sites"] = p["indice2"] - p["indice1"]
    # Get event type
    if p["chr1"] == p["chr2"]:
        p["type"] = sens2type[(p["sens1"], p["sens2"])]
    else:
        p["type"] = "inter"

    return p


def get_thresholds(in_dat, interactive=False):
    """
    Analyses the events in the first million of Hi-C pairs in the library, plots
    the occurrences of each event type according to number of restriction
    fragments, and asks user interactively for the minimum threshold for uncuts
    and loops.

    Parameters
    ----------
    in_dat: str
        Path to the input .dat file containing Hi-C pairs
    interactive: bool
        If True, plots are diplayed and thresholds are require
    Returns
    -------
    dictionary
        dictionary with keys "uncuts" and "loops" where the values are the
        corresponding thresholds entered by the user.
    """
    thr_uncut = None
    thr_loop = None
    max_sites = 500
    n_events = {
        event: np.zeros((1, max_sites))
        for event in ["++ (weirds)", "-- (weirds)", "+- (uncuts)", "-+ (loops)"]
    }
    i = 0
    # Map of sense -> name of event for intrachromosomal pairs.
    sens2type = {
        (0, 0): "++ (weirds)",
        (16, 16): "-- (weirds)",
        (0, 16): "+- (uncuts)",
        (16, 0): "-+ (loops)",
    }
    # open the file for reading (just the first 1 000 000 lines)
    with open(in_dat) as f:
        for line in f:
            i += 1
            if i == 1000000:
                break
            # Process Hi-C pair into a dictionary
            p = process_read_pair(line, sens2type)
            # Type of event and number of restriction site between reads
            etype = p["type"]
            nsite = p["nsites"]
            # Count number of events for intrachrom pairs
            if etype != "inter" and nsites < max_sites:
                n_events[etype][nsites] += 1

    def plot_event(n_events, name):
        plt.xlim([0, 15])
        plot(
            range(1, n_events[name].shape[1] + 1),
            n_events[name],
            "o-",
            label=name,
            linewidth=2.0,
        )

    if interactive:
        # PLot:
        plot_events(n_events, "++ (weirds)")
        plot_events(n_events, "-- (weirds)")
        plot_events(n_events, "+- (uncuts)")
        plot_events(n_events, "-+ (loops)")
        grid()
        plt.xlabel("Number of restriction fragment(s)")
        plt.ylabel("Number of events")
        plt.yscale("log")
        legend()
        # ehavior_file = in_dat.replace(".dat.indices", "_behavior.png")
        # savefig(behavior_file)
        show()

        #  Scanning the alignment file again and count the different events with the determined thresholds:
        thr_uncut = int(input("Enter threshold for the uncuts events (+-):"))
        thr_loop = int(input("Enter threshold for the loops events (-+):"))
    else:
        # Iterate over sites, from furthest to closest
        for site in range(max_sites)[::-1]:
            expected = (
                n_events["++ (weirds)"][site] + n_events["-- (weirds)"][site]
            ) / 2
            # Rough estimation: Keep the last (closest) site where the number
            # occurrences of loop/uncut differs from the exp. value by <= 0.1%
            if abs(n_events["+- (uncuts)"] - expected) <= expected / 1000:
                thr_uncut = site
            if abs(n_events["-+ (loops)"] - expected) <= expected / 1000:
                thr_loop = site
            plot_events(n_events, "-+ (loops)")
        if thr_uncut is None or thr_loop is None:
            raise ValueError(
                "The threshold for loops or uncut could not be estimated. "
                "Please try running with -i to investigate the problem."
            )
    return thr_uncut, thr_loop


def filter_events(in_dat, out_filtered, thr_uncut, thr_loop, plot_events=False):
    """
    Filter out spurious intrachromosomal Hi-C pairs from input file. +- pairs
    that do not exceed the uncut threshold and -+ pairs that do not exceed the
    loop thresholds are excluded from the ouput file. All others are written.
    Parameters
    ----------
    in_dat : str
        Path to the input "dat" file containing Hi-C pairs.
    out_filtered : str
        Path to the output filtered "dat" file.
    thr_uncut : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal +- pair.
    thr_loop : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal -+ pair.
    plot_events : bool
        If True, a plot summarising the proportion of each type of event will be
        shown after filtering.
    """
    n_uncuts = 0
    n_loops = 0
    n_weirds = 0
    n_int = 0
    lrange_intra = 0
    lrange_inter = 0
    n_mito = 0

    # Map of sense -> name of event for intrachromosomal pairs.
    sens2type = {
        (0, 0): "++ (weirds)",
        (16, 16): "-- (weirds)",
        (0, 16): "+- (uncuts)",
        (16, 0): "-+ (loops)",
    }
    i = 0
    # open the files for reading and writing
    with open(in_dat, "r") as inf, open(out_filtered, "w") as outf:
        for line in inf:  # iterate over each line
            p = process_read_pair(line, sens2type)
            nsites = indice2 - indice1
            if chr1 == chr2:
                if indice1 == indice2 and sens1 == sens2:
                    n_weirds += 1
                elif nsites <= thr_loop and (sens1 == 16 and sens2 == 0):
                    n_loops += 1
                elif nsites <= thr_uncut and (sens1 == 0 and sens2 == 16):
                    n_uncuts += 1
                else:
                    lrange_intra += 1
                    outf.write(
                        str(chr1)
                        + "\t"
                        + str(p["pos1"])
                        + "\t"
                        + str(p["sens1"])
                        + "\t"
                        + str(p["indice1"])
                        + "\t"
                        + str(p["chr2"])
                        + "\t"
                        + str(p["pos2"])
                        + "\t"
                        + str(p["sens2"])
                        + "\t"
                        + str(p["indice2"])
                        + "\n"
                    )
            if chr1 != chr2:
                lrange_inter += 1
                outf.write(
                    str(p["chr1"])
                    + "\t"
                    + str(p["pos1"])
                    + "\t"
                    + str(p["sens1"])
                    + "\t"
                    + str(p["indice1"])
                    + "\t"
                    + str(p["chr2"])
                    + "\t"
                    + str(p["pos2"])
                    + "\t"
                    + str(p["sens2"])
                    + "\t"
                    + str(p["indice2"])
                    + "\n"
                )
            if (
                chr1 in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]
                and chr2 not in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]
            ) or (
                chr2 in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]
                and chr1 not in ["CM007981.1", "NC_001224.1", "chrM", "Mito"]
            ):
                n_mito += 1

    if lrange_inter > 0:
        ratio_inter = float(lrange_inter) / float(lrange_intra + lrange_inter) * 100.0
        ratio_mito = float(n_mito) / float(lrange_inter) * 100.0
    else:
        ratio_inter = 0
        ratio_mito = 0

    # Dump quick summary of operation results into stderr
    print("Proportion of inter contacts: {}%".format(ratio_inter), file=sys.stderr)
    print(
        "Pairs discarded: Loops: {0}, Uncuts: {1}, Weirds: {2}".format(
            n_loops, n_uncuts, n_weirds
        ),
        file=sys.stderr,
    )
    print("Pairs kept: {}".format(lrange_intra + lrange_inter), file=sys.stderr)

    # Visualize summary interactively if requested by user
    if plot_events:

        # Plot: make a square figure and axes to plot a pieChart:
        figure(1, figsize=(6, 6))
        ax = axes([0.1, 0.1, 0.8, 0.8])
        # The slices will be ordered and plotted counter-clockwise.
        labels = "Uncuts", "Loops", "Weirds", "3D intra", "3D inter"
        fracs = [n_uncuts, n_loops, n_weirds, lrange_intra, lrange_inter]
        colors = ["salmon", "lightskyblue", "lightcoral", "palegreen", "plum"]
        pie(
            fracs,
            labels=labels,
            colors=colors,
            autopct="%1.1f%%",
            shadow=True,
            startangle=90,
        )
        title("Distribution of library events", bbox={"facecolor": "1.0", "pad": 5})
        plt.text(
            0.3,
            1.15,
            "Threshold Uncuts =" + str(thr_uncut),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            0.3,
            1.05,
            "Threshold Loops =" + str(thr_loop),
            fontdict=None,
            withdash=False,
        )

        plt.text(
            -1.5,
            -1.2,
            "Total number of reads =" + str(i),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.5,
            -1.3,
            "Ratio inter/(intra+inter) =" + str(ratio_inter) + "%",
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.5,
            -1.4,
            "selected reads = {0}%".format(
                float(lrange_inter + lrange_intra)
                / (n_loops + n_uncuts + n_weirds + n_mito + lrange_inter + lrange_intra)
            ),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.6,
            -1.5,
            "selected reads = {0}%".format(ratio_mito),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.5,
            -1.6,
            "Ratio Mito/inter =" + str(ratio_mito) + "%",
            fontdict=None,
            withdash=False,
        )
        piechart_file = infile.replace(".dat.indices", "_piechart.png")
        show()


if __name__ == "__main__":
    args = parse_args()
    input_file = sys.stdin if args.input_file == "-" else args.input_file
    output_file = args.output_file
    if args.thresholds:
        # Thresholds supplied by user beforehand
        uncut_thr, loop_thr = args.thresholds
    else:
        # Threshold defined at runtime
        uncut_thr, loop_thr = get_thresholds(input_file, interactive=args.interactive)
    # Filter library and write to output file
    filter_events(
        input_file, output_file, uncut_thr, loop_thr, plot_events=args.plot_summary
    )
