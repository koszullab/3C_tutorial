# -*- coding: utf-8 -*-
"""
This code assigns a restriction fragment to a locus (chrm-position).
@author: Remi Montagne (reimplementation of Axel Cournac's code)
"""
import sys
import glob
import os.path as op
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction
from Bio.Restriction import *
import argparse


##############################
#        PARSE CL arguments
##############################


def _test():
    import doctest

    doctest.testmod()


def parse_args():
    """ Gets the arguments from the command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input file containing the coordinates of Hi-C "
        "interacting pairs (dat format). Use '-' to read from stdin.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        type=str,
        help="Path to the output file (dat.indices format). Defaults to stdout.",
    )
    parser.add_argument(
        "-e",
        "--enzyme",
        required=True,
        help="The restriction enzyme used in the Hi-C protocol. E.g. DpnII.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=lambda s: map(str, s.split(",")),
        required=True,
        help="Path to the reference genome in FASTA format. Can be multiple "
        " files separated by commas.",
    )
    return parser.parse_args()


##############################
#        FUNCTIONS
##############################


def get_restriction_table(seq, enz):
    """
    Get the restriction table for a genomic sequence.

    Parameters
    ----------
    seq : Seq object
        A biopython Seq object representing a chromosomes or contig.
    enz : str
        The name of the restriction enzyme used.
    Returns
    -------
    list
        List of restriction fragment boundary positions for the input sequence.

    >>> get_restriction_table(Seq("AAGATCGATCGG"),"DpnII")
    [0, 3, 7, 12]
    >>> get_restriction_table(Seq("AA"),"DpnII")
    [0, 2]
    >>> get_restriction_table(Seq("AA"),"aeiou1")
    Traceback (most recent call last):
        ...
    ValueError: aeiou1 is not a valid restriction enzyme.
    >>> get_restriction_table("AA","DpnII")
    Traceback (most recent call last):
        ...
    TypeError: Expected Seq or MutableSeq instance, got <class 'str'> instead

    """
    # Restriction batch containing the restriction enzyme
    try:
        rb = RestrictionBatch([enz])
    except ValueError:
        raise ValueError("{} is not a valid restriction enzyme.".format(enz))
    # Conversion from string type to restriction type
    enz = rb.get(enz)
    map_restriction = rb.search(seq)
    map_restriction = map_restriction[enz]
    map_restriction = [0] + map_restriction
    map_restriction.append(len(seq))

    return map_restriction


def find_frag(pos, rl):
    """
    Use binary search to find the index of a chromosome restriction fragment
    corresponding to an input genomic position.
    Parameters
    ----------
    pos : int
        Genomic position, in base pairs.
    rl : list
        List of genomic positions corresponding to restriction sites.
    Returns
    -------
    int
        The 0-based index of the restriction fragment to which the position belongs.

    >>> find_frag(15, [0, 20, 30])
    0
    >>> find_frag(15, [10, 20, 30])
    Traceback (most recent call last):
        ...
    ValueError: The first position in the restriction table is not 0.
    >>> find_frag(31, [0, 20, 30])
    Traceback (most recent call last):
        ...
    ValueError: Read position is larger than last entry in restriction table.

    """
    if rl[0] != 0:
        raise ValueError("The first position in the restriction table is not 0.")
    if pos > rl[-1]:
        raise ValueError(
            "Read position is larger than last entry in restriction table."
        )
    bottom = 0  # indice of the first restriction fragment
    top = len(rl) - 1  # indice of the last restriction fragment
    found = False
    # binary search for the index of the read
    while not found:
        index = int((top + bottom) / 2)  # $index <- middle of the chr
        # if the position of the read is smaller than the position of the current index
        if top == bottom:
            found = True
        if pos < rl[index]:
            top = index
        elif pos >= rl[index + 1]:
            bottom = index
        else:
            found = True
    return index


def write_indices_file(infile, outfile, restriction_table):
    """
    Writes the "dat.indices" file, which has two more columns than the input
    file corresponding to the restriction fragment index of each read.

    Parameters
    ----------
    infile: file object
        File handle for the input .dat file.
    outfile: file object
        File handle for the output .dat.indices file.
    restriction_table: dict
        Dictionary with chromosome identifiers (str) as keys and list of
        positions (int) of restriction sites as values.
    """
    # eye candy progress animation
    load_anim = ["\\", "|", "/", "-"]
    prog = 0
    # Open files for reading and writing
    i = 0
    for line in infile:  # iterate over each line
        i += 1
        if i % 5000000 == 0:
            prog += 1
            print(
                "Attribution of index to reads: {0}\r".format(load_anim[prog % 4]),
                end="",
                file=sys.stderr,
            )
        chr1, pos1, sens1, chr2, pos2, sens2 = line.split()  # split it by whitespace
        pos1 = int(pos1)
        sens1 = sens1
        pos2 = int(pos2)
        sens2 = sens2

        indice1 = find_frag(pos1, restriction_table[chr1])
        indice2 = find_frag(pos2, restriction_table[chr2])

        outfile.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                chr1, pos1, sens1, indice1, chr2, pos2, sens2, indice2
            )
        )


if __name__ == "__main__":
    print(
        "Assigns restriction fragment ID to Hi-C records in a .dat file.",
        file=sys.stderr,
    )
    # Get arguments
    args = parse_args()
    enz = args.enzyme
    refs = args.reference

    restriction_table = {}
    for ref in refs:
        with open(ref, "rU") as ref_handle:
            for rec in SeqIO.parse(ref_handle, "fasta"):
                chr_name = rec.id
                # Get the lists of restriction sites of each chromosome into
                # a restriction_table
                restriction_table[chr_name] = get_restriction_table(rec.seq, enz)
                print(
                    "{0}: {1} restriction sites".format(
                        chr_name, len(restriction_table[chr_name])
                    ),
                    file=sys.stderr,
                )

    # Associate each read from .dat file with its corresponding restriction
    # fragment index and writes the resulting line
    input_handle = sys.stdin if args.input_file == "-" else open(args.input_file, "r")
    output_handle = sys.stdout if not args.output_file else open(args.output_file, "w")
    write_indices_file(input_handle, output_handle, restriction_table)
