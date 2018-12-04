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
        type=argparse.FileType("r"),
        help="Path to the input file containing the coordinates of Hi-C "
        "interacting pairs (dat format). Use '-' to read from stdin.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        default=sys.stdout,
        help="Path to the output file (dat.indices format). Defaults to stdout.",
    )
    parser.add_argument(
        "-e",
        "--enzyme",
        required=True,
        help="The restriction enzyme used in the Hi-C protocol. E.g. DpnII",
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Path to the reference genome " "in FASTA format.",
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


def find_frag(chrm, pos, restriction_table):
    """
    Find the index of a chromosome restriction fragment
    Parameters
    ----------
    chrm : str
        The identifier of the chromosome
    pos: int
        A genomic position, in base pairs
    Returns
    -------
    int
        The index of the restriction fragment to which the position belongs
    """
    T = restriction_table[chrm]
    bottom = 0  # indice of the first restriction fragment
    top = len(T) - 1  # indice of the last restriction fragment
    found = False
    # binary search for the index of the read
    while not found:
        index = int((top + bottom) / 2)  # $index <- middle of the chr
        # if the position of the read is smaller than the position of the current index
        if top == bottom:
            found = True
        if pos < T[index]:
            top = index
        elif pos >= T[index + 1]:
            bottom = index
        else:
            found = True
    return index


def write_indices_file(infile, outfile, restriction_table):
    """
    Writes the output file which corresponds to the input file
    with 2 more columns, for the restriction fragment index of each read.

    Parameters
    ----------
    infile: str
        Path to the input file
    outfile: str
        Path to the output file
    restriction_table: dict
        Dictionary with chromosome identifiers (str) as keys and list of
        positions (int) of restriction sites as values.
    """
    print("Attribution of index to reads.\n")
    with open(infile, "r") as inf, open(
        outfile, "w"
    ) as outf:  # open the file for reading
        i = 0
        for line in inf:  # iterate over each line
            i += 1
            if i % 5000000 == 0:
                print(
                    "{0} million lines have been attributed restriction fragments.".format(
                        i / 1000000
                    )
                )
            chr1, pos1, sens1, chr2, pos2, sens2 = (
                line.split()
            )  # split it by whitespace
            pos1 = int(pos1)
            sens1 = int(sens1)
            pos2 = int(pos2)
            sens2 = int(sens2)

            indice1 = find_frag(chr1, pos1, restriction_table)
            indice2 = find_frag(chr2, pos2, restriction_table)

            outf.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    chr1, pos1, sens1, indice1, chr2, pos2, sens2, indice2
                )
            )


if __name__ == "__main__":
    _test()
    print("Assigns restriction fragment ID to Hi-C records in a .dat file.")
    args = parse_args()
    # Get arguments
    input_file = sys.stdin if args.input_file == "-" else args.input_file
    output_file = args.output_file
    enz = args.enzyme
    ref = args.reference

    restriction_table = {}
    with open(ref, "rU") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            chr_name = rec.id
            # Get the lists of restriction sites of each chromosome into
            # a restriction_table
            restriction_table[chr_name] = get_restriction_table(rec.seq, enz)
            print(
                "{0}: {1} restriction sites".format(
                    chr_name, len(restriction_table[chr_name])
                )
            )

    # Associate each read from .dat file with its corresponding restriction
    # fragment index and writes the resulting line
    write_indices_file(infile, outfile, restriction_table)
