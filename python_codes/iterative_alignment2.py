# -*- coding: utf-8 -*-
"""
Created on Jan 26 2018
@author: Remi Montagne
aligns iteratively reads from a 3C fastq file
"""

import argparse
import os
import subprocess as sp
import shutil as st
from random import randint
from threading import Thread
import compressed_utils as ct

##############################
#        ARGUMENTS
##############################


def parse_args():
    """ Gets the arguments from the command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument('infile1',
                        help='the end1 mate with extension fastq')
    parser.add_argument('infile2',
                        help='the end2 mate with extension fastq')
    parser.add_argument('-i', '--index', required=True,
                        help='the prefix of the bowtie-built index of '
                        'reference genome')
    parser.add_argument('-p', '--nb_processors', default=1, type=int,
                        help='number of CPUs used for alignment.')
    parser.add_argument('-o1', '--outfile1', default='alignement_1',
                        help='name of the outputfile, without extension')
    parser.add_argument('-o2', '--outfile2', default='alignement_2',
                        help='name of the outputfile, without extension')
    return parser.parse_args()


##############################
#        FUNCTIONS
##############################
def generate_temp_dir(infile):
    """Generates a temporary file named temp by default.
    If the directory already exists, ask another name."""
    directory = infile + 'temporary'
    for i in range(5):
        directory = directory + str(randint(0, 9))
    directory = directory + '/'
    if os.path.exists(directory):
        raise SystemExit('The directory already exists. Please try again.')
    else:
        os.makedirs(directory)
    return directory


def iterative_align(infile, temp_directory, index,
                    nb_processors, outfile):
    """Aligns iteratively the reads of infile with bowtie2."""
    # length of the fragments to align
    n = 20
    # prepare the final output file where the alignments will be reported
    outfile = outfile + '.sam.0'
    f = open(outfile, 'w')
    f.close()
    # set with the name of the unaligned reads :
    my_set = set()
    total_reads = 0

    # Bowtie only accepts uncompressed fastq: uncompress it into a temp file
    if ct.is_compressed(infile):
        uncomp_path = os.path.join(temp_directory,
                                   os.path.basename(infile) + ".tmp")
        with ct.read_compressed(infile) as inf:
            with open(uncomp_path, 'w') as uncomp:
                st.copyfileobj(inf, uncomp)
    else:
        uncomp_path = infile

    # Counting reads
    with open(uncomp_path, 'r') as inf:
        for line in inf:
            total_reads += 1
    total_reads /= 4

    # use first read to guess read length. Stripping newline.
    with open(uncomp_path, 'r') as inf:
        size = inf.readline()
        size = len(inf.readline().rstrip())

    print('{0} reads to parse'.format(total_reads))

    # iterative alignment per se
    while n <= size:
        print('\n' + '-' * 10 + '\nn = {0}'.format(n))

        # Generate a temporary input fastq file with the n first nucleotids
        # of the reads.
        print('Generating truncated reads')
        truncated_reads = truncate_reads(temp_directory, uncomp_path, my_set, n)

        # Align the truncated reads on reference genome
        print('Aligning reads')
        temp_alignment = '{0}/temp_alignment.sam'.format(temp_directory)
        sp.call('bowtie2 -x {0} -p {1} --rdg 500,3 --rfg 500,3 --sam-no-hd --sam-no-sq --quiet --very-sensitive -S {2} {3}'.format(index, nb_processors, temp_alignment, truncated_reads), shell=True)

        # Sort the reads: the reads whose truncated end was aligned are writen
        # in the output file.
        # The reads whose truncated end was not aligned are going to be kept
        # for the next round.
        print('Reporting aligned reads')
        my_set = set()
        my_set = sort_samfile(temp_alignment, outfile, my_set)

        # Writes the unaligned reads in a temporary fastq file.
        # This file will become the input file for the next round.
        # print('Preparing temporary fastq files of unmapped reads')
        # infile = get_next_infile(infile, temp_directory, dico, n)
        n += 20
    
    # one last round
    print('\n' + '-' * 10 + '\nn = {0}'.format(size))
    print('Generating truncated reads')
    truncated_reads = truncate_reads(temp_directory, uncomp_path, my_set, size)
    print('Aligning reads')
    sp.call('bowtie2 -x {0} -p {1} --rdg 500,3 --rfg 500,3 --sam-no-hd --sam-no-sq --quiet --very-sensitive -S {2} {3}'.format(index, nb_processors, temp_alignment, truncated_reads), shell=True)
    print('Reporting aligned reads')
    my_set = set()
    my_set = sort_samfile(temp_alignment, outfile, my_set)

    remaining = len(my_set)
    # Report the unaligned reads as well
    with open(temp_alignment, 'r') as inf:
        with open(outfile, 'a') as outf:
            for line in inf:
                line = line.split()
                name = line[0]
                if name in my_set:
                    status = int(line[1])
                    quality = int(line[4])
                    chromosome = line[2]
                    position = int(line[3])
                    outf.write('{0} {1} {2} {3} {4}\n'.format(name, chromosome, status, position, quality))

    print('{0} reads aligned / {1} total reads.'.format(total_reads - remaining, total_reads))


    return 0


def prepare_my_set(infile):
    """Returns a dictionary with the name of all the reads in the infile
    and the indication unmapped (U)"""
    my_set = set()
    with open(infile, 'r') as f:
        line = f.readline()
        i = 0
        while line:
            if i % 4 == 0:
                my_set.add(line.split()[0].replace('@', ''))
            line = f.readline()
            i += 1
    print('{0} reads to map'.format(len(my_set)))
    return my_set


def truncate_reads(temp_directory, infile, my_set, n):
    """Gets the reads from infile and writes the n first
    nucleotids of each sequence in an auxialiary file in
    the temporary folder 'directory'."""
    outfile = '{0}/truncated.fastq'.format(temp_directory)
    with open(infile, 'r') as inf:
        # the reads whose name is in my_set, i. e. the reads that still have 
        # to be aligned unless n = 25: in this case, report all the reads
        with open(outfile, 'w') as outf:
            row = 0
            for line in inf:
                row += 1
                if row % 4 == 1:
                    name = line.split()[0].replace('@', '')
                    report = True if name in my_set or n == 20 else False
                if report:
                    if row % 2 != 0:
                        outf.write(line)
                    else:
                        outf.write(line[0:n]+'\n')

    return outfile


def sort_samfile(temp_alignment, outfile, my_set):
    """
    Reads all the reads in the global input file (infile).
    If they are aligned with a good quality, write them in the global output
    file (outfile). Else, add their name in the set my_set for the next round
    of alignment.
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the final output file.
    # Keep all fragment name and mapping status in a set for next step
    # (see iterative_align())

    with open(temp_alignment, 'r') as tempf:
        with open(outfile, 'a') as outf:
            for line in tempf:
                line = line.split()
                name = line[0]
                quality = int(line[4])
                status = int(line[1])
                if status in [0, 16] and quality >= 30 and quality < 255:
                    chromosome = line[2]
                    position = int(line[3])
                    outf.write('{0} {1} {2} {3} {4}\n'.format(name, chromosome, status, position, quality))
                else:
                    my_set.add(name)
        print('{0} reads left to map.'.format(len(my_set)))

    return my_set


def get_next_infile(infile, temp_directory, my_set, n):
    """
    Writes an intermediary input file (fastq format) where all the aligned
    reads are stocked for the next round of aligment.
    """
    temp_fastq = '{0}/temporary{1}_1.fastq'.format(temp_directory, n)

    with open(infile, 'r') as inf:
        with open(temp_fastq, 'w') as outf:
            line = inf.readline()
            while line:
                name = line.split()[0].replace('@', '')
                if name in my_set:
                    report = True
                for i in range(4):
                    if report:
                        outf.write(line)
                    line = inf.readline()

    return temp_fastq

##############################
#            MAIN
##############################


def main():
    print('\n' + '-' * 40 + '\nBegining of program\n')
    # Get arguments
    args = parse_args()
    infile1 = args.infile1
    infile2 = args.infile2
    index = args.index
    nb_processors = args.nb_processors
    outfile1 = args.outfile1
    outfile2 = args.outfile2

    # Generates temporary folder
    temp_directory1 = generate_temp_dir(infile1)
    temp_directory2 = generate_temp_dir(infile2)

    # Aligns iteritatively the fastq file
    print('\n' + '-' * 20 + '\nComputing iterative alignment\n')
    t1 = Thread(target=iterative_align, args=(infile1, temp_directory1, index, nb_processors, outfile1))
    t2 = Thread(target=iterative_align, args=(infile2, temp_directory2, index, nb_processors, outfile2))
    t1.start()
    t2.start()
    t1.join()
    t2.join()

    # Deletes the temporary folder
    st.rmtree(temp_directory1)
    st.rmtree(temp_directory2)
    print('\nEnd of program\n' + '-' * 40 + '\n')

if __name__ == '__main__':
    main()
