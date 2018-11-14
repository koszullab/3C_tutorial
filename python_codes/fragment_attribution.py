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

##############################
#        FUNCTIONS         
##############################
def get_ref_genome_list(ref_genome_directory):
    """Gets all the files in the indicated directory and checks that they are in fasta format."""
    if not ref_genome_directory.endswith('/'): ref_genome_directory = ref_genome_directory + '/'
    ref_genome_directory = ref_genome_directory + '*'
    ref_genome_list =  glob.glob(ref_genome_directory)
    ref_genome_list =  [f for f in ref_genome_list if f.split('.')[-1] in ['fa', 'fasta', 'fna', 'fsa']]
    print('Reference genome directory: {0}.\nContains: {1}'.format(ref_genome_directory, ref_genome_list))
    
    return ref_genome_list

def get_restriction_table(ref_genome_list, enz):
    """Get the list of restriction sites for all fasta files in ref_genome_list.
    Returns it as a dictionnary with the chromosome name as a key and 
    a 2D array (1row for starts, 1 row for ends) as values."""
    print('\nRestriction enzyme: '+ enz)
    rb = RestrictionBatch([enz])    #Restriction batch containing the restriction enzyme
    enzyme = rb.get(enz)  #conversion from string type to restriction type
    
    restriction_table={} #dictionary {chromosome_name: [list_of_restriction_sites]}
    for chr_file in ref_genome_list :
        handle = open(chr_file, "rU")
        for rec in SeqIO.parse(handle, "fasta"):
            chr_name=rec.id
            map_restriction = rb.search(rec.seq)
            map_restriction = map_restriction[enzyme]
            map_restriction = [0] + map_restriction
            map_restriction.append(len(rec.seq))
            restriction_table[chr_name] = map_restriction
            print('{0}: {1} restriction sites'.format(chr_name, len(map_restriction)))
        handle.close() 
    return restriction_table

def find_frag(chrm, pos, restriction_table) :
    """Find the index of a chromosome restriction fragment"""
    T=restriction_table[chrm];
    bottom = 0            # indice of the first restriction fragment 
    top = len(T) - 1    # indice of the last restriction fragment 
    found = False
    # search the index of the read
    while not found:
        index = int((top + bottom)/2) # $index <- middle of the chr        
        # if the position of the read is smaller than the position of the current index        
        if top == bottom: found = True
        if pos < T[index]:
            top = index
        elif pos >= T[index + 1]:
            bottom = index
        else:
            found = True
    return index

def write_indices_file(infile, outfile, restriction_table):
    """Writes the output file which corresponds to the input file
    with 2 more columns, for the restriction fragment index of each read."""
    print('Attribution of index to reads.\n')
    with open(infile, 'r') as inf, open(outfile, 'w') as outf: # open the file for reading
        i=0    
        for line in inf: # iterate over each line
            i += 1
            if i % 5000000 == 0:
                print('{0} million lines have been attributed restriction fragments.'.format(i/1000000))
            chr1, pos1, sens1, chr2, pos2, sens2 = line.split() # split it by whitespace
            pos1=int(pos1);sens1=int(sens1)
            pos2=int(pos2);sens2=int(sens2)
    
            indice1=find_frag(chr1, pos1, restriction_table)
            indice2=find_frag(chr2, pos2, restriction_table)
    
            outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(chr1, pos1, sens1, indice1, chr2, pos2, sens2, indice2))


##############################
#           MAIN         
##############################
def main():
    print('Usage : python fragment_attribution.py <infile(.dat format)> <restriction_enzyme> <path/to/reference/genome>')
    
    # Get arguments
    infile = sys.argv[1]
    outfile = infile + '.indices'
    enz = sys.argv[2]
    ref_genome_directory = sys.argv[3]
    
    # Get all the fasta files of individual chromosomes of reference genome(s)
    try:
        if op.isdir(ref_genome_directory):
            print('ref genome directory')
            ref_genome_list = get_ref_genome_list(ref_genome_directory)
        else:
            print('ref genome file')
            ref_genome_list = [ref_genome_directory]
            print('ref genomes list: {0}'.format(ref_genome_list))
    except:
        print('Couldn\'t get the reference genome. Did you enter the right directory or file name ?')

    # Get the lists of restriction sites of each chromosome in restriction_table
    restriction_table = get_restriction_table(ref_genome_list, enz)
        
    # Associate each read from .dat file with its corresponding restriction fragment index and writes the resulting line
    write_indices_file(infile, outfile, restriction_table)
    
    return 0

if __name__ == '__main__':
    main()
