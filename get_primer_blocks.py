#!/usr/bin/env python

'''
A program to read the output of the primer design tool and
produce a TSV file containing the coordinates of the primer
blocks. Each primer block in the output is written as a triplet:
chr,start,end (just like the BED file format).

Example usage:
./get_primer_blocks.py --primers data/PALB2_heeled.idt.log --out blocks.tsv

Authors: Bernie Pope (bjpope@unimelb.edu.au)
Date: 22 September 2013
'''

from argparse import (ArgumentParser)
import sys
import csv

def parse_args():
    parser = ArgumentParser(description="Retrieve the primer block coordinates from the output of the primer design tool")
    parser.add_argument(
        '--primers', type=str, required=True,
        help='primer coordinates')
    parser.add_argument( '--out', metavar='FILE', type=str,
        help='save primer coordintates in FILE, defaults to stdout')
    return parser.parse_args() 

def get_block_coords(primers_file):
    windows = {}
    with open(primers_file) as primers:
        for line in primers:
            # we are starting a new exon
            if line.startswith('chrom:'):
                fields = line.split()
                chr = fields[1].strip()
            elif line.startswith('exon:'):
                fields = line.split()
                exon = fields[1].strip()
            elif line.startswith('Best window:'):
                fields = line.split()
                window_start = fields[2][0:-1] # drop the comma at the end
                best_blocks = windows[window_start]
                for start, end in best_blocks:
                    yield (chr, start, end)
                windows = {}
            elif line.startswith('Scoring window'):
                fields = line.split()
                window_start = fields[4]
                windows[window_start] = []
            elif line.startswith('block start:'):
                fields = line.split()
                block_start = fields[2].strip()
            elif line.startswith('block end:'):
                fields = line.split()
                block_end = fields[2].strip()
                windows[window_start].append((block_start, block_end))


def main():
    args = parse_args()
    if args.out is None:
        outfile = sys.stdout
    else:
        outfile = open(args.out, 'w')
    with outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for coord in get_block_coords(args.primers):
            writer.writerow(coord)

if __name__ == '__main__':
    main()
