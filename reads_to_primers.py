#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
import logging
import sys

def parse_args():
    "Consider mapped reads to amplicon sites"

    parser = ArgumentParser(description="Consider mapped reads to amplicon sites")
    parser.add_argument(
        '--primers', type=str, required=True,
        help='primer coordinates')
    parser.add_argument( '--log', metavar='FILE', type=str,
        help='log progress in FILENAME, defaults to stdout')
    return parser.parse_args() 

def get_amplicon_coords(primers_file):
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
                #print(exon)
                for start, end in best_blocks:
                    print("{} {} {}".format(chr, start, end))
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
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
    logging.basicConfig(
        filename=args.log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))
    get_amplicon_coords(args.primers)


if __name__ == '__main__':
    main()
