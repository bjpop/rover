#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
import logging
import sys
import pysam
import re
import os
from operator import itemgetter
import csv
from version import rover_version

# proportion of read which must overlap region of interest
default_minimum_read_overlap = 0.9
default_proportion_threshold = 0.05
default_absolute_threshold = 2

def parse_args():
    "Consider mapped reads to amplicon sites"

    parser = ArgumentParser(description="Consider mapped reads to amplicon sites")
    parser.add_argument(
    '--version', action='version', version='%(prog)s ' + rover_version)
    parser.add_argument(
        '--primers', type=str, required=True,
        help='File name of primer coordinates in TSV format.')
    parser.add_argument(
        '--overlap', type=float, default=default_minimum_read_overlap,
        help='Minimum fraction overlap of read to block region. '
             'Defaults to {}.'.format(default_minimum_read_overlap))
    parser.add_argument(
        'bams', nargs='+', type=str, help='bam files containing mapped reads')
    parser.add_argument( '--log', metavar='FILE', type=str,
        help='Log progress in FILENAME, defaults to stdout.')
    parser.add_argument('--out', metavar='FILE', type=str,
        required=True, help='Name of output file containing called variants.')
    parser.add_argument('--proportionthresh', metavar='N', type=float,
        default=default_proportion_threshold,
        help='Keep variants which appear in this proportion of the read pairs for '
             'a given target region, and bin otherwise. '
             'Defaults to {}.'.format(default_proportion_threshold))
    parser.add_argument('--absthresh', metavar='N', type=int,
        default=default_absolute_threshold,
        help='Only keep variants which appear in at least this many read pairs. '
             'Defaults to {}.'.format(default_absolute_threshold))
    parser.add_argument('--coverdir',
        required=False,
        help='Directory to write coverage files, defaults to current working directory.')
    return parser.parse_args() 


def get_block_coords(primers_file):
    with open(primers_file) as primers:
        return list(csv.reader(primers, delimiter='\t'))


def lookup_reads(min_overlap, bam, chr, start_col, end_col):
    # arguments are in zero-based indices
    total_reads = 0
    overlapping_reads = 0
    read_pairs = {}
    for read in bam.fetch(chr, start_col, end_col+1):
        total_reads += 1
        # only keep reads which overlap with the block region by a certain proportion
        overlap = proportion_overlap(start_col, end_col, read) 
        if overlap > min_overlap:
            overlapping_reads += 1
            if read.qname not in read_pairs:
                read_pairs[read.qname] = [read]
            else:
                read_pairs[read.qname].append(read)
    logging.info("number of reads intersecting block: {}".format(total_reads))
    logging.info("number of reads sufficiently overlapping block: {}".format(overlapping_reads))
    return read_pairs

def get_MD(read):
    for tag, val in read.tags:
        if tag == 'MD':
            return val
    return None

#M   BAM_CMATCH  0
#I   BAM_CINS    1
#D   BAM_CDEL    2
#N   BAM_CREF_SKIP   3
#S   BAM_CSOFT_CLIP  4
#H   BAM_CHARD_CLIP  5
#P   BAM_CPAD    6
#=   BAM_CEQUAL  7
#X   BAM_CDIFF   8

# find all the variants in a single read (SNVs, Insertions, Deletions)
def read_variants(name, chr, pos, aligned_bases, cigar, md):
    cigar_orig = cigar
    md_orig = md
    seq_index = 0
    result = []

    while cigar and md:
        cigar_code, cigar_segment_extent = cigar[0]
        next_md = md[0]

        if cigar_code == 0:
            # Cigar Match
            if isinstance(next_md, MD_match):
                # MD match 
                if next_md.size >= cigar_segment_extent:
                    next_md.size -= cigar_segment_extent
                    if next_md.size == 0:
                        md = md[1:]
                    cigar = cigar[1:]
                    pos += cigar_segment_extent
                    seq_index += cigar_segment_extent
                else:
                    # next_md.size < cigar_segment_extent
                    cigar = [(cigar_code, cigar_segment_extent - next_md.size)] + cigar[1:]
                    md = md[1:]
                    pos += next_md.size
                    seq_index += next_md.size
            elif isinstance(next_md, MD_mismatch):
                 # MD mismatch
                 seq_base = aligned_bases[seq_index] 
                 result.append(SNV(chr, pos, next_md.ref_base, seq_base))
                 cigar = [(cigar_code, cigar_segment_extent - 1)] + cigar[1:]
                 md = md[1:]
                 pos += 1
                 seq_index += 1
            elif isinstance(next_md, MD_deletion):
                # MD deletion, should not happen in Cigar match
                logging.info("MD del in cigar match {} {}".format(md_orig, cigar_orig))
                exit()
            else:
                logging.info("unexpected MD code {}".format(md_orig))
                exit()
        elif cigar_code == 1: 
            # Insertion
            seq_bases = aligned_bases[seq_index:seq_index + cigar_segment_extent]
            result.append(Insertion(chr, pos, seq_bases))
            cigar = cigar[1:]
            seq_index += cigar_segment_extent
            # pos does not change
        elif cigar_code == 2:
            # Deletion
            if isinstance(next_md, MD_deletion):
                result.append(Deletion(chr, pos, next_md.ref_bases))
                md = md[1:]
                cigar = cigar[1:]
                pos += cigar_segment_extent
                # seq_index does not change
            else:
                logging.info("Non del MD in Del Cigar".format(md_orig, cigar_orig))
                exit()
        else:
            logging.info("unexpected cigar code {}".format(cigar_orig))
            exit()
    return result

class SNV(object):
    def __init__(self, chr, pos, ref_base, seq_base):
        self.chr = chr
        self.pos = pos
        self.ref_base = ref_base
        self.seq_base = seq_base
    def __str__(self):
        return "S: {} {} {} {}".format(self.chr, self.pos, self.ref_base, self.seq_base)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chr, self.pos, self.ref_base, self.seq_base)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return self.ref_base
    def alt(self):
        return self.seq_base


class Insertion(object):
    def __init__(self, chr, pos, inserted_bases):
        self.chr = chr
        self.pos = pos
        self.inserted_bases = inserted_bases
    def __str__(self):
        return "I: {} {} {}".format(self.chr, self.pos, self.inserted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chr, self.pos, self.inserted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return '-'
    def alt(self):
        return self.inserted_bases

class Deletion(object):
    def __init__(self, chr, pos, deleted_bases):
        self.chr = chr
        self.pos = pos
        self.deleted_bases = deleted_bases
    def __str__(self):
        return "D: {} {} {}".format(self.chr, self.pos, self.deleted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chr, self.pos, self.deleted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return self.deleted_bases
    def alt(self):
        return '-'

class MD_match(object):
    def __init__(self, size):
        self.size = size
    def __str__(self):
        return str(self.size)
    def __repr__(self):
        return self.__str__()

class MD_mismatch(object):
    def __init__(self, ref_base):
        self.ref_base = ref_base
    def __str__(self):
        return self.ref_base
    def __repr__(self):
        return self.__str__()

class MD_deletion(object):
    def __init__(self, ref_bases):
        self.ref_bases = ref_bases
    def __str__(self):
        return "^" + self.ref_bases
    def __repr__(self):
        return self.__str__()

# [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
def parse_md(md, result):
    if md:
        number_match = re.match('([0-9]+)(.*)', md)
        if number_match is not None:
            number_groups = number_match.groups()
            number = int(number_groups[0])
            md = number_groups[1]
            return parse_md_snv(md, result + [MD_match(number)])
    return result

def parse_md_snv(md, result):
    if md:
        snv_match = re.match('([A-Z])(.*)', md)
        if snv_match is not None:
            snv_groups = snv_match.groups()
            ref_base = snv_groups[0]
            md = snv_groups[1]
            return parse_md(md, result + [MD_mismatch(ref_base)])
        else:
            return parse_md_del(md, result)
    return result

def parse_md_del(md, result):
    if md:
        del_match = re.match('(\^[A-Z]+)(.*)', md)
        if del_match is not None:
            del_groups = del_match.groups()
            ref_bases = del_groups[0][1:]
            md = del_groups[1]
            return parse_md(md, result + [MD_deletion(ref_bases)])
    return result

def proportion_overlap(block_start, block_end, read):
    read_end = read.pos + read.rlen - 1
    if read_end < block_start or read.pos > block_end:
        # they don't overlap
        return 0.0
    else:
        overlap_start = max(block_start, read.pos)
        overlap_end = min(block_end, read_end)
        overlap_size = overlap_end - overlap_start + 1
        return float(overlap_size) / read.rlen

def write_variant(file, variant, sample):
    file.write(
        '\t'.join([variant.chr, str(variant.pos), str(variant.pos), variant.ref(),
                   variant.alt(), "comments: " + sample]) + '\n')

def process_blocks(args, kept_variants_file, binned_variants_file, bam, sample, block_coords):
    coverage_info = []
    for chr, start, end in block_coords:
        start = int(start)
        end = int(end)
        logging.info("processing block chr: {}, start: {}, end: {}".format(chr, start, end))
        # process all the reads in one block
        block_vars = {}
        num_pairs = 0
        # use 0 based coordinates to lookup reads from bam file
        read_pairs = lookup_reads(args.overlap, bam, chr, start - 1, end - 1)
        for read_name, reads in read_pairs.items():
            if len(reads) == 1:
                logging.warning("read {} with no pair".format(read_name))
            elif len(reads) == 2:
                num_pairs += 1
                read1, read2 = reads
                variants1 = read_variants(read1.qname, chr, read1.pos + 1, read1.query, read1.cigar, parse_md(get_MD(read1), []))
                variants2 = read_variants(read2.qname, chr, read2.pos + 1, read2.query, read2.cigar, parse_md(get_MD(read2), []))
                set_variants1 = set(variants1)
                set_variants2 = set(variants2)
                # find the variants each read in the pair share in common
                same_variants = set_variants1.intersection(set_variants2)
                for var in same_variants:
                    # only consider variants within the bounds of the block
                    if var.pos >= start and var.pos <= end:
                        if var in block_vars:
                            block_vars[var] += 1
                        else:
                            block_vars[var] = 1
            else:
                logging.warning("read {} with more than 2".format(read_name))
        logging.info("number of read pairs in block: {}".format(num_pairs))
        logging.info("number of variants found in block: {}".format(len(block_vars)))
        for var in block_vars:
            num_vars = block_vars[var]
            proportion = float(num_vars) / num_pairs
            proportion_str = "{:.2f}".format(proportion)
            if num_vars >= args.absthresh and proportion >= args.proportionthresh:
                write_variant(kept_variants_file, var, sample)
            else:
                write_variant(binned_variants_file, var, sample)
        coverage_info.append((chr, start, end, num_pairs))
    coverage_filename = sample + '.coverage'
    if args.coverdir is not None:
        coverage_filename = os.path.join(args.coverdir, coverage_filename)
    with open(coverage_filename, 'w') as coverage_file:
        coverage_file.write('chr\tblock_start\tblock_end\tnum_pairs\n')
        for chr, start, end, num_pairs in sorted(coverage_info, key=itemgetter(3)):
            coverage_file.write('{}\t{}\t{}\t{}\n'.format(chr, start, end, num_pairs))



# output_header = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "NUM_PAIRS_WITH_VAR", "NUM_PAIRS_AT_POS", "PERCENT"])

def process_bams(args):
    block_coords = get_block_coords(args.primers)
    with open(args.out, "w") as kept_variants_file, \
         open(args.out + '.binned', "w") as binned_variants_file:
        #kept_variants_file.write(output_header + '\n')
        #binned_variants_file.write(output_header + '\n')
        for bam_filename in args.bams:
            base = os.path.basename(bam_filename)
            sample = base.split('.')
            if len(sample) > 0:
                sample = sample[0]
            else:
                exit('Cannot deduce sample name from bam filename {}'.format(bam_filename))
            with pysam.Samfile(bam_filename, "rb") as bam:
                logging.info("processing bam file {}".format(bam_filename))
                process_blocks(args, kept_variants_file,
                               binned_variants_file, bam, sample, block_coords)

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
    process_bams(args)



if __name__ == '__main__':
    main()
