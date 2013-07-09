#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
import logging
import sys
import pysam
import re

def parse_args():
    "Consider mapped reads to amplicon sites"

    parser = ArgumentParser(description="Consider mapped reads to amplicon sites")
    parser.add_argument(
        '--primers', type=str, required=True,
        help='primer coordinates')
    parser.add_argument(
        '--overlap', type=float, default=0.5,
        help='minimum percentage overlap of read to amplicon region')
    parser.add_argument(
        'bam', type=str, help='bam file containing mapped reads')
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

def lookup_reads(min_overlap, bam, chr, start_col, end_col):
    # arguments are in zero-based indices
    read_pairs = {}
    for read in bam.fetch(chr, start_col, end_col+1):
        # only keep reads which overlap with the amplicon region by a certain percentage
        overlap = percent_overlap(start_col, end_col, read) 
        if overlap > min_overlap:
            if read.qname not in read_pairs:
                read_pairs[read.qname] = [read]
            else:
                read_pairs[read.qname].append(read)
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
            seq_bases = aligned_bases[seq_index:seq_index + cigar_segment_extent + 1]
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

class Insertion(object):
    def __init__(self, chr, pos, inserted_bases):
        self.chr = chr
        self.pos = pos
        self.inserted_bases = inserted_bases
    def __str__(self):
        return "I: {} {} {}".format(self.chr, self.pos, self.inserted_bases)
    def __repr__(self):
        return str(self)

class Deletion(object):
    def __init__(self, chr, pos, deleted_bases):
        self.chr = chr
        self.pos = pos
        self.deleted_bases = deleted_bases
    def __str__(self):
        return "D: {} {} {}".format(self.chr, self.pos, self.deleted_bases)
    def __repr__(self):
        return str(self)


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

def percent_overlap(amplicon_start, amplicon_end, read):
    # XXX not sure about indels
    read_end = read.pos + read.rlen - 1
    if read_end < amplicon_start or read.pos > amplicon_end:
        # they don't overlap
        return 0.0
    else:
        overlap_start = max(amplicon_start, read.pos)
        overlap_end = min(amplicon_end, read_end)
        overlap_size = overlap_end - overlap_start + 1
        return float(overlap_size) / read.rlen

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
    with pysam.Samfile(args.bam, "rb") as bam:
        for chr, start, end in get_amplicon_coords(args.primers):
            read_pairs = lookup_reads(args.overlap, bam, chr, int(start), int(end))
            for read_name, reads in read_pairs.items():
                if len(reads) == 1:
                    logging.info("read {} with no pair".format(read_name))
                elif len(reads) == 2:
                    read1, read2 = reads
                    variants1 = read_variants(read1.qname, chr, read1.pos + 1, read1.query, read1.cigar, parse_md(get_MD(read1), []))
                    variants2 = read_variants(read2.qname, chr, read2.pos + 1, read2.query, read2.cigar, parse_md(get_MD(read2), []))
                    print("read {}".format(read_name))
                    print("vars1: {}".format(variants1))
                    print("vars2: {}".format(variants2))
                else:
                    logging.info("read {} with more than 2".format(read_name))


if __name__ == '__main__':
    main()
