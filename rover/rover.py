#!/usr/bin/env python

"""Rover - ROVER-PCR Variant Caller: read-pair overlap considerate
variant-calling software for PCR-based massively parallel sequencing datasets"""

from argparse import ArgumentParser
import datetime
import logging
import sys
import pysam
import re
import os
import vcf
from operator import itemgetter
import csv
from version import rover_version
from itertools import (izip, chain, repeat)
from Bio import pairwise2

# proportion of block which must be overlapped by read
DEFAULT_MINIMUM_READ_OVERLAP_BLOCK = 0.9
DEFAULT_PROPORTION_THRESHOLD = 0.05
DEFAULT_ABSOLUTE_THRESHOLD = 2
DEFAULT_PRIMER_THRESHOLD = 5
DEFAULT_GAP_PENALTY = 2.0

def parse_args():
    """ Consider mapped reads to amplicon sites."""
    parser = ArgumentParser(description="Consider mapped reads to amplicon \
        sites")
    parser.add_argument('--version', action='version', version='%(prog)s ' \
        + rover_version)
    parser.add_argument(
        '--primers', type=str, required=True,
        help='File name of primer coordinates in TSV format.')
    parser.add_argument(
        '--overlap', type=float, default=DEFAULT_MINIMUM_READ_OVERLAP_BLOCK,
        help='Minimum proportion of block which must be overlapped by a read. '
             'Defaults to {}.'.format(DEFAULT_MINIMUM_READ_OVERLAP_BLOCK))
    parser.add_argument(
        'bams', nargs='+', type=str, help='BAM files containing mapped reads.')
    parser.add_argument('--log', metavar='FILE', type=str, \
        help='Log progress in FILENAME, defaults to stdout.')
    parser.add_argument('--out', metavar='FILE', type=str, \
        required=True, help='Name of output file containing called variants.')
    parser.add_argument('--proportionthresh', metavar='N', type=float, \
        default=DEFAULT_PROPORTION_THRESHOLD, \
        help='Keep variants which appear in this proportion of the read pairs '
             'for a given target region, and bin otherwise. '
             'Defaults to {}.'.format(DEFAULT_PROPORTION_THRESHOLD))
    parser.add_argument('--absthresh', metavar='N', type=int, \
        default=DEFAULT_ABSOLUTE_THRESHOLD, \
        help='Only keep variants which appear in at least this many \
        read pairs. ' \
        'Defaults to {}.'.format(DEFAULT_ABSOLUTE_THRESHOLD))
    parser.add_argument('--qualthresh', metavar='N', type=int, \
        help='Minimum base quality score (phred).')
    parser.add_argument('--primercheck', metavar='FILE', type=str, \
    help='Expected base sequences and locations of primers as determined by a \
    primer generating program.')
    parser.add_argument('--primerthresh', metavar='N', type=int, \
        default=DEFAULT_PRIMER_THRESHOLD, \
        help='Maximum allowed variance in base sequence of primers. Defaults \
        to {}.'.format(DEFAULT_PRIMER_THRESHOLD))
    parser.add_argument('--gap_penalty', metavar='N', type=float, \
        default=DEFAULT_GAP_PENALTY, \
        help='Score deduction on gap in pairwise2 global alignment. Deduction \
        of gap_penalty * 2 for opening a gap and gap_penalty for extending \
        a gap. Defaults to {}.'.format(DEFAULT_GAP_PENALTY))
    parser.add_argument('--id_info', type=str, \
    help='File containing rs ID information.')
    parser.add_argument('--coverdir', required=False, \
        help='Directory to write coverage files, defaults to current working \
        directory.')
    return parser.parse_args()

def get_block_coords(primers_file):
    with open(primers_file) as primers:
        return list(csv.reader(primers, delimiter='\t'))

def get_primer_sequence(primers_coords_file):
    with open(primers_coords_file) as primer_coords:
        return list(csv.reader(primer_coords, delimiter='\t'))

def lookup_reads(min_overlap, bam, chrsm, start_col, end_col):
    total_reads = 0
    overlapping_reads = 0
    read_pairs = {}
    for read in bam.fetch(chrsm, start_col, end_col+1):
        total_reads += 1
        # only keep reads which overlap with the block region by a
        # certain proportion
        overlap = proportion_overlap(start_col, end_col, read)
        if overlap > min_overlap:
            overlapping_reads += 1
            if read.qname not in read_pairs:
                read_pairs[read.qname] = [read]
            else:
                read_pairs[read.qname].append(read)
    logging.info("Number of reads intersecting block: {}".format(total_reads))
    logging.info("Number of reads sufficiently overlapping block: {}" \
        .format(overlapping_reads))
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
def read_variants(args, name, chrsm, pos, aligned_bases, cigar, md):
    cigar_orig = cigar
    md_orig = md
    seq_index = 0
    result = []
    context = None

    while cigar and md:
        cigar_code, cigar_segment_extent = cigar[0]
        next_md = md[0]
        if cigar_code == 0:
            if isinstance(next_md, MD_match):
                # MD match
                if next_md.size >= cigar_segment_extent:
                    next_md.size -= cigar_segment_extent
                    if next_md.size == 0:
                        md = md[1:]
                    context = aligned_bases[seq_index + \
                    cigar_segment_extent - 1].base
                    cigar = cigar[1:]
                    pos += cigar_segment_extent
                    seq_index += cigar_segment_extent
                else:
                    # next_md.size < cigar_segment_extent
                    cigar = [(cigar_code, cigar_segment_extent - \
                        next_md.size)] + cigar[1:]
                    context = aligned_bases[seq_index + next_md.size - 1].base
                    md = md[1:]
                    pos += next_md.size
                    seq_index += next_md.size
            elif isinstance(next_md, MD_mismatch):
                # MD mismatch
                seq_base_qual = aligned_bases[seq_index]
                # check if the read base is above the minimum quality score
                if (args.qualthresh is None) or (seq_base_qual.qual >= \
                    args.qualthresh):
                    seq_base = seq_base_qual.base
                    result.append(SNV(chrsm, pos, next_md.ref_base, seq_base, \
                        '.', None))
                else:
                    seq_base = seq_base_qual.base
                    result.append(SNV(chrsm, pos, next_md.ref_base, seq_base, \
                        '.', ";qlt"))
                cigar = [(cigar_code, cigar_segment_extent - 1)] + cigar[1:]
                context = next_md.ref_base
                md = md[1:]
                pos += 1
                seq_index += 1
            elif isinstance(next_md, MD_deletion):
                # MD deletion, should not happen in Cigar match
                logging.info("MD del in cigar match {} {}".format(md_orig, \
                    cigar_orig))
                exit()
            else:
                logging.info("Unexpected MD code {}".format(md_orig))
                exit()
        elif cigar_code == 1:
            # Insertion
            seq_bases_quals = aligned_bases[seq_index:seq_index + \
            cigar_segment_extent]
            seq_bases = ''.join([b.base for b in seq_bases_quals])
            # check that all the bases are above the minimum quality threshold
            if (args.qualthresh is None) or all([b.qual >= args.qualthresh for \
                b in seq_bases_quals]):
                result.append(Insertion(chrsm, pos, seq_bases, '.', None, \
                    context))
            else:
                result.append(Insertion(chrsm, pos, seq_bases, '.', ";qlt", \
                    context))
            cigar = cigar[1:]
            seq_index += cigar_segment_extent
            # pos does not change
        elif cigar_code == 2:
            # Deletion
            if isinstance(next_md, MD_deletion):
                seq_base = aligned_bases[seq_index]
                if seq_base.qual >= args.qualthresh:
                    result.append(Deletion(chrsm, pos, next_md.ref_bases, '.', \
                        None, context))
                else:
                    result.append(Deletion(chrsm, pos, next_md.ref_bases, '.', \
                        ";qlt", context))
                context = next_md.ref_bases[-1]
                md = md[1:]
                cigar = cigar[1:]
                pos += cigar_segment_extent
                # seq_index does not change
            else:
                logging.info("Non del MD in Del Cigar".format(md_orig, \
                    cigar_orig))
                exit()
        elif cigar_code == 4:
            # soft clipping
            context = 'S'
            md = md[1:]
            cigar = cigar[1:]
            seq_index += cigar_segment_extent
        elif cigar_code == 5:
            # hard clipping
            context = 'H'
            md = md[1:]
            cigar = cigar[1:]
        else:
            logging.info("Unexpected cigar code {}".format(cigar_orig))
            exit()
    return result

# SAM/BAM files store the quality score of a base as a byte (ascii character)
# in "Qual plus 33 format". So we subtract off 33 from the ascii code
# to get the actual score
# See: http://samtools.sourceforge.net/SAMv1.pdf
# ASCII codes 32 an above are the so-called printable characters, but 32
# is a whitespace character, so SAM uses 33 and above.
def ascii_to_phred(ascii):
    return ord(ascii) - 33

def make_base_seq(name, bases, qualities):
    """ Take a list of DNA bases and a corresponding list of quality scores
    and return a list of Base objects where the base and score are
    paired together."""
    num_bases = len(bases)
    num_qualities = len(qualities)
    if num_bases <= num_qualities:
        return [Base(b, ascii_to_phred(q)) for (b, q) in izip(bases, qualities)]
    else:
        logging.warning("In read {} fewer quality scores {} than bases {}" \
            .format(name, num_qualities, num_bases))
        # we have fewer quality scores than bases
        # pad the end with 0 scores (which is ord('!') - 33)
        return [Base(b, ascii_to_phred(q)) \
            for (b, q) in izip(bases, chain(qualities, repeat('!')))]

class Base(object):
    """ A DNA base paired with its quality score."""
    def __init__(self, base, qual):
        self.base = base # a string
        self.qual = qual # an int
    def as_tuple(self):
        return (self.base, self.qual)
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __str__(self):
        return str(self.as_tuple())
    def __repr__(self):
        return str(self)
    def __hash__(self):
        return hash(self.as_tuple)

class SNV(object):
    # bases are represented just as DNA strings
    def __init__(self, chrsm, pos, ref_base, seq_base, qual, filter_reason):
        self.chrsm = chrsm
        self.pos = pos
        self.ref_base = ref_base
        self.seq_base = seq_base
        self.qual = qual
        self.filter_reason = filter_reason
        self.info = []
    def __str__(self):
        return "S: {} {} {} {}".format(self.chrsm, self.pos, self.ref_base, \
            self.seq_base)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chrsm, self.pos, self.ref_base, self.seq_base)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return self.ref_base
    def alt(self):
        return self.seq_base
    def fil(self):
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        return self.pos
    def quality(self):
        return self.qual

class Insertion(object):
    # bases are represented just as DNA strings
    def __init__(self, chrsm, pos, inserted_bases, qual, filter_reason, \
        context):
        self.chrsm = chrsm
        self.pos = pos
        self.inserted_bases = inserted_bases
        self.qual = qual
        self.filter_reason = filter_reason
        self.info = []
        self.context = context
        if self.context == None:
            self.info.append("BS=T")
            self.context = '-'
        elif self.context == 'S':
            self.info.append("SC=T")
            self.context = '-'
        elif self.context == 'H':
            self.info.append("HC=T")
            self.context = '-'
    def __str__(self):
        return "I: {} {} {}".format(self.chrsm, self.pos, self.inserted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chrsm, self.pos, self.inserted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return self.context
    def alt(self):
        return self.context + self.inserted_bases
    def fil(self):
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        return self.pos - 1
    def quality(self):
        return self.qual

class Deletion(object):
    # bases are represented just as DNA strings
    def __init__(self, chrsm, pos, deleted_bases, qual, filter_reason, context):
        self.chrsm = chrsm
        self.pos = pos
        self.deleted_bases = deleted_bases
        self.qual = qual
        self.filter_reason = filter_reason
        self.info = []
        self.context = context
        if self.context == None:
            self.info.append("BS=T")
            self.context = '-'
        elif self.context == 'S':
            self.info.append("SC=T")
            self.context = '-'
        elif self.context == 'H':
            self.info.append("HC=T")
            self.context = '-'
    def __str__(self):
        return "D: {} {} {}".format(self.chrsm, self.pos, self.deleted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        return (self.chrsm, self.pos, self.deleted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        return self.context + self.deleted_bases
    def alt(self):
        return self.context
    def fil(self):
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        return self.pos - 1
    def quality(self):
        return self.qual

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
    """ Compute the proportion of the block that is overlapped by the read

          block_start               block_end
               |-------------------------|

        ^---------------------------^
    read.pos                      read_end

               |--------------------|
         overlap_start        overlap_end

    """
    read_end = read.pos + read.rlen - 1
    if read.rlen <= 0:
        # read is degenerate, zero length
        # treat it as no overla$p
        logging.warn("Degenerate read: {}, length: {}".format(read.qname, \
            read.rlen))
        return 0.0
    if read_end < block_start or read.pos > block_end:
        # they don't overlap
        return 0.0
    else:
        overlap_start = max(block_start, read.pos)
        overlap_end = min(block_end, read_end)
        overlap_size = overlap_end - overlap_start + 1
        block_size = block_end - block_start + 1
        return float(overlap_size) / block_size

def write_variant(vcf_file, variant, id_info, args):
    """ Writes variant to vcf_file, while also finding the relevant rs number
    from dbsnp if applicable."""
    info = 0
    record_info = None
    # If the variant is deemed a "PASS", find the relevant rs number from dbsnp.
    if variant.fil() == "PASS" and args.id_info:
        for record in id_info.fetch(variant.chrsm, variant.position(), variant.\
            position() + max(len(variant.ref()), len(variant.alt())) + 1):
            if record.POS == variant.position() and record.REF == \
            variant.ref() and (variant.alt() in record.ALT):
                info = 1
                record_info = record
    if info == 1:
        vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
str(record_info.ID), variant.ref(), variant.alt(), variant.quality(), variant\
.fil(), ';'.join(variant.info)]) + '\n')
    else:
        vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
'.', variant.ref(), variant.alt(), variant.quality(), variant.fil(), ';'\
.join(variant.info)]) + '\n')

def nts(none_string):
    """ Turns None into an empty string."""
    if none_string is None:
        return ''
    return str(none_string)

def reverse_complement(sequence):
    """ Return the reverse complement of a DNA string."""
    complementary_bases = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    rc_bases = []
    for base in sequence:
        rc_bases.append(complementary_bases[str(base)])
    rc_seq = "".join([b for b in rc_bases])
    return rc_seq[::-1]

def possible_primer(primer_sequence, block_info, bases, direction):
    """ Generates the sequence of bases at the location where we expect
    the primer to be located"""
    forward_primer_length = len(primer_sequence[block_info[3]])
    reverse_primer_length = len(primer_sequence[block_info[4]])

    if direction == -1:
        primer_bases = []
        for primer_base in bases[:forward_primer_length]:
            primer_bases.append(primer_base.base)
        return "".join([b for b in primer_bases])

    if direction == 1:
        primer_bases = []
        for primer_base in bases[-1 * reverse_primer_length:]:
            primer_bases.append(primer_base.base)
        return "".join([b for b in primer_bases])

def primer_diff(primer1, primer2, gap_penalty):
    """ Compares two primers (in string representation), and assigns a score
    using the global alignment algorithm from the pairwise2 module.

    Score is calculated as follows: 1 point for a match, 0 points for a
    mismatch, -2 * gap_penalty for opening a gap and -1 * gap_penalty for
    extending a gap. Returned as the length of the primer minus the score,
    which means a return value of 0 corresponds to an exact match."""
    score = pairwise2.align.globalxs(primer2, primer1, -2 * gap_penalty, -1 * \
        gap_penalty, score_only=1)
    if isinstance(score, float):
        return len(primer1) - score
    else:
        return len(primer1)

def check_forward_primer(primer_sequence, block_info, bases, gap_penalty):
    """ Returns the score for the forward primer."""
    ref_primer_forward = primer_sequence[block_info[3]]
    forward_primer_region = possible_primer(primer_sequence, block_info, \
        bases, -1)
    return primer_diff(ref_primer_forward, forward_primer_region, gap_penalty)

def check_reverse_primer(primer_sequence, block_info, bases, gap_penalty):
    """ Returns the score for the reverse primer."""
    ref_primer_reverse = primer_sequence[block_info[4]]
    reverse_primer_region = possible_primer(primer_sequence, block_info, \
        bases, 1)
    return primer_diff(ref_primer_reverse, reverse_complement\
        (reverse_primer_region), gap_penalty)

def process_blocks(args, kept_variants_file, bam, sample, block_coords, \
    primer_sequence, data, data2, id_info):
    coverage_info = []
    total_scores = {}
    data.write('\t'.join(["Primer name", "0.0", "1.0", "2.0", "3.0", "4.0", \
        "5.0", "6.0", "7.0", "8.0", "9.0"]))
    data.write('\n')
    for block_info in block_coords:
        # process all the reads in one block
        chrsm, start, end = block_info[:3]
        start = int(start)
        end = int(end)
        logging.info("Processing block chrsm: {}, start: {}, end: {}"\
            .format(chrsm, start, end))
        block_vars = {}
        num_pairs = 0
        num_discards = 0
        num_unexpected = 0
        forward_scores = {}
        reverse_scores = {}
        # use 0 based coordinates to lookup reads from bam file
        read_pairs = lookup_reads(args.overlap, bam, chrsm, start - 1, end - 1)
        for read_name, reads in read_pairs.items():
            if len(reads) == 1:
                logging.warning("Read {} with no pair".format(read_name))
            elif len(reads) == 2:
                discard = 0
                num_pairs += 1
                read1, read2 = reads
                read1_bases = make_base_seq(read1.qname, read1.query, \
                    read1.qqual)
                read2_bases = make_base_seq(read2.qname, read2.query, \
                    read2.qqual)
                variants1 = read_variants(args, read1.qname, chrsm, read1.pos \
                    + 1, read1_bases, read1.cigar, parse_md(get_md(read1), []))
                variants2 = read_variants(args, read2.qname, chrsm, read2.pos \
                    + 1, read2_bases, read2.cigar, parse_md(get_md(read2), []))
                set_variants1 = set(variants1)
                set_variants2 = set(variants2)

                # find the variants each read in the pair share in common
                same_variants = set_variants1.intersection(set_variants2)

                if args.primercheck:
                    # the read we want for the forward primer is usually the
                    # second read, but this is not always the case, so we check
                    # the starting position of the second read to see if it
                    # matches the expected primer start position
                    if start - len(primer_sequence[block_info[3]]) == \
                    read2.pos + 1:
                        reverse_check = check_reverse_primer(primer_sequence, \
                            block_info, read1_bases, args.gap_penalty)
                        forward_check = check_forward_primer(primer_sequence, \
                            block_info, read2_bases, args.gap_penalty)
                    # if that didn't match, we also check the first read
                    elif start - len(primer_sequence[block_info[3]]) == \
                    read1.pos + 1:
                        forward_check = check_forward_primer(primer_sequence, \
                            block_info, read1_bases, args.gap_penalty)
                        reverse_check = check_reverse_primer(primer_sequence, \
                            block_info, read2_bases, args.gap_penalty)
                    # if neither matched, we discard the read pair as at least
                    # one of the reads should begin at the expected primer start
                    # position
                    else:
                        logging.warning("Read {} discarded due to unexpected \
start position".format(read_name))
                        discard = 2
                        num_unexpected += 1
                        num_pairs -= 1

                    if discard == 0:
                        forward_score = forward_check
                        reverse_score = reverse_check

                        scores = record_scores(forward_score, reverse_score, \
                            forward_scores, reverse_scores, total_scores)

                        forward_scores, reverse_scores, total_scores = \
                        scores[0], scores[1], scores[2]

                    if forward_score > args.primerthresh or reverse_score \
                    > args.primerthresh:
                        discard = 1

                if discard == 0:
                    for var in same_variants:
                        # only consider variants within the bounds of the block
                        if var.pos >= start and var.pos <= end:
                            if var in block_vars:
                                block_vars[var] += 1
                            else:
                                block_vars[var] = 1
                elif discard == 1:
                    # one of the reads had a primer sequence which was too
                    # different from what we expected (based on argument
                    # primerthresh)
                    # both reads were discarded as a result
                    logging.warning("Read {} discarded due to variant primer \
sequence".format(read_name))
                    num_pairs -= 1
                    num_discards += 1
            else:
                logging.warning("Read {} with more than 2".format(read_name))
        if args.primercheck:
            logging.warning("Number of read pairs discarded due to unexpected \
start position: {}".format(num_unexpected))
            logging.warning("Number of read pairs discarded due to unexpected \
primer sequence: {}".format(num_discards))
            logging.info("Number of acceptable read pairs remaining in block: \
{}".format(num_pairs))
        else:
            logging.info("Number of read pairs in block: {}".format(num_pairs))
        logging.info("Number of variants found in block: {}"\
            .format(len(block_vars)))

        if args.primercheck:
            write_block_data(block_info, forward_scores, reverse_scores, data)

        for var in block_vars:
            num_vars = block_vars[var]
            proportion = float(num_vars) / num_pairs
            var.info.append("Sample=" + str(sample))
            var.info.append("NV=" + str(num_vars))
            var.info.append("NP=" + str(num_pairs))
            var.info.append("PCT=" + str('{:.2%}'.format(proportion)))
            if num_vars < args.absthresh:
                var.filter_reason = ''.join([nts(var.filter_reason), ";at"])
            if proportion < args.proportionthresh:
                var.filter_reason = ''.join([nts(var.filter_reason), ";pt"])
            write_variant(kept_variants_file, var, id_info, args)
        if args.primercheck:
            coverage_info.append((chrsm, start, end, num_pairs, num_discards \
                + num_unexpected))
        else:
            coverage_info.append((chrsm, start, end, num_pairs))
    coverage_filename = sample + '.coverage'

    if args.primercheck:
        write_total_data(sample, total_scores, data2)

    if args.coverdir is not None:
        coverage_filename = os.path.join(args.coverdir, coverage_filename)
    with open(coverage_filename, 'w') as coverage_file:
        if args.primercheck:
            coverage_file.write('chrsm\tblock_start\tblock_end\tnum_pairs\t\
num_pairs_discarded\n')
            for chrsm, start, end, num_pairs, num_discards in sorted(\
                coverage_info, key=itemgetter(3)):
                coverage_file.write('{}\t{}\t{}\t{}\t{}\n'.format(chrsm, \
                    start, end, num_pairs, num_discards))
        else:
            coverage_file.write('chrsm\tblock_start\tblock_end\tnum_pairs\n')
            for chrsm, start, end, num_pairs in sorted(coverage_info, \
                key=itemgetter(3)):
                coverage_file.write('{}\t{}\t{}\t{}\n'.format(chrsm, start, \
                    end, num_pairs))

def record_scores(forward_score, reverse_score, forward_scores, \
    reverse_scores, total_scores):
    """ Records the score data which can be written to .dat files for further
    analysis."""
    if forward_score in forward_scores:
        forward_scores[forward_score] += 1
    else:
        forward_scores[forward_score] = 1
    if forward_score in total_scores:
        total_scores[forward_score] += 1
    else:
        total_scores[forward_score] = 1
    if reverse_score in reverse_scores:
        reverse_scores[reverse_score] += 1
    else:
        reverse_scores[reverse_score] = 1
    if reverse_score in total_scores:
        total_scores[reverse_score] += 1
    else:
        total_scores[reverse_score] = 1
    return [forward_scores, reverse_scores, total_scores]

def write_block_data(block_info, forward_scores, reverse_scores, data):
    """ Write data to .dat file for forward and reverse primers."""
    data.write(block_info[3] + '\t')
    forward_total = sum(forward_scores.values())
    for mismatch in range(0, 10):
        if float(mismatch) in forward_scores.keys():
            data.write("{:.2%}".format(forward_scores[mismatch]/\
                float(forward_total)) + '\t')
        else:
            data.write("-" + '\t')
    data.write("\n")

    data.write(block_info[4] + '\t')
    reverse_total = sum(reverse_scores.values())
    for mismatch in range(0, 10):
        if float(mismatch) in reverse_scores.keys():
            data.write("{:.2%}".format(reverse_scores[mismatch]/\
                float(reverse_total)) + '\t')
        else:
            data.write("-" + '\t')
    data.write("\n")

def write_total_data(sample, total_scores, data2):
    """ Write data to .dat file for total primers."""
    data2.write("# " + sample + '\n')
    total2 = sum(total_scores.values())
    for mismatch in sorted(total_scores):
        if mismatch < 10:
            data2.write(str(mismatch) + '\t' + "{:.2%}"\
                .format(total_scores[mismatch]/float(total2)) + '\n')

def write_metadata(args, vcf_file):
    """ Write the opening lines of metadata to the vcf file."""
    vcf_file.write("##fileformat=VCFv4.2" + '\n')
    today = datetime.date.today()
    vcf_file.write("##fileDate=" + str(today)[:4] + str(today)[5:7] + \
        str(today)[8:] + '\n')
    vcf_file.write("##source=ROVER-PCR Variant Caller" + '\n')
    vcf_file.write("##INFO=<ID=Sample,Number=1,Type=String,Description=\
\"Sample Name\">" + '\n')
    vcf_file.write("##INFO=<ID=NV,Number=1,Type=Float,Description=\
\"Number of read pairs with variant\">" + '\n')
    vcf_file.write("##INFO=<ID=NP,Number=1,Type=Float,Description=\
\"Number of read pairs at POS\">" + '\n')
    vcf_file.write("##INFO=<ID=PCT,Number=1,Type=Float,Description=\
\"Percentage of read pairs at POS with variant\">" + '\n')
    vcf_file.write("##INFO=<ID=BS,Number=1,Type=String,Description=\
\"Context base cannot be determined as indel is located near a region \
not covered by the MD string\">" + '\n')
    vcf_file.write("##INFO=<ID=HC,Number=1,Type=String,Description=\
\"Context base cannot be determined due to \
hard clipping on the aligned sequence prior to indel event\">" + '\n')
    vcf_file.write("##INFO=<ID=SC,Number=1,Type=String,Description=\"Context \
base cannot be determined due to \
soft clipping on the aligned sequence prior to indel event\">" + '\n')
    vcf_file.write("##INFO=<ID=IP,Number=1,Type=String,Description=\
\"Misaligned or incorrect base sequence in primer region\">" + '\n')
    if args.qualthresh:
        vcf_file.write("##FILTER=<ID=qlt,Description=\"Variant has phred \
quality score below " + str(args.qualthresh) + "\">" + '\n')
    if args.absthresh:
        vcf_file.write("##FILTER=<ID=at,Description=\"Variant does not appear \
in at least " + str(args.absthresh) + " read pairs\">" + '\n')
    if args.proportionthresh:
        vcf_file.write("##FILTER=<ID=pt,Descroption=\"Variant does not appear \
in at least " + str(args.proportionthresh*100) \
+ "% of read pairs for the given region\">" + '\n')

# Proper tab separated column headings
OUTPUT_HEADER = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", \
    "FILTER", "INFO"])

def process_bams(args):
    block_coords = get_block_coords(args.primers)
    primer_sequence = {}
    vcf_reader = 0
    # a dictionary of primers and their sequences
    if args.primercheck:
        primer_info = get_primer_sequence(args.primercheck)
        for primer in primer_info:
            primer_sequence[primer[0]] = primer[1]
    with open(args.out, "w") as kept_variants_file:
        graph_data = open("data.dat", "w")
        graph_total_data = open("data2.dat", "w")
        write_metadata(args, kept_variants_file)
        if args.id_info:
            vcf_reader = vcf.Reader(filename=args.id_info)
        kept_variants_file.write(OUTPUT_HEADER + '\n')
        for bam_filename in args.bams:
            base = os.path.basename(bam_filename)
            sample = base.split('.')
            if len(sample) > 0:
                sample = sample[0]
            else:
                exit('Cannot deduce sample name from bam filename {}'\
                    .format(bam_filename))
            with pysam.Samfile(bam_filename, "rb") as bam:
                logging.info("Processing bam file {}".format(bam_filename))
                process_blocks(args, kept_variants_file, bam, sample, \
                    block_coords, primer_sequence, graph_data, \
                    graph_total_data, vcf_reader)

def main():
    args = parse_args()
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('Program started')
    logging.info('Command line: {0}'.format(' '.join(sys.argv)))
    process_bams(args)

if __name__ == '__main__':
    main()
