#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
# from pyfaidx import Fasta
import datetime
import logging
import sys
import pysam
import re
import os
from operator import itemgetter
import csv
from version import rover_version
from itertools import (izip, chain, repeat)

# proportion of block which must be overlapped by read 
default_minimum_read_overlap_block = 0.9
default_proportion_threshold = 0.05
default_absolute_threshold = 2

def parse_args():
    "Consider mapped reads to amplicon sites"

    parser = ArgumentParser(description="Consider mapped reads to amplicon sites")
    parser.add_argument('--reference', type=str, 
	 help='File name of reference DNA sequence in FASTA format.')
    parser.add_argument(
    '--version', action='version', version='%(prog)s ' + rover_version)
    parser.add_argument(
        '--primers', type=str, required=True,
        help='File name of primer coordinates in TSV format.')
    parser.add_argument(
        '--overlap', type=float, default=default_minimum_read_overlap_block,
        help='Minimum proportion of block which must be overlapped by a read. '
             'Defaults to {}.'.format(default_minimum_read_overlap_block))
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
    parser.add_argument('--qualthresh', metavar='N', type=int,
        help='Minimum base quality score (phred).')
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
def read_variants(args, name, chr, pos, aligned_bases, cigar, md):
    cigar_orig = cigar
    md_orig = md
    seq_index = 0
    result = []
    context = None    

#    while cigar and seq_index < len(aligned_bases):
#        cigar_code, cigar_segment_extent = cigar[0]
#	if cigar_code == 0:
#	    # Cigar Match
#	    if fasta[ref + 1].upper() == aligned_bases[seq_index].base:
#	    	pos += cigar_segment_extent
#		ref += cigar_segment_extent
#		seq_index += cigar_segment_extent
#		cigar = cigar[1:]
#	    else:
#		seq_base_qual = aligned_bases[seq_index]
#		seq_base = seq_base_qual.base
#		if (args.qualthresh is None) or (seq_base_qual.qual >= args.qualthresh):
#	            result.append(SNV(chr, pos, fasta[ref + 1], seq_base, seq_base_qual.qual, None))
#		else:
#		    result.append(SNV(chr, pos, fasta[ref + 1], seq_base, seq_base_qual.qual, ";qlt"))
#		cigar = [(cigar_code, cigar_segment_extent - 1)] + cigar[1:]
#		seq_index += 1
#		ref += 1
#		pos += 1
#	elif cigar_code == 1:
#	    extra_bases_quals = aligned_bases[(seq_index):(seq_index + cigar_segment_extent)]
#	    extra_bases = ''.join([b.base for b in extra_bases_quals])
#	    context = fasta[ref - 1]
#	    if (args.qualthresh is None) or all([b.qual >= args.qualthresh for b in extra_bases_quals]):
#	        result.append(Insertion(chr, pos, extra_bases, 15, None, context))
#	    else:
#		result.append(Insertion(chr, pos, extra_bases, 15, ";qlt", context))
#	    cigar = cigar[1:]
#	    seq_index += cigar_segment_extent
#	elif cigar_code == 2:
#	    deleted_bases = fasta[ref:(ref + cigar_segment_extent)]
#	    context = fasta[ref - 1]
#	    seq_base = aligned_bases[seq_index]
#	    if seq_base.qual >= args.qualthresh:
#               result.append(Deletion(chr, pos, deleted_bases, 15, None, context))
#	    else:
#		result.append(Deletion(chr, pos, deleted_bases, 15, ";qlt", context))
#	    pos += cigar_segment_extent
#	    ref += cigar_segment_extent
#	    cigar = cigar[1:]
#	else:
#	    logging.info("unexpected cigar code {}".format(cigar_orig))
#	    exit()
#    return result
    
    while cigar and md:
	cigar_code, cigar_segment_extent = cigar[0]
	next_md = md[0]
	if cigar_code == 0:
	    if isinstance(next_md, MD_match):
                # MD match
		if next_md.size >= cigar_segment_extent:
                    next_md.size -= cigar_segment_extent
                    if next_md.size  == 0:
		        md = md[1:]
                    context = aligned_bases[seq_index + cigar_segment_extent - 1].base
		    cigar = cigar[1:]
                    pos += cigar_segment_extent
                    seq_index += cigar_segment_extent
                else:
                    # next_md.size < cigar_segment_extent
                    cigar = [(cigar_code, cigar_segment_extent - next_md.size)] + cigar[1:]
		    context = aligned_bases[seq_index + next_md.size - 1].base
		    md = md[1:]
                    pos += next_md.size
                    seq_index += next_md.size
            elif isinstance(next_md, MD_mismatch):
		# MD mismatch
                seq_base_qual = aligned_bases[seq_index]
                # check if the read base is above the minimum quality score
                if (args.qualthresh is None) or (seq_base_qual.qual >= args.qualthresh):
                    seq_base = seq_base_qual.base
                    result.append(SNV(chr, pos, next_md.ref_base, seq_base, seq_base_qual.qual, None))
		else:
		    seq_base = seq_base_qual.base
		    result.append(SNV(chr, pos, next_md.ref_base, seq_base, seq_base_qual.qual, ";qlt"))
                cigar = [(cigar_code, cigar_segment_extent - 1)] + cigar[1:]
                context = next_md.ref_base
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
	    seq_bases_quals = aligned_bases[seq_index:seq_index + cigar_segment_extent]
            seq_bases = ''.join([b.base for b in seq_bases_quals])
            # check that all the bases are above the minimum quality threshold
            if (args.qualthresh is None) or all([b.qual >= args.qualthresh for b in seq_bases_quals]):
                result.append(Insertion(chr, pos, seq_bases, '-', None, context))
	    else:
	        result.append(Insertion(chr, pos, seq_bases, '-', ";qlt", context))
	    cigar = cigar[1:]
            seq_index += cigar_segment_extent
	    # pos does not change
        elif cigar_code == 2:
            # Deletion
            if isinstance(next_md, MD_deletion):
                seq_base = aligned_bases[seq_index]
		if seq_base.qual >= args.qualthresh:
		    result.append(Deletion(chr, pos, next_md.ref_bases, '-', None, context))
                else:
		    result.append(Deletion(chr, pos, next_md.ref_bases, '-', ";qlt", context))
		context = next_md.ref_bases[-1]
		md = md[1:]
                cigar = cigar[1:]
                pos += cigar_segment_extent
                # seq_index does not change
            else:
                logging.info("Non del MD in Del Cigar".format(md_orig, cigar_orig))
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
	    logging.info("unexpected cigar code {}".format(cigar_orig))
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
    '''Take a list of DNA bases and a corresponding list of quality scores
    and return a list of Base objects where the base and score are
    paired together.'''
    num_bases = len(bases)
    num_qualities = len(qualities)
    if num_bases <= num_qualities:
        return [Base(b, ascii_to_phred(q)) for (b, q) in izip(bases, qualities)]
    else:
        logging.warning("In read {} fewer quality scores {} than bases {}"
            .format(name, num_qualities, num_bases))
        # we have fewer quality scores than bases
        # pad the end with 0 scores (which is ord('!') - 33)
        return [Base(b, ascii_to_phred(q))
            for (b, q) in izip(bases, chain(qualities, repeat('!')))]

# a DNA base paired with its quality score
class Base(object):
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
    def __init__(self, chr, pos, ref_base, seq_base, qual, filter):
        self.chr = chr
        self.pos = pos
        self.ref_base = ref_base
        self.seq_base = seq_base
	self.qual = qual
	self.filter = filter
	self.info = []
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
    def fil(self):
	if self.filter is None:
	    return "PASS"
	else:
	    return self.filter[1:]
    def position(self):
	return self.pos
    def quality(self):
	if self.qual == '-':
	    return '.'
	else:
	    return str('{:.2%}'.format(self.qual/100.0))

class Insertion(object):
    # bases are represented just as DNA strings
    def __init__(self, chr, pos, inserted_bases, qual, filter, context):
        self.chr = chr
        self.pos = pos
        self.inserted_bases = inserted_bases
	self.qual = qual
	self.filter = filter
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
        return self.context
    def alt(self):
	return self.context + self.inserted_bases
    def fil(self):
        if self.filter is None:
	    return "PASS"
	else:
	    return self.filter[1:]
    def position(self):
	return self.pos - 1
    def quality(self):
	if self.qual == '-':
	    return '.'
	else:
	    return self.qual

class Deletion(object):
    # bases are represented just as DNA strings
    def __init__(self, chr, pos, deleted_bases, qual, filter, context):
        self.chr = chr
        self.pos = pos
        self.deleted_bases = deleted_bases
	self.qual = qual
	self.filter = filter
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
	return self.context + self.deleted_bases
    def alt(self):
	return self.context
    def fil(self):
	if self.filter is None:
	    return "PASS"
	else:
	    return self.filter[1:]
    def position(self):
	return self.pos - 1
    def quality(self):
	if self.qual == '-':
	    return '.'
	else:
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
    '''Compute the proportion of the block that is overlapped by the read

          block_start               block_end
               |-------------------------|

        ^---------------------------^
    read.pos                      read_end

               |--------------------|
         overlap_start        overlap_end

    '''
    read_end = read.pos + read.rlen - 1
    if read.rlen <= 0:
        # read is degenerate, zero length
        # treat it as no overlap
        logging.warn("Degenerate read: {}, length: {}".format(read.qname, read.rlen))
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

def write_variant(file, variant, sample, args):
    """
    global repeats
    global line
    global kept_variants_file
    if (sample, variant.pos) in repeats.keys():
	print "repeated " + sample + '' + str(variant.pos)
	kept_variants_file.close()
	kept_variants_file = open(args.out, "r+")
	line_offset = []
	offset = 0
	for row in kept_variants_file:
	    line_offset.append(offset)
	    offset += len(row)
	kept_variants_file.seek(0)
	# file.seek(line_offset[repeats[(sample, variant.pos, variant.ref(), variant.alt())]])
	kept_variants_file.seek(line_offset[repeats[(sample, variant.pos)]] - 1)
	kept_variants_file.write("Repeat would go here." + '\n')
	kept_variants_file.close()
	kept_variants_file = open(args.out, "w")
    else:
	# repeats[(sample, variant.pos, variant.ref(), variant.alt())] = line
	repeats[(sample, variant.pos)] = line
    """
    file.write('\t'.join([variant.chr[3:], str(variant.position()), \
'.', variant.ref(), variant.alt(), variant.quality(), variant.fil(), ';'.join(variant.info)]) + '\n')

def nts(s):
    # Turns None into an empty string
    if s is None:
	return ''
    return str(s)

def process_blocks(args, kept_variants_file, bam, sample, block_coords):
    coverage_info = []
    for block_info in block_coords:
        chr, start, end = block_info[:3]
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
                #print(read1.query)
                #print([ord(x) - 33 for x in read1.qqual])
                #print(read2.query)
                #print([ord(x) - 33 for x in read2.qqual])
                #exit()
                read1_bases = make_base_seq(read1.qname, read1.query, read1.qqual)
                read2_bases = make_base_seq(read2.qname, read2.query, read2.qqual)
		variants1 = read_variants(args, read1.qname, chr, read1.pos + 1, read1_bases, read1.cigar, parse_md(get_MD(read1), []))
                variants2 = read_variants(args, read2.qname, chr, read2.pos + 1, read2_bases, read2.cigar, parse_md(get_MD(read2), []))
                set_variants1 = set(variants1)
                set_variants2 = set(variants2)
                # find the variants each read in the pair share in common
                same_variants = set_variants1.intersection(set_variants2)
                for var in same_variants:
                    # only consider variants within the bounds of the block
                    if var.pos >= start and var.pos <= end:
                        if var in block_vars:
                            if var.qual == '-':
				block_vars[var] = (block_vars[var][0] + 1, var.qual)
			    else:
				block_vars[var] = tuple(map(sum, zip(block_vars[var], (1, var.qual))))
                        else:
                            block_vars[var] = (1, var.qual)
            else:
                logging.warning("read {} with more than 2".format(read_name))
        logging.info("number of read pairs in block: {}".format(num_pairs))
        logging.info("number of variants found in block: {}".format(len(block_vars)))
        for var in block_vars:
            num_vars = block_vars[var][0]
            proportion = float(num_vars) / num_pairs
            proportion_str = "{:.2f}".format(proportion)
  	    var.info.append("Sample=" + str(sample))
	    var.info.append("NV=" + str(num_vars))
	    var.info.append("NP=" + str(num_pairs))
	    var.info.append("PCT=" + str('{:.2%}'.format(proportion)))
	    if num_vars < args.absthresh:
		var.filter = ''.join([nts(var.filter), ";at"])
	    if proportion < args.proportionthresh:
		var.filter = ''.join([nts(var.filter), ";pt"])
	    if block_vars[var][1] == '-':
		var.qual = '-'
	    else:
		var.qual = (block_vars[var][1])/float(num_vars)
	    write_variant(kept_variants_file, var, sample, args)
        coverage_info.append((chr, start, end, num_pairs))
    coverage_filename = sample + '.coverage'
    if args.coverdir is not None:
        coverage_filename = os.path.join(args.coverdir, coverage_filename)
    with open(coverage_filename, 'w') as coverage_file:
        coverage_file.write('chr\tblock_start\tblock_end\tnum_pairs\n')
        for chr, start, end, num_pairs in sorted(coverage_info, key=itemgetter(3)):
            coverage_file.write('{}\t{}\t{}\t{}\n'.format(chr, start, end, num_pairs))

def write_metadata(args, file):
    file.write("##fileformat=VCFv4.2" + '\n')
    today = datetime.date.today()
    file.write("##fileDate=" + str(today)[:4] + str(today)[5:7] + str(today)[8:] + '\n')
    file.write("##source=ROVER-PCR Variant Caller" + '\n')
    if args.reference:
	file.write("##reference=file:///" + str(args.reference) + '\n')
    # file.write("##contig=" + '\n')
    # file.write("##phasing=" + '\n')
    file.write("##INFO=<ID=Sample,Number=1,Type=String,Description=\"Sample Name\">" + '\n')
    file.write("##INFO=<ID=NV,Number=1,Type=Float,Description=\"Number of read pairs with variant\">" + '\n')
    file.write("##INFO=<ID=NP,Number=1,Type=Float,Description=\"Number of read pairs at POS\">" + '\n')
    file.write("##INFO=<ID=PCT,Number=1,Type=Float,Description=\"Percentage of read pairs at POS with variant\">" + '\n')
    file.write("##INFO=<ID=BS,Number=1,Type=String,Description=\"Context base cannot be determined as indel is located near start of block\">" + '\n')
    file.write("##INFO=<ID=HC,Number=1,Type=String,Description=\"Context base cannot be determined due to \
hard clipping on the aligned sequence prior to indel event\">" + '\n')
    file.write("##INFO=<ID=SC,Number=1,Type=String,Description=\"Context base cannot be determined due to \
soft clipping on the aligned sequence prior to indel event\">" + '\n')
    # file.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples with Data\">" + '\n')
    # file.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" + '\n')
    if args.qualthresh: 
        file.write("##FILTER=<ID=qlt,Description=\"Variant has phred quality score below " + str(args.qualthresh) + "\">" + '\n')
    if args.absthresh:
	file.write("##FILTER=<ID=at,Description=\"Variant does not appear in at least " + str(args.absthresh) + " read pairs\">" + '\n')
    if args.proportionthresh:
	file.write("##FILTER=<ID=pt,Descroption=\"Variant does not appear in at least " + str(args.proportionthresh*100) \
		+ "% of read pairs for the given region\">" + '\n')

# Extra formatting applied to column headings so that everything lines up
output_header = '\t'.join(["#CHROM", "POS", '', "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

# Proper tab separated column headings
# output_header = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

def process_bams(args):
    block_coords = get_block_coords(args.primers)
    # with open(args.out, "w") as kept_variants_file, \
    #      open(args.out + '.binned', "w") as binned_variants_file:
    with open(args.out, "w") as kept_variants_file:
	write_metadata(args, kept_variants_file)
	# write_metadata(args, binned_variants_file)
	kept_variants_file.write(output_header + '\n')
        # binned_variants_file.write(output_header + '\n')
	# ref_dict = Fasta(args.reference)
	for bam_filename in args.bams:
            base = os.path.basename(bam_filename)
            sample = base.split('.')
            if len(sample) > 0:
                sample = sample[0]
            else:
                exit('Cannot deduce sample name from bam filename {}'.format(bam_filename))
            with pysam.Samfile(bam_filename, "rb") as bam:
                logging.info("processing bam file {}".format(bam_filename))
                process_blocks(args, kept_variants_file, bam, sample, block_coords)

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
